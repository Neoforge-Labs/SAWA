"""Pure data aggregation for Fig 2 — resistance landscape heatmap.

Given the DR results parquet (Sierra output) and the curated silver
parquet (metadata: country + subtype + year), produce a long-format
table with one row per (country, year, drug_class, subtype) giving:

  * n_sequences: how many sequences fell in that cell
  * n_resistant: how many had at least one drug in that class scored
    at Sierra level `high` or `intermediate`
  * frac_resistant: n_resistant / n_sequences
  * drug_class

That table is the direct input to the heatmap renderer in the notebook.
Kept in a separate module so the aggregation is unit-testable without
the matplotlib import cost.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass

import polars as pl

logger = logging.getLogger(__name__)

# Sierra level values that count as "resistant" for prevalence stats.
# Matches the binary target used in the DR feature builder and Sierra's
# clinical-decision boundary.
RESISTANT_LEVELS: frozenset[str] = frozenset({"high", "intermediate"})

# Minimum n per cell before we report prevalence. Cells below this get
# `frac_resistant=None` so the heatmap can render them as an explicit
# "not enough data" grey.
MIN_CELL_N: int = 5


@dataclass(frozen=True)
class AggregationCounts:
    """High-level counts for the "Fig 2 built OK" summary log."""

    n_dr_rows: int
    n_curated_rows: int
    n_joined_rows: int
    n_output_cells: int
    n_cells_above_min: int


def aggregate_resistance_landscape(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    min_cell_n: int = MIN_CELL_N,
) -> tuple[pl.DataFrame, AggregationCounts]:
    """Return (long-format landscape frame, aggregation counts).

    Inputs (only the columns we consume are listed; others ignored):
      dr_results: [accession, drug_class, drug, level]
      curated:    [accession, country, collection_year, subtype]

    Output columns:
      [country, collection_year, drug_class, subtype,
       n_sequences, n_resistant, frac_resistant]

    Sparse cells (below `min_cell_n`) keep the counts but get
    `frac_resistant=None`.
    """
    for name, df, cols in (
        ("dr_results", dr_results, {"accession", "drug_class", "level"}),
        ("curated", curated, {"accession", "country", "collection_year", "subtype"}),
    ):
        missing = cols - set(df.columns)
        if missing:
            raise ValueError(f"{name} missing required columns: {sorted(missing)}")

    n_dr = dr_results.height
    n_curated = curated.height

    # Reduce Sierra results to per-(accession, drug_class) worst call.
    # "Worst" defined as: is the sequence resistant to *any* drug in
    # this class? That's the clinically-relevant unit — a sequence is
    # counted as NNRTI-resistant if resistant to any NNRTI drug.
    per_seq_class = (
        dr_results.with_columns(
            pl.col("level").is_in(list(RESISTANT_LEVELS)).alias("_is_resistant"),
        )
        .group_by(["accession", "drug_class"])
        .agg(pl.col("_is_resistant").any().alias("is_resistant"))
    )

    # Join country + year + subtype from curated. Inner join drops
    # accessions we can't stratify.
    meta_cols = ["accession", "country", "collection_year", "subtype"]
    joined = per_seq_class.join(
        curated.select(meta_cols),
        on="accession",
        how="inner",
    ).filter(
        pl.col("country").is_not_null()
        & pl.col("collection_year").is_not_null()
        & pl.col("subtype").is_not_null()
    )

    logger.info(
        "resistance landscape: %d DR rows × %d curated rows → %d joined per-sequence-per-class rows",
        n_dr,
        n_curated,
        joined.height,
    )

    landscape = (
        joined.group_by(["country", "collection_year", "drug_class", "subtype"])
        .agg(
            n_sequences=pl.len(),
            n_resistant=pl.col("is_resistant").sum(),
        )
        .with_columns(
            frac_resistant=pl.when(pl.col("n_sequences") >= min_cell_n)
            .then(pl.col("n_resistant") / pl.col("n_sequences"))
            .otherwise(None),
        )
        .sort(["country", "collection_year", "drug_class", "subtype"])
    )

    counts = AggregationCounts(
        n_dr_rows=n_dr,
        n_curated_rows=n_curated,
        n_joined_rows=joined.height,
        n_output_cells=landscape.height,
        n_cells_above_min=int(landscape["frac_resistant"].is_not_null().sum()),
    )
    logger.info(
        "landscape: %d cells (%d above min_cell_n=%d)",
        counts.n_output_cells,
        counts.n_cells_above_min,
        min_cell_n,
    )
    return landscape, counts


@dataclass(frozen=True)
class OmittedSubtypes:
    """Observed subtypes a figure's display order does not cover."""

    labels: tuple[str, ...]
    n_cells: int
    # One sequence contributes once per drug class, so this counts
    # sequence-by-drug-class observations, not sequences.
    n_observations: int

    def __bool__(self) -> bool:
        return bool(self.labels)

    def describe(self) -> str:
        if not self.labels:
            return "no observed subtype is omitted by the display order"
        return (
            f"{len(self.labels)} observed subtype(s) omitted by the display "
            f"order: {', '.join(self.labels)} "
            f"({self.n_cells} cells, {self.n_observations} sequence-class observations)"
        )


def subtypes_omitted_by(landscape: pl.DataFrame, reported: Iterable[str]) -> OmittedSubtypes:
    """Which observed subtypes `reported` leaves out of the figure.

    The aggregation above groups by every subtype present in the data, but
    the heatmap renders one column per entry in a fixed display order and
    skips any label not in it. That was invisible while the only omitted
    labels were `unknown` and a few singletons. It stopped being harmless
    once ingest started resolving CRFs and parent clades correctly: `A` and
    `CRF01_AE` together are ~3,900 curated rows that no display order
    currently has a slot for.

    Call this before rendering, so a stratum cannot leave the figure
    without anyone noticing. Whether to widen the display order, roll
    parent clades into their sub-clades for display, or give CRFs their own
    column is a presentation decision — this only makes the omission
    visible.
    """
    if "subtype" not in landscape.columns:
        raise ValueError("landscape frame missing 'subtype' column")

    reported_set = set(reported)
    omitted = landscape.filter(~pl.col("subtype").is_in(list(reported_set)))
    labels = tuple(sorted(omitted["subtype"].unique().to_list()))
    n_observations = (
        int(omitted["n_sequences"].sum()) if "n_sequences" in omitted.columns else 0
    )
    result = OmittedSubtypes(
        labels=labels, n_cells=omitted.height, n_observations=n_observations
    )
    if result:
        logger.warning("fig02: %s", result.describe())
    return result
