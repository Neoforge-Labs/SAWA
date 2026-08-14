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
        n_dr, n_curated, joined.height,
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
        counts.n_output_cells, counts.n_cells_above_min, min_cell_n,
    )
    return landscape, counts
