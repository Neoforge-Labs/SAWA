"""Country x year x drug-class resistance view for the dashboard.

The aggregation itself is `paper1.fig02_aggregate.aggregate_resistance_landscape`,
reused so the page and the manuscript figure cannot disagree about what a cell is. This
module adds the three things the page needs on top of it and the manuscript figure does
not:

* a roll-up over subtype, because four dimensions do not fit one page figure;
* a restriction to the targeted region, the countries the ingest search was aimed at;
* an accounting of every scored accession that is *not* in the result, by reason, so the
  page can state what it left out instead of letting the omission pass unnoticed.

The upstream aggregation drops an accession with a null country, year or subtype without
saying so. The accounting here makes each of those drops countable, and the per-class
check below refuses a result in which the cells and the accounting disagree.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass

import polars as pl

from hivdrwatch.dashboard.export.aggregate import UnknownDrugClassError
from hivdrwatch.ingest.schema import EAST_AFRICA_COUNTRIES
from hivdrwatch.paper1.fig02_aggregate import MIN_CELL_N, aggregate_resistance_landscape
from hivdrwatch.paper1.style import DRUG_CLASS_ORDER

logger = logging.getLogger(__name__)

# The countries the ingest search targeted, in the search's own order. Taken from the ingest
# module so "in region" cannot drift from what was searched for.
TARGETED_COUNTRIES: tuple[str, ...] = tuple(EAST_AFRICA_COUNTRIES)

_CLASS_RANK: dict[str, int] = {name: i for i, name in enumerate(DRUG_CLASS_ORDER)}

CURATED_COLUMNS: tuple[str, ...] = ("accession", "country", "collection_year", "subtype")


class DuplicateCuratedAccessionError(ValueError):
    """An accession appears twice in the curated table, so the join would double its counts."""


class LandscapeAccountingError(ValueError):
    """The landscape cells and the accession accounting disagree about who is in the figure."""


@dataclass(frozen=True)
class LandscapeAccounting:
    """Every scored accession falls in exactly one of these, in this order of precedence."""

    n_in_heatmap: int
    n_out_of_region: int
    n_no_country: int
    n_no_year: int
    # Not in the curated table at all, or no subtype: the upstream aggregation needs one.
    n_other: int

    @property
    def n_total(self) -> int:
        return (
            self.n_in_heatmap
            + self.n_out_of_region
            + self.n_no_country
            + self.n_no_year
            + self.n_other
        )


def _require_columns(curated: pl.DataFrame) -> None:
    missing = set(CURATED_COLUMNS) - set(curated.columns)
    if missing:
        raise ValueError(f"curated missing required columns: {sorted(missing)}")


def _placement(
    dr_results: pl.DataFrame, curated: pl.DataFrame, countries: Sequence[str]
) -> pl.DataFrame:
    """One row per scored accession with the reason it is in or out of the heatmap."""
    _require_columns(curated)
    duplicated = curated["accession"].is_duplicated()
    if duplicated.any():
        raise DuplicateCuratedAccessionError(
            f"{int(duplicated.sum())} curated rows share an accession with another row"
        )
    scored = dr_results.select("accession").unique()
    in_curated = curated.select(*CURATED_COLUMNS).with_columns(_in_curated=pl.lit(True))
    return (
        scored.join(in_curated, on="accession", how="left")
        .with_columns(
            placement=pl.when(pl.col("_in_curated").is_null())
            .then(pl.lit("other"))
            .when(pl.col("country").is_null())
            .then(pl.lit("no_country"))
            .when(~pl.col("country").is_in(list(countries)))
            .then(pl.lit("out_of_region"))
            .when(pl.col("collection_year").is_null())
            .then(pl.lit("no_year"))
            .when(pl.col("subtype").is_null())
            .then(pl.lit("other"))
            .otherwise(pl.lit("in_heatmap"))
        )
        .select("accession", "placement")
    )


def landscape_accounting(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    countries: Sequence[str] = TARGETED_COUNTRIES,
) -> LandscapeAccounting:
    """Count the scored accessions in the heatmap and, by reason, those that are not."""
    return _account(_placement(dr_results, curated, countries))


def _account(placement: pl.DataFrame) -> LandscapeAccounting:
    counts = dict(placement["placement"].value_counts().iter_rows())
    return LandscapeAccounting(
        n_in_heatmap=counts.get("in_heatmap", 0),
        n_out_of_region=counts.get("out_of_region", 0),
        n_no_country=counts.get("no_country", 0),
        n_no_year=counts.get("no_year", 0),
        n_other=counts.get("other", 0),
    )


def _check_drug_classes(landscape: pl.DataFrame) -> None:
    unknown = sorted(
        set(landscape["drug_class"].drop_nulls().unique()) - set(DRUG_CLASS_ORDER)
    )
    if unknown:
        raise UnknownDrugClassError(
            f"landscape: drug classes {unknown} are not in the canonical order {DRUG_CLASS_ORDER}"
        )


def rollup_over_subtype(
    landscape: pl.DataFrame,
    *,
    countries: Sequence[str] = TARGETED_COUNTRIES,
    min_cell_n: int = MIN_CELL_N,
) -> pl.DataFrame:
    """Sum the subtype-stratified landscape within each (country, year, drug class).

    Output: [country, collection_year, drug_class, n_sequences, n_resistant, frac_resistant],
    restricted to `countries`. The floor is applied again to the summed count, not inherited:
    two subtype cells that are each below it can sum to a cell above it, and `frac_resistant`
    of a sum is not a sum of fractions. Counts are kept below the floor, as upstream does.
    """
    _check_drug_classes(landscape)
    country_rank = {code: i for i, code in enumerate(countries)}
    return (
        landscape.filter(pl.col("country").is_in(list(countries)))
        .group_by(["country", "collection_year", "drug_class"])
        .agg(
            n_sequences=pl.col("n_sequences").sum(),
            n_resistant=pl.col("n_resistant").sum(),
        )
        .with_columns(
            frac_resistant=pl.when(pl.col("n_sequences") >= min_cell_n)
            .then(pl.col("n_resistant") / pl.col("n_sequences"))
            .otherwise(None)
        )
        .sort(
            pl.col("drug_class").replace_strict(_CLASS_RANK, return_dtype=pl.UInt8),
            pl.col("country").replace_strict(country_rank, return_dtype=pl.UInt8),
            "collection_year",
        )
        .select(
            "country",
            "collection_year",
            "drug_class",
            "n_sequences",
            "n_resistant",
            "frac_resistant",
        )
    )


def _check_cells_match_accounting(
    cells: pl.DataFrame, placement: pl.DataFrame, dr_results: pl.DataFrame
) -> None:
    """Per class, the cells must hold exactly the in-heatmap accessions scored for it."""
    included = placement.filter(pl.col("placement") == "in_heatmap").select("accession")
    expected = (
        dr_results.select("accession", "drug_class")
        .unique()
        .join(included, on="accession", how="inner")
        .group_by("drug_class")
        .agg(n=pl.len())
    )
    actual = cells.group_by("drug_class").agg(n=pl.col("n_sequences").sum())
    want = dict(expected.iter_rows())
    got = dict(actual.iter_rows())
    if want != got:
        raise LandscapeAccountingError(
            f"landscape cells hold {got} sequences per class but the accounting expects {want}"
        )


def resistance_landscape(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    countries: Sequence[str] = TARGETED_COUNTRIES,
    min_cell_n: int = MIN_CELL_N,
) -> tuple[pl.DataFrame, LandscapeAccounting]:
    """The page's landscape cells and the accounting of who is not in them.

    `dr_results` must already be deduplicated (one row per accession and drug).
    """
    landscape, _ = aggregate_resistance_landscape(dr_results, curated, min_cell_n=min_cell_n)
    cells = rollup_over_subtype(landscape, countries=countries, min_cell_n=min_cell_n)
    placement = _placement(dr_results, curated, countries)
    _check_cells_match_accounting(cells, placement, dr_results)
    accounting = _account(placement)
    logger.info(
        "landscape: %d cells over %d accessions; excluded %d out of region, %d without a "
        "country, %d without a collection year, %d other",
        cells.height,
        accounting.n_in_heatmap,
        accounting.n_out_of_region,
        accounting.n_no_country,
        accounting.n_no_year,
        accounting.n_other,
    )
    return cells, accounting
