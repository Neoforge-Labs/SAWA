"""Tests for the dashboard's country x year x drug-class view.

The aggregation is reused from the manuscript figure, so these tests are about what the
dashboard adds: the roll-up over subtype, the region restriction, and the accounting of who
is left out. Each test names the break it would catch. Frames are built in code with a small
`min_cell_n` where the floor matters, so a cell can be above it without hundreds of rows.
"""

from __future__ import annotations

import subprocess
import sys

import polars as pl
import pytest

from hivdrwatch.dashboard.export import landscape
from hivdrwatch.dashboard.export.aggregate import UnknownDrugClassError
from hivdrwatch.dashboard.export.landscape import (
    TARGETED_COUNTRIES,
    DuplicateCuratedAccessionError,
    LandscapeAccountingError,
    landscape_accounting,
    resistance_landscape,
    rollup_over_subtype,
)
from hivdrwatch.paper1 import fig02_aggregate

_LANDSCAPE_SCHEMA = {
    "country": pl.String,
    "collection_year": pl.Int64,
    "drug_class": pl.String,
    "subtype": pl.String,
    "n_sequences": pl.UInt32,
    "n_resistant": pl.UInt32,
    "frac_resistant": pl.Float64,
}


def _cells(*rows: tuple[str, int, str, str, int, int]) -> pl.DataFrame:
    """Subtype-stratified cells as the manuscript aggregation emits them (floor not applied)."""
    return pl.DataFrame(
        [(*r, None) for r in rows],
        schema=_LANDSCAPE_SCHEMA,
        orient="row",
    )


def _roll(*rows: tuple[str, int, str, str, int, int], floor: int = 5) -> pl.DataFrame:
    return rollup_over_subtype(_cells(*rows), min_cell_n=floor)


def _cell(df: pl.DataFrame, country: str, year: int, cls: str) -> dict[str, object]:
    (row,) = df.filter(
        (pl.col("country") == country)
        & (pl.col("collection_year") == year)
        & (pl.col("drug_class") == cls)
    ).to_dicts()
    return row


# --- roll-up ----------------------------------------------------------------------------


def test_rollup_sums_counts_across_subtypes() -> None:
    out = _roll(("KE", 2010, "PI", "A1", 7, 1), ("KE", 2010, "PI", "D", 5, 2))
    assert _cell(out, "KE", 2010, "PI") == {
        "country": "KE",
        "collection_year": 2010,
        "drug_class": "PI",
        "n_sequences": 12,
        "n_resistant": 3,
        "frac_resistant": 0.25,
    }


def test_rollup_recomputes_the_fraction_instead_of_averaging_the_subtype_fractions() -> None:
    # 10/10 and 0/90: the mean of the two fractions is 0.5, the share of the 100 is 0.1.
    out = _roll(("KE", 2010, "PI", "A1", 10, 10), ("KE", 2010, "PI", "D", 90, 0))
    assert _cell(out, "KE", 2010, "PI")["frac_resistant"] == pytest.approx(0.1)


def test_rollup_applies_the_floor_to_the_summed_count_so_two_sparse_subtypes_can_clear_it() -> (
    None
):
    # 3 + 4 sequences: each subtype cell is below a floor of 5, their sum is not. A roll-up
    # that inherited the subtype cells' null fraction would leave this cell grey.
    out = _roll(("KE", 2010, "PI", "A1", 3, 1), ("KE", 2010, "PI", "D", 4, 2))
    cell = _cell(out, "KE", 2010, "PI")
    assert cell["n_sequences"] == 7
    assert cell["frac_resistant"] == pytest.approx(3 / 7)


def test_rollup_keeps_the_counts_but_not_the_fraction_when_the_sum_is_still_below_the_floor() -> (
    None
):
    out = _roll(("KE", 2010, "PI", "A1", 2, 1), ("KE", 2010, "PI", "D", 2, 1))
    cell = _cell(out, "KE", 2010, "PI")
    assert (cell["n_sequences"], cell["n_resistant"]) == (4, 2)
    assert cell["frac_resistant"] is None


def test_the_floor_is_inclusive_at_the_minimum() -> None:
    out = _roll(("KE", 2010, "PI", "A1", 5, 1), ("KE", 2011, "PI", "A1", 4, 1))
    assert _cell(out, "KE", 2010, "PI")["frac_resistant"] == pytest.approx(0.2)
    assert _cell(out, "KE", 2011, "PI")["frac_resistant"] is None


def test_the_default_floor_is_the_manuscript_figures_floor() -> None:
    at_floor = fig02_aggregate.MIN_CELL_N
    cells = _cells(
        ("KE", 2010, "PI", "A1", at_floor, 1),
        ("KE", 2011, "PI", "A1", at_floor - 1, 1),
    )
    out = rollup_over_subtype(cells)
    assert _cell(out, "KE", 2010, "PI")["frac_resistant"] is not None
    assert _cell(out, "KE", 2011, "PI")["frac_resistant"] is None


def test_rollup_keeps_only_the_targeted_countries() -> None:
    out = _roll(
        ("KE", 2010, "PI", "A1", 9, 1),
        ("TH", 2010, "PI", "A1", 9, 1),
        ("ZA", 2010, "PI", "A1", 9, 1),
    )
    assert out["country"].to_list() == ["KE"]


def test_rollup_orders_by_drug_class_then_search_order_of_country_then_year() -> None:
    out = _roll(
        ("UG", 2012, "INSTI", "A1", 9, 1),
        ("UG", 2011, "NNRTI", "A1", 9, 1),
        ("TZ", 2011, "NNRTI", "A1", 9, 1),
        ("KE", 2012, "NNRTI", "A1", 9, 1),
        ("KE", 2011, "NNRTI", "A1", 9, 1),
        ("ET", 2011, "PI", "A1", 9, 1),
    )
    assert list(out.select("drug_class", "country", "collection_year").iter_rows()) == [
        ("PI", "ET", 2011),
        ("NNRTI", "KE", 2011),
        ("NNRTI", "KE", 2012),
        ("NNRTI", "UG", 2011),
        ("NNRTI", "TZ", 2011),
        ("INSTI", "UG", 2012),
    ]


def test_rollup_output_does_not_depend_on_input_row_order() -> None:
    cells = _cells(
        ("KE", 2010, "PI", "A1", 7, 1),
        ("KE", 2010, "PI", "D", 5, 2),
        ("UG", 2011, "NNRTI", "C", 9, 4),
        ("TZ", 2012, "NRTI", "A1", 6, 0),
    )
    reordered = cells.reverse()
    assert reordered["country"].to_list() != cells["country"].to_list()
    assert rollup_over_subtype(reordered).equals(rollup_over_subtype(cells))


def test_rollup_refuses_a_drug_class_with_no_canonical_position() -> None:
    with pytest.raises(UnknownDrugClassError):
        _roll(("KE", 2010, "CAI", "A1", 9, 1))


def test_the_targeted_region_is_the_search_set_in_search_order() -> None:
    assert TARGETED_COUNTRIES == ("KE", "UG", "TZ", "RW", "BI", "ET", "SS", "SO")


def test_the_aggregation_is_the_manuscript_figures_not_a_copy() -> None:
    assert (
        landscape.aggregate_resistance_landscape
        is fig02_aggregate.aggregate_resistance_landscape
    )


# --- accounting -------------------------------------------------------------------------


def _dr(accessions: list[str]) -> pl.DataFrame:
    """Every accession scored for PI and NNRTI, one non-resistant call each."""
    rows = [
        {"accession": acc, "drug_class": cls, "drug": drug, "level": "susceptible"}
        for acc in accessions
        for cls, drug in (("PI", "ATV/r"), ("NNRTI", "EFV"))
    ]
    return pl.DataFrame(rows)


def _curated(*rows: tuple[str, str | None, int | None, str | None]) -> pl.DataFrame:
    return pl.DataFrame(
        rows,
        schema={
            "accession": pl.String,
            "country": pl.String,
            "collection_year": pl.Int64,
            "subtype": pl.String,
        },
        orient="row",
    )


def test_every_scored_accession_is_counted_once_under_the_reason_it_is_left_out() -> None:
    curated = _curated(
        ("A1", "KE", 2010, "A1"),  # drawn
        ("A2", "UG", 2011, "D"),  # drawn
        ("B1", "TH", 2010, "C"),  # outside the region
        ("B2", "BR", 2010, "B"),  # outside the region
        ("C1", None, 2010, "A1"),  # states no country
        ("D1", "KE", None, "A1"),  # in region, no collection year
        ("D2", "TZ", None, "D"),
        ("D3", "RW", None, "D"),
        ("E1", "KE", 2010, None),  # no subtype: the manuscript aggregation cannot place it
    )
    scored = _dr(
        ["A1", "A2", "B1", "B2", "C1", "D1", "D2", "D3", "E1", "F1"]
    )  # F1 not curated
    got = landscape_accounting(scored, curated)
    assert got.n_in_heatmap == 2
    assert got.n_out_of_region == 2
    assert got.n_no_country == 1
    assert got.n_no_year == 3
    assert got.n_other == 2
    assert got.n_total == scored["accession"].n_unique() == 10


def test_a_record_with_no_country_is_counted_as_that_and_not_as_out_of_region() -> None:
    # Null is not "outside the set"; an `is_in` test alone would file it under the wrong reason.
    got = landscape_accounting(_dr(["X"]), _curated(("X", None, 2010, "A1")))
    assert (got.n_no_country, got.n_out_of_region) == (1, 0)


def test_country_takes_precedence_over_a_missing_year() -> None:
    curated = _curated(("X", "TH", None, "C"), ("Y", None, None, "C"))
    got = landscape_accounting(_dr(["X", "Y"]), curated)
    assert (got.n_out_of_region, got.n_no_country, got.n_no_year) == (1, 1, 0)


def test_curated_accessions_that_were_never_scored_are_not_counted() -> None:
    curated = _curated(
        ("A", "KE", 2010, "A1"), ("Z", "TH", 2010, "C"), ("Y", None, None, None)
    )
    got = landscape_accounting(_dr(["A"]), curated)
    assert got.n_total == 1 and got.n_in_heatmap == 1


def test_a_duplicated_curated_accession_is_refused_rather_than_double_counted() -> None:
    curated = _curated(("A", "KE", 2010, "A1"), ("A", "KE", 2011, "A1"))
    with pytest.raises(DuplicateCuratedAccessionError):
        landscape_accounting(_dr(["A"]), curated)


def test_curated_missing_a_required_column_is_refused() -> None:
    with pytest.raises(ValueError, match="collection_year"):
        landscape_accounting(
            _dr(["A"]),
            pl.DataFrame({"accession": ["A"], "country": ["KE"], "subtype": ["A1"]}),
        )


# --- the whole view ---------------------------------------------------------------------


def _scored_with_outcomes(outcomes: dict[str, tuple[str, str]]) -> pl.DataFrame:
    """accession -> (PI level, NNRTI level)."""
    rows = [
        {"accession": acc, "drug_class": cls, "drug": drug, "level": level}
        for acc, (pi_level, nnrti_level) in outcomes.items()
        for cls, drug, level in (("PI", "ATV/r", pi_level), ("NNRTI", "EFV", nnrti_level))
    ]
    return pl.DataFrame(rows)


def test_resistance_landscape_counts_resistant_sequences_per_cell_and_excludes_the_region_rest() -> (
    None
):
    outcomes = {f"K{i}": ("susceptible", "high" if i < 2 else "susceptible") for i in range(6)}
    outcomes |= {f"T{i}": ("high", "high") for i in range(6)}  # Thai: must not be drawn
    curated = _curated(
        *[(f"K{i}", "KE", 2010, "A1" if i % 2 else "D") for i in range(6)],
        *[(f"T{i}", "TH", 2010, "C") for i in range(6)],
    )
    cells, accounting = resistance_landscape(_scored_with_outcomes(outcomes), curated)
    assert cells["country"].unique().to_list() == ["KE"]
    assert _cell(cells, "KE", 2010, "NNRTI") == {
        "country": "KE",
        "collection_year": 2010,
        "drug_class": "NNRTI",
        "n_sequences": 6,
        "n_resistant": 2,
        "frac_resistant": pytest.approx(2 / 6),
    }
    assert _cell(cells, "KE", 2010, "PI")["n_resistant"] == 0
    assert (accounting.n_in_heatmap, accounting.n_out_of_region) == (6, 6)


def test_resistance_landscape_refuses_cells_that_disagree_with_the_accounting(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # If the manuscript aggregation ever changed who it counts (a different join, a new filter),
    # the page would publish cells whose sequences the exclusion note does not account for.
    real = landscape.aggregate_resistance_landscape

    def inflated(*args: object, **kwargs: object) -> tuple[pl.DataFrame, object]:
        cells, counts = real(*args, **kwargs)  # type: ignore[arg-type]
        return cells.with_columns(n_sequences=pl.col("n_sequences") + 1), counts

    monkeypatch.setattr(landscape, "aggregate_resistance_landscape", inflated)
    curated = _curated(("A", "KE", 2010, "A1"), ("B", "KE", 2010, "D"))
    with pytest.raises(LandscapeAccountingError):
        resistance_landscape(_dr(["A", "B"]), curated)


def test_the_cells_must_agree_with_the_accounting() -> None:
    curated = _curated(("A", "KE", 2010, "A1"), ("B", "KE", 2010, "D"))
    scored = _dr(["A", "B"])
    cells, _ = resistance_landscape(scored, curated)
    placement = landscape._placement(scored, curated, TARGETED_COUNTRIES)
    landscape._check_cells_match_accounting(cells, placement, scored)

    tampered = cells.with_columns(n_sequences=pl.col("n_sequences") + 1)
    with pytest.raises(LandscapeAccountingError):
        landscape._check_cells_match_accounting(tampered, placement, scored)


def test_importing_the_export_does_not_load_biopython_or_the_genbank_ingester() -> None:
    # The drift job imports `summary` to compare two md5 strings; the dashboard parses no
    # sequences, so it must not pay for the ingester's dependencies to read a country list.
    code = (
        "import sys, hivdrwatch.dashboard.export.summary;"
        "print(sorted(m for m in ('Bio', 'hivdrwatch.ingest.genbank') if m in sys.modules))"
    )
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert done.stdout.strip() == "[]"


def test_the_ingest_package_still_exposes_its_ingesters_and_one_country_table() -> None:
    from hivdrwatch import ingest
    from hivdrwatch.ingest import genbank, schema

    assert ingest.GenBankIngester is genbank.GenBankIngester
    assert ingest.LANLIngester is not None and callable(ingest.merge_sources)
    assert (
        ingest.EAST_AFRICA_COUNTRIES
        is genbank.EAST_AFRICA_COUNTRIES
        is schema.EAST_AFRICA_COUNTRIES
    )
    with pytest.raises(AttributeError):
        ingest.NoSuchThing  # noqa: B018
