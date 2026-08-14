"""Unit tests for the Fig 2 aggregator (pure Python; no matplotlib).

The heatmap rendering itself is exercised by running the notebook
against a real parquet; that's an integration test. This file locks in
the aggregation semantics: worst-per-class, join shape, MIN_CELL_N gate.
"""

from __future__ import annotations

import polars as pl
import pytest

from hivdrwatch.paper1.fig02_aggregate import (
    MIN_CELL_N,
    RESISTANT_LEVELS,
    aggregate_resistance_landscape,
)


def _dr(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


def _curated(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


class TestAggregation:
    def test_shape_and_columns(self) -> None:
        dr = _dr(
            [
                {"accession": "A1", "drug_class": "NRTI", "level": "high"},
                {"accession": "A1", "drug_class": "NNRTI", "level": "susceptible"},
                {"accession": "A2", "drug_class": "NRTI", "level": "susceptible"},
            ]
        )
        curated = _curated(
            [
                {"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "A2", "country": "KE", "collection_year": 2020, "subtype": "A1"},
            ]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        assert set(landscape.columns) >= {
            "country",
            "collection_year",
            "drug_class",
            "subtype",
            "n_sequences",
            "n_resistant",
            "frac_resistant",
        }
        # 2 sequences × 2 drug classes (NRTI, NNRTI) = at most 2 rows
        # after per-seq-per-class dedup and one-country-year-subtype group
        assert landscape.height >= 1
        nrti_row = landscape.filter(pl.col("drug_class") == "NRTI").to_dicts()[0]
        assert nrti_row["n_sequences"] == 2
        assert nrti_row["n_resistant"] == 1
        assert nrti_row["frac_resistant"] == pytest.approx(0.5)

    def test_worst_per_class(self) -> None:
        """A sequence resistant to ANY drug in the class counts as resistant."""
        dr = _dr(
            [
                # A1: susceptible to 3TC, high to FTC — should count as NRTI-resistant
                {"accession": "A1", "drug_class": "NRTI", "level": "susceptible"},
                {"accession": "A1", "drug_class": "NRTI", "level": "high"},
            ]
        )
        curated = _curated(
            [
                {"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"},
            ]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        assert landscape.height == 1
        assert landscape["n_resistant"][0] == 1

    def test_intermediate_counts_as_resistant(self) -> None:
        dr = _dr([{"accession": "A1", "drug_class": "NRTI", "level": "intermediate"}])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        assert landscape["n_resistant"][0] == 1
        assert landscape["frac_resistant"][0] == 1.0

    def test_low_low_dont_count(self) -> None:
        dr = _dr(
            [
                {"accession": "A1", "drug_class": "NRTI", "level": "low"},
                {"accession": "A1", "drug_class": "NRTI", "level": "potential_low"},
                {"accession": "A1", "drug_class": "NRTI", "level": "susceptible"},
            ]
        )
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        assert landscape["n_resistant"][0] == 0

    def test_min_cell_n_masks_frac_but_keeps_counts(self) -> None:
        # 3 sequences × 1 class → n=3, below MIN_CELL_N=5
        dr = _dr(
            [{"accession": f"A{i}", "drug_class": "NRTI", "level": "high"} for i in range(3)]
        )
        curated = _curated(
            [
                {
                    "accession": f"A{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "A1",
                }
                for i in range(3)
            ]
        )
        landscape, counts = aggregate_resistance_landscape(dr, curated, min_cell_n=5)
        assert landscape.height == 1
        assert landscape["n_sequences"][0] == 3
        assert landscape["n_resistant"][0] == 3
        assert landscape["frac_resistant"][0] is None  # below min
        assert counts.n_cells_above_min == 0

    def test_missing_country_or_year_dropped(self) -> None:
        dr = _dr(
            [
                {"accession": "A1", "drug_class": "NRTI", "level": "high"},
                {"accession": "A2", "drug_class": "NRTI", "level": "high"},
            ]
        )
        curated = _curated(
            [
                {"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "A2", "country": None, "collection_year": 2020, "subtype": "A1"},
            ]
        )
        _landscape, counts = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        assert counts.n_joined_rows == 1  # A2 dropped for missing country

    def test_missing_required_columns_raises(self) -> None:
        with pytest.raises(ValueError, match="missing required columns"):
            aggregate_resistance_landscape(
                pl.DataFrame({"accession": ["A"]}),
                pl.DataFrame({"accession": ["A"]}),
            )

    def test_stratifies_by_subtype(self) -> None:
        dr = _dr(
            [
                {"accession": "A1", "drug_class": "NRTI", "level": "high"},
                {"accession": "B1", "drug_class": "NRTI", "level": "susceptible"},
            ]
        )
        curated = _curated(
            [
                {"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "B1", "country": "KE", "collection_year": 2020, "subtype": "C"},
            ]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        # One row per subtype
        assert set(landscape["subtype"].to_list()) == {"A1", "C"}
        a1_row = landscape.filter(pl.col("subtype") == "A1").to_dicts()[0]
        c_row = landscape.filter(pl.col("subtype") == "C").to_dicts()[0]
        assert a1_row["n_resistant"] == 1
        assert c_row["n_resistant"] == 0


class TestConstants:
    def test_resistant_levels_pinned(self) -> None:
        assert frozenset({"high", "intermediate"}) == RESISTANT_LEVELS

    def test_min_cell_n_is_five(self) -> None:
        assert MIN_CELL_N == 5
