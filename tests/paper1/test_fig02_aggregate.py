"""Unit tests for the Fig 2 aggregator (pure Python; no matplotlib).

The heatmap rendering itself is exercised by running the notebook
against a real parquet; that's an integration test. This file locks in
the aggregation semantics: worst-per-class, join shape, MIN_CELL_N gate.
"""

from __future__ import annotations

import polars as pl
import pytest

from hivdrwatch.paper1 import SUBTYPE_ORDER
from hivdrwatch.paper1.fig02_aggregate import (
    MIN_CELL_N,
    RESISTANT_LEVELS,
    aggregate_resistance_landscape,
    subtypes_omitted_by,
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


class TestSubtypesOmittedByDisplayOrder:
    """The heatmap draws one column per SUBTYPE_ORDER entry and skips the
    rest, so the omission has to be visible. `A` and `CRF01_AE` are the
    labels ingest 0.2.0 started resolving that no display order covers.
    """

    def _landscape(self, *subtypes: str) -> pl.DataFrame:
        dr = _dr(
            [
                {"accession": f"S{i}", "drug_class": "NRTI", "level": "high"}
                for i, _ in enumerate(subtypes)
            ]
        )
        curated = _curated(
            [
                {
                    "accession": f"S{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": s,
                }
                for i, s in enumerate(subtypes)
            ]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        return landscape

    def test_reports_labels_with_no_column(self) -> None:
        landscape = self._landscape("A1", "C", "A", "CRF01_AE")
        omitted = subtypes_omitted_by(landscape, SUBTYPE_ORDER)
        assert omitted.labels == ("A", "CRF01_AE")
        assert omitted.n_cells == 2
        assert omitted.n_observations == 2
        assert bool(omitted) is True
        assert "CRF01_AE" in omitted.describe()

    def test_empty_when_the_order_covers_everything(self) -> None:
        landscape = self._landscape("A1", "C", "D")
        omitted = subtypes_omitted_by(landscape, SUBTYPE_ORDER)
        assert omitted.labels == ()
        assert omitted.n_cells == 0
        assert bool(omitted) is False
        assert "no observed subtype is omitted" in omitted.describe()

    def test_counts_every_drug_class_cell(self) -> None:
        dr = _dr(
            [
                {"accession": "S0", "drug_class": "NRTI", "level": "high"},
                {"accession": "S0", "drug_class": "NNRTI", "level": "high"},
                {"accession": "S0", "drug_class": "INSTI", "level": "susceptible"},
            ]
        )
        curated = _curated(
            [
                {
                    "accession": "S0",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "CRF01_AE",
                }
            ]
        )
        landscape, _ = aggregate_resistance_landscape(dr, curated, min_cell_n=1)
        omitted = subtypes_omitted_by(landscape, SUBTYPE_ORDER)
        assert omitted.labels == ("CRF01_AE",)
        # One sequence, three drug classes: three cells, three observations.
        assert omitted.n_cells == 3
        assert omitted.n_observations == 3

    def test_warns_so_a_dropped_stratum_cannot_pass_silently(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        landscape = self._landscape("A1", "CRF01_AE")
        with caplog.at_level("WARNING"):
            subtypes_omitted_by(landscape, SUBTYPE_ORDER)
        assert any("CRF01_AE" in r.getMessage() for r in caplog.records)

    def test_requires_a_subtype_column(self) -> None:
        with pytest.raises(ValueError, match="missing 'subtype' column"):
            subtypes_omitted_by(pl.DataFrame({"country": ["KE"]}), SUBTYPE_ORDER)

    def test_current_display_order_has_no_slot_for_parent_clades_or_crfs(self) -> None:
        """Pins the gap itself: if SUBTYPE_ORDER is widened, this fails and
        the figure-2 decision in docs/evaluations/ingest-provenance.md has
        been made."""
        assert "A" not in SUBTYPE_ORDER
        assert not any(s.startswith("CRF") for s in SUBTYPE_ORDER)
