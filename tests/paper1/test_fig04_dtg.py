"""Unit tests for the Fig 4 DTG case study aggregator.

Locks in the preregistered aggregation semantics from
`paper1_tools/preregistration_fig04.md`. Any test that fails here
means either the code drifted from the plan (fix the code) or the
plan changed (amend the preregistration).
"""

from __future__ import annotations

import polars as pl
import pytest

from hivdrwatch.paper1.fig04_dtg import (
    INSTI_CLASS,
    MIN_EPOCH_N,
    POST_YEARS,
    PRE_YEARS,
    RESISTANT_LEVELS,
    ROLLOUT_YEAR,
    aggregate_dtg_landscape,
)


def _dr(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


def _curated(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


def _mk_stratum(
    country: str,
    subtype: str,
    pre_year: int,
    post_year: int,
    n_pre: int,
    n_pre_r: int,
    n_post: int,
    n_post_r: int,
) -> tuple[list[dict], list[dict]]:
    """Emit dr + curated rows for one stratum in one call."""
    dr_rows: list[dict] = []
    curated_rows: list[dict] = []

    def _add(prefix: str, year: int, count: int, r_count: int) -> None:
        for i in range(count):
            acc = f"{prefix}{country}{subtype}{i:04d}"
            level = "high" if i < r_count else "susceptible"
            dr_rows.append({"accession": acc, "drug_class": INSTI_CLASS, "level": level})
            curated_rows.append(
                {"accession": acc, "country": country, "collection_year": year, "subtype": subtype}
            )

    _add("PRE", pre_year, n_pre, n_pre_r)
    _add("POST", post_year, n_post, n_post_r)
    return dr_rows, curated_rows


class TestConstants:
    def test_resistant_levels(self) -> None:
        assert RESISTANT_LEVELS == frozenset({"high", "intermediate"})

    def test_insti_class(self) -> None:
        assert INSTI_CLASS == "INSTI"

    def test_epochs(self) -> None:
        assert PRE_YEARS == (2010, 2018)
        assert POST_YEARS[0] == 2020
        assert ROLLOUT_YEAR == 2019

    def test_min_epoch_n(self) -> None:
        assert MIN_EPOCH_N == 30


class TestEpochAssignment:
    def test_rollout_year_2019_dropped(self) -> None:
        dr = _dr([{"accession": "A1", "drug_class": INSTI_CLASS, "level": "high"}])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2019, "subtype": "A1"}]
        )
        per_epoch, delta, summary = aggregate_dtg_landscape(dr, curated, min_epoch_n=1)
        # 2019 is buffer — sequence should not appear in either epoch.
        assert summary.n_pre_seqs == 0
        assert summary.n_post_seqs == 0
        assert per_epoch.height == 0
        assert delta.height == 0

    def test_pre_boundary_2018_included(self) -> None:
        dr = _dr([{"accession": "A1", "drug_class": INSTI_CLASS, "level": "high"}])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2018, "subtype": "A1"}]
        )
        per_epoch, _, summary = aggregate_dtg_landscape(dr, curated, min_epoch_n=1)
        assert summary.n_pre_seqs == 1
        assert per_epoch.filter(pl.col("epoch") == "pre").height == 1

    def test_post_boundary_2020_included(self) -> None:
        dr = _dr([{"accession": "A1", "drug_class": INSTI_CLASS, "level": "high"}])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        per_epoch, _, summary = aggregate_dtg_landscape(dr, curated, min_epoch_n=1)
        assert summary.n_post_seqs == 1
        assert per_epoch.filter(pl.col("epoch") == "post").height == 1


class TestResistanceCounting:
    def test_intermediate_counts(self) -> None:
        dr_rows, curated_rows = _mk_stratum(
            country="KE", subtype="A1",
            pre_year=2015, post_year=2021,
            n_pre=1, n_pre_r=0,
            n_post=1, n_post_r=0,
        )
        dr_rows[1] = {"accession": dr_rows[1]["accession"], "drug_class": INSTI_CLASS, "level": "intermediate"}
        per_epoch, _, _ = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n=1
        )
        post = per_epoch.filter(pl.col("epoch") == "post").to_dicts()[0]
        assert post["n_insti_resistant"] == 1

    def test_low_and_susceptible_dont_count(self) -> None:
        dr = _dr(
            [
                {"accession": "A1", "drug_class": INSTI_CLASS, "level": "low"},
                {"accession": "A2", "drug_class": INSTI_CLASS, "level": "potential_low"},
                {"accession": "A3", "drug_class": INSTI_CLASS, "level": "susceptible"},
            ]
        )
        curated = _curated(
            [
                {"accession": a, "country": "KE", "collection_year": 2015, "subtype": "A1"}
                for a in ("A1", "A2", "A3")
            ]
        )
        per_epoch, _, _ = aggregate_dtg_landscape(dr, curated, min_epoch_n=1)
        pre = per_epoch.filter(pl.col("epoch") == "pre").to_dicts()[0]
        assert pre["n_insti_resistant"] == 0


class TestMinEpochN:
    def test_undersized_stratum_omitted_from_delta(self) -> None:
        # 5 pre + 5 post, min_epoch_n=30 → should be omitted from delta.
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021, n_pre=5, n_pre_r=0, n_post=5, n_post_r=2
        )
        per_epoch, delta, _ = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n=30
        )
        # Per-epoch frame still shows counts; delta is empty.
        assert per_epoch.height == 2
        assert delta.height == 0

    def test_signal_flagged_when_ci_excludes_zero(self) -> None:
        # Big jump with enough N: pre 0/100, post 20/100 → CI clearly > 0.
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021,
            n_pre=100, n_pre_r=0,
            n_post=100, n_post_r=20,
        )
        _, delta, summary = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n=30
        )
        assert delta.height == 1
        row = delta.to_dicts()[0]
        assert row["signal"] is True
        assert row["delta_prev"] > 0
        assert row["delta_ci_low"] > 0
        assert summary.n_strata_with_signal == 1

    def test_no_signal_when_ci_spans_zero(self) -> None:
        # Small change with enough N: pre 5/100, post 6/100 → CI spans 0.
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021,
            n_pre=100, n_pre_r=5,
            n_post=100, n_post_r=6,
        )
        _, delta, summary = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n=30
        )
        assert delta.height == 1
        row = delta.to_dicts()[0]
        assert row["signal"] is False
        assert summary.n_strata_with_signal == 0


class TestSubtypeStratification:
    def test_separate_deltas_per_subtype(self) -> None:
        dr_rows: list[dict] = []
        curated_rows: list[dict] = []
        for subtype in ("A1", "C"):
            d, c = _mk_stratum(
                "KE", subtype, 2015, 2021,
                n_pre=50, n_pre_r=0, n_post=50, n_post_r=10,
            )
            dr_rows.extend(d)
            curated_rows.extend(c)
        _, delta, _ = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n=30
        )
        assert set(delta["subtype"].to_list()) == {"A1", "C"}


class TestNonInstiIgnored:
    def test_nrti_rows_dont_leak_in(self) -> None:
        dr_rows = [
            {"accession": "A1", "drug_class": "NRTI", "level": "high"},
            {"accession": "A1", "drug_class": INSTI_CLASS, "level": "susceptible"},
        ]
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2021, "subtype": "A1"}]
        )
        per_epoch, _, _ = aggregate_dtg_landscape(
            _dr(dr_rows), curated, min_epoch_n=1
        )
        post = per_epoch.filter(pl.col("epoch") == "post").to_dicts()[0]
        # NRTI row must not count as INSTI resistance.
        assert post["n_insti_resistant"] == 0


class TestValidation:
    def test_missing_columns_raises(self) -> None:
        with pytest.raises(ValueError, match="missing required columns"):
            aggregate_dtg_landscape(
                pl.DataFrame({"accession": ["A"]}),
                pl.DataFrame({"accession": ["A"]}),
            )
