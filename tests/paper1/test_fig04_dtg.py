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
    FDR_Q,
    INSTI_CLASS,
    MIN_EPOCH_N_PRIMARY,
    MIN_EPOCH_N_SECONDARY,
    POST_YEARS,
    PRE_YEARS,
    RESISTANT_LEVELS,
    ROLLOUT_YEAR,
    _benjamini_hochberg,
    _clopper_pearson_ci,
    _fisher_exact,
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
                {
                    "accession": acc,
                    "country": country,
                    "collection_year": year,
                    "subtype": subtype,
                }
            )

    _add("PRE", pre_year, n_pre, n_pre_r)
    _add("POST", post_year, n_post, n_post_r)
    return dr_rows, curated_rows


class TestConstants:
    def test_resistant_levels(self) -> None:
        assert frozenset({"high", "intermediate"}) == RESISTANT_LEVELS

    def test_insti_class(self) -> None:
        assert INSTI_CLASS == "INSTI"

    def test_epochs(self) -> None:
        assert PRE_YEARS == (2010, 2018)
        assert POST_YEARS[0] == 2020
        assert ROLLOUT_YEAR == 2019

    def test_min_epoch_n_tiers(self) -> None:
        assert MIN_EPOCH_N_PRIMARY == 30
        assert MIN_EPOCH_N_SECONDARY == 5

    def test_fdr_q(self) -> None:
        assert FDR_Q == 0.05


class TestClopperPearson:
    def test_zero_events_lower_bound_is_zero(self) -> None:
        p, lo, hi = _clopper_pearson_ci(0, 30)
        assert lo == 0.0 and p == 0.0 and 0 < hi < 1

    def test_all_events_upper_bound_is_one(self) -> None:
        p, lo, hi = _clopper_pearson_ci(30, 30)
        assert p == 1.0 and hi == 1.0 and 0 < lo < 1

    def test_bracket_always_holds(self) -> None:
        for k, n in [(0, 5), (1, 5), (5, 5), (50, 100), (0, 1000)]:
            p, lo, hi = _clopper_pearson_ci(k, n)
            assert 0 <= lo <= p <= hi <= 1


class TestFisherExact:
    def test_no_difference_gives_p_near_one(self) -> None:
        _, p = _fisher_exact(5, 100, 5, 100)
        assert p > 0.5

    def test_clear_signal_gives_small_p(self) -> None:
        _, p = _fisher_exact(30, 100, 3, 100)
        assert p < 0.001

    def test_zero_n_returns_p_one(self) -> None:
        _, p = _fisher_exact(0, 0, 3, 30)
        assert p == 1.0


class TestBH:
    def test_all_significant_pass(self) -> None:
        qs = _benjamini_hochberg([0.001, 0.002, 0.003], q=0.05)
        assert all(q <= 0.05 for q in qs)

    def test_monotone_ordering_preserved(self) -> None:
        pvals = [0.9, 0.001, 0.5, 0.02]
        qs = _benjamini_hochberg(pvals, q=0.05)
        # Smaller p → smaller q
        rank_by_p = sorted(range(len(pvals)), key=lambda i: pvals[i])
        prev_q = -1.0
        for i in rank_by_p:
            assert qs[i] >= prev_q
            prev_q = qs[i]

    def test_empty(self) -> None:
        assert _benjamini_hochberg([], q=0.05) == []


class TestEpochAssignment:
    def test_rollout_year_2019_dropped(self) -> None:
        dr = _dr([{"accession": "A1", "drug_class": INSTI_CLASS, "level": "high"}])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2019, "subtype": "A1"}]
        )
        per_epoch, delta, summary = aggregate_dtg_landscape(
            dr, curated, min_epoch_n_secondary=1
        )
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
        per_epoch, _, summary = aggregate_dtg_landscape(dr, curated, min_epoch_n_secondary=1)
        assert summary.n_pre_seqs == 1
        assert per_epoch.filter(pl.col("epoch") == "pre").height == 1

    def test_post_boundary_2020_included(self) -> None:
        dr = _dr([{"accession": "A1", "drug_class": INSTI_CLASS, "level": "high"}])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        per_epoch, _, summary = aggregate_dtg_landscape(dr, curated, min_epoch_n_secondary=1)
        assert summary.n_post_seqs == 1
        assert per_epoch.filter(pl.col("epoch") == "post").height == 1


class TestResistanceCounting:
    def test_intermediate_counts(self) -> None:
        dr_rows, curated_rows = _mk_stratum(
            country="KE",
            subtype="A1",
            pre_year=2015,
            post_year=2021,
            n_pre=1,
            n_pre_r=0,
            n_post=1,
            n_post_r=0,
        )
        dr_rows[1] = {
            "accession": dr_rows[1]["accession"],
            "drug_class": INSTI_CLASS,
            "level": "intermediate",
        }
        per_epoch, _, _ = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n_secondary=1
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
        per_epoch, _, _ = aggregate_dtg_landscape(dr, curated, min_epoch_n_secondary=1)
        pre = per_epoch.filter(pl.col("epoch") == "pre").to_dicts()[0]
        assert pre["n_insti_resistant"] == 0


class TestTiering:
    def test_below_secondary_threshold_omitted_from_delta(self) -> None:
        # 3 pre + 3 post — below secondary threshold (5) — omitted entirely
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021, n_pre=3, n_pre_r=0, n_post=3, n_post_r=1
        )
        per_epoch, delta, _ = aggregate_dtg_landscape(
            _dr(dr_rows),
            _curated(curated_rows),
            min_epoch_n_secondary=5,
        )
        assert per_epoch.height == 2
        assert delta.height == 0

    def test_secondary_tier_scored_but_not_headline(self) -> None:
        # 10 pre + 10 post — passes secondary (≥5) but below primary (≥30)
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021, n_pre=10, n_pre_r=0, n_post=10, n_post_r=5
        )
        _, delta, summary = aggregate_dtg_landscape(_dr(dr_rows), _curated(curated_rows))
        assert delta.height == 1
        row = delta.to_dicts()[0]
        assert row["tier"] == "secondary"
        # Even if BH-significant, secondary strata don't count toward
        # the preregistered headline.
        assert summary.n_strata_secondary == 1
        assert summary.n_strata_primary == 0

    def test_primary_tier_with_bh_signal(self) -> None:
        # Big jump with primary-tier N: pre 0/100, post 25/100
        # Fisher p is extremely small; BH q also small → signal True.
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021, n_pre=100, n_pre_r=0, n_post=100, n_post_r=25
        )
        _, delta, summary = aggregate_dtg_landscape(_dr(dr_rows), _curated(curated_rows))
        assert delta.height == 1
        row = delta.to_dicts()[0]
        assert row["tier"] == "primary"
        assert row["signal"] is True
        assert row["signal_direction"] == "up"
        assert row["bh_q"] < 0.05
        assert row["fisher_p"] < 0.001
        assert summary.n_headline_after_fdr == 1

    def test_primary_no_signal_when_fisher_non_significant(self) -> None:
        # Small change with primary-tier N: pre 5/100, post 6/100 → Fisher p large
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021, n_pre=100, n_pre_r=5, n_post=100, n_post_r=6
        )
        _, delta, summary = aggregate_dtg_landscape(_dr(dr_rows), _curated(curated_rows))
        assert delta.height == 1
        row = delta.to_dicts()[0]
        assert row["tier"] == "primary"
        assert row["signal"] is False
        assert summary.n_headline_after_fdr == 0

    def test_downward_direction_never_flagged_as_signal(self) -> None:
        # Post lower than pre — signal must be False regardless of Fisher p
        dr_rows, curated_rows = _mk_stratum(
            "KE", "A1", 2015, 2021, n_pre=100, n_pre_r=30, n_post=100, n_post_r=5
        )
        _, delta, _ = aggregate_dtg_landscape(_dr(dr_rows), _curated(curated_rows))
        row = delta.to_dicts()[0]
        assert row["signal_direction"] == "down"
        assert row["signal"] is False


class TestSubtypeStratification:
    def test_separate_deltas_per_subtype(self) -> None:
        dr_rows: list[dict] = []
        curated_rows: list[dict] = []
        for subtype in ("A1", "C"):
            d, c = _mk_stratum(
                "KE",
                subtype,
                2015,
                2021,
                n_pre=50,
                n_pre_r=0,
                n_post=50,
                n_post_r=10,
            )
            dr_rows.extend(d)
            curated_rows.extend(c)
        _, delta, _ = aggregate_dtg_landscape(
            _dr(dr_rows), _curated(curated_rows), min_epoch_n_primary=30
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
            _dr(dr_rows), curated, min_epoch_n_secondary=1
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
