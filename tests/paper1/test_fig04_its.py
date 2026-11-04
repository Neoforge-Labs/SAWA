"""Unit tests for the Fig 4 secondary analysis — segmented regression (ITS).

Preregistration §6 "Secondary". The OLS is hand-rolled on numpy rather
than pulling in statsmodels (CLAUDE.md dependency rule 2 — scipy/numpy
are already in the tree), so these tests check it against an
independent implementation wherever one exists: `_ols` reduced to a
single predictor must agree with `scipy.stats.linregress` on every
quantity it reports.
"""

from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest
from scipy import stats

from hivdrwatch.paper1.fig04_dtg import INSTI_CLASS, ROLLOUT_YEAR
from hivdrwatch.paper1.fig04_its import (
    MIN_OBSERVATIONS,
    ITSNotIdentifiableError,
    _ols,
    build_yearly_prevalence,
    fit_country_its,
    fit_its,
)


def _rows(country: str, per_year: dict[int, tuple[int, int]]) -> tuple[list[dict], list[dict]]:
    """(n_sequences, n_resistant) per year -> dr + curated rows."""
    dr_rows: list[dict] = []
    curated_rows: list[dict] = []
    for year, (n, n_r) in per_year.items():
        for i in range(n):
            acc = f"{country}{year}{i:04d}"
            dr_rows.append(
                {
                    "accession": acc,
                    "drug_class": INSTI_CLASS,
                    "level": "high" if i < n_r else "susceptible",
                }
            )
            curated_rows.append(
                {
                    "accession": acc,
                    "country": country,
                    "collection_year": year,
                    "subtype": "A1",
                }
            )
    return dr_rows, curated_rows


class TestOlsAgainstScipy:
    """`_ols` must match scipy on the case scipy can also fit."""

    def test_matches_linregress_on_single_predictor(self) -> None:
        rng = np.random.default_rng(20260922)
        x = np.arange(12, dtype=float)
        y = 0.4 + 0.07 * x + rng.normal(0, 0.05, size=12)

        design = np.column_stack([np.ones_like(x), x])
        fit = _ols(design, y)
        ref = stats.linregress(x, y)

        assert fit.betas[0] == pytest.approx(ref.intercept, rel=1e-10)
        assert fit.betas[1] == pytest.approx(ref.slope, rel=1e-10)
        assert fit.std_errors[1] == pytest.approx(ref.stderr, rel=1e-10)
        assert fit.std_errors[0] == pytest.approx(ref.intercept_stderr, rel=1e-10)
        assert fit.p_values[1] == pytest.approx(ref.pvalue, rel=1e-10)
        assert fit.r_squared == pytest.approx(ref.rvalue**2, rel=1e-10)
        assert fit.residual_df == 10

    def test_recovers_a_noise_free_generating_process_exactly(self) -> None:
        # y = 2 + 3*a - 1.5*b, fit exactly -> r2 == 1, zero residual
        a = np.array([0.0, 1, 2, 3, 4, 5])
        b = np.array([0.0, 0, 0, 1, 1, 1])
        y = 2 + 3 * a - 1.5 * b
        fit = _ols(np.column_stack([np.ones_like(a), a, b]), y)
        assert fit.betas == pytest.approx([2.0, 3.0, -1.5], abs=1e-9)
        assert fit.r_squared == pytest.approx(1.0, abs=1e-12)

    def test_raises_when_underdetermined(self) -> None:
        # 3 observations, 4 parameters -> no residual df
        x = np.arange(3, dtype=float)
        design = np.column_stack([np.ones_like(x), x, x**2, x**3])
        with pytest.raises(ITSNotIdentifiableError):
            _ols(design, np.array([1.0, 2.0, 3.0]))


class TestYearlyPrevalence:
    def test_one_row_per_country_year_with_prevalence(self) -> None:
        dr, cur = _rows("KE", {2015: (10, 2), 2016: (20, 5)})
        yearly = build_yearly_prevalence(pl.from_dicts(dr), pl.from_dicts(cur))
        assert yearly.height == 2
        row = yearly.filter(pl.col("year") == 2016).to_dicts()[0]
        assert row["n_sequences"] == 20
        assert row["n_resistant"] == 5
        assert row["prevalence"] == pytest.approx(0.25)

    def test_rollout_year_excluded_as_buffer(self) -> None:
        dr, cur = _rows("KE", {2018: (5, 1), ROLLOUT_YEAR: (5, 5), 2020: (5, 2)})
        yearly = build_yearly_prevalence(pl.from_dicts(dr), pl.from_dicts(cur))
        assert ROLLOUT_YEAR not in yearly["year"].to_list()

    def test_post_flag_marks_years_after_the_buffer(self) -> None:
        dr, cur = _rows("KE", {2015: (5, 0), 2021: (5, 1)})
        yearly = build_yearly_prevalence(pl.from_dicts(dr), pl.from_dicts(cur))
        flags = dict(zip(yearly["year"], yearly["post"], strict=True))
        assert flags[2015] == 0
        assert flags[2021] == 1

    def test_time_is_centred_on_the_intervention(self) -> None:
        dr, cur = _rows("KE", {2015: (5, 0), 2021: (5, 1)})
        yearly = build_yearly_prevalence(pl.from_dicts(dr), pl.from_dicts(cur))
        t = dict(zip(yearly["year"], yearly["t"], strict=True))
        assert t[2015] == 2015 - ROLLOUT_YEAR
        assert t[2021] == 2021 - ROLLOUT_YEAR


class TestCountryFit:
    def _yearly(self, per_year: dict[int, tuple[int, int]]) -> pl.DataFrame:
        dr, cur = _rows("KE", per_year)
        return build_yearly_prevalence(pl.from_dicts(dr), pl.from_dicts(cur))

    def test_recovers_a_planted_level_and_slope_change(self) -> None:
        # Construct prevalence exactly linear in t with a level jump at
        # the interruption, then check the four coefficients come back.
        # pre:  prev = 0.02 + 0.001*t
        # post: prev = 0.02 + 0.001*t + 0.05 + 0.004*t
        per_year: dict[int, tuple[int, int]] = {}
        n = 1000
        for year in list(range(2011, 2019)) + list(range(2020, 2025)):
            t = year - ROLLOUT_YEAR
            post = 1 if year > ROLLOUT_YEAR else 0
            prev = 0.02 + 0.001 * t + post * (0.05 + 0.004 * t)
            per_year[year] = (n, round(prev * n))

        fit = fit_country_its(self._yearly(per_year))
        est = {k: v.estimate for k, v in fit.coefficients.items()}
        assert est["intercept"] == pytest.approx(0.02, abs=2e-3)
        assert est["t"] == pytest.approx(0.001, abs=2e-4)
        assert est["post"] == pytest.approx(0.05, abs=2e-3)
        assert est["t_post"] == pytest.approx(0.004, abs=2e-4)

    def test_flat_series_gives_no_significant_change(self) -> None:
        per_year = {y: (500, 25) for y in list(range(2011, 2019)) + list(range(2020, 2025))}
        fit = fit_country_its(self._yearly(per_year))
        assert fit.coefficients["post"].p_value > 0.05
        assert fit.coefficients["t_post"].p_value > 0.05

    def test_coefficient_reporting_invariants(self) -> None:
        per_year: dict[int, tuple[int, int]] = {}
        for i, year in enumerate(list(range(2011, 2019)) + list(range(2020, 2025))):
            per_year[year] = (300, 10 + i)
        fit = fit_country_its(self._yearly(per_year))
        for coef in fit.coefficients.values():
            assert coef.std_error > 0
            assert 0.0 <= coef.p_value <= 1.0
            assert coef.ci_low < coef.estimate < coef.ci_high
            assert math.isfinite(coef.t_statistic)

    def test_requires_both_segments(self) -> None:
        # Pre-period only — `post` and the interaction are all zero, so
        # the design matrix is rank-deficient and the model is not
        # identifiable. Must say so rather than return garbage.
        per_year = {y: (100, 5) for y in range(2011, 2019)}
        with pytest.raises(ITSNotIdentifiableError, match="post"):
            fit_country_its(self._yearly(per_year))

    def test_requires_enough_observations(self) -> None:
        per_year = {2017: (100, 5), 2018: (100, 6), 2020: (100, 9)}
        with pytest.raises(ITSNotIdentifiableError):
            fit_country_its(self._yearly(per_year))

    def test_min_observations_leaves_residual_df(self) -> None:
        assert MIN_OBSERVATIONS >= 5

    def test_rejects_a_multi_country_frame(self) -> None:
        # fit_country_its fits one country; handing it a pooled frame
        # would silently fit the wrong model, so it must refuse.
        dr_ke, cur_ke = _rows("KE", {y: (100, 5) for y in range(2011, 2019)})
        dr_ug, cur_ug = _rows("UG", {y: (100, 5) for y in range(2020, 2025)})
        yearly = build_yearly_prevalence(
            pl.from_dicts(dr_ke + dr_ug), pl.from_dicts(cur_ke + cur_ug)
        )
        with pytest.raises(ValueError, match="exactly one country"):
            fit_country_its(yearly)


class TestFitIts:
    def test_one_row_per_country_coefficient(self) -> None:
        dr_all: list[dict] = []
        cur_all: list[dict] = []
        for country in ("KE", "UG"):
            per_year = {
                y: (200, 5 + i)
                for i, y in enumerate(list(range(2011, 2019)) + list(range(2020, 2025)))
            }
            dr, cur = _rows(country, per_year)
            dr_all += dr
            cur_all += cur

        yearly = build_yearly_prevalence(pl.from_dicts(dr_all), pl.from_dicts(cur_all))
        coefs, fits = fit_its(yearly)
        assert {f.country for f in fits} == {"KE", "UG"}
        assert coefs.height == 2 * 4
        assert set(coefs["term"].unique()) == {"intercept", "t", "post", "t_post"}

    def test_unfittable_country_is_skipped_not_fatal(self) -> None:
        # KE is fittable; UG has pre-period only and must be dropped
        # without taking the whole run down.
        dr_ke, cur_ke = _rows(
            "KE",
            {
                y: (200, 5 + i)
                for i, y in enumerate(list(range(2011, 2019)) + list(range(2020, 2025)))
            },
        )
        dr_ug, cur_ug = _rows("UG", {y: (100, 5) for y in range(2011, 2019)})
        yearly = build_yearly_prevalence(
            pl.from_dicts(dr_ke + dr_ug), pl.from_dicts(cur_ke + cur_ug)
        )
        coefs, fits = fit_its(yearly)
        assert {f.country for f in fits} == {"KE"}
        assert set(coefs["country"].unique()) == {"KE"}
