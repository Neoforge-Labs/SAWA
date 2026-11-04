"""Fig 4 secondary analysis — segmented regression (interrupted time series).

Preregistration §6 "Secondary": per-country segmented regression of
yearly INSTI-resistance prevalence on

    prevalence ~ t + post + t:post

where `t` is the year centred on the intervention and `post` marks the
post-rollout segment. The `post` coefficient is the level change at the
interruption; `t_post` is the change in slope. Those two are the effect
estimates of interest.

**Operationalisation.** The rollout year itself (`ROLLOUT_YEAR`, 2019)
is dropped, exactly as the primary analysis drops it — WHO guidance
landed mid-2019 and national rollout was staggered through 2020, so the
year belongs to neither segment. `t = year - ROLLOUT_YEAR`, so t = 0 is
the interruption and `post`'s coefficient reads as the discontinuity
there rather than as an extrapolation back to year zero. `post = 1` for
years after the buffer. Keeping the same cohort and the same outcome
definition as the primary analysis (`fig04_dtg.build_insti_indicator`)
means a reviewer never has to ask why the two analyses saw different
data.

**Why hand-rolled OLS.** statsmodels would be the obvious import, but
CLAUDE.md's dependency rule asks whether something already in the tree
does the job, and numpy/scipy do: this is a four-parameter linear model
solved by least squares, with standard errors from the usual
(XᵀX)⁻¹σ̂² and a t reference distribution. `_ols` is generic in its
design matrix precisely so the tests can reduce it to one predictor and
check every reported quantity against `scipy.stats.linregress`.

**Known limitation.** This is OLS on a proportion, which is what the
preregistration locked and what the standard segmented-regression
literature does, but the residual variance of a proportion depends on
both the prevalence and the per-year denominator, so the homoscedastic
OLS standard errors are approximate. Per-year `n_sequences` is carried
on the output frame so a reader can see where that bites. The primary
Fisher analysis carries the confirmatory claim; this is secondary and
reported as such.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import polars as pl
from scipy import stats

from .fig04_dtg import ROLLOUT_YEAR, build_insti_indicator

logger = logging.getLogger(__name__)

# Four parameters plus at least one residual degree of freedom.
MIN_OBSERVATIONS: int = 5

# Earliest collection year in the Fig 4 cohort (preregistration §4).
MIN_YEAR: int = 2010

TERMS: tuple[str, ...] = ("intercept", "t", "post", "t_post")


class ITSNotIdentifiableError(RuntimeError):
    """The segmented model cannot be fit from the data supplied."""


@dataclass(frozen=True)
class OLSFit:
    """Raw least-squares output for an arbitrary design matrix."""

    betas: list[float]
    std_errors: list[float]
    t_statistics: list[float]
    p_values: list[float]
    residual_df: int
    r_squared: float


@dataclass(frozen=True)
class ITSCoefficient:
    estimate: float
    std_error: float
    t_statistic: float
    p_value: float
    ci_low: float
    ci_high: float


@dataclass(frozen=True)
class ITSFit:
    """One country's segmented regression."""

    country: str
    n_years: int
    n_years_pre: int
    n_years_post: int
    n_sequences: int
    residual_df: int
    r_squared: float
    coefficients: dict[str, ITSCoefficient]


def _ols(design: np.ndarray, y: np.ndarray) -> OLSFit:
    """Least squares with t-based inference.

    Generic in the design matrix so the tests can check it against
    `scipy.stats.linregress` on the single-predictor case.
    """
    n_obs, n_params = design.shape
    residual_df = n_obs - n_params
    if residual_df < 1:
        raise ITSNotIdentifiableError(
            f"{n_obs} observations cannot fit {n_params} parameters with any "
            f"residual degrees of freedom"
        )
    if np.linalg.matrix_rank(design) < n_params:
        raise ITSNotIdentifiableError(
            "design matrix is rank-deficient — a term is constant or collinear "
            "(usually one segment is missing, so `post` is all zeros)"
        )

    betas, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ betas
    ss_resid = float(resid @ resid)
    sigma2 = ss_resid / residual_df
    xtx_inv = np.linalg.inv(design.T @ design)
    std_errors = np.sqrt(np.clip(np.diag(xtx_inv) * sigma2, 0.0, None))

    with np.errstate(divide="ignore", invalid="ignore"):
        t_stats = np.where(std_errors > 0, betas / std_errors, np.inf)
    p_values = 2.0 * stats.t.sf(np.abs(t_stats), residual_df)

    centred = y - y.mean()
    ss_total = float(centred @ centred)
    r_squared = 1.0 - ss_resid / ss_total if ss_total > 0 else 1.0

    return OLSFit(
        betas=[float(b) for b in betas],
        std_errors=[float(s) for s in std_errors],
        t_statistics=[float(t) for t in t_stats],
        p_values=[float(p) for p in p_values],
        residual_df=int(residual_df),
        r_squared=float(r_squared),
    )


def build_yearly_prevalence(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    rollout_year: int = ROLLOUT_YEAR,
    min_year: int = MIN_YEAR,
) -> pl.DataFrame:
    """Country × year prevalence series for the segmented regression.

    Columns: country, year, t, post, n_sequences, n_resistant, prevalence.
    The rollout year is excluded as an implementation buffer.
    """
    joined = build_insti_indicator(dr_results, curated)
    return (
        joined.filter(
            (pl.col("collection_year") >= min_year)
            & (pl.col("collection_year") != rollout_year)
        )
        .group_by(["country", "collection_year"])
        .agg(
            pl.len().alias("n_sequences"),
            pl.col("is_insti_resistant").sum().alias("n_resistant"),
        )
        .rename({"collection_year": "year"})
        .with_columns(
            (pl.col("n_resistant") / pl.col("n_sequences")).alias("prevalence"),
            (pl.col("year") - rollout_year).alias("t"),
            (pl.col("year") > rollout_year).cast(pl.Int64).alias("post"),
        )
        .sort(["country", "year"])
    )


def fit_country_its(
    yearly: pl.DataFrame, *, alpha: float = 0.05, min_observations: int = MIN_OBSERVATIONS
) -> ITSFit:
    """Fit `prevalence ~ t + post + t:post` for a single country."""
    countries = yearly["country"].unique().to_list()
    if len(countries) != 1:
        raise ValueError(f"expected exactly one country, got {sorted(countries)}")
    country = str(countries[0])

    frame = yearly.sort("year")
    if frame.height < min_observations:
        raise ITSNotIdentifiableError(
            f"{country}: {frame.height} year(s) of data, need at least "
            f"{min_observations} to fit four parameters with residual df"
        )

    t = frame["t"].to_numpy().astype(float)
    post = frame["post"].to_numpy().astype(float)
    y = frame["prevalence"].to_numpy().astype(float)
    design = np.column_stack([np.ones_like(t), t, post, t * post])

    fit = _ols(design, y)
    crit = float(stats.t.ppf(1.0 - alpha / 2.0, fit.residual_df))
    coefficients = {
        term: ITSCoefficient(
            estimate=fit.betas[i],
            std_error=fit.std_errors[i],
            t_statistic=fit.t_statistics[i],
            p_value=fit.p_values[i],
            ci_low=fit.betas[i] - crit * fit.std_errors[i],
            ci_high=fit.betas[i] + crit * fit.std_errors[i],
        )
        for i, term in enumerate(TERMS)
    }

    return ITSFit(
        country=country,
        n_years=frame.height,
        n_years_pre=int((post == 0).sum()),
        n_years_post=int((post == 1).sum()),
        n_sequences=int(frame["n_sequences"].sum()),
        residual_df=fit.residual_df,
        r_squared=fit.r_squared,
        coefficients=coefficients,
    )


def fit_its(yearly: pl.DataFrame, *, alpha: float = 0.05) -> tuple[pl.DataFrame, list[ITSFit]]:
    """Fit every country, skipping those the model cannot identify.

    A country with only one segment, or too few years, is dropped with a
    WARNING rather than failing the run — the Fig 4 cohort has countries
    that entered late and cannot support a segmented fit.
    """
    rows: list[dict[str, object]] = []
    fits: list[ITSFit] = []
    for country in sorted(yearly["country"].unique().to_list()):
        subset = yearly.filter(pl.col("country") == country)
        try:
            fit = fit_country_its(subset, alpha=alpha)
        except ITSNotIdentifiableError as exc:
            logger.warning("ITS skipped for %s: %s", country, exc)
            continue
        fits.append(fit)
        for term, coef in fit.coefficients.items():
            rows.append(
                {
                    "country": fit.country,
                    "term": term,
                    "estimate": coef.estimate,
                    "std_error": coef.std_error,
                    "t_statistic": coef.t_statistic,
                    "p_value": coef.p_value,
                    "ci_low": coef.ci_low,
                    "ci_high": coef.ci_high,
                    "n_years": fit.n_years,
                    "n_years_pre": fit.n_years_pre,
                    "n_years_post": fit.n_years_post,
                    "n_sequences": fit.n_sequences,
                    "residual_df": fit.residual_df,
                    "r_squared": fit.r_squared,
                }
            )

    schema = {
        "country": pl.Utf8,
        "term": pl.Utf8,
        "estimate": pl.Float64,
        "std_error": pl.Float64,
        "t_statistic": pl.Float64,
        "p_value": pl.Float64,
        "ci_low": pl.Float64,
        "ci_high": pl.Float64,
        "n_years": pl.Int64,
        "n_years_pre": pl.Int64,
        "n_years_post": pl.Int64,
        "n_sequences": pl.Int64,
        "residual_df": pl.Int64,
        "r_squared": pl.Float64,
    }
    coefs = pl.from_dicts(rows, schema=schema) if rows else pl.DataFrame(schema=schema)
    logger.info(
        "ITS: fit %d of %d countries (%d skipped as unidentifiable)",
        len(fits),
        yearly["country"].n_unique(),
        yearly["country"].n_unique() - len(fits),
    )
    return coefs.sort(["country", "term"]), fits
