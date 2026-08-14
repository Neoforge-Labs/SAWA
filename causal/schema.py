"""Canonical schemas for causal analysis outputs.

One `CausalEffect` row per (treatment, outcome, adjustment_set) triple.
Refutation results ride along in the same row so a downstream reader can
tell "did this estimate survive placebo / random-common-cause / subset
challenges" without a second join.

We deliberately don't collapse effects across DAGs (e.g. per-drug for DR,
per-country×year for outbreak). Every row states its data scope in
`data_scope` so the "feature X matters for drug Y in subtype Z" question
has a natural row identity.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum

from pydantic import BaseModel, Field

CAUSAL_SCHEMA_VERSION = "0.1.0"


class EstimationMethod(StrEnum):
    LINEAR_REGRESSION = "backdoor.linear_regression"
    LOGISTIC_REGRESSION = "backdoor.generalized_linear_model"
    PROPENSITY_MATCHING = "backdoor.propensity_score_matching"
    PROPENSITY_STRATIFICATION = "backdoor.propensity_score_stratification"


class RefutationVerdict(StrEnum):
    """Whether a refuter's result is consistent with a real effect.

    dowhy returns numeric p-values; we bucket for readability.
    """

    PASSES = "passes"          # estimate survives the challenge (real signal)
    FAILS = "fails"            # estimate changes meaningfully → spurious
    INCONCLUSIVE = "inconclusive"  # refuter couldn't decide (small n, degenerate)


class RefutationResult(BaseModel):
    method: str                                   # "placebo_treatment_refuter", etc.
    original_estimate: float
    new_estimate: float
    p_value: float | None = None
    verdict: RefutationVerdict


class CausalEffect(BaseModel):
    """One backdoor-adjusted causal effect estimate."""

    # --- what's being estimated ---
    treatment: str
    outcome: str
    adjustment_set: list[str] = Field(default_factory=list)
    data_scope: str = Field(
        ...,
        description="Free-text scope, e.g. 'DR/3TC/subtype=A1' or 'outbreak/country=KE'",
    )

    # --- estimate ---
    n_observations: int
    method: EstimationMethod
    estimate: float
    ci_low: float | None = None
    ci_high: float | None = None
    p_value: float | None = None

    # --- identification ---
    identified: bool
    identification_criterion: str | None = None   # "backdoor" / "iv" / "frontdoor"

    # --- refutation ---
    refutations: list[RefutationResult] = Field(default_factory=list)

    # --- provenance ---
    dag_version: str
    estimated_at: date
    schema_version: str = CAUSAL_SCHEMA_VERSION
