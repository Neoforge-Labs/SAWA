"""Estimate causal effects with dowhy.

For each (treatment, outcome) pair we care about:
  1. Identify the estimand from the DAG (backdoor / IV / frontdoor).
  2. Estimate under the identified estimand with an appropriate model
     (linear regression for continuous, GLM logistic for binary).
  3. Refute the estimate — placebo treatment, random common cause,
     subset validity. An estimate that fails refutation is not
     "wrong" but is meaningfully less trustworthy.

Everything here is thin glue over dowhy. The value we add on top:
  * Consistent CausalEffect output schema, ready to write to parquet.
  * Type-clean adjustment set derivation from our own DAG modules.
  * Graceful degradation when dowhy can't identify (no crash — the
    row is returned with `identified=False` so the report shows the gap).
"""

from __future__ import annotations

import logging
from datetime import date
from typing import Any

import polars as pl

from .dag import DAG_VERSION, DRDag, OutbreakDag, to_dowhy_gml
from .schema import (
    CausalEffect,
    EstimationMethod,
    RefutationResult,
    RefutationVerdict,
)

logger = logging.getLogger(__name__)


REFUTATION_TOLERANCE = 0.20  # 20% shift → refutation fails

DEFAULT_BOOTSTRAP_SAMPLES = 200
DEFAULT_BOOTSTRAP_ALPHA = 0.05  # 95% CI


# ------------------------------------------------------------------
# public API
# ------------------------------------------------------------------


def estimate_dr_effect(
    df: pl.DataFrame,
    mutation: str,
    *,
    outcome: str = "label_binary",
    data_scope: str = "DR/all",
    bootstrap_samples: int = 0,
    bootstrap_alpha: float = DEFAULT_BOOTSTRAP_ALPHA,
    random_state: int = 0,
) -> CausalEffect:
    """Estimate mutation → resistance causal effect on the DR gold frame.

    Expects `df` to contain the mutation column, `outcome` (binary 0/1),
    and the DR-DAG covariates (`subtype`, `country`, `collection_year`).

    Set `bootstrap_samples > 0` (e.g. 200) to compute percentile
    bootstrap CIs on the point estimate. Off by default because
    bootstrap ~doubles runtime and small candidate sets don't need it.
    """
    dag = DRDag()
    graph = dag.build_graph()
    adjustment = dag.minimum_adjustment_for_mutation(mutation)
    return _run_dowhy(
        df=df,
        treatment=mutation,
        outcome=outcome,
        graph=graph,
        adjustment=adjustment,
        method=EstimationMethod.LOGISTIC_REGRESSION,
        data_scope=data_scope,
        bootstrap_samples=bootstrap_samples,
        bootstrap_alpha=bootstrap_alpha,
        random_state=random_state,
    )


def estimate_outbreak_effect(
    df: pl.DataFrame,
    feature: str,
    *,
    outcome: str = "label_incidence",
    data_scope: str = "outbreak/all",
    bootstrap_samples: int = 0,
    bootstrap_alpha: float = DEFAULT_BOOTSTRAP_ALPHA,
    random_state: int = 0,
) -> CausalEffect:
    """Estimate feature → incidence causal effect on the outbreak gold frame."""
    dag = OutbreakDag()
    graph = dag.build_graph()
    adjustment = dag.minimum_adjustment_for_feature(feature)
    return _run_dowhy(
        df=df,
        treatment=feature,
        outcome=outcome,
        graph=graph,
        adjustment=adjustment,
        method=EstimationMethod.LINEAR_REGRESSION,
        data_scope=data_scope,
        bootstrap_samples=bootstrap_samples,
        bootstrap_alpha=bootstrap_alpha,
        random_state=random_state,
    )


def estimate_dr_effects_per_subtype(
    df: pl.DataFrame,
    mutation: str,
    *,
    outcome: str = "label_binary",
    min_samples: int = 100,
    subtype_col: str = "subtype",
    bootstrap_samples: int = 0,
    bootstrap_alpha: float = DEFAULT_BOOTSTRAP_ALPHA,
) -> list[CausalEffect]:
    """Per-subtype stratified estimation — the cheap interaction proxy.

    A full mutation × subtype interaction term in the DAG would require
    rewriting the graph per drug (subtype-specific edges). Stratifying
    the data by subtype and running the estimator on each slice gives
    us the same interpretive power without the DAG rewrite: an effect
    that's 0.9 in A1 and 0.1 in C tells you the effect is subtype-modified.

    Returns one CausalEffect per subtype with at least `min_samples` rows.
    Small subtypes are skipped rather than reported at low power.
    """
    if subtype_col not in df.columns:
        logger.warning(
            "%s not in dataframe; can't stratify — returning single-scope estimate",
            subtype_col,
        )
        return [estimate_dr_effect(df, mutation, outcome=outcome)]

    # Sorted: the caller appends these effects to one list that becomes the
    # effects Parquet and the HTML report, so an unsorted `unique()` reorders
    # both on every run. Sorting is stable across vintages too, which
    # `maintain_order=True` would not be.
    out: list[CausalEffect] = []
    for subtype in sorted(s for s in df[subtype_col].unique().to_list() if s is not None):
        sub_df = df.filter(pl.col(subtype_col) == subtype)
        if sub_df.height < min_samples:
            logger.info("skipping subtype=%s (n=%d < %d)", subtype, sub_df.height, min_samples)
            continue
        eff = estimate_dr_effect(
            sub_df,
            mutation,
            outcome=outcome,
            data_scope=f"DR/subtype={subtype}",
            bootstrap_samples=bootstrap_samples,
            bootstrap_alpha=bootstrap_alpha,
        )
        out.append(eff)
    return out


# ------------------------------------------------------------------
# internals
# ------------------------------------------------------------------


def _run_dowhy(
    *,
    df: pl.DataFrame,
    treatment: str,
    outcome: str,
    graph,
    adjustment: set[str],
    method: EstimationMethod,
    data_scope: str,
    bootstrap_samples: int = 0,
    bootstrap_alpha: float = DEFAULT_BOOTSTRAP_ALPHA,
    random_state: int = 0,
) -> CausalEffect:
    """Run identify → estimate → refute. Returns CausalEffect either way."""
    if treatment not in df.columns or outcome not in df.columns:
        return _unidentified_row(
            treatment,
            outcome,
            adjustment,
            data_scope,
            method,
            reason=f"missing column ({treatment} or {outcome})",
        )

    pandas_df = (
        df.select([treatment, outcome] + [c for c in adjustment if c in df.columns])
        .drop_nulls()
        .to_pandas()
    )

    if len(pandas_df) < 30:
        return _unidentified_row(
            treatment,
            outcome,
            adjustment,
            data_scope,
            method,
            reason=f"insufficient rows: {len(pandas_df)}",
        )

    try:
        from dowhy import CausalModel  # type: ignore[import-not-found]
    except ImportError:
        logger.warning("dowhy not installed; skipping causal estimation")
        return _unidentified_row(
            treatment,
            outcome,
            adjustment,
            data_scope,
            method,
            reason="dowhy not installed",
        )

    gml = to_dowhy_gml(graph)

    try:
        model = CausalModel(
            data=pandas_df,
            treatment=treatment,
            outcome=outcome,
            graph=gml,
        )
        identified = model.identify_effect(proceed_when_unidentifiable=False)
    except Exception as e:
        logger.warning("identify_effect failed for %s → %s: %s", treatment, outcome, e)
        return _unidentified_row(
            treatment,
            outcome,
            adjustment,
            data_scope,
            method,
            reason=f"unidentifiable: {e}",
        )

    try:
        estimate = model.estimate_effect(
            identified,
            method_name=method.value,
        )
    except Exception as e:
        logger.warning("estimate_effect failed for %s → %s: %s", treatment, outcome, e)
        return _unidentified_row(
            treatment,
            outcome,
            adjustment,
            data_scope,
            method,
            reason=f"estimation error: {e}",
            identified=True,
        )

    point = float(getattr(estimate, "value", 0.0) or 0.0)
    ci_low, ci_high, p_val = _bootstrap_ci(
        pandas_df=pandas_df,
        treatment=treatment,
        outcome=outcome,
        gml=gml,
        method=method,
        point=point,
        n_samples=bootstrap_samples,
        alpha=bootstrap_alpha,
        random_state=random_state,
    )
    refutations = _run_refutations(model, identified, estimate, point)

    return CausalEffect(
        treatment=treatment,
        outcome=outcome,
        adjustment_set=sorted(adjustment),
        data_scope=data_scope,
        n_observations=len(pandas_df),
        method=method,
        estimate=point,
        ci_low=ci_low,
        ci_high=ci_high,
        p_value=p_val,
        identified=True,
        identification_criterion=_identification_criterion(identified),
        refutations=refutations,
        dag_version=DAG_VERSION,
        estimated_at=date.today(),
    )


def _bootstrap_ci(
    *,
    pandas_df,
    treatment: str,
    outcome: str,
    gml: str,
    method: EstimationMethod,
    point: float,
    n_samples: int,
    alpha: float,
    random_state: int,
) -> tuple[float | None, float | None, float | None]:
    """Percentile-bootstrap CI + naive two-sided p-value about 0.

    Returns (None, None, None) when `n_samples == 0` — we don't
    silently substitute a degenerate estimate.

    The p-value is a rough sign test: fraction of bootstrap estimates
    that flip sign vs the point estimate, times 2 (two-sided). It's a
    non-standard construction but good enough for "does the CI include
    zero" style reasoning; report intervals, not p-values, in the paper.
    """
    if n_samples <= 0:
        return (None, None, None)

    try:
        import numpy as np
        from dowhy import CausalModel  # type: ignore[import-not-found]
    except ImportError:
        return (None, None, None)

    rng = np.random.default_rng(random_state)
    n = len(pandas_df)
    estimates: list[float] = []
    for _ in range(n_samples):
        idx = rng.integers(0, n, size=n)
        boot_df = pandas_df.iloc[idx].reset_index(drop=True)
        # Guard: bootstrap sample might have a single class in `treatment`
        # or `outcome`; skip that resample (do not append).
        if boot_df[treatment].nunique() < 2 or boot_df[outcome].nunique() < 2:
            continue
        try:
            m = CausalModel(data=boot_df, treatment=treatment, outcome=outcome, graph=gml)
            ident = m.identify_effect(proceed_when_unidentifiable=True)
            est = m.estimate_effect(ident, method_name=method.value)
            v = getattr(est, "value", None)
            if v is not None:
                estimates.append(float(v))
        except Exception as e:
            logger.debug("bootstrap resample failed: %s", e)

    if len(estimates) < max(20, n_samples // 5):
        # Too few successful resamples to draw an honest interval.
        logger.warning(
            "bootstrap yielded only %d/%d valid resamples; skipping CI",
            len(estimates),
            n_samples,
        )
        return (None, None, None)

    arr = np.array(estimates)
    lo = float(np.quantile(arr, alpha / 2))
    hi = float(np.quantile(arr, 1 - alpha / 2))

    # Sign-flip p-value (two-sided) relative to the point estimate direction.
    if point > 0:
        p = 2.0 * float((arr <= 0).mean())
    elif point < 0:
        p = 2.0 * float((arr >= 0).mean())
    else:
        p = 1.0
    p = min(p, 1.0)

    return (lo, hi, p)


def _run_refutations(
    model, identified, estimate, original_estimate: float
) -> list[RefutationResult]:
    """Best-effort refutation — placebo, random common cause, subset.

    Refuters can be slow (bootstrap-y) so we cap at three and swallow
    failures per-refuter rather than aborting the whole row.
    """
    refuters = [
        ("placebo_treatment_refuter", "placebo_treatment_refuter"),
        ("random_common_cause", "random_common_cause"),
        ("data_subset_refuter", "data_subset_refuter"),
    ]
    out: list[RefutationResult] = []
    for label, method in refuters:
        try:
            r = model.refute_estimate(identified, estimate, method_name=method)
            new_val = float(getattr(r, "new_effect", original_estimate) or 0.0)
            p = getattr(r, "refutation_result", None)
            p_val: float | None = None
            if isinstance(p, dict):
                p_val = p.get("p_value")
            verdict = _score_refutation(original=original_estimate, new=new_val, method=label)
            out.append(
                RefutationResult(
                    method=label,
                    original_estimate=original_estimate,
                    new_estimate=new_val,
                    p_value=float(p_val) if p_val is not None else None,
                    verdict=verdict,
                )
            )
        except Exception as e:
            logger.debug("refuter %s failed: %s", label, e)
            out.append(
                RefutationResult(
                    method=label,
                    original_estimate=original_estimate,
                    new_estimate=original_estimate,
                    p_value=None,
                    verdict=RefutationVerdict.INCONCLUSIVE,
                )
            )
    return out


def _score_refutation(*, original: float, new: float, method: str) -> RefutationVerdict:
    """Bucket a refuter's numeric result into pass/fail/inconclusive.

    Rules:
      - `placebo_treatment_refuter`: new estimate should be ~0 → passes
        iff |new| ≤ tolerance * |original|.
      - `random_common_cause` and `data_subset_refuter`: new estimate
        should be close to original → passes iff |Δ| ≤ tolerance * |original|.
    """
    if original == 0:
        # Zero estimate — placebo passes trivially; others inconclusive.
        if method == "placebo_treatment_refuter":
            return (
                RefutationVerdict.PASSES if abs(new) < 1e-6 else RefutationVerdict.INCONCLUSIVE
            )
        return RefutationVerdict.INCONCLUSIVE

    if method == "placebo_treatment_refuter":
        return (
            RefutationVerdict.PASSES
            if abs(new) <= REFUTATION_TOLERANCE * abs(original)
            else RefutationVerdict.FAILS
        )
    # Sensitivity-style refuters: original vs new should be close.
    rel_change = abs(new - original) / abs(original)
    return (
        RefutationVerdict.PASSES
        if rel_change <= REFUTATION_TOLERANCE
        else RefutationVerdict.FAILS
    )


def _identification_criterion(identified: Any) -> str | None:
    """Extract dowhy's chosen criterion (backdoor / iv / frontdoor)."""
    for attr in ("estimand_type", "identifier_method", "criterion"):
        v = getattr(identified, attr, None)
        if v:
            return str(v)
    return None


def _unidentified_row(
    treatment: str,
    outcome: str,
    adjustment: set[str],
    data_scope: str,
    method: EstimationMethod,
    *,
    reason: str,
    identified: bool = False,
) -> CausalEffect:
    logger.info("unidentified/skipped %s → %s: %s", treatment, outcome, reason)
    return CausalEffect(
        treatment=treatment,
        outcome=outcome,
        adjustment_set=sorted(adjustment),
        data_scope=f"{data_scope} [skipped: {reason}]",
        n_observations=0,
        method=method,
        estimate=0.0,
        identified=identified,
        identification_criterion=None,
        refutations=[],
        dag_version=DAG_VERSION,
        estimated_at=date.today(),
    )


# ------------------------------------------------------------------
# batch drivers — used by the flow
# ------------------------------------------------------------------


def estimate_all_dr_effects(
    df: pl.DataFrame,
    mutations: list[str] | None = None,
    *,
    drug: str | None = None,
    bootstrap_samples: int = 0,
    per_subtype: bool = False,
    subtype_min_samples: int = 100,
) -> list[CausalEffect]:
    """For every mutation in the DR DAG, estimate its causal effect on
    binary resistance in `df`. Optionally filtered to a single drug.

    * `bootstrap_samples > 0` → attach percentile-bootstrap 95% CIs.
    * `per_subtype=True` → additionally emit one row per subtype (data
      scope `DR/drug=X/subtype=Y`), giving you the interaction proxy.
    """
    if drug:
        df = df.filter(pl.col("drug") == drug)
        scope_suffix = f"/drug={drug}"
    else:
        scope_suffix = ""

    # Derive binary target from Sierra 5-level label (mirrors models/dr).
    if "label_binary" not in df.columns and "label_level" in df.columns:
        df = df.with_columns(
            label_binary=pl.col("label_level").is_in(["high", "intermediate"]).cast(pl.Int8)
        )

    dag = DRDag()
    to_check = mutations or list(dag.mutation_features)
    effects: list[CausalEffect] = []
    for mutation in to_check:
        effects.append(
            estimate_dr_effect(
                df,
                mutation,
                data_scope=f"DR{scope_suffix}",
                bootstrap_samples=bootstrap_samples,
            )
        )
        if per_subtype:
            per_st = estimate_dr_effects_per_subtype(
                df,
                mutation,
                min_samples=subtype_min_samples,
                bootstrap_samples=bootstrap_samples,
            )
            # Prefix drug scope for cross-drug disambiguation.
            for eff in per_st:
                if drug and not eff.data_scope.startswith(f"DR/drug={drug}"):
                    eff = eff.model_copy(
                        update={
                            "data_scope": eff.data_scope.replace("DR/", f"DR/drug={drug}/")
                        }
                    )
                effects.append(eff)
    return effects


def estimate_all_outbreak_effects(
    df: pl.DataFrame,
    features: list[str] | None = None,
    *,
    bootstrap_samples: int = 0,
) -> list[CausalEffect]:
    """Estimate each outbreak feature's causal effect on incidence."""
    candidate_features = features or [
        "n_new_cluster_members",
        "n_active_clusters",
        "top_cluster_size",
        "mean_growth_rate",
        "mean_subtype_entropy",
        "subtype_a1_fraction",
        "subtype_c_fraction",
        "subtype_d_fraction",
        # deposit_rate is included so its own causal effect on incidence
        # (should be ~0 if properly modelled) is visible in the report.
        "deposit_rate",
    ]
    return [
        estimate_outbreak_effect(
            df, f, data_scope="outbreak", bootstrap_samples=bootstrap_samples
        )
        for f in candidate_features
    ]
