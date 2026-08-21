"""Pure aggregation for Fig 4 — DTG-era INSTI resistance case study.

Given the DR results parquet (Sierra output) and the curated silver
parquet (metadata: country + subtype + year), produce a stratum-level
frame with:

  * epoch: "pre" (2010-2018) or "post" (2020-latest); 2019 excluded
  * country, subtype
  * n_sequences, n_insti_resistant
  * prev, prev_ci_low, prev_ci_high  (**Clopper–Pearson exact 95 % CI**)

plus a paired-comparison frame with:

  * (country, subtype,
     delta_prev, delta_ci_low, delta_ci_high,   ← Newcombe hybrid, kept
                                                  for the delta CI plot
     fisher_p, fisher_or,                       ← **Fisher exact test**
     bh_q,                                       ← **BH-FDR adjusted p**
     tier,                                       ← "primary" (n ≥ 30
                                                  per epoch) or
                                                  "secondary" (n ≥ 5)
     signal,                                     ← retained in headline
     signal_direction)

**Statistical approach — Option C (preregistered rationale).**

The primary hypothesis test is **Fisher's exact test** on the 2 × 2
pre-vs-post × resistant-vs-susceptible contingency table. Fisher's
exact is preferred over the two-proportion z-test because it makes no
normal-approximation assumption and remains valid down to expected
counts of 1 per cell — critical for the small non-Ugandan strata
(Rwanda, Burundi, Ethiopia) that entered the cohort at v0.2.0.

Per-epoch prevalence CIs use **Clopper–Pearson exact binomial
intervals**, which are conservative but always valid regardless of n.
The Wilson score interval, previously used, is well-behaved down to
n ≈ 5 for prevalence in [0.05, 0.95] but can be anti-conservative at
extreme p̂. Clopper–Pearson trades tightness for correctness — a
deliberate choice for a paper whose headline is a public-health
claim.

**Tiered inclusion.**

* **Primary tier**: n_pre ≥ 30 AND n_post ≥ 30. This matches the
  preregistered threshold (Fig 4 preregistration §6). Signals in this
  tier count toward the confirmatory headline.
* **Secondary tier**: n_pre ≥ 5 AND n_post ≥ 5, but not both ≥ 30.
  Signals in this tier are reported as *exploratory*, do not count
  toward the preregistered FDR budget, and are labelled accordingly
  in the manuscript.

**Multiple testing.** BH-FDR at q = 0.05 is applied within the
primary-tier set. Secondary-tier p-values are BH-adjusted separately
and reported for context.

Preregistered analysis plan lives in
`paper1_tools/preregistration_fig04.md`. Any deviation from that
plan is exploratory and must be labelled as such in the manuscript.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Literal

import polars as pl

logger = logging.getLogger(__name__)

# Sierra levels that count as "resistant" — matches Fig 2 convention.
RESISTANT_LEVELS: frozenset[str] = frozenset({"high", "intermediate"})

# Drug class we're interested in for this case study.
INSTI_CLASS: str = "INSTI"

# Preregistered epoch boundaries. 2019 is deliberately excluded as an
# implementation-year buffer (WHO DTG guidance issued mid-2019 and
# rollout was staggered across countries through 2020).
PRE_YEARS: tuple[int, int] = (2010, 2018)
POST_YEARS: tuple[int, int] = (2020, 9999)  # "latest" — no upper bound in data
ROLLOUT_YEAR: int = 2019

# Preregistered threshold for the primary hypothesis test.
MIN_EPOCH_N_PRIMARY: int = 30
# Exploratory secondary threshold — reported but does not count toward
# the preregistered FDR budget.
MIN_EPOCH_N_SECONDARY: int = 5

# BH-FDR level — preregistered.
FDR_Q: float = 0.05

Tier = Literal["primary", "secondary"]


@dataclass(frozen=True)
class DTGSummary:
    """Compact "Fig 4 built OK" summary."""

    n_pre_seqs: int
    n_post_seqs: int
    n_strata_primary: int
    n_strata_secondary: int
    n_signal_primary: int
    n_signal_secondary: int
    n_headline_after_fdr: int


def _epoch(year: int) -> str | None:
    if PRE_YEARS[0] <= year <= PRE_YEARS[1]:
        return "pre"
    if POST_YEARS[0] <= year <= POST_YEARS[1]:
        return "post"
    return None


# ------------------------------------------------------------------
# Exact / Newcombe interval + Fisher's exact test
# ------------------------------------------------------------------


def _clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float, float]:
    """Clopper–Pearson exact binomial 100·(1−α)% CI.

    Returns (point, lo, hi). Uses the Beta-distribution characterisation:
        lo = BetaInv(α/2; k, n−k+1)   (0 when k = 0)
        hi = BetaInv(1−α/2; k+1, n−k) (1 when k = n)

    Handles n = 0 → (0, 0, 0).
    """
    if n == 0:
        return 0.0, 0.0, 0.0
    phat = k / n
    try:
        from scipy.stats import beta
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy required for Clopper–Pearson CI; add scipy to project deps"
        ) from e
    lo = float(beta.ppf(alpha / 2.0, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)) if k < n else 1.0
    lo = min(lo, phat)  # Guarantee lo ≤ phat ≤ hi against fp noise.
    hi = max(hi, phat)
    return phat, lo, hi


def _newcombe_delta_ci(
    k_post: int, n_post: int, k_pre: int, n_pre: int, alpha: float = 0.05
) -> tuple[float, float, float]:
    """Newcombe hybrid-score CI for the difference of two proportions.

    Uses Wilson individual CIs combined via Newcombe (1998) Method 10.
    Kept for the delta-plot CI; the hypothesis test itself uses
    Fisher's exact (below).
    """
    if n_post == 0 or n_pre == 0:
        return 0.0, 0.0, 0.0
    z = _z_from_alpha(alpha)
    p1, l1, u1 = _wilson_ci(k_post, n_post, z=z)
    p2, l2, u2 = _wilson_ci(k_pre, n_pre, z=z)
    delta = p1 - p2
    lower = delta - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    upper = delta + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return delta, lower, upper


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    """Wilson score CI — kept as internal helper for Newcombe combination."""
    if n == 0:
        return 0.0, 0.0, 0.0
    phat = k / n
    denom = 1 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = (z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))) / denom
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    lo = min(lo, phat)
    hi = max(hi, phat)
    return phat, lo, hi


def _fisher_exact(k_post: int, n_post: int, k_pre: int, n_pre: int) -> tuple[float, float]:
    """Fisher's exact test on the 2×2 (epoch × resistant) contingency.

    Returns (odds_ratio, two-sided_p). The odds ratio is
    (k_post × (n_pre − k_pre)) / (k_pre × (n_post − k_post))
    with the standard 0.5 correction when any cell is 0, so the OR
    remains finite for reporting.

    Returns (nan, 1.0) if either n is 0.
    """
    if n_post == 0 or n_pre == 0:
        return float("nan"), 1.0
    try:
        from scipy.stats import fisher_exact
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy required for Fisher's exact test; add scipy to project deps"
        ) from e
    table = [[k_post, n_post - k_post], [k_pre, n_pre - k_pre]]
    odds, p = fisher_exact(table, alternative="two-sided")
    if not math.isfinite(odds):
        # Add 0.5 to every cell (Haldane–Anscombe) to keep OR finite for
        # reporting; the p-value from Fisher exact is unaffected.
        odds = ((k_post + 0.5) * (n_pre - k_pre + 0.5)) / (
            (k_pre + 0.5) * (n_post - k_post + 0.5)
        )
    return float(odds), float(p)


def _z_from_alpha(alpha: float) -> float:
    """Two-sided z for confidence level 1 − α, without invoking scipy for α=0.05."""
    if abs(alpha - 0.05) < 1e-12:
        return 1.959963984540054  # scipy.stats.norm.ppf(0.975)
    from scipy.stats import norm

    return float(norm.ppf(1 - alpha / 2.0))


def _benjamini_hochberg(pvals: list[float], q: float) -> list[float]:
    """Return BH-adjusted q-values (one per input p, same order).

    Standard Benjamini–Hochberg step-up. Returns q-values so that the
    caller can decide the significance threshold; a stratum is
    BH-significant iff its q-value ≤ q.

    Empty input returns an empty list.
    """
    if not pvals:
        return []
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    ranked = sorted(pvals)
    # Compute the monotone-adjusted q-value in reverse rank order.
    q_sorted = [0.0] * m
    prev = 1.0
    for i in range(m - 1, -1, -1):
        raw = ranked[i] * m / (i + 1)
        prev = min(prev, raw)
        q_sorted[i] = min(prev, 1.0)
    # Un-sort back to original order.
    q_by_orig = [0.0] * m
    for rank_i, orig_i in enumerate(order):
        q_by_orig[orig_i] = q_sorted[rank_i]
    return q_by_orig


# ------------------------------------------------------------------
# Aggregator
# ------------------------------------------------------------------


def aggregate_dtg_landscape(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    min_epoch_n_primary: int = MIN_EPOCH_N_PRIMARY,
    min_epoch_n_secondary: int = MIN_EPOCH_N_SECONDARY,
    fdr_q: float = FDR_Q,
) -> tuple[pl.DataFrame, pl.DataFrame, DTGSummary]:
    """Return (per-epoch frame, delta frame, summary).

    Delta frame columns:
      country, subtype, n_pre, n_post, prev_pre, prev_post,
      delta_prev, delta_ci_low, delta_ci_high,
      fisher_or, fisher_p, bh_q,
      tier, signal, signal_direction

    Strata with n_pre < min_epoch_n_secondary OR n_post < min_epoch_n_secondary
    are omitted from the delta frame entirely (still present in the
    per-epoch frame).
    """
    for name, df, cols in (
        ("dr_results", dr_results, {"accession", "drug_class", "level"}),
        ("curated", curated, {"accession", "country", "collection_year", "subtype"}),
    ):
        missing = cols - set(df.columns)
        if missing:
            raise ValueError(f"{name} missing required columns: {sorted(missing)}")

    # Per-sequence "any INSTI resistant" indicator.
    per_seq = (
        dr_results.filter(pl.col("drug_class") == INSTI_CLASS)
        .with_columns(pl.col("level").is_in(list(RESISTANT_LEVELS)).alias("_is_r"))
        .group_by("accession")
        .agg(pl.col("_is_r").any().alias("is_insti_resistant"))
    )

    joined = per_seq.join(
        curated.select(["accession", "country", "collection_year", "subtype"]),
        on="accession",
        how="inner",
    ).filter(
        pl.col("country").is_not_null()
        & pl.col("collection_year").is_not_null()
        & pl.col("subtype").is_not_null()
    )

    joined = joined.with_columns(
        pl.col("collection_year").map_elements(_epoch, return_dtype=pl.Utf8).alias("epoch")
    ).filter(pl.col("epoch").is_not_null())

    per_epoch = (
        joined.group_by(["country", "subtype", "epoch"])
        .agg(
            n_sequences=pl.len(),
            n_insti_resistant=pl.col("is_insti_resistant").sum(),
        )
        .sort(["country", "subtype", "epoch"])
    )

    ci_rows = [
        _clopper_pearson_ci(int(r["n_insti_resistant"]), int(r["n_sequences"]))
        for r in per_epoch.iter_rows(named=True)
    ]
    per_epoch = per_epoch.with_columns(
        pl.Series("prev", [r[0] for r in ci_rows]),
        pl.Series("prev_ci_low", [r[1] for r in ci_rows]),
        pl.Series("prev_ci_high", [r[2] for r in ci_rows]),
    )

    # Pair pre/post per (country, subtype) — inline pivot without on-conflicts.
    wide = per_epoch.pivot(
        values=["n_sequences", "n_insti_resistant"],
        index=["country", "subtype"],
        on="epoch",
    )
    delta_rows: list[dict] = []
    for row in wide.iter_rows(named=True):
        n_pre = int(row.get("n_sequences_pre") or 0)
        n_post = int(row.get("n_sequences_post") or 0)
        if n_pre < min_epoch_n_secondary or n_post < min_epoch_n_secondary:
            continue
        k_pre = int(row.get("n_insti_resistant_pre") or 0)
        k_post = int(row.get("n_insti_resistant_post") or 0)

        d, lo, hi = _newcombe_delta_ci(k_post, n_post, k_pre, n_pre)
        odds, p = _fisher_exact(k_post, n_post, k_pre, n_pre)
        tier: Tier = (
            "primary"
            if (n_pre >= min_epoch_n_primary and n_post >= min_epoch_n_primary)
            else "secondary"
        )
        prev_pre = k_pre / n_pre if n_pre else 0.0
        prev_post = k_post / n_post if n_post else 0.0
        direction = (
            "up" if prev_post > prev_pre else ("down" if prev_post < prev_pre else "flat")
        )

        delta_rows.append(
            {
                "country": row["country"],
                "subtype": row["subtype"],
                "n_pre": n_pre,
                "n_post": n_post,
                "prev_pre": prev_pre,
                "prev_post": prev_post,
                "delta_prev": d,
                "delta_ci_low": lo,
                "delta_ci_high": hi,
                "fisher_or": odds,
                "fisher_p": p,
                "tier": tier,
                "signal_direction": direction,
            }
        )

    if not delta_rows:
        delta = pl.DataFrame(
            schema={
                "country": pl.Utf8,
                "subtype": pl.Utf8,
                "n_pre": pl.Int64,
                "n_post": pl.Int64,
                "prev_pre": pl.Float64,
                "prev_post": pl.Float64,
                "delta_prev": pl.Float64,
                "delta_ci_low": pl.Float64,
                "delta_ci_high": pl.Float64,
                "fisher_or": pl.Float64,
                "fisher_p": pl.Float64,
                "bh_q": pl.Float64,
                "tier": pl.Utf8,
                "signal": pl.Boolean,
                "signal_direction": pl.Utf8,
            }
        )
    else:
        # BH-FDR is applied *within each tier* — the preregistered FDR
        # budget belongs to primary; secondary q-values are for context.
        primary_idx = [i for i, r in enumerate(delta_rows) if r["tier"] == "primary"]
        secondary_idx = [i for i, r in enumerate(delta_rows) if r["tier"] == "secondary"]
        primary_qs = _benjamini_hochberg(
            [delta_rows[i]["fisher_p"] for i in primary_idx], q=fdr_q
        )
        secondary_qs = _benjamini_hochberg(
            [delta_rows[i]["fisher_p"] for i in secondary_idx], q=fdr_q
        )
        for i, q in zip(primary_idx, primary_qs, strict=True):
            delta_rows[i]["bh_q"] = q
        for i, q in zip(secondary_idx, secondary_qs, strict=True):
            delta_rows[i]["bh_q"] = q

        # Signal criterion: post > pre AND BH-adjusted q ≤ fdr_q.
        for r in delta_rows:
            r["signal"] = bool((r["signal_direction"] == "up") and (r["bh_q"] <= fdr_q))

        delta = pl.from_dicts(delta_rows).sort(["tier", "country", "subtype"])

    n_pre_seqs = int(joined.filter(pl.col("epoch") == "pre").height)
    n_post_seqs = int(joined.filter(pl.col("epoch") == "post").height)
    n_primary = int(delta.filter(pl.col("tier") == "primary").height) if delta.height else 0
    n_secondary = (
        int(delta.filter(pl.col("tier") == "secondary").height) if delta.height else 0
    )
    n_sig_primary = (
        int(delta.filter((pl.col("tier") == "primary") & pl.col("signal")).height)
        if delta.height
        else 0
    )
    n_sig_secondary = (
        int(delta.filter((pl.col("tier") == "secondary") & pl.col("signal")).height)
        if delta.height
        else 0
    )

    summary = DTGSummary(
        n_pre_seqs=n_pre_seqs,
        n_post_seqs=n_post_seqs,
        n_strata_primary=n_primary,
        n_strata_secondary=n_secondary,
        n_signal_primary=n_sig_primary,
        n_signal_secondary=n_sig_secondary,
        n_headline_after_fdr=n_sig_primary,
    )
    logger.info(
        "DTG landscape (Option C): pre=%d, post=%d — %d primary (%d signal, BH-FDR q<=%.2g), "
        "%d secondary (%d signal, exploratory)",
        n_pre_seqs,
        n_post_seqs,
        n_primary,
        n_sig_primary,
        fdr_q,
        n_secondary,
        n_sig_secondary,
    )
    return per_epoch, delta, summary


__all__ = [
    "FDR_Q",
    "INSTI_CLASS",
    "MIN_EPOCH_N_PRIMARY",
    "MIN_EPOCH_N_SECONDARY",
    "POST_YEARS",
    "PRE_YEARS",
    "RESISTANT_LEVELS",
    "ROLLOUT_YEAR",
    "DTGSummary",
    "aggregate_dtg_landscape",
]
