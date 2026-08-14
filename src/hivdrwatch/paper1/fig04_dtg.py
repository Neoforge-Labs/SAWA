"""Pure aggregation for Fig 4 — DTG-era INSTI resistance case study.

Given the DR results parquet (Sierra output) and the curated silver
parquet (metadata: country + subtype + year), produce a stratum-level
frame with:

  * epoch: "pre" (2010-2018) or "post" (2020-latest); 2019 excluded
  * country, subtype
  * n_sequences, n_insti_resistant
  * prev, prev_ci_low, prev_ci_high  (Wilson score interval, 95%)

plus a paired-comparison frame with:
  * (country, subtype, delta_prev, delta_ci_low, delta_ci_high, signal)

`signal` = True iff post > pre AND the 95% CI of the delta excludes 0.

**Preregistered analysis plan lives in
`paper1_tools/preregistration_fig04.md`.** Any deviation from that
plan is exploratory and must be labeled as such in the manuscript.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

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
POST_YEARS: tuple[int, int] = (2020, 9999)   # "latest" — no upper bound in data
ROLLOUT_YEAR: int = 2019

# Minimum N per stratum before we score prevalence. Preregistered:
# stratum-level tests need ≥ 30 sequences per epoch to be admissible.
MIN_EPOCH_N: int = 30


@dataclass(frozen=True)
class DTGSummary:
    """Compact "Fig 4 built OK" summary."""

    n_pre_seqs: int
    n_post_seqs: int
    n_strata_scored: int
    n_strata_with_signal: int


def _epoch(year: int) -> str | None:
    if PRE_YEARS[0] <= year <= PRE_YEARS[1]:
        return "pre"
    if POST_YEARS[0] <= year <= POST_YEARS[1]:
        return "post"
    return None


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    """Wilson score interval for a binomial proportion.

    Returns (point, ci_low, ci_high). Handles n=0 by returning (0, 0, 0)
    — no data, no signal, no interval.
    """
    if n == 0:
        return 0.0, 0.0, 0.0
    phat = k / n
    denom = 1 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = (z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))) / denom
    return phat, max(0.0, center - half), min(1.0, center + half)


def _delta_ci(
    k_post: int, n_post: int, k_pre: int, n_pre: int, z: float = 1.96
) -> tuple[float, float, float]:
    """Newcombe hybrid score interval for the difference of two proportions.

    We use the Wilson individual intervals and combine them via the
    method from Newcombe (1998, Method 10). Handles n=0 by returning
    (0, 0, 0).
    """
    if n_post == 0 or n_pre == 0:
        return 0.0, 0.0, 0.0
    p1, l1, u1 = _wilson_ci(k_post, n_post, z)
    p2, l2, u2 = _wilson_ci(k_pre, n_pre, z)
    delta = p1 - p2
    # Newcombe hybrid: lower = delta - sqrt((p1-l1)^2 + (u2-p2)^2)
    lower = delta - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    upper = delta + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return delta, lower, upper


def aggregate_dtg_landscape(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    min_epoch_n: int = MIN_EPOCH_N,
) -> tuple[pl.DataFrame, pl.DataFrame, DTGSummary]:
    """Return (per-epoch frame, delta frame, summary).

    Inputs (only consumed columns listed):
      dr_results: [accession, drug_class, level]
      curated:    [accession, country, collection_year, subtype]

    Per-epoch frame columns:
      [country, subtype, epoch, n_sequences, n_insti_resistant,
       prev, prev_ci_low, prev_ci_high]

    Delta frame columns:
      [country, subtype, delta_prev, delta_ci_low, delta_ci_high, signal,
       n_pre, n_post]

    Strata with n_pre < min_epoch_n or n_post < min_epoch_n are omitted
    from the delta frame (still present in the per-epoch frame with
    prev but null delta).
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
        .with_columns(
            pl.col("level").is_in(list(RESISTANT_LEVELS)).alias("_is_r"),
        )
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

    # Assign epoch.
    joined = joined.with_columns(
        pl.col("collection_year")
        .map_elements(_epoch, return_dtype=pl.Utf8)
        .alias("epoch")
    ).filter(pl.col("epoch").is_not_null())

    # Per-(country, subtype, epoch) aggregation.
    per_epoch = (
        joined.group_by(["country", "subtype", "epoch"])
        .agg(
            n_sequences=pl.len(),
            n_insti_resistant=pl.col("is_insti_resistant").sum(),
        )
        .sort(["country", "subtype", "epoch"])
    )

    # Wilson CI per row.
    def _wilson_row(row: dict) -> dict:
        p, lo, hi = _wilson_ci(int(row["n_insti_resistant"]), int(row["n_sequences"]))
        return {"prev": p, "prev_ci_low": lo, "prev_ci_high": hi}

    ci_rows = [_wilson_row(row) for row in per_epoch.iter_rows(named=True)]
    per_epoch = per_epoch.with_columns(
        pl.Series("prev", [r["prev"] for r in ci_rows]),
        pl.Series("prev_ci_low", [r["prev_ci_low"] for r in ci_rows]),
        pl.Series("prev_ci_high", [r["prev_ci_high"] for r in ci_rows]),
    )

    # Delta frame: pair pre/post per (country, subtype).
    wide = per_epoch.pivot(
        values=["n_sequences", "n_insti_resistant"],
        index=["country", "subtype"],
        on="epoch",
    )
    delta_rows: list[dict] = []
    for row in wide.iter_rows(named=True):
        n_pre = int(row.get("n_sequences_pre") or 0)
        n_post = int(row.get("n_sequences_post") or 0)
        if n_pre < min_epoch_n or n_post < min_epoch_n:
            continue
        k_pre = int(row.get("n_insti_resistant_pre") or 0)
        k_post = int(row.get("n_insti_resistant_post") or 0)
        d, lo, hi = _delta_ci(k_post, n_post, k_pre, n_pre)
        signal = (d > 0) and (lo > 0)
        delta_rows.append(
            {
                "country": row["country"],
                "subtype": row["subtype"],
                "n_pre": n_pre,
                "n_post": n_post,
                "prev_pre": k_pre / n_pre if n_pre else 0.0,
                "prev_post": k_post / n_post if n_post else 0.0,
                "delta_prev": d,
                "delta_ci_low": lo,
                "delta_ci_high": hi,
                "signal": signal,
            }
        )
    delta = pl.from_dicts(delta_rows) if delta_rows else pl.DataFrame(
        schema={
            "country": pl.Utf8, "subtype": pl.Utf8,
            "n_pre": pl.Int64, "n_post": pl.Int64,
            "prev_pre": pl.Float64, "prev_post": pl.Float64,
            "delta_prev": pl.Float64,
            "delta_ci_low": pl.Float64, "delta_ci_high": pl.Float64,
            "signal": pl.Boolean,
        }
    )

    n_pre_seqs = int(
        joined.filter(pl.col("epoch") == "pre").height
    )
    n_post_seqs = int(
        joined.filter(pl.col("epoch") == "post").height
    )
    n_signal = int(delta.filter(pl.col("signal")).height) if delta.height else 0

    summary = DTGSummary(
        n_pre_seqs=n_pre_seqs,
        n_post_seqs=n_post_seqs,
        n_strata_scored=delta.height,
        n_strata_with_signal=n_signal,
    )
    logger.info(
        "DTG landscape: pre=%d, post=%d, %d strata scored, %d with signal",
        n_pre_seqs, n_post_seqs, delta.height, n_signal,
    )
    return per_epoch, delta, summary


__all__ = [
    "INSTI_CLASS",
    "MIN_EPOCH_N",
    "POST_YEARS",
    "PRE_YEARS",
    "RESISTANT_LEVELS",
    "ROLLOUT_YEAR",
    "DTGSummary",
    "aggregate_dtg_landscape",
]
