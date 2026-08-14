"""Automated feature pruning from causal effect estimates.

Consumes the parquet emitted by `flows/causal.py` and emits a JSON
config the trainers can honor:

    {
      "dag_version": "0.1.0",
      "generated_at": "2026-08-12",
      "prune_reasons": {
        "PR_82_A": "CI includes 0 in DR/drug=3TC",
        "n_active_clusters": "fails 2/3 refutations in outbreak"
      },
      "keep": ["RT_184_V", "RT_103_N", ...],
      "drop": ["PR_82_A", ...]
    }

The pruning rules are conservative on purpose. A feature is dropped
only if BOTH:
  1. It's identified in the DAG (unidentified rows are ignored, not
     pruned — we don't know their effect either way).
  2. It has either (a) a CI that spans 0 with the point estimate small,
     or (b) failed at least 2 of 3 refutations.

Everything else is kept. Trainer code sees the config as a hint, not a
mandate — the model card records which features were dropped and why.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# A CI is "consistent with zero" if it spans 0 AND the point estimate
# is smaller than SMALL_EFFECT_THRESHOLD * (ci_high - ci_low). That is,
# the estimate is <<half the CI width — a noisy zero, not a real effect.
SMALL_EFFECT_THRESHOLD = 0.5

# Minimum failed refutations to trigger pruning on refutation grounds.
MIN_FAILED_REFUTATIONS = 2


@dataclass(frozen=True)
class PruneConfig:
    dag_version: str
    generated_at: date
    keep: tuple[str, ...]
    drop: tuple[str, ...]
    prune_reasons: dict[str, str]

    def to_json(self) -> str:
        payload = {
            "dag_version": self.dag_version,
            "generated_at": self.generated_at.isoformat(),
            "keep": list(self.keep),
            "drop": list(self.drop),
            "prune_reasons": dict(self.prune_reasons),
        }
        return json.dumps(payload, indent=2, sort_keys=True)

    def write(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json(), encoding="utf-8")
        return path


def load_prune_config(path: Path) -> PruneConfig | None:
    """Load a prune config emitted by `derive_prune_config`.

    Returns None (not an error) if the file doesn't exist — trainers
    should treat "no config" as "keep everything."
    """
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    return PruneConfig(
        dag_version=payload["dag_version"],
        generated_at=date.fromisoformat(payload["generated_at"]),
        keep=tuple(payload.get("keep", [])),
        drop=tuple(payload.get("drop", [])),
        prune_reasons=dict(payload.get("prune_reasons", {})),
    )


def derive_prune_config(
    effects_df: pl.DataFrame,
    *,
    dag_version: str = "0.1.0",
    scope_prefix: str | None = None,
) -> PruneConfig:
    """Turn a `causal_effects.parquet` frame into a prune config.

    `scope_prefix` restricts to rows whose `data_scope` starts with the
    prefix (e.g. `"DR"` to prune only DR features, ignoring outbreak).
    """
    if effects_df.height == 0:
        return PruneConfig(
            dag_version=dag_version,
            generated_at=date.today(),
            keep=(),
            drop=(),
            prune_reasons={},
        )

    df = effects_df
    if scope_prefix:
        df = df.filter(pl.col("data_scope").str.starts_with(scope_prefix))
    df = df.filter(pl.col("identified") == True)  # noqa: E712

    keep: set[str] = set()
    drop: dict[str, str] = {}

    for row in df.iter_rows(named=True):
        treatment = row["treatment"]
        # Rule 1: CI spans zero AND point estimate is small vs interval width.
        ci_low = row.get("ci_low")
        ci_high = row.get("ci_high")
        estimate = row.get("estimate") or 0.0
        ci_reason = None
        if ci_low is not None and ci_high is not None:
            spans_zero = ci_low <= 0 <= ci_high
            width = ci_high - ci_low
            if spans_zero and width > 0 and abs(estimate) < SMALL_EFFECT_THRESHOLD * width:
                ci_reason = (
                    f"CI [{ci_low:.3f}, {ci_high:.3f}] spans 0 with small "
                    f"effect {estimate:.3f} in {row['data_scope']}"
                )

        # Rule 2: fails ≥ MIN_FAILED_REFUTATIONS of 3 refutations.
        n_failed = int(row.get("ref_fail", 0) or 0)
        ref_reason = None
        if n_failed >= MIN_FAILED_REFUTATIONS:
            ref_reason = (
                f"fails {n_failed}/3 refutations in {row['data_scope']}"
            )

        if ci_reason or ref_reason:
            # First-observed reason wins to keep the config readable.
            if treatment not in drop:
                drop[treatment] = ci_reason or ref_reason  # type: ignore[assignment]
        else:
            keep.add(treatment)

    # A feature dropped in one scope but kept in another → keep it (a
    # per-scope pruning table can come later; treatment-level v0.1 stays simple).
    for feat in list(drop):
        if feat in keep:
            drop.pop(feat)

    logger.info(
        "prune config: %d keep, %d drop (scope_prefix=%r)",
        len(keep),
        len(drop),
        scope_prefix,
    )
    return PruneConfig(
        dag_version=dag_version,
        generated_at=date.today(),
        keep=tuple(sorted(keep)),
        drop=tuple(sorted(drop)),
        prune_reasons=dict(drop),
    )


def apply_prune_config(
    feature_columns: list[str], config: PruneConfig | None
) -> tuple[list[str], list[str]]:
    """Filter feature columns against a prune config.

    Returns `(kept_cols, dropped_cols)`. Columns not mentioned in the
    config are kept — this is a subtractive filter, not an allowlist.
    """
    if config is None or not config.drop:
        return (list(feature_columns), [])
    drop_set = set(config.drop)
    kept: list[str] = []
    dropped: list[str] = []
    for c in feature_columns:
        # Match either exact feature name or any column whose prefix is
        # the treatment (e.g. "RT_184" pruning kills RT_184_V, RT_184_I,
        # RT_184_covered together).
        if c in drop_set or any(c.startswith(f"{d}_") or c == d for d in drop_set):
            dropped.append(c)
        else:
            kept.append(c)
    return (kept, dropped)
