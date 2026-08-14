"""Pure aggregation for Fig 3 — cluster-detection parity.

Reads the JSON emitted by `pipeline/cluster_validate.py --with-hivtrace`
and produces two Polars frames:

  * metrics_frame — one row per (clusterer, metric) with the value.
    Used for the grouped-bar comparison.
  * confusion_frame — one row per (clusterer, {tp, fp, fn, tn}).

Also exposes `parity_agreement(payload)` returning the pair-agreement
fraction between the two clusterers (or None if HIV-TRACE wasn't run).

Kept separate from the notebook so shape and semantics are unit-testable.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# Metrics the Fig 3 bars show, in the order they should render.
FIG3_METRICS: tuple[str, ...] = ("precision", "recall", "f1", "adjusted_rand_index")

# Canonical clusterer names (used for consistent color mapping).
CLUSTERER_ORDER: tuple[str, ...] = (
    "tn93+connected_components",
    "hivtrace-hivnetworkcsv",
)


@dataclass(frozen=True)
class ParitySummary:
    """High-level parity numbers for the notebook's caption text."""

    n_clusterers: int
    n_sequences: int
    n_truth_clusters: int
    pair_agreement: float | None
    ours_f1: float
    hivtrace_f1: float | None


def load_parity_payload(json_path: Path) -> dict:
    """Read the cluster_validate JSON, tolerating the pre-parity format.

    Pre-parity JSONs only have top-level `clusterer`/`metrics`. We
    normalize into the new `results` list shape so downstream aggregators
    can assume a single shape.
    """
    payload = json.loads(json_path.read_text())
    if "results" in payload:
        return payload

    # Backfill from the pre-parity single-clusterer format.
    payload["results"] = [
        {
            "clusterer": payload["clusterer"],
            "threshold": payload.get("threshold"),
            "n_predicted_clusters": payload.get("metrics", {}).get(
                "n_predicted_clusters", 0
            ),
            "metrics": payload["metrics"],
        }
    ]
    payload.setdefault("parity", None)
    return payload


def metrics_frame(payload: dict) -> pl.DataFrame:
    """Long frame with one row per (clusterer, metric, value)."""
    rows: list[dict] = []
    for r in payload.get("results", []):
        clusterer = r["clusterer"]
        m = r.get("metrics", {})
        for metric in FIG3_METRICS:
            v = m.get(metric)
            if v is None:
                continue
            rows.append({"clusterer": clusterer, "metric": metric, "value": float(v)})
    if not rows:
        raise ValueError("parity payload contains no scored metrics")
    return pl.from_dicts(rows)


def confusion_frame(payload: dict) -> pl.DataFrame:
    """One row per (clusterer, tp/fp/fn/tn)."""
    rows: list[dict] = []
    for r in payload.get("results", []):
        m = r.get("metrics", {})
        rows.append(
            {
                "clusterer": r["clusterer"],
                "tp": int(m.get("tp", 0)),
                "fp": int(m.get("fp", 0)),
                "fn": int(m.get("fn", 0)),
                "tn": int(m.get("tn", 0)),
            }
        )
    return pl.from_dicts(rows) if rows else pl.DataFrame()


def summarize(payload: dict) -> ParitySummary:
    """Compact caption-shaped summary."""
    results = payload.get("results", [])
    sim = payload.get("sim", {})
    parity = payload.get("parity")

    ours_f1 = 0.0
    hivtrace_f1: float | None = None
    for r in results:
        m = r.get("metrics", {})
        if r["clusterer"].startswith("tn93"):
            ours_f1 = float(m.get("f1", 0.0))
        elif r["clusterer"].startswith("hivtrace"):
            hivtrace_f1 = float(m.get("f1", 0.0))

    return ParitySummary(
        n_clusterers=len(results),
        n_sequences=int(sim.get("n_sequences", 0)),
        n_truth_clusters=int(sim.get("n_clusters", 0)),
        pair_agreement=float(parity["fraction_agree"]) if parity else None,
        ours_f1=ours_f1,
        hivtrace_f1=hivtrace_f1,
    )


__all__ = [
    "CLUSTERER_ORDER",
    "FIG3_METRICS",
    "ParitySummary",
    "confusion_frame",
    "load_parity_payload",
    "metrics_frame",
    "summarize",
]
