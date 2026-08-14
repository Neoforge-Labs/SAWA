"""Pure aggregation for Fig 5 — pipeline scalability (wall time vs N).

Reads the JSON results emitted by `benchmarks/scale_*.py` and produces
a long-format Polars frame the notebook plots on a log-log scale.

Kept separate from the notebook so it's unit-testable without matplotlib.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# Canonical ordering for the legend + stage colors on the figure.
STAGE_ORDER: tuple[str, ...] = ("ingest", "align", "cluster")


@dataclass(frozen=True)
class ScalabilitySummary:
    """High-level stats for the summary log below the plot."""

    n_stages: int
    n_points: int
    max_n: int
    min_wall_seconds: float
    max_wall_seconds: float


def load_scalability_frame(results_dir: Path) -> tuple[pl.DataFrame, ScalabilitySummary]:
    """Load every `scale_*.json` under `results_dir` into a long frame.

    Returned columns:
        stage, n_sequences, wall_seconds, peak_mb, throughput_per_sec,
        generated_at, seed, host_platform

    Missing files are silently skipped; benchmarks the user hasn't
    run yet just don't appear in the plot. The notebook shows a
    "run `make benchmark-<stage>` to populate this" hint if a stage
    is missing.
    """
    if not results_dir.exists():
        raise FileNotFoundError(f"results dir not found: {results_dir}")

    rows: list[dict] = []
    files = sorted(results_dir.glob("scale_*.json"))
    for path in files:
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError as e:
            logger.warning("skipping malformed benchmark JSON %s: %s", path, e)
            continue

        stage = payload.get("stage")
        if not stage:
            logger.warning("skipping benchmark JSON missing `stage`: %s", path)
            continue

        for m in payload.get("measurements", []):
            rows.append(
                {
                    "stage": stage,
                    "n_sequences": int(m["n_sequences"]),
                    "wall_seconds": float(m["wall_seconds"]),
                    "peak_mb": float(m["peak_mb"]),
                    "throughput_per_sec": float(m["throughput_per_sec"]),
                    "generated_at": payload.get("generated_at"),
                    "seed": int(payload.get("seed", 0)),
                    "host_platform": payload.get("host", {}).get("platform", "unknown"),
                }
            )

    if not rows:
        raise ValueError(
            f"no benchmark measurements found under {results_dir}. "
            "run `make benchmarks` to populate."
        )

    frame = pl.from_dicts(rows).sort(["stage", "n_sequences"])

    summary = ScalabilitySummary(
        n_stages=frame["stage"].n_unique(),
        n_points=frame.height,
        max_n=int(frame["n_sequences"].max()),
        min_wall_seconds=float(frame["wall_seconds"].min()),
        max_wall_seconds=float(frame["wall_seconds"].max()),
    )
    logger.info(
        "loaded scalability: %d stages, %d points, N up to %d",
        summary.n_stages,
        summary.n_points,
        summary.max_n,
    )
    return frame, summary


def fit_power_law(frame: pl.DataFrame, stage: str) -> tuple[float, float] | None:
    """Fit wall = a * N^b to a stage's measurements; return (a, b).

    Uses ordinary least squares in log-log space. Returns None if fewer
    than 2 measurements are available for the stage.
    """
    sub = frame.filter(pl.col("stage") == stage).sort("n_sequences")
    if sub.height < 2:
        return None

    import math

    ns = [math.log(n) for n in sub["n_sequences"].to_list()]
    ws = [math.log(w) for w in sub["wall_seconds"].to_list()]
    n = len(ns)
    mean_x = sum(ns) / n
    mean_y = sum(ws) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(ns, ws, strict=True))
    den = sum((x - mean_x) ** 2 for x in ns)
    if den == 0:
        return None
    b = num / den
    log_a = mean_y - b * mean_x
    return math.exp(log_a), b


__all__ = [
    "STAGE_ORDER",
    "ScalabilitySummary",
    "fit_power_law",
    "load_scalability_frame",
]
