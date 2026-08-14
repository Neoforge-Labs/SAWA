"""Shared plumbing for benchmark scripts.

Every scale_* script generates a `BenchmarkResult` and serializes it to
JSON with `write_benchmark_result`. The schema is stable — the Fig 5
notebook consumes it, and CI diffs successive results to detect
regressions.

Measurement uses `resource.getrusage` for peak RSS (portable across
Linux + macOS) and `time.perf_counter` for wall time. We deliberately
don't reach for pyperf / airspeed-velocity — a JSON file with a wall
time and a peak MB is exactly what a reviewer wants to see.
"""

from __future__ import annotations

import json
import logging
import platform
import resource
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

logger = logging.getLogger(__name__)

BENCHMARK_SCHEMA_VERSION = "0.1.0"


@dataclass(frozen=True)
class HostInfo:
    """Bare-minimum host info for reproducibility annotation.

    Deliberately no GPU / cloud metadata — the scale story for Paper 1
    is a single-machine story ("runs on a laptop"), and cross-cloud
    numbers would obscure that.
    """

    platform: str
    machine: str
    python: str
    cpu_count: int | None

    @classmethod
    def collect(cls) -> HostInfo:
        return cls(
            platform=platform.platform(),
            machine=platform.machine(),
            python=sys.version.split()[0],
            cpu_count=_cpu_count(),
        )


def _cpu_count() -> int | None:
    try:
        import os

        return os.cpu_count()
    except Exception:  # pragma: no cover - trivially defensive
        return None


@dataclass(frozen=True)
class BenchmarkMeasurement:
    """One (n, wall, mem) triple for a single input size."""

    n_sequences: int
    wall_seconds: float
    peak_mb: float
    throughput_per_sec: float
    extras: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class BenchmarkResult:
    """A full benchmark run — one stage, multiple input sizes."""

    stage: str
    seed: int
    generated_at: str  # ISO date
    host: HostInfo
    measurements: tuple[BenchmarkMeasurement, ...]
    schema_version: str = BENCHMARK_SCHEMA_VERSION
    notes: str = ""

    def to_dict(self) -> dict:
        return {
            "stage": self.stage,
            "seed": self.seed,
            "generated_at": self.generated_at,
            "host": asdict(self.host),
            "measurements": [asdict(m) for m in self.measurements],
            "schema_version": self.schema_version,
            "notes": self.notes,
        }


@contextmanager
def measure_block(n_sequences: int) -> Iterator[dict]:
    """Time and RSS-measure a code block.

    Yields a mutable dict; on exit populates `wall_seconds`, `peak_mb`,
    `throughput_per_sec`, `n_sequences`. Any keys set inside the block
    (e.g. `extras`) are preserved and merged into the final measurement.

    Usage::

        with measure_block(n_sequences=1000) as m:
            do_work()
            m["extras"] = {"n_pairs": 42}

        # m is now a fully populated dict shaped like BenchmarkMeasurement
    """
    start_wall = time.perf_counter()
    start_rss = _peak_rss_mb()
    payload: dict = {}
    try:
        yield payload
    finally:
        wall = time.perf_counter() - start_wall
        peak = max(_peak_rss_mb(), start_rss)
        payload["n_sequences"] = n_sequences
        payload["wall_seconds"] = float(wall)
        payload["peak_mb"] = float(peak)
        payload["throughput_per_sec"] = float(n_sequences) / wall if wall > 0 else float("inf")
        payload.setdefault("extras", {})


def _peak_rss_mb() -> float:
    """Return peak RSS in MB, portable across Linux + macOS.

    Linux getrusage returns kB; macOS returns bytes. We normalize.
    """
    ru = resource.getrusage(resource.RUSAGE_SELF)
    raw = ru.ru_maxrss
    if sys.platform == "darwin":
        return raw / (1024 * 1024)
    return raw / 1024


def write_benchmark_result(result: BenchmarkResult, out_path: Path) -> Path:
    """Serialize a `BenchmarkResult` to JSON with 2-space indent."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=False) + "\n")
    logger.info("wrote benchmark result → %s", out_path)
    return out_path


def today_iso() -> str:
    """UTC ISO date for the `generated_at` field.

    Kept as a function so tests can monkeypatch.
    """
    return datetime.now(UTC).date().isoformat()


def parse_sizes(raw: list[str] | list[int]) -> tuple[int, ...]:
    """Parse `--sizes 100 200 500` into a sorted unique int tuple."""
    parsed = sorted({int(x) for x in raw})
    if not parsed or any(x <= 0 for x in parsed):
        raise ValueError(f"sizes must be positive; got {raw}")
    return tuple(parsed)


__all__ = [
    "BENCHMARK_SCHEMA_VERSION",
    "BenchmarkMeasurement",
    "BenchmarkResult",
    "HostInfo",
    "measure_block",
    "parse_sizes",
    "today_iso",
    "write_benchmark_result",
]
