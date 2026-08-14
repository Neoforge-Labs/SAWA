"""Unit tests for the Fig 5 scalability aggregator + power-law fit."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from hivdrwatch.paper1.fig05_scalability import (
    STAGE_ORDER,
    fit_power_law,
    load_scalability_frame,
)


def _write_bench(path: Path, stage: str, measurements: list[dict]) -> Path:
    payload = {
        "stage": stage,
        "seed": 42,
        "generated_at": "2026-08-13",
        "host": {"platform": "test", "machine": "x86", "python": "3.12", "cpu_count": 4},
        "measurements": measurements,
        "schema_version": "0.1.0",
        "notes": "test",
    }
    path.write_text(json.dumps(payload))
    return path


class TestLoadFrame:
    def test_loads_multi_stage(self, tmp_path: Path) -> None:
        _write_bench(
            tmp_path / "scale_align.json",
            "align",
            [
                {"n_sequences": 100, "wall_seconds": 0.5, "peak_mb": 30.0, "throughput_per_sec": 200.0},
                {"n_sequences": 200, "wall_seconds": 1.0, "peak_mb": 40.0, "throughput_per_sec": 200.0},
            ],
        )
        _write_bench(
            tmp_path / "scale_cluster.json",
            "cluster",
            [
                {"n_sequences": 100, "wall_seconds": 2.0, "peak_mb": 50.0, "throughput_per_sec": 50.0},
            ],
        )
        frame, summary = load_scalability_frame(tmp_path)
        assert summary.n_stages == 2
        assert summary.n_points == 3
        assert summary.max_n == 200
        assert set(frame["stage"].unique().to_list()) == {"align", "cluster"}

    def test_empty_dir_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="no benchmark measurements"):
            load_scalability_frame(tmp_path)

    def test_missing_dir_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_scalability_frame(tmp_path / "nope")

    def test_malformed_json_is_skipped(self, tmp_path: Path) -> None:
        (tmp_path / "scale_bogus.json").write_text("not json at all {{{")
        _write_bench(
            tmp_path / "scale_align.json",
            "align",
            [{"n_sequences": 100, "wall_seconds": 0.5, "peak_mb": 30.0, "throughput_per_sec": 200.0}],
        )
        frame, _ = load_scalability_frame(tmp_path)
        assert frame.height == 1

    def test_missing_stage_field_skipped(self, tmp_path: Path) -> None:
        (tmp_path / "scale_bad.json").write_text(json.dumps({"measurements": []}))
        _write_bench(
            tmp_path / "scale_align.json",
            "align",
            [{"n_sequences": 100, "wall_seconds": 0.5, "peak_mb": 30.0, "throughput_per_sec": 200.0}],
        )
        frame, summary = load_scalability_frame(tmp_path)
        assert summary.n_stages == 1


class TestPowerLawFit:
    def test_linear_recovers_exponent_1(self, tmp_path: Path) -> None:
        # wall = 0.01 * N
        _write_bench(
            tmp_path / "scale_align.json",
            "align",
            [
                {"n_sequences": 100, "wall_seconds": 1.0, "peak_mb": 10, "throughput_per_sec": 100},
                {"n_sequences": 200, "wall_seconds": 2.0, "peak_mb": 10, "throughput_per_sec": 100},
                {"n_sequences": 400, "wall_seconds": 4.0, "peak_mb": 10, "throughput_per_sec": 100},
            ],
        )
        frame, _ = load_scalability_frame(tmp_path)
        fit = fit_power_law(frame, "align")
        assert fit is not None
        _a, b = fit
        assert math.isclose(b, 1.0, abs_tol=1e-6)

    def test_quadratic_recovers_exponent_2(self, tmp_path: Path) -> None:
        # wall = N^2
        _write_bench(
            tmp_path / "scale_cluster.json",
            "cluster",
            [
                {"n_sequences": 100, "wall_seconds": 1e4, "peak_mb": 10, "throughput_per_sec": 100},
                {"n_sequences": 200, "wall_seconds": 4e4, "peak_mb": 10, "throughput_per_sec": 100},
                {"n_sequences": 400, "wall_seconds": 16e4, "peak_mb": 10, "throughput_per_sec": 100},
            ],
        )
        frame, _ = load_scalability_frame(tmp_path)
        fit = fit_power_law(frame, "cluster")
        assert fit is not None
        _a, b = fit
        assert math.isclose(b, 2.0, abs_tol=1e-6)

    def test_single_point_returns_none(self, tmp_path: Path) -> None:
        _write_bench(
            tmp_path / "scale_ingest.json",
            "ingest",
            [{"n_sequences": 100, "wall_seconds": 0.5, "peak_mb": 10, "throughput_per_sec": 200}],
        )
        frame, _ = load_scalability_frame(tmp_path)
        assert fit_power_law(frame, "ingest") is None


class TestConstants:
    def test_stage_order_locked(self) -> None:
        assert STAGE_ORDER == ("ingest", "align", "cluster")
