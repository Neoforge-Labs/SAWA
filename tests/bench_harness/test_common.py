"""Unit tests for the benchmark harness common module."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from benchmarks.common import (
    BENCHMARK_SCHEMA_VERSION,
    BenchmarkMeasurement,
    BenchmarkResult,
    HostInfo,
    measure_block,
    parse_sizes,
    today_iso,
    write_benchmark_result,
)


class TestMeasureBlock:
    def test_populates_expected_keys(self) -> None:
        with measure_block(n_sequences=100) as m:
            time.sleep(0.01)
        assert m["n_sequences"] == 100
        assert m["wall_seconds"] >= 0.01
        assert m["peak_mb"] > 0
        assert m["throughput_per_sec"] > 0
        assert "extras" in m

    def test_preserves_user_set_extras(self) -> None:
        with measure_block(n_sequences=50) as m:
            m["extras"] = {"n_pairs": 42.0}
        assert m["extras"] == {"n_pairs": 42.0}

    def test_can_build_measurement_from_dict(self) -> None:
        with measure_block(n_sequences=100) as m:
            pass
        bm = BenchmarkMeasurement(**m)
        assert bm.n_sequences == 100


class TestParseSizes:
    def test_string_input(self) -> None:
        assert parse_sizes(["100", "500", "1000"]) == (100, 500, 1000)

    def test_dedupe_and_sort(self) -> None:
        assert parse_sizes(["1000", "100", "500", "100"]) == (100, 500, 1000)

    def test_rejects_zero_and_negative(self) -> None:
        with pytest.raises(ValueError):
            parse_sizes(["0", "100"])
        with pytest.raises(ValueError):
            parse_sizes(["-1"])

    def test_rejects_empty(self) -> None:
        with pytest.raises(ValueError):
            parse_sizes([])


class TestSerialization:
    def test_roundtrip(self, tmp_path: Path) -> None:
        result = BenchmarkResult(
            stage="align",
            seed=42,
            generated_at=today_iso(),
            host=HostInfo.collect(),
            measurements=(
                BenchmarkMeasurement(
                    n_sequences=100,
                    wall_seconds=0.5,
                    peak_mb=30.0,
                    throughput_per_sec=200.0,
                ),
            ),
            notes="test",
        )
        out = tmp_path / "scale_align.json"
        write_benchmark_result(result, out)
        payload = json.loads(out.read_text())
        assert payload["stage"] == "align"
        assert payload["schema_version"] == BENCHMARK_SCHEMA_VERSION
        assert payload["measurements"][0]["n_sequences"] == 100
        assert payload["measurements"][0]["throughput_per_sec"] == 200.0


class TestHostInfo:
    def test_collects_fields(self) -> None:
        host = HostInfo.collect()
        assert host.platform
        assert host.machine
        assert host.python.count(".") >= 1  # semantic-versionish
