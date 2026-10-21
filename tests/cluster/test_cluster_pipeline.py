"""Wiring tests for the cluster flow: the two new gates must reach the engine.

A default that never arrives is the defect this change exists to fix,
so the CLI defaults and their forwarding are pinned here. Nothing in
this module runs Prefect, Docker or a real cohort.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from prefect.logging import disable_run_logger

from hivdrwatch.cluster.distance import DEFAULT_MIN_OVERLAP, DEFAULT_THRESHOLD
from hivdrwatch.phylo.build import MIN_ALIGNED_COVERAGE
from hivdrwatch.pipeline import cluster as pipeline_cluster

ALIGNMENT_WIDTH = 3012


def _window(seed: str, n: int) -> str:
    return (seed * ((n // len(seed)) + 1))[:n]


def _padded(prefix: str) -> str:
    return prefix + "-" * (ALIGNMENT_WIDTH - len(prefix))


class TestRunDistances:
    def test_python_engine_honours_min_overlap(self, tmp_path: Path) -> None:
        fasta = tmp_path / "aligned.fasta"
        shared = _window("ACGT", 150)
        fasta.write_text(f">A\n{_padded(shared)}\n>B\n{_padded(shared)}\n")

        with disable_run_logger():
            emitted = pipeline_cluster.run_distances.fn(
                fasta, tmp_path / "work", DEFAULT_THRESHOLD, min_overlap=100
            )
            refused = pipeline_cluster.run_distances.fn(
                fasta, tmp_path / "work2", DEFAULT_THRESHOLD, min_overlap=200
            )

        assert emitted.read_text().splitlines()[1:] == ["A,B,0.000000"]
        assert refused.read_text().splitlines()[1:] == []

    def test_docker_engine_dispatches_to_run_tn93(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fasta = tmp_path / "aligned.fasta"
        fasta.write_text(f">A\n{_padded(_window('ACGT', 2000))}\n")
        seen: dict[str, object] = {}

        def fake_run_tn93(aligned: Path, out: Path, **kwargs: object) -> None:
            seen.update({"aligned": aligned, "out": out, **kwargs})
            out.write_text("ID1,ID2,Distance\n")

        monkeypatch.setattr(pipeline_cluster, "run_tn93", fake_run_tn93)
        with disable_run_logger():
            pipeline_cluster.run_distances.fn(fasta, tmp_path / "work", 0.015, engine="docker")
        assert seen["threshold"] == 0.015


class TestClusterTask:
    def test_assigns_every_curated_accession(self, tmp_path: Path) -> None:
        distances = tmp_path / "distances.csv"
        distances.write_text("ID1,ID2,Distance\nA,B,0.001\n")
        curated = pl.DataFrame({"accession": ["A", "B", "C"]})
        with disable_run_logger():
            out = pipeline_cluster.cluster.fn(distances, curated, 0.015)
        by_acc = {r["accession"]: r["cluster_size"] for r in out.iter_rows(named=True)}
        assert by_acc == {"A": 2, "B": 2, "C": 1}


class TestCli:
    @pytest.fixture
    def captured(self, monkeypatch: pytest.MonkeyPatch) -> dict:
        seen: dict = {}

        def fake_flow(**kwargs: object) -> tuple[Path, Path]:
            seen.update(kwargs)
            return Path("cluster.parquet"), Path("features.parquet")

        monkeypatch.setattr(pipeline_cluster, "cluster_flow", fake_flow)
        return seen

    def test_defaults_carry_both_gates(self, captured: dict) -> None:
        assert pipeline_cluster.main([]) == 0
        assert captured["min_overlap"] == DEFAULT_MIN_OVERLAP
        assert captured["min_coverage"] == MIN_ALIGNED_COVERAGE
        assert captured["threshold"] == DEFAULT_THRESHOLD
        assert captured["distance_engine"] == "python"

    def test_both_gates_are_overridable(self, captured: dict) -> None:
        assert pipeline_cluster.main(["--min-overlap", "500", "--min-coverage", "0.8"]) == 0
        assert captured["min_overlap"] == 500
        assert captured["min_coverage"] == 0.8
