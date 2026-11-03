"""Each flow must surface the storage-tail outcome, not discard it.

`storage_tail` never raises, so before it returned a status the flows
called it and threw the answer away: a database write that broke was
indistinguishable, in the flow's own log, from one nobody requested.
These tests run the real flow bodies with their heavy tasks stubbed and
the real `storage_tail`, steering only `require_postgres`.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

import hivdrwatch.pipeline._storage_tail as tail_mod
from hivdrwatch.pipeline import cluster, dr_ingest, ingest

FlowRunner = Callable[[Path, bool], None]


@pytest.fixture(autouse=True)
def _postgres_down(monkeypatch: pytest.MonkeyPatch) -> None:
    """With loading requested, the database is unreachable: a FAILED tail."""
    monkeypatch.setattr(tail_mod, "require_postgres", lambda: False)


def _flow_logger(monkeypatch: pytest.MonkeyPatch, module: Any) -> None:
    monkeypatch.setattr(module, "get_run_logger", lambda: logging.getLogger("test.flow"))


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x")
    return path


class _FakeRun:
    def __init__(self) -> None:
        self.tags: dict[str, str] = {}

    def set_tags(self, tags: dict[str, str]) -> None:
        self.tags.update(tags)

    def log_params(self, params: dict[str, Any]) -> None: ...

    def log_metrics(self, metrics: dict[str, float]) -> None: ...

    def log_artifact(self, path: Path) -> None: ...


def _run_ingest(
    tmp_path: Path, to_postgres: bool, monkeypatch: pytest.MonkeyPatch
) -> _FakeRun:
    fake = _FakeRun()

    class _Ctx:
        def __enter__(self) -> _FakeRun:
            return fake

        def __exit__(self, *exc: object) -> None:
            return None

    out = _touch(tmp_path / "bronze.parquet")
    monkeypatch.setattr(ingest, "load_lanl", lambda *a, **k: SimpleNamespace(height=0))
    monkeypatch.setattr(ingest, "merge_and_write", lambda *a, **k: out)
    monkeypatch.setattr(ingest, "tracked_run", lambda *a, **k: _Ctx())
    _flow_logger(monkeypatch, ingest)
    ingest.ingest_flow.fn(
        None,
        [],
        lanl_dir=tmp_path,
        genbank_dir=tmp_path,
        out_path=out,
        skip_genbank=True,
        to_postgres=to_postgres,
    )
    return fake


def _run_dr_ingest(tmp_path: Path, to_postgres: bool, monkeypatch: pytest.MonkeyPatch) -> None:
    dr = _touch(tmp_path / "dr.parquet")
    st = _touch(tmp_path / "st.parquet")
    monkeypatch.setattr(dr_ingest, "load_merged_sequences", lambda *a, **k: None)
    monkeypatch.setattr(
        dr_ingest, "call_sierra", lambda *a, **k: ([], [], SimpleNamespace(version="t"))
    )
    monkeypatch.setattr(dr_ingest, "write_dr_parquet", lambda *a, **k: (dr, st))
    _flow_logger(monkeypatch, dr_ingest)
    dr_ingest.dr_ingest_flow.fn(tmp_path, tmp_path, tmp_path, to_postgres=to_postgres)


def _run_cluster(tmp_path: Path, to_postgres: bool, monkeypatch: pytest.MonkeyPatch) -> None:
    gold = _touch(tmp_path / "cluster_assignments.parquet")
    feats = _touch(tmp_path / "outbreak_features.parquet")
    for name in ("load_curated", "materialize_aligned_fasta", "run_distances", "cluster"):
        monkeypatch.setattr(cluster, name, lambda *a, **k: None)
    monkeypatch.setattr(cluster, "outbreak_rollup", lambda *a, **k: None)
    monkeypatch.setattr(cluster, "write_gold", lambda *a, **k: (gold, feats))
    _flow_logger(monkeypatch, cluster)
    cluster.cluster_flow.fn(tmp_path, tmp_path, tmp_path, tmp_path, to_postgres=to_postgres)


def _run_ingest_discarding_run(
    tmp_path: Path, to_postgres: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    _run_ingest(tmp_path, to_postgres, monkeypatch)


FLOWS = {
    "ingest": _run_ingest_discarding_run,
    "dr_ingest": _run_dr_ingest,
    "cluster": _run_cluster,
}


class TestFlowLogs:
    @pytest.mark.parametrize("flow_name", sorted(FLOWS))
    def test_failed_tail_is_an_error_in_the_flow_log(
        self,
        flow_name: str,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        with caplog.at_level(logging.INFO, logger="test.flow"):
            FLOWS[flow_name](tmp_path, True, monkeypatch)

        flow_records = [r for r in caplog.records if r.name == "test.flow"]
        assert [r.levelno for r in flow_records if "storage tail" in r.getMessage()] == [
            logging.ERROR
        ]
        assert "postgres_unavailable" in flow_records[-1].getMessage()

    @pytest.mark.parametrize("flow_name", sorted(FLOWS))
    def test_unrequested_tail_is_not_an_error(
        self,
        flow_name: str,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        with caplog.at_level(logging.INFO, logger="test.flow"):
            FLOWS[flow_name](tmp_path, False, monkeypatch)

        levels = [
            r.levelno
            for r in caplog.records
            if r.name == "test.flow" and "storage tail" in r.getMessage()
        ]
        assert levels == [logging.INFO]


class TestIngestMlflowTag:
    def test_failed_tail_is_tagged_on_the_run(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        run = _run_ingest(tmp_path, True, monkeypatch)
        assert run.tags["storage_tail"] == "failed"

    def test_unrequested_tail_is_tagged_skipped(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        run = _run_ingest(tmp_path, False, monkeypatch)
        assert run.tags["storage_tail"] == "skipped"
