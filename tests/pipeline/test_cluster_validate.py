"""The results JSON must say whether HIV-TRACE parity was requested and what became of it.

A requested-but-skipped run used to serialise identically to a deliberate
tn93-only run, so a Figure 3 parity claim could rest on a degraded run
with nothing on the artifact to show it. These tests drive the real
flow body end to end on a tiny simulation; only the HIV-TRACE container
call is replaced. The tasks run as plain functions so no Prefect server
is started.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from prefect.logging import disable_run_logger

from hivdrwatch.pipeline import cluster_validate
from hivdrwatch.pipeline.cluster_validate import (
    HIVTRACE_SKIP_IMAGE_UNAVAILABLE,
    cluster_validate_flow,
)

_TINY = {
    "n_clusters": 2,
    "cluster_size_lo": 3,
    "cluster_size_hi": 4,
    "n_singletons": 2,
    "seq_length": 300,
    "seed": 7,
}


@pytest.fixture(autouse=True)
def _run_tasks_in_process(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in (
        "simulate",
        "compute_distances",
        "cluster_ours",
        "cluster_hivtrace",
        "write_results",
    ):
        monkeypatch.setattr(cluster_validate, name, getattr(cluster_validate, name).fn)


def _run(tmp_path: Path, *, with_hivtrace: bool) -> dict[str, object]:
    out = tmp_path / "results.json"
    with disable_run_logger():
        cluster_validate_flow.fn(
            out_path=out, work_dir=tmp_path / "work", with_hivtrace=with_hivtrace, **_TINY
        )
    payload: dict[str, object] = json.loads(out.read_text())
    return payload


def _hivtrace_unavailable(*_args: object, **_kwargs: object) -> None:
    raise cluster_validate.HIVTraceNotAvailable("image not present")


def _fake_hivtrace(
    _csv: Path, _net: Path, accessions: list[str], **_kwargs: object
) -> tuple[dict[str, str], None]:
    # Every sequence its own cluster: enough to produce a scored second result.
    return {a: f"SINGLETON_{a}" for a in accessions}, None


class TestHIVTraceRecord:
    def test_requested_but_unavailable_is_recorded_in_the_artifact(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(cluster_validate, "run_hivtrace", _hivtrace_unavailable)

        payload = _run(tmp_path, with_hivtrace=True)

        assert payload["with_hivtrace"] is True
        assert payload["hivtrace_skipped"] == HIVTRACE_SKIP_IMAGE_UNAVAILABLE
        assert payload["parity"] is None
        assert len(payload["results"]) == 1  # type: ignore[arg-type]

    def test_deliberate_tn93_only_run_is_not_marked_skipped(self, tmp_path: Path) -> None:
        payload = _run(tmp_path, with_hivtrace=False)

        assert payload["with_hivtrace"] is False
        assert payload["hivtrace_skipped"] is None
        assert payload["parity"] is None

    def test_skipped_and_deliberate_runs_differ_on_the_artifact(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        deliberate = _run(tmp_path / "a", with_hivtrace=False)
        monkeypatch.setattr(cluster_validate, "run_hivtrace", _hivtrace_unavailable)
        skipped = _run(tmp_path / "b", with_hivtrace=True)

        assert deliberate["parity"] == skipped["parity"] is None
        assert deliberate["results"] == skipped["results"]
        assert deliberate != skipped

    def test_successful_parity_run_records_request_and_no_skip(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(cluster_validate, "run_hivtrace", _fake_hivtrace)

        payload = _run(tmp_path, with_hivtrace=True)

        assert payload["with_hivtrace"] is True
        assert payload["hivtrace_skipped"] is None
        assert payload["parity"] is not None
        assert len(payload["results"]) == 2  # type: ignore[arg-type]
