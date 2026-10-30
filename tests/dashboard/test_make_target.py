"""`make dashboard` must not supply an ingested count the operator did not.

The export refuses to run without `--ingested-count` because the parquet cannot
tell it how many accessions entered ingest. A default in the Makefile would put that
default back one layer up and write a stale number into a published provenance claim.
These tests run the real target with a stand-in `python`, so what they observe is the
command line `make` would have launched.
"""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.skipif(shutil.which("make") is None, reason="make is not installed")


def _make_dashboard(
    tmp_path: Path, *args: str, env_count: str | None = None
) -> subprocess.CompletedProcess[str]:
    shim_dir = tmp_path / "shim"
    shim_dir.mkdir()
    shim = shim_dir / "python"
    shim.write_text('#!/bin/sh\necho "PYTHON-RAN $*"\n')
    shim.chmod(shim.stat().st_mode | stat.S_IXUSR)
    env = {k: v for k, v in os.environ.items() if k != "INGESTED_COUNT"}
    env["PATH"] = f"{shim_dir}{os.pathsep}{env['PATH']}"
    if env_count is not None:
        env["INGESTED_COUNT"] = env_count
    return subprocess.run(
        ["make", "-C", str(REPO), "dashboard", *args],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def test_make_dashboard_without_a_count_stops_before_the_export_runs(tmp_path: Path) -> None:
    done = _make_dashboard(tmp_path)
    assert done.returncode != 0
    assert "INGESTED_COUNT is required" in done.stderr
    assert "PYTHON-RAN" not in done.stdout


def test_make_dashboard_treats_an_empty_count_as_missing(tmp_path: Path) -> None:
    done = _make_dashboard(tmp_path, env_count="")
    assert done.returncode != 0
    assert "PYTHON-RAN" not in done.stdout


@pytest.mark.parametrize("how", ["argument", "environment"])
def test_make_dashboard_passes_the_operators_count_through_unchanged(
    tmp_path: Path, how: str
) -> None:
    done = (
        _make_dashboard(tmp_path, "INGESTED_COUNT=41234")
        if how == "argument"
        else _make_dashboard(tmp_path, env_count="41234")
    )
    assert done.returncode == 0, done.stderr
    assert "--ingested-count 41234" in done.stdout


def test_make_dashboard_points_the_export_at_the_curated_table_under_the_interim_dir(
    tmp_path: Path,
) -> None:
    done = _make_dashboard(tmp_path, "INGESTED_COUNT=41234", "INTERIM_DIR=elsewhere/interim")
    assert done.returncode == 0, done.stderr
    assert "--input elsewhere/interim/dr_results.parquet" in done.stdout
    assert "--curated elsewhere/interim/curated/aligned.parquet" in done.stdout
