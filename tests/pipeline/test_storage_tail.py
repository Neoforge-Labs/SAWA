"""Unit tests for the storage tail — Parquet → Postgres load.

Covers:
  * enable/disable via CLI flag and env var
  * skip path when Parquet is missing
  * skip path when DB is unreachable (mocked)
  * happy path calls the correct loader with the right args
  * unknown stage doesn't crash
  * every stage in _LOADERS is one of the declared Stage values

Never touches a real database — session + require_postgres + loaders
are all monkeypatched. That's the whole point of the helper: it's the
thin dispatch layer between flows and storage.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

# --- Helpers ---


def _make_parquet(tmp_path: Path, name: str) -> Path:
    p = tmp_path / name
    p.write_bytes(b"parquet-magic-placeholder")  # loader is mocked; content doesn't matter
    return p


@pytest.fixture
def patch_storage_tail(monkeypatch: pytest.MonkeyPatch):
    """Yield a namespace holding the mock objects the test can inspect.

    Replaces `require_postgres`, `get_session`, `ensure_schema`, and
    the three loaders inside `pipeline._storage_tail`. Tests set
    return values on the returned mocks.
    """
    import hivdrwatch.pipeline._storage_tail as tail

    class _Ns:
        require_postgres = MagicMock(return_value=True)
        ensure_schema = MagicMock()
        # session context manager returning a MagicMock session
        session = MagicMock()

        # loader mocks — each returns a LoadCounts-shaped MagicMock
        load_sequences = MagicMock()
        load_dr_calls = MagicMock()
        load_cluster_memberships = MagicMock()

    ns = _Ns()

    # Configure get_session as a context manager
    ctx = MagicMock()
    ctx.__enter__ = MagicMock(return_value=ns.session)
    ctx.__exit__ = MagicMock(return_value=False)
    get_session_mock = MagicMock(return_value=ctx)

    # Configure loader return values
    for loader_mock in (ns.load_sequences, ns.load_dr_calls, ns.load_cluster_memberships):
        loader_mock.return_value = MagicMock(inserted=42, skipped=3, pipeline_run_id=7)

    monkeypatch.setattr(tail, "require_postgres", ns.require_postgres)
    monkeypatch.setattr(tail, "get_session", get_session_mock)
    monkeypatch.setattr(tail, "ensure_schema", ns.ensure_schema)
    monkeypatch.setattr(
        tail,
        "_LOADERS",
        {
            "ingest": ns.load_sequences,
            "sierra": ns.load_dr_calls,
            "cluster": ns.load_cluster_memberships,
        },
    )
    # Reset env var so tests control it explicitly.
    monkeypatch.delenv("HIVDRWATCH_LOAD_TO_POSTGRES", raising=False)
    return ns


# --- Enable/disable logic ---


class TestEnableGate:
    def test_disabled_by_default(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "x.parquet")
        result = storage_tail("ingest", p)  # no flag, no env
        assert result is None
        assert patch_storage_tail.require_postgres.call_count == 0
        assert patch_storage_tail.load_sequences.call_count == 0

    def test_enabled_via_arg(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "x.parquet")
        result = storage_tail("ingest", p, enabled=True)
        assert result is not None
        assert patch_storage_tail.load_sequences.call_count == 1

    def test_enabled_via_env_var(
        self, patch_storage_tail, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from hivdrwatch.pipeline._storage_tail import ENV_ENABLE, storage_tail

        monkeypatch.setenv(ENV_ENABLE, "1")
        p = _make_parquet(tmp_path, "x.parquet")
        result = storage_tail("ingest", p)  # no explicit enabled arg
        assert result is not None
        assert patch_storage_tail.load_sequences.call_count == 1

    def test_arg_overrides_env(
        self, patch_storage_tail, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from hivdrwatch.pipeline._storage_tail import ENV_ENABLE, storage_tail

        monkeypatch.setenv(ENV_ENABLE, "1")
        p = _make_parquet(tmp_path, "x.parquet")
        # enabled=False should still enable because we only fall through
        # to env when the arg is *None* (the "unspecified" sentinel);
        # False is a explicit disable.
        result = storage_tail("ingest", p, enabled=False)
        assert result is None
        assert patch_storage_tail.load_sequences.call_count == 0


class TestSkipPaths:
    def test_missing_parquet_skips_cleanly(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        result = storage_tail("ingest", tmp_path / "does_not_exist.parquet", enabled=True)
        assert result is None
        assert patch_storage_tail.require_postgres.call_count == 0

    def test_db_unreachable_skips_cleanly(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        patch_storage_tail.require_postgres.return_value = False
        p = _make_parquet(tmp_path, "x.parquet")
        result = storage_tail("ingest", p, enabled=True)
        assert result is None
        assert patch_storage_tail.load_sequences.call_count == 0

    def test_loader_exception_never_crashes_flow(
        self, patch_storage_tail, tmp_path: Path
    ) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        patch_storage_tail.load_sequences.side_effect = RuntimeError("db went away mid-load")
        p = _make_parquet(tmp_path, "x.parquet")
        # Must not raise.
        result = storage_tail("ingest", p, enabled=True)
        assert result is None

    def test_unknown_stage_skips_cleanly(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "x.parquet")
        # Type-check would catch this, but the helper should still
        # degrade gracefully at runtime.
        result = storage_tail("nope", p, enabled=True)  # type: ignore[arg-type]
        assert result is None


class TestLoaderDispatch:
    def test_ingest_calls_load_sequences(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "seqs.parquet")
        storage_tail("ingest", p, enabled=True, git_sha="abc123")

        patch_storage_tail.load_sequences.assert_called_once()
        (args, kwargs) = (
            patch_storage_tail.load_sequences.call_args.args,
            patch_storage_tail.load_sequences.call_args.kwargs,
        )
        assert args[0] == p
        assert kwargs["session"] is patch_storage_tail.session
        assert kwargs["git_sha"] == "abc123"

    def test_sierra_calls_load_dr_calls(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "dr.parquet")
        storage_tail("sierra", p, enabled=True)
        patch_storage_tail.load_dr_calls.assert_called_once()

    def test_cluster_calls_load_cluster_memberships(
        self, patch_storage_tail, tmp_path: Path
    ) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "cm.parquet")
        storage_tail("cluster", p, enabled=True)
        patch_storage_tail.load_cluster_memberships.assert_called_once()


class TestReturnedCounts:
    def test_forwards_loader_result_shape(self, patch_storage_tail, tmp_path: Path) -> None:
        from hivdrwatch.pipeline._storage_tail import storage_tail

        p = _make_parquet(tmp_path, "x.parquet")
        counts = storage_tail("ingest", p, enabled=True)
        assert counts is not None
        assert counts.inserted == 42
        assert counts.skipped == 3
        assert counts.pipeline_run_id == 7


class TestContract:
    def test_every_stage_has_a_loader(self) -> None:
        from hivdrwatch.pipeline._storage_tail import _LOADERS

        assert set(_LOADERS) == {"ingest", "sierra", "cluster"}
