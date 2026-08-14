"""Unit tests for storage.session.require_postgres()."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest


def test_returns_true_on_successful_ping(monkeypatch: pytest.MonkeyPatch) -> None:
    import hivdrwatch.storage.session as sess_mod

    fake_conn = MagicMock()
    ctx = MagicMock()
    ctx.__enter__ = MagicMock(return_value=fake_conn)
    ctx.__exit__ = MagicMock(return_value=False)

    fake_engine = MagicMock()
    fake_engine.connect = MagicMock(return_value=ctx)

    monkeypatch.setattr(sess_mod, "get_engine", lambda url=None: fake_engine)

    assert sess_mod.require_postgres() is True
    fake_conn.execute.assert_called_once()


def test_returns_false_when_connect_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    import hivdrwatch.storage.session as sess_mod

    fake_engine = MagicMock()
    fake_engine.connect.side_effect = ConnectionError("no route to host")

    monkeypatch.setattr(sess_mod, "get_engine", lambda url=None: fake_engine)

    assert sess_mod.require_postgres() is False


def test_returns_false_when_ping_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    import hivdrwatch.storage.session as sess_mod

    fake_conn = MagicMock()
    fake_conn.execute.side_effect = RuntimeError("connection lost after handshake")
    ctx = MagicMock()
    ctx.__enter__ = MagicMock(return_value=fake_conn)
    ctx.__exit__ = MagicMock(return_value=False)

    fake_engine = MagicMock()
    fake_engine.connect = MagicMock(return_value=ctx)

    monkeypatch.setattr(sess_mod, "get_engine", lambda url=None: fake_engine)

    assert sess_mod.require_postgres() is False


def test_never_raises_even_on_completely_broken_engine(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import hivdrwatch.storage.session as sess_mod

    def _broken(url=None):
        raise ValueError("bad URL")

    monkeypatch.setattr(sess_mod, "get_engine", _broken)
    # Must not raise — pipeline flows depend on graceful degradation.
    assert sess_mod.require_postgres() is False
