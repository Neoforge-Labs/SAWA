"""`_insert_many`'s generic fallback must skip duplicates and nothing else.

PostgreSQL and SQLite take dedicated `ON CONFLICT` branches, so the
fallback is unreachable on the dialects this project runs. To exercise it
anyway these tests run a real SQLite engine whose dialect reports an
unrecognised name: the statements and the integrity errors are genuine
driver output, only the dialect dispatch is steered. The fallback has no
production caller to test through.
"""

from __future__ import annotations

from collections.abc import Iterator
from datetime import UTC, datetime
from types import SimpleNamespace
from typing import Any

import pytest
from sqlalchemy import create_engine, func, select, text
from sqlalchemy.exc import IntegrityError, OperationalError
from sqlalchemy.orm import Session

from hivdrwatch.storage.loader import _insert_many
from hivdrwatch.storage.schema import Base, DRCall


def _call(accession: str, **overrides: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "accession": accession,
        "drug_class": "NNRTI",
        "drug": "EFV",
        "level": "high",
        "mutations": None,
        "algorithm_version": "9.6.1",
        "pipeline_run_id": 1,
        "called_at": datetime(2026, 1, 1, tzinfo=UTC),
    }
    row.update(overrides)
    return row


def _count(session: Session) -> int:
    return session.execute(select(func.count()).select_from(DRCall)).scalar_one()


@pytest.fixture
def session(monkeypatch: pytest.MonkeyPatch) -> Iterator[Session]:
    engine = create_engine("sqlite:///:memory:", future=True)
    Base.metadata.create_all(engine)
    # Anything other than "postgresql"/"sqlite" selects the fallback.
    monkeypatch.setattr(engine.dialect, "name", "fallback-under-test")
    with Session(engine, expire_on_commit=False) as sess:
        yield sess


class TestFallbackPath:
    def test_is_the_path_under_test(self, session: Session) -> None:
        assert session.get_bind().dialect.name == "fallback-under-test"

    def test_skips_a_true_duplicate_and_keeps_going(self, session: Session) -> None:
        inserted, skipped = _insert_many(session, DRCall, [_call("A1")])
        assert (inserted, skipped) == (1, 0)

        inserted, skipped = _insert_many(session, DRCall, [_call("A1"), _call("A2")])

        assert (inserted, skipped) == (1, 1)
        assert _count(session) == 2

    def test_non_unique_integrity_error_propagates(self, session: Session) -> None:
        # NOT NULL on `level` is an IntegrityError, but not "row already present".
        with pytest.raises(IntegrityError):
            _insert_many(session, DRCall, [_call("A1", level=None)])
        assert _count(session) == 0

    def test_non_integrity_error_propagates(self, session: Session) -> None:
        session.execute(text("DROP TABLE dr_calls"))

        with pytest.raises(OperationalError):
            _insert_many(session, DRCall, [_call("A1")])


class TestDuplicateClassification:
    @pytest.mark.parametrize(
        ("orig", "is_duplicate"),
        [
            (SimpleNamespace(sqlstate="23505"), True),
            (SimpleNamespace(pgcode="23505"), True),
            (SimpleNamespace(sqlstate="23502"), False),  # not_null_violation
            (SimpleNamespace(pgcode="23503"), False),  # foreign_key_violation
            (SimpleNamespace(sqlite_errorname="SQLITE_CONSTRAINT_UNIQUE"), True),
            (SimpleNamespace(sqlite_errorname="SQLITE_CONSTRAINT_PRIMARYKEY"), True),
            (SimpleNamespace(sqlite_errorname="SQLITE_CONSTRAINT_NOTNULL"), False),
            (SimpleNamespace(args=(1062, "Duplicate entry")), True),
            (SimpleNamespace(args=(1048, "Column cannot be null")), False),
            (SimpleNamespace(), True),  # driver exposes no code
        ],
    )
    def test_driver_codes(self, orig: object, is_duplicate: bool) -> None:
        from hivdrwatch.storage.loader import _is_duplicate_row

        exc = IntegrityError("INSERT", {}, orig)  # type: ignore[arg-type]
        assert _is_duplicate_row(exc) is is_duplicate
