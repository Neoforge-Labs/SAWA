"""FastAPI dependency-injection wiring.

Session dependency yields a SQLAlchemy `Session` bound to the same
engine as the rest of the application. Tests override it via
`app.dependency_overrides[get_db]` to point at an in-memory SQLite.
"""

from __future__ import annotations

from collections.abc import Iterator

from sqlalchemy.orm import Session

from ..storage.session import get_session


def get_db() -> Iterator[Session]:
    """FastAPI dependency — yields a session, commits on success."""
    with get_session() as sess:
        yield sess


__all__ = ["get_db"]
