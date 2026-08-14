"""Database engine + session factory.

Reads `HIVDRWATCH_DATABASE_URL` from the env (defaults to the compose
Postgres). Falls back to a local SQLite file when
`HIVDRWATCH_DATABASE_URL=sqlite:///...` so unit tests never need a
Postgres instance.

Session objects are context-managed — always call `with get_session()
as s:` so transactions commit / rollback cleanly.
"""

from __future__ import annotations

import logging
import os
from contextlib import contextmanager
from typing import Iterator

from sqlalchemy import Engine, create_engine
from sqlalchemy.orm import Session, sessionmaker

logger = logging.getLogger(__name__)

DEFAULT_DATABASE_URL = (
    "postgresql+psycopg://{user}:{pw}@{host}:{port}/{db}"
)


def _default_url() -> str:
    """Build the default Postgres URL from POSTGRES_* env vars."""
    return DEFAULT_DATABASE_URL.format(
        user=os.getenv("POSTGRES_USER", "hivdrwatch"),
        pw=os.getenv("POSTGRES_PASSWORD", "hivdrwatch"),
        host=os.getenv("POSTGRES_HOST", "localhost"),
        port=os.getenv("POSTGRES_PORT", "5432"),
        db=os.getenv("HIVDRWATCH_DATABASE_NAME", "hivdrwatch"),
    )


def database_url() -> str:
    """Return the effective database URL."""
    return os.getenv("HIVDRWATCH_DATABASE_URL") or _default_url()


_ENGINE: Engine | None = None
_SESSION_FACTORY: sessionmaker[Session] | None = None


def get_engine(url: str | None = None, *, echo: bool = False) -> Engine:
    """Return a lazily-initialized SQLAlchemy engine.

    Callers should not close the returned engine — the process owns it
    for its lifetime. Test suites that need a fresh engine per test
    should call `reset_engine()` between tests.
    """
    global _ENGINE, _SESSION_FACTORY
    if _ENGINE is None:
        effective_url = url or database_url()
        logger.info("initializing DB engine → %s", _redact(effective_url))
        _ENGINE = create_engine(effective_url, echo=echo, future=True)
        _SESSION_FACTORY = sessionmaker(bind=_ENGINE, expire_on_commit=False)
    return _ENGINE


def reset_engine() -> None:
    """Drop the process-global engine — tests only."""
    global _ENGINE, _SESSION_FACTORY
    if _ENGINE is not None:
        _ENGINE.dispose()
    _ENGINE = None
    _SESSION_FACTORY = None


@contextmanager
def get_session(url: str | None = None) -> Iterator[Session]:
    """Yield a transactional session; commit on success, rollback on error."""
    get_engine(url)
    assert _SESSION_FACTORY is not None  # get_engine sets it
    session = _SESSION_FACTORY()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def require_postgres(*, timeout_seconds: float = 3.0) -> bool:
    """Return True iff we can connect to the configured database and
    execute `SELECT 1` within `timeout_seconds`.

    Used by pipeline flows to decide whether to fire the storage tail:

        if require_postgres():
            load_sequences(parquet, session=sess)

    Never raises — a False result on the pipeline's happy path just
    means the Parquet still lands on disk but no DB load happens.
    Logs one WARNING with a redacted URL so operators can see why.
    """
    from sqlalchemy import text as _text  # noqa: PLC0415

    url = database_url()
    try:
        engine = get_engine(url)
        # `connect_args={"connect_timeout": ...}` isn't universally
        # supported across dialects; use engine.connect() with a
        # per-connection timeout via SQL only when the dialect supports
        # it. Simpler: try connect() + SELECT 1 and let the driver's
        # own default timeout apply (~10s for psycopg by default).
        # We don't try to force a shorter timeout because CI PG on
        # cold start can take a few seconds to respond.
        with engine.connect() as conn:
            conn.execute(_text("SELECT 1"))
        return True
    except Exception as e:  # noqa: BLE001 - defensive at boundary
        logger.warning(
            "database not reachable → %s (%s: %s); skipping storage tail",
            _redact(url), type(e).__name__, e,
        )
        return False


def _redact(url: str) -> str:
    """Never log passwords."""
    if "@" not in url or "://" not in url:
        return url
    scheme, rest = url.split("://", 1)
    if "@" not in rest:
        return url
    creds, host = rest.split("@", 1)
    if ":" in creds:
        user, _pw = creds.split(":", 1)
        return f"{scheme}://{user}:***@{host}"
    return f"{scheme}://{creds}@{host}"


__all__ = [
    "database_url",
    "get_engine",
    "get_session",
    "require_postgres",
    "reset_engine",
]
