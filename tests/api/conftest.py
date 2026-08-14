"""Shared API test fixtures.

Overrides the `get_db` dependency to hand out an in-memory SQLite
session with the schema created, so tests never touch a live Postgres.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Iterator

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker
from sqlalchemy.pool import StaticPool

from hivdrwatch.api.dependencies import get_db
from hivdrwatch.api.main import create_app
from hivdrwatch.storage.schema import (
    Base,
    ClusterMembership,
    DRCall,
    PipelineRun,
    Sequence,
)


@pytest.fixture
def engine():
    # sqlite:///:memory: gives each connection a fresh DB by default; the
    # API's session dependency opens its own connection separate from the
    # seeder, so tables would appear "not found." StaticPool + a shared
    # cache pins the DB to a single connection reused everywhere.
    engine = create_engine(
        "sqlite:///:memory:",
        future=True,
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(engine)
    return engine


@pytest.fixture
def session_factory(engine):
    return sessionmaker(bind=engine, expire_on_commit=False)


@pytest.fixture
def db_session(session_factory) -> Iterator[Session]:
    sess = session_factory()
    try:
        yield sess
    finally:
        sess.close()


@pytest.fixture
def seeded_session(db_session: Session) -> Session:
    """Populate the DB with a few rows of each table for route tests."""
    run = PipelineRun()
    run.stage = "ingest"
    run.schema_version = "0.1.0"
    run.storage_schema_version = "0.1.0"
    db_session.add(run)
    db_session.flush()

    for i in range(5):
        s = Sequence()
        s.accession = f"AB{i:06d}"
        s.schema_version = "0.1.0"
        s.source = "genbank" if i % 2 == 0 else "lanl"
        s.sequence = "ACGT" * 100
        s.length = 400
        s.country = "KE" if i < 3 else "UG"
        s.subtype = "A1" if i < 3 else "C"
        s.collection_year = 2020 + i
        s.pipeline_run_id = run.id
        s.raw_record_hash = "a" * 64
        db_session.add(s)

    for i in range(5):
        d = DRCall()
        d.accession = f"AB{i:06d}"
        d.drug_class = "NNRTI"
        d.drug = "EFV"
        d.level = "high" if i % 2 == 0 else "susceptible"
        d.algorithm_version = "9.6.1"
        d.pipeline_run_id = run.id
        d.called_at = datetime.now(timezone.utc)
        db_session.add(d)

    for i in range(4):
        c = ClusterMembership()
        c.accession = f"AB{i:06d}"
        c.cluster_id = "C-0001" if i < 2 else "C-0002"
        c.cluster_size = 2
        c.tn93_threshold = 0.015
        c.pipeline_run_id = run.id
        c.computed_at = datetime.now(timezone.utc)
        db_session.add(c)

    db_session.commit()
    return db_session


@pytest.fixture
def app(session_factory) -> FastAPI:
    app = create_app()

    def _override_db() -> Iterator[Session]:
        sess = session_factory()
        try:
            yield sess
            sess.commit()
        except Exception:
            sess.rollback()
            raise
        finally:
            sess.close()

    app.dependency_overrides[get_db] = _override_db
    return app


@pytest.fixture
def client(app: FastAPI, seeded_session: Session) -> TestClient:
    return TestClient(app)
