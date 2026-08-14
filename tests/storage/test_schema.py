"""Schema smoke tests — run against SQLite so tests never need Postgres."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from hivdrwatch.storage.schema import (
    STORAGE_SCHEMA_VERSION,
    Base,
    ClusterMembership,
    DRCall,
    PipelineRun,
    Sequence,
)


@pytest.fixture
def session() -> Session:
    engine = create_engine("sqlite:///:memory:", future=True)
    Base.metadata.create_all(engine)
    with Session(engine, expire_on_commit=False) as sess:
        yield sess


def _mk_run(sess: Session, stage: str = "ingest") -> PipelineRun:
    run = PipelineRun()
    run.stage = stage
    run.schema_version = "0.1.0"
    run.storage_schema_version = STORAGE_SCHEMA_VERSION
    sess.add(run)
    sess.flush()
    return run


class TestSchemaVersion:
    def test_pinned(self) -> None:
        assert STORAGE_SCHEMA_VERSION == "0.1.0"


class TestPipelineRun:
    def test_insert_and_close(self, session: Session) -> None:
        run = _mk_run(session, stage="ingest")
        session.commit()
        assert run.id > 0
        assert run.started_at is not None
        assert run.completed_at is None

    def test_close(self, session: Session) -> None:
        run = _mk_run(session)
        run.completed_at = datetime.now(timezone.utc)
        run.row_count = 12345
        session.commit()
        loaded = session.get(PipelineRun, run.id)
        assert loaded is not None
        assert loaded.row_count == 12345


class TestSequenceUniqueness:
    def test_same_accession_same_schema_rejects(self, session: Session) -> None:
        run = _mk_run(session)
        s1 = Sequence()
        s1.accession = "AB123456"
        s1.schema_version = "0.1.0"
        s1.source = "genbank"
        s1.sequence = "ACGT" * 100
        s1.length = 400
        s1.pipeline_run_id = run.id
        s1.raw_record_hash = "a" * 64
        session.add(s1)
        session.commit()

        s2 = Sequence()
        s2.accession = "AB123456"
        s2.schema_version = "0.1.0"
        s2.source = "lanl"
        s2.sequence = "ACGT" * 100
        s2.length = 400
        s2.pipeline_run_id = run.id
        s2.raw_record_hash = "b" * 64
        session.add(s2)
        with pytest.raises(IntegrityError):
            session.commit()

    def test_same_accession_different_schema_ok(self, session: Session) -> None:
        run = _mk_run(session)
        for schema in ("0.1.0", "0.2.0"):
            s = Sequence()
            s.accession = "AB123456"
            s.schema_version = schema
            s.source = "genbank"
            s.sequence = "ACGT" * 100
            s.length = 400
            s.pipeline_run_id = run.id
            s.raw_record_hash = "a" * 64
            session.add(s)
        session.commit()
        count = session.execute(
            select(Sequence).where(Sequence.accession == "AB123456")
        ).scalars().all()
        assert len(count) == 2


class TestDRCallUniqueness:
    def test_uniq_on_seq_drug_algo(self, session: Session) -> None:
        run = _mk_run(session, stage="sierra")
        call = DRCall()
        call.accession = "AB123456"
        call.drug_class = "NNRTI"
        call.drug = "EFV"
        call.level = "high"
        call.algorithm_version = "9.6.1"
        call.pipeline_run_id = run.id
        session.add(call)
        session.commit()

        dup = DRCall()
        dup.accession = "AB123456"
        dup.drug_class = "NNRTI"
        dup.drug = "EFV"
        dup.level = "intermediate"
        dup.algorithm_version = "9.6.1"
        dup.pipeline_run_id = run.id
        session.add(dup)
        with pytest.raises(IntegrityError):
            session.commit()

    def test_different_algorithm_version_ok(self, session: Session) -> None:
        run = _mk_run(session, stage="sierra")
        for algo in ("9.6.0", "9.6.1"):
            call = DRCall()
            call.accession = "AB123456"
            call.drug_class = "NNRTI"
            call.drug = "EFV"
            call.level = "high"
            call.algorithm_version = algo
            call.pipeline_run_id = run.id
            session.add(call)
        session.commit()


class TestClusterMembership:
    def test_uniq_on_seq_run(self, session: Session) -> None:
        run = _mk_run(session, stage="cluster")
        m = ClusterMembership()
        m.accession = "AB123456"
        m.cluster_id = "C-0001"
        m.cluster_size = 12
        m.tn93_threshold = 0.015
        m.pipeline_run_id = run.id
        session.add(m)
        session.commit()

        dup = ClusterMembership()
        dup.accession = "AB123456"
        dup.cluster_id = "C-0002"
        dup.cluster_size = 5
        dup.tn93_threshold = 0.015
        dup.pipeline_run_id = run.id
        session.add(dup)
        with pytest.raises(IntegrityError):
            session.commit()

    def test_same_seq_different_run_ok(self, session: Session) -> None:
        r1 = _mk_run(session, stage="cluster")
        r2 = _mk_run(session, stage="cluster")
        for run in (r1, r2):
            m = ClusterMembership()
            m.accession = "AB123456"
            m.cluster_id = "C-0001"
            m.cluster_size = 12
            m.tn93_threshold = 0.015
            m.pipeline_run_id = run.id
            session.add(m)
        session.commit()
