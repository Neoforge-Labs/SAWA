"""Loader tests — Parquet → SQLite via the same code path Postgres uses."""

from __future__ import annotations

from datetime import date
from pathlib import Path

import polars as pl
import pytest
from sqlalchemy import create_engine, func, select
from sqlalchemy.orm import Session

from hivdrwatch.storage.loader import (
    load_cluster_memberships,
    load_dr_calls,
    load_sequences,
)
from hivdrwatch.storage.schema import Base, ClusterMembership, DRCall, Sequence


@pytest.fixture
def session() -> Session:
    engine = create_engine("sqlite:///:memory:", future=True)
    Base.metadata.create_all(engine)
    with Session(engine, expire_on_commit=False) as sess:
        yield sess


def _write_sequences_parquet(path: Path, n: int) -> Path:
    df = pl.DataFrame(
        {
            "accession": [f"AB{i:06d}" for i in range(n)],
            "source": ["genbank"] * n,
            "sequence": ["ACGT" * 100] * n,
            "length": [400] * n,
            "country": ["KE"] * n,
            "collection_year": [2020] * n,
            "subtype": ["A1"] * n,
            "raw_record_hash": ["a" * 64] * n,
            "ingested_at": [date(2026, 8, 12)] * n,
            "schema_version": ["0.1.0"] * n,
        }
    )
    df.write_parquet(path)
    return path


def _write_dr_parquet(path: Path, n: int) -> Path:
    df = pl.DataFrame(
        {
            "accession": [f"AB{i:06d}" for i in range(n)],
            "drug_class": ["NNRTI"] * n,
            "drug": ["EFV"] * n,
            "level": ["high"] * n,
            "algorithm_version": ["9.6.1"] * n,
        }
    )
    df.write_parquet(path)
    return path


def _write_cluster_parquet(path: Path, n: int) -> Path:
    df = pl.DataFrame(
        {
            "accession": [f"AB{i:06d}" for i in range(n)],
            "cluster_id": [f"C-{i // 5:04d}" for i in range(n)],
            "cluster_size": [5] * n,
            "tn93_threshold": [0.015] * n,
            "schema_version": ["0.1.0"] * n,
        }
    )
    df.write_parquet(path)
    return path


class TestLoadSequences:
    def test_first_load_inserts_all(self, session: Session, tmp_path: Path) -> None:
        p = _write_sequences_parquet(tmp_path / "s.parquet", n=10)
        counts = load_sequences(p, session=session)
        assert counts.inserted == 10
        assert counts.skipped == 0
        assert session.execute(select(func.count()).select_from(Sequence)).scalar() == 10

    def test_reload_skips_duplicates(self, session: Session, tmp_path: Path) -> None:
        p = _write_sequences_parquet(tmp_path / "s.parquet", n=10)
        load_sequences(p, session=session)
        second = load_sequences(p, session=session)
        assert second.inserted == 0
        assert second.skipped == 10
        # Still only 10 rows total.
        assert session.execute(select(func.count()).select_from(Sequence)).scalar() == 10

    def test_pipeline_run_id_recorded(self, session: Session, tmp_path: Path) -> None:
        p = _write_sequences_parquet(tmp_path / "s.parquet", n=3)
        counts = load_sequences(p, session=session)
        assert counts.pipeline_run_id > 0
        rows = session.execute(select(Sequence)).scalars().all()
        assert all(r.pipeline_run_id == counts.pipeline_run_id for r in rows)


class TestLoadDRCalls:
    def test_inserts(self, session: Session, tmp_path: Path) -> None:
        p = _write_dr_parquet(tmp_path / "dr.parquet", n=5)
        counts = load_dr_calls(p, session=session)
        assert counts.inserted == 5
        assert session.execute(select(func.count()).select_from(DRCall)).scalar() == 5

    def test_dedup_same_algo(self, session: Session, tmp_path: Path) -> None:
        p = _write_dr_parquet(tmp_path / "dr.parquet", n=5)
        load_dr_calls(p, session=session)
        second = load_dr_calls(p, session=session)
        assert second.inserted == 0
        assert second.skipped == 5


class TestLoadClusters:
    def test_inserts(self, session: Session, tmp_path: Path) -> None:
        p = _write_cluster_parquet(tmp_path / "c.parquet", n=10)
        counts = load_cluster_memberships(p, session=session)
        assert counts.inserted == 10
        assert (
            session.execute(select(func.count()).select_from(ClusterMembership)).scalar() == 10
        )

    def test_second_run_ok_same_accession(self, session: Session, tmp_path: Path) -> None:
        # Two separate runs on the same accession set → both should
        # insert (uniqueness is per-run).
        p = _write_cluster_parquet(tmp_path / "c.parquet", n=5)
        first = load_cluster_memberships(p, session=session)
        second = load_cluster_memberships(p, session=session)
        assert first.inserted == 5
        assert second.inserted == 5
        assert first.pipeline_run_id != second.pipeline_run_id


class TestEmptyParquet:
    def test_empty_returns_zero(self, session: Session, tmp_path: Path) -> None:
        empty = pl.DataFrame(
            schema={
                "accession": pl.Utf8,
                "source": pl.Utf8,
                "sequence": pl.Utf8,
                "length": pl.Int64,
                "raw_record_hash": pl.Utf8,
                "schema_version": pl.Utf8,
            }
        )
        p = tmp_path / "empty.parquet"
        empty.write_parquet(p)
        counts = load_sequences(p, session=session)
        assert counts.inserted == 0
        assert counts.skipped == 0
