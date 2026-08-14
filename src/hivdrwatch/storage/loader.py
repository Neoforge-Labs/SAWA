"""Bulk-load Parquet outputs into Postgres.

Every pipeline stage lands its result in a Parquet file first
(bronze/silver/gold layers). This module loads those Parquets into the
`sequences`, `dr_calls`, and `cluster_memberships` tables, opening a
`PipelineRun` row to tie every batch back to a canonical run.

Idempotency: `sequences` is unique on `(schema_version, accession)`,
`dr_calls` on `(accession, drug_class, drug, algorithm_version)`,
`cluster_memberships` on `(accession, pipeline_run_id)`. Reruns
upsert-or-skip via `ON CONFLICT DO NOTHING`.

We prefer SQLAlchemy Core `insert(...).on_conflict_do_nothing(...)`
over psycopg-level `COPY` — Core is portable across SQLite (tests)
and Postgres, whereas `COPY` is Postgres-only.
"""

from __future__ import annotations

import hashlib
import json
import logging
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import UTC, date, datetime
from pathlib import Path

import polars as pl
from sqlalchemy import inspect
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.orm import Session

from .schema import STORAGE_SCHEMA_VERSION, ClusterMembership, DRCall, PipelineRun, Sequence

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LoadCounts:
    """How many rows we inserted vs skipped as duplicates."""

    inserted: int
    skipped: int
    pipeline_run_id: int


def _iso(v: date | datetime | None) -> datetime | None:
    if v is None:
        return None
    if isinstance(v, datetime):
        return v if v.tzinfo else v.replace(tzinfo=UTC)
    return datetime(v.year, v.month, v.day, tzinfo=UTC)


def _config_hash(config: dict) -> str:
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode("utf-8")).hexdigest()


def _open_run(
    session: Session,
    *,
    stage: str,
    schema_version: str,
    git_sha: str | None,
    config: dict | None,
    notes: str | None = None,
) -> PipelineRun:
    run = PipelineRun()
    run.stage = stage
    run.schema_version = schema_version
    run.storage_schema_version = STORAGE_SCHEMA_VERSION
    run.git_sha = git_sha
    run.config_hash = _config_hash(config) if config else None
    run.notes = notes
    session.add(run)
    session.flush()  # populate `id`
    logger.info("pipeline run %d opened for stage=%s", run.id, stage)
    return run


def _close_run(session: Session, run: PipelineRun, row_count: int) -> None:
    run.completed_at = datetime.now(UTC)
    run.row_count = row_count
    session.add(run)
    session.flush()


def _insert_many(session: Session, table, rows: Iterable[dict]) -> tuple[int, int]:
    """Insert with dialect-specific idempotency; return (inserted, skipped)."""
    rows = list(rows)
    if not rows:
        return 0, 0

    dialect = session.get_bind().dialect.name
    if dialect == "postgresql":
        stmt = pg_insert(table).values(rows).on_conflict_do_nothing()
        result = session.execute(stmt)
        inserted = result.rowcount if result.rowcount is not None else len(rows)
        return inserted, len(rows) - inserted

    if dialect == "sqlite":
        # `INSERT OR IGNORE` is a single-statement idempotent insert on
        # SQLite. Batching per-row with a savepoint is O(n) round-trips
        # and, more importantly, `session.rollback()` on IntegrityError
        # rolls back the current transaction (including the pipeline_run
        # row we just inserted) — so we can't fall back to try/except
        # per row without breaking the run linkage.
        from sqlalchemy.dialects.sqlite import insert as sqlite_insert

        stmt = sqlite_insert(table).values(rows).prefix_with("OR IGNORE")
        result = session.execute(stmt)
        inserted = result.rowcount if result.rowcount is not None else len(rows)
        return inserted, len(rows) - inserted

    # Generic fallback — one savepoint per row.
    from sqlalchemy import insert

    inserted = 0
    for row in rows:
        savepoint = session.begin_nested()
        try:
            session.execute(insert(table).values(**row))
            savepoint.commit()
            inserted += 1
        except Exception:  # duplicate
            savepoint.rollback()
    return inserted, len(rows) - inserted


def load_sequences(
    parquet_path: Path,
    *,
    session: Session,
    git_sha: str | None = None,
    stage: str = "ingest",
) -> LoadCounts:
    """Load a bronze sequences Parquet into `sequences` + open a run row."""
    df = pl.read_parquet(parquet_path)
    if df.is_empty():
        logger.warning("empty sequences parquet: %s", parquet_path)
        return LoadCounts(0, 0, -1)

    schema_version = (
        str(df["schema_version"][0]) if "schema_version" in df.columns else "unknown"
    )
    run = _open_run(
        session,
        stage=stage,
        schema_version=schema_version,
        git_sha=git_sha,
        config={"parquet": str(parquet_path)},
    )
    now = datetime.now(UTC)
    rows: list[dict] = []
    for row in df.iter_rows(named=True):
        rows.append(
            {
                "accession": row["accession"],
                "schema_version": row.get("schema_version") or schema_version,
                "source": str(row.get("source", "unknown")),
                "sequence": row.get("sequence", ""),
                "length": int(row.get("length") or len(row.get("sequence") or "")),
                "gene_region": row.get("gene_region") or "pol",
                "country": row.get("country"),
                "collection_year": row.get("collection_year"),
                "subtype": str(row.get("subtype") or "unknown"),
                "subtype_source": row.get("subtype_source"),
                "is_drug_naive": row.get("is_drug_naive"),
                "risk_group": row.get("risk_group"),
                "patient_id_hashed": row.get("patient_id_hashed"),
                "pubmed_id": row.get("pubmed_id"),
                "study_accession": row.get("study_accession"),
                "pipeline_run_id": run.id,
                "raw_record_hash": row.get("raw_record_hash") or "",
                "ingested_at": _iso(row.get("ingested_at")) or now,
            }
        )

    inserted, skipped = _insert_many(session, Sequence.__table__, rows)
    _close_run(session, run, row_count=inserted)
    logger.info("sequences loaded: %d inserted, %d skipped", inserted, skipped)
    return LoadCounts(inserted, skipped, run.id)


def load_dr_calls(
    parquet_path: Path,
    *,
    session: Session,
    git_sha: str | None = None,
    stage: str = "sierra",
) -> LoadCounts:
    """Load a Sierra dr_calls Parquet into `dr_calls`."""
    df = pl.read_parquet(parquet_path)
    if df.is_empty():
        logger.warning("empty dr_calls parquet: %s", parquet_path)
        return LoadCounts(0, 0, -1)

    schema_version = (
        str(df["algorithm_version"][0]) if "algorithm_version" in df.columns else "unknown"
    )
    run = _open_run(
        session,
        stage=stage,
        schema_version=schema_version,
        git_sha=git_sha,
        config={"parquet": str(parquet_path)},
    )
    now = datetime.now(UTC)
    rows: list[dict] = []
    for row in df.iter_rows(named=True):
        rows.append(
            {
                "accession": row["accession"],
                "drug_class": str(row.get("drug_class", "unknown")),
                "drug": str(row.get("drug", "unknown")),
                "level": str(row.get("level", "unknown")),
                "mutations": row.get("mutations"),
                "algorithm_version": row.get("algorithm_version") or schema_version,
                "pipeline_run_id": run.id,
                "called_at": _iso(row.get("called_at")) or now,
            }
        )

    inserted, skipped = _insert_many(session, DRCall.__table__, rows)
    _close_run(session, run, row_count=inserted)
    logger.info("dr_calls loaded: %d inserted, %d skipped", inserted, skipped)
    return LoadCounts(inserted, skipped, run.id)


def load_cluster_memberships(
    parquet_path: Path,
    *,
    session: Session,
    git_sha: str | None = None,
    stage: str = "cluster",
) -> LoadCounts:
    """Load a cluster-assignments Parquet into `cluster_memberships`."""
    df = pl.read_parquet(parquet_path)
    if df.is_empty():
        return LoadCounts(0, 0, -1)

    schema_version = (
        str(df["schema_version"][0]) if "schema_version" in df.columns else "unknown"
    )
    run = _open_run(
        session,
        stage=stage,
        schema_version=schema_version,
        git_sha=git_sha,
        config={"parquet": str(parquet_path)},
    )
    now = datetime.now(UTC)
    rows: list[dict] = []
    for row in df.iter_rows(named=True):
        rows.append(
            {
                "accession": row["accession"],
                "cluster_id": row["cluster_id"],
                "cluster_size": int(row.get("cluster_size", 0)),
                "tn93_threshold": float(row.get("tn93_threshold", 0.015)),
                "pipeline_run_id": run.id,
                "computed_at": _iso(row.get("computed_at")) or now,
            }
        )

    inserted, skipped = _insert_many(session, ClusterMembership.__table__, rows)
    _close_run(session, run, row_count=inserted)
    logger.info("cluster_memberships loaded: %d inserted, %d skipped", inserted, skipped)
    return LoadCounts(inserted, skipped, run.id)


def ensure_schema(session: Session) -> None:
    """Create every table if it doesn't exist — dev convenience.

    Production uses Alembic migrations instead. Only call this in
    tests or a first-time local setup.
    """
    from .schema import Base

    engine = session.get_bind()
    inspector = inspect(engine)
    existing = set(inspector.get_table_names())
    missing = set(Base.metadata.tables) - existing
    if missing:
        logger.warning("creating tables: %s (dev-only path)", sorted(missing))
        Base.metadata.create_all(engine)


__all__ = [
    "LoadCounts",
    "ensure_schema",
    "load_cluster_memberships",
    "load_dr_calls",
    "load_sequences",
]
