"""SQLModel tables for the HIV-DR-Watch persistence layer.

Tables:
  * `pipeline_runs`      — one row per canonical pipeline run.
  * `sequences`          — one row per (schema_version, accession).
                           Country centroid stored as PostGIS POINT
                           when the extension is available; otherwise
                           just lat/lon columns.
  * `dr_calls`            — one row per (sequence, drug_class, drug).
                           Sierra algorithm version pinned in-row.
  * `cluster_memberships` — one row per (sequence, cluster_id, run).

Every downstream row carries `pipeline_run_id` — you can always trace
a resistance call back to the ingest run that produced its sequence.

## PostGIS handling

Country centroids use a plain `Column(Geometry("POINT", srid=4326))`
when PostGIS is installed. For local dev without PostGIS, we fall
back to `latitude`/`longitude` FLOAT columns (checked by `session.py`
on first connect).
"""

from __future__ import annotations

from datetime import UTC, datetime

from sqlalchemy import (
    JSON,
    DateTime,
    Float,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, MappedAsDataclass, mapped_column


def _utcnow() -> datetime:
    """Portable timestamp default — works on Postgres, SQLite, DuckDB."""
    return datetime.now(UTC)


STORAGE_SCHEMA_VERSION = "0.1.0"


class Base(MappedAsDataclass, DeclarativeBase, init=False):
    """Shared SQLAlchemy declarative base."""


class PipelineRun(Base):
    """One canonical pipeline run — every downstream row cites its id."""

    __tablename__ = "pipeline_runs"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    started_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )
    completed_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    stage: Mapped[str] = mapped_column(
        String(32), nullable=False
    )  # ingest, sierra, curate, ...
    schema_version: Mapped[str] = mapped_column(String(16), nullable=False)
    storage_schema_version: Mapped[str] = mapped_column(
        String(16), nullable=False, default=STORAGE_SCHEMA_VERSION
    )
    git_sha: Mapped[str | None] = mapped_column(String(40), nullable=True)
    config_hash: Mapped[str | None] = mapped_column(
        String(64), nullable=True
    )  # sha256 of run config
    row_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    notes: Mapped[str | None] = mapped_column(Text, nullable=True)


class Sequence(Base):
    """One canonical HIV-1 pol sequence with metadata + ingest provenance."""

    __tablename__ = "sequences"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(String(20), nullable=False, index=True)
    schema_version: Mapped[str] = mapped_column(String(16), nullable=False)
    source: Mapped[str] = mapped_column(String(16), nullable=False)  # genbank, lanl, kemri

    # Sequence + coord
    sequence: Mapped[str] = mapped_column(Text, nullable=False)
    length: Mapped[int] = mapped_column(Integer, nullable=False)
    gene_region: Mapped[str] = mapped_column(String(8), nullable=False, default="pol")

    # Provenance
    country: Mapped[str | None] = mapped_column(String(2), nullable=True, index=True)
    latitude: Mapped[float | None] = mapped_column(Float, nullable=True)
    longitude: Mapped[float | None] = mapped_column(Float, nullable=True)
    collection_year: Mapped[int | None] = mapped_column(Integer, nullable=True, index=True)
    subtype: Mapped[str] = mapped_column(
        String(16), nullable=False, default="unknown", index=True
    )
    subtype_source: Mapped[str | None] = mapped_column(String(16), nullable=True)

    # Clinical (rarely available)
    is_drug_naive: Mapped[bool | None] = mapped_column(nullable=True)
    risk_group: Mapped[str | None] = mapped_column(String(32), nullable=True)
    patient_id_hashed: Mapped[str | None] = mapped_column(String(64), nullable=True)

    # Publication
    pubmed_id: Mapped[str | None] = mapped_column(String(16), nullable=True)
    study_accession: Mapped[str | None] = mapped_column(String(32), nullable=True)

    # Ingestion trail
    pipeline_run_id: Mapped[int] = mapped_column(Integer, nullable=False, index=True)
    raw_record_hash: Mapped[str] = mapped_column(String(64), nullable=False)
    ingested_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    __table_args__ = (
        UniqueConstraint("schema_version", "accession", name="uq_sequences_schema_accession"),
    )


class DRCall(Base):
    """One Stanford HIVdb Sierra resistance call for one (sequence, drug)."""

    __tablename__ = "dr_calls"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(String(20), nullable=False, index=True)
    drug_class: Mapped[str] = mapped_column(
        String(8), nullable=False, index=True
    )  # PI/NRTI/NNRTI/INSTI
    drug: Mapped[str] = mapped_column(String(16), nullable=False)
    level: Mapped[str] = mapped_column(
        String(16), nullable=False, index=True
    )  # susceptible..high
    mutations: Mapped[dict | None] = mapped_column(JSON, nullable=True)  # {position: mutation}
    algorithm_version: Mapped[str] = mapped_column(String(16), nullable=False)

    pipeline_run_id: Mapped[int] = mapped_column(Integer, nullable=False, index=True)
    called_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    __table_args__ = (
        UniqueConstraint(
            "accession",
            "drug_class",
            "drug",
            "algorithm_version",
            name="uq_dr_calls_seq_drug_algo",
        ),
    )


class ClusterMembership(Base):
    """One sequence's cluster assignment under one clustering run."""

    __tablename__ = "cluster_memberships"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(String(20), nullable=False, index=True)
    cluster_id: Mapped[str] = mapped_column(String(32), nullable=False, index=True)
    cluster_size: Mapped[int] = mapped_column(Integer, nullable=False)
    tn93_threshold: Mapped[float] = mapped_column(Float, nullable=False)

    pipeline_run_id: Mapped[int] = mapped_column(Integer, nullable=False, index=True)
    computed_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=_utcnow
    )

    __table_args__ = (
        UniqueConstraint(
            "accession",
            "pipeline_run_id",
            name="uq_cluster_memberships_seq_run",
        ),
    )


# Re-export the SQLAlchemy metadata so Alembic env.py can pick it up
# without importing every table individually.
metadata = Base.metadata


__all__ = [
    "STORAGE_SCHEMA_VERSION",
    "Base",
    "ClusterMembership",
    "DRCall",
    "PipelineRun",
    "Sequence",
    "metadata",
]
