"""Persistence layer for HIV-DR-Watch.

Two backends, one contract:
  * Postgres (with PostGIS for country centroids) is the system of
    record for the live surveillance service.
  * DuckDB is the read-side analytics engine — reads Parquet directly
    for ad-hoc queries and reads Postgres via the `postgres_scan`
    extension when hosted.

Both speak SQLModel/SQLAlchemy tables defined in `schema.py`. Migrations
are Alembic-managed under `migrations/`. The `loader` module handles
bulk Parquet → Postgres loads via `COPY FROM STDIN`.

Design decisions:
  * **One `pipeline_runs` table** ties every downstream row back to a
    reproducible run (schema version, git SHA, config hash).
  * **`(schema_version, accession)` uniqueness** so schema-breaking
    changes don't collide with prior data.
  * **All timestamps UTC** stored as TIMESTAMPTZ.
"""

from .schema import (
    STORAGE_SCHEMA_VERSION,
    ClusterMembership,
    DRCall,
    PipelineRun,
    Sequence,
    metadata,
)
from .session import get_engine, get_session

__all__ = [
    "STORAGE_SCHEMA_VERSION",
    "ClusterMembership",
    "DRCall",
    "PipelineRun",
    "Sequence",
    "get_engine",
    "get_session",
    "metadata",
]
