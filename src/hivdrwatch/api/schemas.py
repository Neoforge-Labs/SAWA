"""Pydantic response models — the API's public contract.

Kept distinct from the SQLAlchemy tables in `hivdrwatch.storage.schema`
so we can evolve the DB shape independently of the wire format (and
never leak internal columns like `pipeline_run_id`).
"""

from __future__ import annotations

from datetime import datetime

from pydantic import BaseModel, Field


class HealthResponse(BaseModel):
    status: str = Field(default="ok")
    checked_at: datetime


class VersionResponse(BaseModel):
    hivdrwatch: str
    ingest_schema: str
    cluster_schema: str
    storage_schema: str
    api_schema: str = Field(default="0.1.0")


class Sequence(BaseModel):
    accession: str
    source: str
    country: str | None = None
    collection_year: int | None = None
    subtype: str
    length: int
    gene_region: str
    is_drug_naive: bool | None = None


class SequenceListResponse(BaseModel):
    items: list[Sequence]
    total: int
    limit: int
    offset: int


class DRCall(BaseModel):
    accession: str
    drug_class: str
    drug: str
    level: str
    algorithm_version: str
    called_at: datetime


class DRCallListResponse(BaseModel):
    items: list[DRCall]
    total: int
    limit: int
    offset: int


class ClusterMember(BaseModel):
    accession: str
    cluster_id: str
    cluster_size: int
    tn93_threshold: float
    computed_at: datetime


class ClusterListResponse(BaseModel):
    items: list[ClusterMember]
    total: int
    limit: int
    offset: int


class ErrorResponse(BaseModel):
    detail: str
    request_id: str | None = None


__all__ = [
    "ClusterListResponse",
    "ClusterMember",
    "DRCall",
    "DRCallListResponse",
    "ErrorResponse",
    "HealthResponse",
    "Sequence",
    "SequenceListResponse",
    "VersionResponse",
]
