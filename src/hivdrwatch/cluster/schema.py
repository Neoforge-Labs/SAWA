"""Canonical schemas for phylogenetic clustering.

Two shapes:
- `PairwiseDistance` — one row per (accession1, accession2) below the
  clustering threshold. Sparse — above-threshold pairs are dropped.
- `ClusterMembership` — one row per (accession, cluster_id). A sequence
  can belong to at most one cluster in the connected-components model.
"""

from __future__ import annotations

from datetime import date
from pydantic import BaseModel, Field

from ..ingest.schema import AccessionStr

CLUSTER_SCHEMA_VERSION = "0.1.0"


class PairwiseDistance(BaseModel):
    """One TN93 distance below threshold."""

    accession_a: AccessionStr
    accession_b: AccessionStr
    distance: float = Field(..., ge=0.0, le=1.0)


class ClusterMembership(BaseModel):
    """One sequence's cluster assignment."""

    accession: AccessionStr
    cluster_id: str                     # deterministic, e.g. "C-0001"
    cluster_size: int
    tn93_threshold: float               # threshold used to define the cluster

    # --- provenance ---
    computed_at: date
    schema_version: str = CLUSTER_SCHEMA_VERSION
