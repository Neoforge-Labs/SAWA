"""GET /clusters — transmission-cluster membership listing."""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, Depends, Query
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ...storage.schema import ClusterMembership as ClusterRow
from ..dependencies import get_db
from ..schemas import ClusterListResponse, ClusterMember

router = APIRouter(prefix="/clusters", tags=["clusters"])


@router.get("", response_model=ClusterListResponse, summary="List cluster memberships")
def list_clusters(
    cluster_id: Optional[str] = Query(None, max_length=32),
    accession: Optional[str] = Query(None, max_length=20),
    min_cluster_size: Optional[int] = Query(None, ge=1),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
) -> ClusterListResponse:
    stmt = select(ClusterRow)
    if cluster_id:
        stmt = stmt.where(ClusterRow.cluster_id == cluster_id)
    if accession:
        stmt = stmt.where(ClusterRow.accession == accession.upper())
    if min_cluster_size is not None:
        stmt = stmt.where(ClusterRow.cluster_size >= min_cluster_size)

    total = db.execute(select(func.count()).select_from(stmt.subquery())).scalar_one()
    rows = db.execute(
        stmt.order_by(ClusterRow.cluster_id, ClusterRow.accession)
        .limit(limit)
        .offset(offset)
    ).scalars().all()

    return ClusterListResponse(
        items=[
            ClusterMember(
                accession=r.accession,
                cluster_id=r.cluster_id,
                cluster_size=r.cluster_size,
                tn93_threshold=r.tn93_threshold,
                computed_at=r.computed_at,
            )
            for r in rows
        ],
        total=int(total),
        limit=limit,
        offset=offset,
    )
