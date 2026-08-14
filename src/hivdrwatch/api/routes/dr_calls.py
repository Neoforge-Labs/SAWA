"""GET /dr-calls — Sierra resistance-call listing."""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, Depends, Query
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ...storage.schema import DRCall as DRCallRow
from ..dependencies import get_db
from ..schemas import DRCall as DRCallResp, DRCallListResponse

router = APIRouter(prefix="/dr-calls", tags=["dr-calls"])


@router.get("", response_model=DRCallListResponse, summary="List Sierra resistance calls")
def list_dr_calls(
    accession: Optional[str] = Query(None, max_length=20),
    drug_class: Optional[str] = Query(None, max_length=8),
    level: Optional[str] = Query(None, max_length=16),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
) -> DRCallListResponse:
    stmt = select(DRCallRow)
    if accession:
        stmt = stmt.where(DRCallRow.accession == accession.upper())
    if drug_class:
        stmt = stmt.where(DRCallRow.drug_class == drug_class)
    if level:
        stmt = stmt.where(DRCallRow.level == level)

    total = db.execute(select(func.count()).select_from(stmt.subquery())).scalar_one()
    rows = db.execute(
        stmt.order_by(DRCallRow.accession, DRCallRow.drug_class, DRCallRow.drug)
        .limit(limit)
        .offset(offset)
    ).scalars().all()

    return DRCallListResponse(
        items=[
            DRCallResp(
                accession=r.accession,
                drug_class=r.drug_class,
                drug=r.drug,
                level=r.level,
                algorithm_version=r.algorithm_version,
                called_at=r.called_at,
            )
            for r in rows
        ],
        total=int(total),
        limit=limit,
        offset=offset,
    )
