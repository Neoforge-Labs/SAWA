"""GET /sequences — filtered sequence listing."""

from __future__ import annotations

from fastapi import APIRouter, Depends, Query
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from ...storage.schema import Sequence as SequenceRow
from ..dependencies import get_db
from ..schemas import Sequence as SequenceResp
from ..schemas import SequenceListResponse

router = APIRouter(prefix="/sequences", tags=["sequences"])


@router.get("", response_model=SequenceListResponse, summary="List sequences")
def list_sequences(
    country: str | None = Query(
        None, min_length=2, max_length=2, description="ISO 3166-1 alpha-2"
    ),
    subtype: str | None = Query(None, max_length=16),
    collection_year: int | None = Query(None, ge=1980, le=2100),
    source: str | None = Query(None, max_length=16),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
) -> SequenceListResponse:
    stmt = select(SequenceRow)
    if country:
        stmt = stmt.where(SequenceRow.country == country.upper())
    if subtype:
        stmt = stmt.where(SequenceRow.subtype == subtype)
    if collection_year is not None:
        stmt = stmt.where(SequenceRow.collection_year == collection_year)
    if source:
        stmt = stmt.where(SequenceRow.source == source)

    total = db.execute(select(func.count()).select_from(stmt.subquery())).scalar_one()
    rows = (
        db.execute(stmt.order_by(SequenceRow.accession).limit(limit).offset(offset))
        .scalars()
        .all()
    )

    return SequenceListResponse(
        items=[
            SequenceResp(
                accession=r.accession,
                source=r.source,
                country=r.country,
                collection_year=r.collection_year,
                subtype=r.subtype,
                length=r.length,
                gene_region=r.gene_region,
                is_drug_naive=r.is_drug_naive,
            )
            for r in rows
        ],
        total=int(total),
        limit=limit,
        offset=offset,
    )
