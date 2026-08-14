"""GET /health — liveness probe."""

from __future__ import annotations

from datetime import datetime, timezone

from fastapi import APIRouter

from ..schemas import HealthResponse

router = APIRouter(tags=["health"])


@router.get("/health", response_model=HealthResponse, summary="Liveness probe")
def health() -> HealthResponse:
    return HealthResponse(status="ok", checked_at=datetime.now(timezone.utc))
