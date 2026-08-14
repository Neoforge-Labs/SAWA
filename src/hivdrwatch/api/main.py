"""FastAPI app factory.

Split into `create_app()` so tests can build a fresh app without
inheriting a real DB session or a locked rate-limiter state.
"""

from __future__ import annotations

import logging
from datetime import UTC, datetime

from fastapi import FastAPI, Request, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address

from .config import get_settings
from .middleware import REQUEST_ID_HEADER, RequestIdMiddleware
from .routes import clusters, dr_calls, health, meta, sequences
from .schemas import ErrorResponse

logger = logging.getLogger(__name__)


def create_app() -> FastAPI:
    """Build a fully-wired FastAPI app.

    Called once at import time (`app = create_app()`); called per-test
    from `tests/api/` fixtures with dependency overrides in place.
    """
    settings = get_settings()

    app = FastAPI(
        title="SAWA API",
        description=(
            "SAWA — Sub-Saharan Africa Watch for Antiretroviral resistance. "
            "Open real-time HIV-1 drug-resistance surveillance for "
            "sub-Saharan Africa. Read-only in v1; writes gated behind the "
            "KEMRI Data Sharing Agreement. (Python package: `hivdrwatch`.)"
        ),
        version="0.1.0",
        root_path=settings.root_path,
        openapi_url="/openapi.json" if settings.enable_docs else None,
        docs_url="/docs" if settings.enable_docs else None,
        redoc_url="/redoc" if settings.enable_docs else None,
    )

    # --- Middleware ---
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cors_origins,
        allow_credentials=settings.cors_allow_credentials,
        allow_methods=["GET", "OPTIONS"],  # read-only in v1
        allow_headers=["*"],
        expose_headers=[REQUEST_ID_HEADER, "X-Response-Time-ms"],
    )
    app.add_middleware(RequestIdMiddleware)

    # --- Rate limiter ---
    limiter = Limiter(
        key_func=get_remote_address, default_limits=[settings.rate_limit_default]
    )
    app.state.limiter = limiter
    app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

    # --- Routers ---
    app.include_router(health.router)
    app.include_router(meta.router)
    app.include_router(sequences.router)
    app.include_router(dr_calls.router)
    app.include_router(clusters.router)

    # --- Uniform error envelope ---
    @app.exception_handler(Exception)
    async def _unhandled(request: Request, exc: Exception) -> JSONResponse:
        rid = getattr(request.state, "request_id", None)
        logger.exception("unhandled exception rid=%s", rid)
        return JSONResponse(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            content=ErrorResponse(detail="internal server error", request_id=rid).model_dump(),
        )

    logger.info(
        "hivdrwatch api ready (cors=%s, rate_limit=%s)",
        settings.cors_origins,
        settings.rate_limit_default,
    )
    return app


# Module-level `app` for `uvicorn hivdrwatch.api.main:app`.
app = create_app()


def _now_utc() -> datetime:
    return datetime.now(UTC)


__all__ = ["app", "create_app"]
