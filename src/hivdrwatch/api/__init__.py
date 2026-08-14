"""Public FastAPI service for HIV-DR-Watch.

Read-only in v1. Writes (partner-clinical intake, KEMRI-side ingest)
land in v2 after the DSA is signed — currently blocked by CLAUDE.md's
`require_dsa()` gate.

Design:
  * Handlers are thin; the storage layer is the single source of truth.
  * OpenAPI auto-generated at `/docs`; JSON schema at `/openapi.json`.
  * Every request gets an `X-Request-Id` header propagated from
    middleware so cross-service tracing works out of the box.
  * Rate-limited via `slowapi` (per-IP token bucket).

Bootstrap:
    uvicorn hivdrwatch.api.main:app --host 0.0.0.0 --port 8080
"""

from .main import create_app

__all__ = ["create_app"]
