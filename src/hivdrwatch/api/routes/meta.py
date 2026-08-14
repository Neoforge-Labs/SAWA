"""GET /version — package + schema version stamps.

Schema-version imports are lazy so a missing/unimportable upstream
module (e.g. `ingest` when Biopython isn't installed) doesn't crash
the entire API at startup — `/version` degrades to `unknown` for that
component instead.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import logging

from fastapi import APIRouter

from ..schemas import VersionResponse

logger = logging.getLogger(__name__)

router = APIRouter(tags=["meta"])


def _try_version(module_path: str, attr: str) -> str:
    try:
        module = importlib.import_module(module_path)
        return str(getattr(module, attr))
    except Exception as e:  # noqa: BLE001 - defensive at import boundary
        logger.warning("could not load %s.%s: %s", module_path, attr, e)
        return "unknown"


@router.get("/version", response_model=VersionResponse, summary="Package + schema versions")
def version() -> VersionResponse:
    try:
        pkg_version = importlib.metadata.version("hivdrwatch")
    except importlib.metadata.PackageNotFoundError:
        pkg_version = "0.0.0+dev"

    return VersionResponse(
        hivdrwatch=pkg_version,
        ingest_schema=_try_version("hivdrwatch.ingest.schema", "SCHEMA_VERSION"),
        cluster_schema=_try_version("hivdrwatch.cluster.schema", "CLUSTER_SCHEMA_VERSION"),
        storage_schema=_try_version("hivdrwatch.storage.schema", "STORAGE_SCHEMA_VERSION"),
    )
