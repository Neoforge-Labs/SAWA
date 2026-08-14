"""Shared tail-task helper — Parquet → Postgres load.

Every pipeline flow (`ingest`, `dr_ingest`, `cluster`) ends with a call
to `storage_tail(stage, parquet_path, ...)` that:

  1. Checks whether loading is enabled (CLI flag, env, or explicit arg)
  2. Verifies the DB is reachable (`require_postgres`)
  3. Auto-creates tables on first run (dev/test convenience — production
     uses Alembic migrations)
  4. Calls the right loader for the stage
  5. Logs inserted/skipped counts

**Never crashes the flow.** If any step fails, the Parquet output on
disk is still the source of truth; loading is a best-effort side
effect that surveillance operators can retry.

The `require_dsa()` guard from CLAUDE.md is a HARD block: KEMRI code
paths must never invoke this helper. Enforced upstream in the KEMRI
adapters (Phase 2), not here — this helper doesn't need to know about
data provenance because it only runs on the LANL/GenBank public paths.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Callable
from pathlib import Path
from typing import Literal

from ..storage.loader import (
    LoadCounts,
    ensure_schema,
    load_cluster_memberships,
    load_dr_calls,
    load_sequences,
)
from ..storage.session import get_session, require_postgres

logger = logging.getLogger(__name__)

Stage = Literal["ingest", "sierra", "cluster"]

# Env var alternative to the --to-postgres CLI flag. Set to any truthy
# value to enable loading by default across every flow.
ENV_ENABLE = "HIVDRWATCH_LOAD_TO_POSTGRES"


def _env_enabled() -> bool:
    raw = os.getenv(ENV_ENABLE, "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


_LOADERS: dict[Stage, Callable[..., LoadCounts]] = {
    "ingest": load_sequences,
    "sierra": load_dr_calls,
    "cluster": load_cluster_memberships,
}


def storage_tail(
    stage: Stage,
    parquet_path: Path,
    *,
    enabled: bool | None = None,
    git_sha: str | None = None,
) -> LoadCounts | None:
    """Load `parquet_path` into Postgres for the given `stage`.

    Enable logic (any one is sufficient):
      * `enabled=True` — caller override, wins over env
      * `HIVDRWATCH_LOAD_TO_POSTGRES=1` env var

    Returns the `LoadCounts` on success, or None when skipped or when
    any step fails. Never raises.
    """
    if enabled is None:
        enabled = _env_enabled()
    if not enabled:
        logger.info(
            "storage tail skipped for stage=%s (pass --to-postgres or set %s=1 to enable)",
            stage,
            ENV_ENABLE,
        )
        return None

    if stage not in _LOADERS:
        logger.error("unknown storage-tail stage: %s", stage)
        return None

    if not parquet_path.exists():
        logger.warning("storage tail skipped: parquet not found at %s", parquet_path)
        return None

    if not require_postgres():
        # `require_postgres` already logged the WARNING with reason.
        return None

    loader = _LOADERS[stage]
    try:
        with get_session() as sess:
            ensure_schema(sess)
            counts = loader(parquet_path, session=sess, git_sha=git_sha)
    except Exception as e:
        logger.exception("storage tail failed for stage=%s: %s", stage, e)
        return None

    logger.info(
        "storage tail ok for stage=%s: inserted=%d skipped=%d run_id=%d",
        stage,
        counts.inserted,
        counts.skipped,
        counts.pipeline_run_id,
    )
    return counts


__all__ = ["ENV_ENABLE", "Stage", "storage_tail"]
