"""Data-governance gates — the KEMRI Data Sharing Agreement.

KEMRI clinical data and any PHI is untouchable until the DSA is signed
and Kelyn confirms in writing (CLAUDE.md, "Data and schemas"). This
module is the single place that rule is enforced, so a future KEMRI
adapter cannot quietly ship without it.

Enforcement has two layers:

  * `require_dsa(context)` — an explicit call any KEMRI-reachable entry
    point makes before doing work.
  * `SequenceRecord.source` — a validator in `ingest.schema` that calls
    `require_dsa` for every gated source. This is the choke point every
    ingestion path already crosses, so the gate holds even if a new
    adapter forgets the explicit call.

The gate is env-driven on purpose: flipping `HIVDRWATCH_DSA_SIGNED` is a
conscious, auditable act by an operator, not something a code path can
do for itself. Opening it logs at WARNING so it is visible in any run
log that is later audited.
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

DSA_ENV_VAR = "HIVDRWATCH_DSA_SIGNED"

# `Source` values that may only be processed once the DSA is signed.
# Kept as strings so `ingest.schema` can import this module without a
# circular import back into the schema it guards.
GATED_SOURCES: frozenset[str] = frozenset({"kemri"})

_AFFIRMATIVE: frozenset[str] = frozenset({"1", "true", "yes", "on"})


class DSANotSignedError(RuntimeError):
    """Raised when a DSA-gated code path runs before the DSA is signed."""


def dsa_is_signed() -> bool:
    """True iff the operator has declared the KEMRI DSA signed.

    Anything other than an explicit affirmative reads as unsigned — a
    typo'd env var must fail closed, never open.
    """
    return os.getenv(DSA_ENV_VAR, "").strip().lower() in _AFFIRMATIVE


def require_dsa(context: str) -> None:
    """Gate `context` behind the signed KEMRI DSA.

    Call this at the entry point of anything that could touch KEMRI
    data or PHI. `context` names the operation in both the exception
    and the audit log, so "kemri ingest" beats "guard".
    """
    if not dsa_is_signed():
        raise DSANotSignedError(
            f"{context} is blocked: the KEMRI Data Sharing Agreement is not "
            f"signed. This gate protects partner clinical data and PHI. Set "
            f"{DSA_ENV_VAR}=1 only once the DSA is executed and Kelyn has "
            f"confirmed in writing."
        )
    logger.warning(
        "DSA gate open (%s=1) — running %s against DSA-gated data",
        DSA_ENV_VAR,
        context,
    )


__all__ = [
    "DSA_ENV_VAR",
    "GATED_SOURCES",
    "DSANotSignedError",
    "dsa_is_signed",
    "require_dsa",
]
