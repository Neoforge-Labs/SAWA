"""SHA256-keyed cache for Stanford HIVdb Sierra responses.

Sierra output is a pure function of `(sequence, algorithm_version)`. We key
the cache by `sha256(sequence)` and store the raw GraphQL response as JSON.
Downstream parsers replay from cache without re-hitting Sierra.

Cache layout:
    <cache_dir>/
        <sha256[0:2]>/<sha256>.json   # raw GraphQL response for one sequence
        algorithm_version.json         # last-observed algorithm version + date

Two-char sharding keeps directory sizes reasonable at 100k+ sequences.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def sequence_sha256(sequence: str) -> str:
    """Canonical hash: uppercased, U→T, whitespace stripped."""
    norm = sequence.upper().replace("U", "T").strip()
    return hashlib.sha256(norm.encode("utf-8")).hexdigest()


class SierraCache:
    """Filesystem cache for Sierra GraphQL responses, keyed by sequence hash."""

    def __init__(self, cache_dir: Path) -> None:
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ---- sequence-level cache ----

    def path_for(self, sha: str) -> Path:
        return self.cache_dir / sha[:2] / f"{sha}.json"

    def get(self, sha: str) -> dict[str, Any] | None:
        p = self.path_for(sha)
        if not p.exists():
            return None
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            logger.warning("corrupt cache entry %s; ignoring", p)
            return None

    def put(self, sha: str, payload: dict[str, Any]) -> None:
        p = self.path_for(sha)
        p.parent.mkdir(parents=True, exist_ok=True)
        # Atomic write: temp file + rename so a crash mid-write doesn't corrupt.
        tmp = p.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        tmp.replace(p)

    def has(self, sha: str) -> bool:
        return self.path_for(sha).exists()

    def missing(self, shas: list[str]) -> list[str]:
        return [s for s in shas if not self.has(s)]

    # ---- algorithm-version pin ----

    def _algo_path(self) -> Path:
        return self.cache_dir / "algorithm_version.json"

    def write_algorithm(self, version: str, publish_date: str | None) -> None:
        self._algo_path().write_text(
            json.dumps(
                {"version": version, "publish_date": publish_date},
                separators=(",", ":"),
            ),
            encoding="utf-8",
        )

    def read_algorithm(self) -> dict[str, Any] | None:
        p = self._algo_path()
        if not p.exists():
            return None
        return json.loads(p.read_text(encoding="utf-8"))
