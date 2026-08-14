"""Stanford HIVdb Sierra client — GraphQL, batching, caching, offline replay.

Sierra endpoint:
    https://hivdb.stanford.edu/graphql

Design decisions:
- **Cache is the primary interface.** Every sequence goes through
  `SierraCache` first; only cache misses hit the network. This means CI can
  run `--offline` against a checked-in fixture cache, and re-runs cost zero.
- **Batch up to 100 sequences per GraphQL request.** Sierra accepts larger
  bodies but 100 is the sweet spot for latency + error blast radius.
- **Pin the algorithm version on every call.** Sierra returns
  `currentVersion { text publishDate }` — we record that into `DRResult`
  rows so a Stanford algorithm bump doesn't silently change labels.
- **Retries + rate limiting.** Sierra is polite but not infinite; we
  bound concurrency to 4 and back off on 5xx / rate-limit responses.

Usage:
    from hivdrwatch.resistance.sierra import SierraClient
    from hivdrwatch.resistance.cache import SierraCache

    cache = SierraCache(Path("data/raw/hivdb"))
    client = SierraClient(cache=cache)
    results, subtypes = client.analyze_sequences(
        [("AF004885", "ACGT..."), ("MK123456", "GTCA...")],
    )
"""

from __future__ import annotations

import logging
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import date, datetime
from typing import Any

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from .cache import SierraCache, sequence_sha256
from .schema import DRResult, DrugClass, ResistanceLevel, SubtypeCall

logger = logging.getLogger(__name__)

SIERRA_URL = "https://hivdb.stanford.edu/graphql"
DEFAULT_BATCH_SIZE = 100
DEFAULT_TIMEOUT_S = 60.0

# Minimal GraphQL query — deliberately narrow. Grab only what DRResult needs.
GRAPHQL_QUERY = """
query AnalyzeSequences($sequences: [UnalignedSequenceInput!]!) {
  currentVersion { text publishDate }
  viewer {
    sequenceAnalysis(sequences: $sequences) {
      inputSequence { header }
      bestMatchingSubtype { display }
      drugResistance {
        gene { name }
        drugScores {
          drugClass { name }
          drug { name displayAbbr fullName }
          score
          text
        }
        mutationsByTypes {
          mutationType
          mutations { text }
        }
      }
    }
  }
}
"""


@dataclass(frozen=True)
class AlgorithmVersion:
    version: str
    publish_date: date | None


class SierraError(RuntimeError):
    """Sierra returned an error we can't recover from."""


class SierraClient:
    """Thin, cache-first wrapper around the Sierra GraphQL endpoint."""

    def __init__(
        self,
        cache: SierraCache,
        *,
        endpoint: str = SIERRA_URL,
        batch_size: int = DEFAULT_BATCH_SIZE,
        offline: bool = False,
        client: httpx.Client | None = None,
    ) -> None:
        if batch_size < 1 or batch_size > 100:
            raise ValueError("batch_size must be 1..100")
        self.cache = cache
        self.endpoint = endpoint
        self.batch_size = batch_size
        self.offline = offline
        # Injectable client for tests; owned lifecycle if we create it.
        self._client = client
        self._owns_client = client is None

    def __enter__(self) -> "SierraClient":
        if self._client is None:
            self._client = httpx.Client(timeout=DEFAULT_TIMEOUT_S)
        return self

    def __exit__(self, *exc: object) -> None:
        if self._owns_client and self._client is not None:
            self._client.close()
            self._client = None

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def analyze_sequences(
        self, sequences: Iterable[tuple[str, str]]
    ) -> tuple[list[DRResult], list[SubtypeCall]]:
        """Analyze `(accession, sequence)` pairs. Returns (dr_results, subtypes).

        Cache-first: sequences whose sha256 is already cached don't hit
        the network. In offline mode a cache miss raises SierraError.
        """
        pairs = list(sequences)
        # Deduplicate by sequence hash — many accessions can share a sequence.
        by_sha: dict[str, list[str]] = {}
        seq_of_sha: dict[str, str] = {}
        for accession, seq in pairs:
            sha = sequence_sha256(seq)
            by_sha.setdefault(sha, []).append(accession)
            seq_of_sha.setdefault(sha, seq)

        misses = [sha for sha in by_sha if not self.cache.has(sha)]
        logger.info(
            "Sierra: %d sequences (%d unique), %d cache hits, %d misses",
            len(pairs),
            len(by_sha),
            len(by_sha) - len(misses),
            len(misses),
        )

        if misses:
            if self.offline:
                raise SierraError(
                    f"offline mode but {len(misses)} cache misses; "
                    f"first: {misses[0]}"
                )
            self._fetch_and_cache(misses, seq_of_sha)

        dr_results: list[DRResult] = []
        subtypes: list[SubtypeCall] = []
        today = date.today()
        for sha, accessions in by_sha.items():
            payload = self.cache.get(sha)
            if payload is None:
                raise SierraError(f"expected cache entry after fetch: {sha}")
            algo = _parse_algorithm(payload)
            for accession in accessions:
                dr_results.extend(_parse_dr_results(payload, accession, sha, algo, today))
                st = _parse_subtype(payload, accession, sha, algo, today)
                if st is not None:
                    subtypes.append(st)
        return dr_results, subtypes

    def current_algorithm(self) -> AlgorithmVersion:
        """Fetch the current HIVdb algorithm version. Also writes to cache."""
        if self.offline:
            cached = self.cache.read_algorithm()
            if cached is None:
                raise SierraError("offline mode and no cached algorithm version")
            return AlgorithmVersion(
                version=cached["version"],
                publish_date=_maybe_date(cached.get("publish_date")),
            )
        # Send a trivial query with no sequences to get currentVersion.
        payload = self._post(
            "query { currentVersion { text publishDate } }", variables={}
        )
        cur = payload["data"]["currentVersion"]
        algo = AlgorithmVersion(
            version=cur["text"], publish_date=_maybe_date(cur.get("publishDate"))
        )
        self.cache.write_algorithm(
            algo.version,
            algo.publish_date.isoformat() if algo.publish_date else None,
        )
        return algo

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------

    def _fetch_and_cache(
        self, missing_shas: list[str], seq_of_sha: dict[str, str]
    ) -> None:
        for batch in _chunked(missing_shas, self.batch_size):
            inputs = [
                {"header": sha, "sequence": seq_of_sha[sha]} for sha in batch
            ]
            payload = self._post(GRAPHQL_QUERY, {"sequences": inputs})
            self._split_and_cache(payload, batch)

    def _split_and_cache(
        self, payload: dict[str, Any], sha_order: list[str]
    ) -> None:
        """Sierra returns one analysis object per input, in order. We slice
        the batch response into per-sequence payloads so each cache entry
        stands alone and is replayable."""
        analyses = (
            payload.get("data", {})
            .get("viewer", {})
            .get("sequenceAnalysis", [])
        )
        if len(analyses) != len(sha_order):
            raise SierraError(
                f"Sierra returned {len(analyses)} analyses for "
                f"{len(sha_order)} inputs"
            )
        current_version = payload.get("data", {}).get("currentVersion", {})
        for sha, analysis in zip(sha_order, analyses, strict=True):
            per_seq = {
                "currentVersion": current_version,
                "analysis": analysis,
            }
            self.cache.put(sha, per_seq)

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=60),
        retry=retry_if_exception_type(
            (httpx.HTTPError, httpx.TimeoutException, SierraError)
        ),
        reraise=True,
    )
    def _post(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        if self._client is None:
            self._client = httpx.Client(timeout=DEFAULT_TIMEOUT_S)
            self._owns_client = True
        # Gentle client-side rate limit — Sierra doesn't publish a hard limit
        # but we've observed 429s at aggressive concurrencies.
        time.sleep(0.05)
        resp = self._client.post(
            self.endpoint, json={"query": query, "variables": variables}
        )
        if resp.status_code >= 500 or resp.status_code == 429:
            raise SierraError(f"Sierra transient error: {resp.status_code}")
        if resp.status_code != 200:
            raise SierraError(
                f"Sierra error {resp.status_code}: {resp.text[:400]}"
            )
        payload = resp.json()
        if "errors" in payload:
            raise SierraError(f"Sierra GraphQL errors: {payload['errors']}")
        return payload


# ---- pure parsers (unit-testable without network) ----


def _parse_algorithm(cached_payload: dict[str, Any]) -> AlgorithmVersion:
    cur = cached_payload.get("currentVersion", {}) or {}
    return AlgorithmVersion(
        version=str(cur.get("text", "unknown")),
        publish_date=_maybe_date(cur.get("publishDate")),
    )


def _parse_dr_results(
    cached_payload: dict[str, Any],
    accession: str,
    sha: str,
    algo: AlgorithmVersion,
    called_at: date,
) -> list[DRResult]:
    analysis = cached_payload.get("analysis") or {}
    out: list[DRResult] = []
    for dr in analysis.get("drugResistance", []) or []:
        gene = (dr.get("gene") or {}).get("name") or ""
        mutations = _flatten_mutations(dr.get("mutationsByTypes") or [])
        for score_entry in dr.get("drugScores") or []:
            drug_class_name = (score_entry.get("drugClass") or {}).get("name")
            if drug_class_name not in DrugClass.__members__:
                continue
            drug_info = score_entry.get("drug") or {}
            out.append(
                DRResult(
                    accession=accession,
                    sequence_sha256=sha,
                    algorithm_version=algo.version,
                    algorithm_publish_date=algo.publish_date,
                    drug_class=DrugClass[drug_class_name],
                    drug=drug_info.get("displayAbbr") or drug_info.get("name") or "",
                    drug_display=drug_info.get("fullName"),
                    score=float(score_entry.get("score") or 0.0),
                    level=ResistanceLevel.from_sierra_text(
                        str(score_entry.get("text") or "Susceptible")
                    ),
                    gene=gene,
                    mutations=mutations,
                    called_at=called_at,
                )
            )
    return out


def _parse_subtype(
    cached_payload: dict[str, Any],
    accession: str,
    sha: str,
    algo: AlgorithmVersion,
    called_at: date,
) -> SubtypeCall | None:
    analysis = cached_payload.get("analysis") or {}
    subtype = (analysis.get("bestMatchingSubtype") or {}).get("display")
    if not subtype:
        return None
    return SubtypeCall(
        accession=accession,
        sequence_sha256=sha,
        subtype=str(subtype),
        algorithm_version=algo.version,
        called_at=called_at,
    )


def _flatten_mutations(mutations_by_types: list[dict[str, Any]]) -> list[str]:
    out: list[str] = []
    for entry in mutations_by_types:
        for m in entry.get("mutations") or []:
            text = m.get("text")
            if text:
                out.append(str(text))
    return sorted(set(out))


def _maybe_date(value: str | None) -> date | None:
    if not value:
        return None
    # Sierra publishDate: "2024-08-15" or ISO datetime — accept both.
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).date()
    except ValueError:
        try:
            return datetime.strptime(value, "%Y-%m-%d").date()
        except ValueError:
            return None


def _chunked(items: list[str], size: int) -> Iterator[list[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]
