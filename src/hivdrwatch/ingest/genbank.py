"""Pull HIV-1 pol sequences from GenBank via NCBI E-utilities.

Design notes:
- E-utilities is polite-by-request: 3 req/s without an API key, 10 req/s with.
  Set NCBI_API_KEY in env; the fetcher uses it if present.
- We batch efetch calls (default 200/batch) — larger risks timeouts, smaller
  is chatty.
- Every fetched record is written to raw/ as a GenBank flat file before
  parsing. If parsing changes, we replay from raw/ without re-hitting NCBI.
- The `since` param lets us do incremental pulls: pass the last successful
  ingestion date and only new deposits come back.

Not implemented yet (intentional, day-1 minimum):
- Parallel fetches. E-utilities rate limits make this a foot-gun.
- Async. Add when the surveillance service needs it; the batch ingester
  doesn't.
"""

from __future__ import annotations

import hashlib
import http.client
import logging
import os
import re
import time
from collections.abc import Iterator
from datetime import date
from pathlib import Path
from typing import Any

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from .schema import (
    SequenceRecord,
    Source,
    Subtype,
    country_name_to_code,
    subtype_from_designation,
)

logger = logging.getLogger(__name__)

# East African priority countries (ISO 3166-1 alpha-2 → GenBank country name)
EAST_AFRICA_COUNTRIES: dict[str, str] = {
    "KE": "Kenya",
    "UG": "Uganda",
    "TZ": "Tanzania",
    "RW": "Rwanda",
    "BI": "Burundi",
    "ET": "Ethiopia",
    "SS": "South Sudan",
    "SO": "Somalia",
}

# NCBI E-utilities gets flaky above ~50 records under load. Empirically
# batches of 50 succeed on the first try >95% of the time; larger sizes
# invite `IncompleteRead` even with retries.
DEFAULT_BATCH_SIZE = 50

# When a batch of size N fails all its retries, bisect: try N/2 twice.
# Recurse until batch size hits this floor, at which point failures are
# surfaced (something's actually wrong with the specific record).
MIN_BISECTION_BATCH = 5

# Exceptions we consider transient and worth retrying on inside a batch.
# `IncompleteRead` is the big one — NCBI cuts responses mid-stream when
# it's under load. Keep this list narrow so real bugs still surface.
_TRANSIENT_EXCEPTIONS = (
    OSError,
    RuntimeError,
    http.client.IncompleteRead,
    http.client.HTTPException,
)
MIN_POL_LENGTH = 800  # under this and phylogenetic signal is too weak

# HXB2 pol is 3,012 nt. A deposit comfortably longer than that spans other
# genes too (the gag→nef submissions are 8.6-8.8 kb) and must not be
# labelled "pol"; the 3,012-4,000 band is pol plus flanking sequence.
MAX_POL_LENGTH = 4000


class GenBankIngester:
    """Fetch, parse, and validate HIV pol sequences from NCBI.

    Usage:
        ing = GenBankIngester(email="you@neoforge.io", raw_dir=Path("data/raw/genbank"))
        records = list(ing.fetch_country("KE", since=date(2020, 1, 1)))
    """

    def __init__(
        self,
        email: str,
        raw_dir: Path,
        api_key: str | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> None:
        if not email or "@" not in email:
            raise ValueError("NCBI requires a real email address for E-utilities access")
        Entrez.email = email
        Entrez.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.raw_dir = raw_dir
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.batch_size = batch_size
        # Rate limit: 10/s with key, 3/s without. We stay well under.
        self._min_interval = 0.11 if Entrez.api_key else 0.34
        self._last_call = 0.0

    # ---- public API ----

    def fetch_country(
        self,
        country_code: str,
        since: date | None = None,
        until: date | None = None,
        min_length: int = MIN_POL_LENGTH,
    ) -> Iterator[SequenceRecord]:
        """Yield validated SequenceRecords for one country.

        `since`/`until` filter by GenBank *modification* date, not sample
        collection date. That's deliberate — for incremental ingestion we
        want anything newly deposited, not anything newly collected.
        """
        country_name = EAST_AFRICA_COUNTRIES.get(country_code.upper())
        if country_name is None:
            raise ValueError(f"unknown country code: {country_code}")

        query = self._build_query(country_name, since, until, min_length)
        logger.info("GenBank query: %s", query)

        ids = self._esearch_all(query)
        logger.info("found %d matching accessions for %s", len(ids), country_code)

        n_batches = (len(ids) + self.batch_size - 1) // self.batch_size
        for batch_ix, start in enumerate(range(0, len(ids), self.batch_size)):
            batch = ids[start : start + self.batch_size]
            logger.info(
                "%s batch %d/%d (%d ids)",
                country_code,
                batch_ix + 1,
                n_batches,
                len(batch),
            )
            for raw_record in self._efetch_batch_resilient(batch):
                parsed = self._parse(raw_record, country_code)
                if parsed is not None:
                    yield parsed

    # ---- internals ----

    @staticmethod
    def _build_query(
        country_name: str,
        since: date | None,
        until: date | None,
        min_length: int,
    ) -> str:
        """Build an E-utilities query covering full-pol + subgene submissions.

        The previous query used `pol[Gene]` alone, which matched only records
        where the entire pol polyprotein was annotated. Most HIV drug-resistance
        surveillance labs submit **subregion** sequences (RT-only for NRTIs/NNRTIs,
        integrase-only for INSTIs, protease for PIs). Those get missed unless
        we OR across the subgene names.

        The OR clause matches gene-name variants (`gene[Gene]`, `gene[Product]`,
        `gene[Title]`) because NCBI annotation isn't consistent across depositors.
        Empirically bumps East Africa coverage from ~23k → ~45–60k sequences
        with reasonable PR/RT/IN mix instead of integrase-heavy skew.
        """
        gene_clause = (
            "("
            "pol[Gene] OR pol[Product] "
            "OR protease[Gene] OR protease[Product] "
            'OR "reverse transcriptase"[Product] OR RT[Gene] '
            "OR integrase[Gene] OR integrase[Product]"
            ")"
        )
        parts = [
            '"Human immunodeficiency virus 1"[Organism]',
            gene_clause,
            f'"{country_name}"[Country]',
            f"{min_length}:100000[SLEN]",
        ]
        if since or until:
            lo = since.strftime("%Y/%m/%d") if since else "1980/01/01"
            hi = until.strftime("%Y/%m/%d") if until else "3000/01/01"
            parts.append(f'("{lo}"[MDAT] : "{hi}"[MDAT])')
        return " AND ".join(parts)

    def _rate_limit(self) -> None:
        elapsed = time.monotonic() - self._last_call
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_call = time.monotonic()

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        retry=retry_if_exception_type(_TRANSIENT_EXCEPTIONS),
        reraise=True,
    )
    def _esearch_all(self, query: str) -> list[str]:
        """esearch with WebEnv to page through arbitrarily large result sets."""
        self._rate_limit()
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0)
        head = Entrez.read(handle)
        handle.close()
        total = int(head["Count"])
        if total == 0:
            return []
        webenv = head["WebEnv"]
        query_key = head["QueryKey"]

        ids: list[str] = []
        page = 10_000  # esearch retmax cap
        for start in range(0, total, page):
            self._rate_limit()
            h = Entrez.esearch(
                db="nucleotide",
                term=query,
                retstart=start,
                retmax=page,
                webenv=webenv,
                query_key=query_key,
            )
            chunk = Entrez.read(h)
            h.close()
            ids.extend(chunk["IdList"])
        return ids

    def _efetch_batch_resilient(self, ids: list[str]) -> Iterator[SeqRecord]:
        """Fetch with per-batch retries + bisection fallback.

        If the retry-decorated `_efetch_batch` still fails after all
        attempts, we assume this batch is too large for current network
        conditions (or contains a poison record) and split it in half.
        Recursion base case is `MIN_BISECTION_BATCH` — below that,
        failures surface as real errors.
        """
        try:
            yield from self._efetch_batch(ids)
            return
        except _TRANSIENT_EXCEPTIONS as e:
            # Never format the exception object into log strings directly —
            # some (e.g. http.client.IncompleteRead) have a __repr__ that
            # crashes when their `partial` field isn't bytes.
            err_name = type(e).__name__
            if len(ids) <= MIN_BISECTION_BATCH:
                logger.error(
                    "batch of %d failed after retries; giving up (type=%s) ids=%s",
                    len(ids),
                    err_name,
                    ids[:3] + (["..."] if len(ids) > 3 else []),
                )
                raise
            mid = len(ids) // 2
            logger.warning(
                "batch of %d exhausted retries (%s); bisecting to %d + %d",
                len(ids),
                type(e).__name__,
                mid,
                len(ids) - mid,
            )
            yield from self._efetch_batch_resilient(ids[:mid])
            yield from self._efetch_batch_resilient(ids[mid:])

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        retry=retry_if_exception_type(_TRANSIENT_EXCEPTIONS),
        reraise=True,
    )
    def _efetch_batch(self, ids: list[str]) -> Iterator[SeqRecord]:
        """Fetch a batch as GenBank flat file; also persist raw for replay.

        Resume-safe: if `<raw_dir>/batch_<first>_<last>.gb` already
        exists from a previous run AND passes the terminator check
        (`_looks_valid_gb`), we skip the network call and re-parse the
        cached file. Corrupt caches (partial writes from crashed prior
        runs before atomic-write landed) are automatically re-fetched.
        """
        raw_path = self.raw_dir / f"batch_{ids[0]}_{ids[-1]}.gb"

        if raw_path.exists():
            if _looks_valid_gb(raw_path):
                logger.info(
                    "cache hit: %s (%d bytes)",
                    raw_path.name,
                    raw_path.stat().st_size,
                )
                with raw_path.open() as fh:
                    yield from SeqIO.parse(fh, "genbank")
                return
            logger.warning("cache file %s looks truncated; re-fetching", raw_path.name)
            raw_path.unlink(missing_ok=True)

        self._rate_limit()
        logger.info("fetching batch of %d (first id=%s)", len(ids), ids[0])
        handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="gb", retmode="text")
        raw_text = handle.read()
        handle.close()

        # Validate BEFORE persisting so we don't cache a truncated response.
        if not _text_looks_valid_gb(raw_text):
            raise http.client.IncompleteRead(
                len(raw_text.encode("utf-8")),
                expected=None,
            )

        # Atomic write: temp file + rename so a mid-write crash can't
        # leave a truncated cache entry that we'd mistake for valid.
        tmp_path = raw_path.with_suffix(".gb.tmp")
        tmp_path.write_text(raw_text, encoding="utf-8")
        tmp_path.replace(raw_path)

        # Re-parse from disk so parsing is identical to replay.
        with raw_path.open() as fh:
            yield from SeqIO.parse(fh, "genbank")

    def _parse(self, record: SeqRecord, requested_country: str) -> SequenceRecord | None:
        try:
            metadata = _extract_source_features(record)
            collection_date = _parse_collection_date(metadata.get("collection_date"))

            # The record's own geography, never the query's. The country
            # search over-matches — Thai, Brazilian and South African
            # deposits all come back from East African queries — so
            # inheriting `requested_country` fabricated the one field two
            # papers stratify on.
            country_raw = _stated_country(metadata)
            country = _normalize_country(country_raw)
            if country_raw and country is None:
                logger.warning(
                    "%s: country %r is not in the name map; storing it verbatim "
                    "with no ISO code",
                    record.id,
                    country_raw,
                )
            elif country is not None and country != requested_country.upper():
                logger.info(
                    "%s: states %s but was returned by the %s query",
                    record.id,
                    country,
                    requested_country.upper(),
                )

            subtype = _normalize_subtype(metadata.get("note", "") or "")
            sequence = str(record.seq)
            record_hash = hashlib.sha256(sequence.encode("utf-8")).hexdigest()

            return SequenceRecord(
                accession=record.id,
                source=Source.GENBANK,
                sequence=sequence,
                length=len(sequence),
                gene_region=_infer_gene_region(len(sequence)),
                country=country,
                country_raw=country_raw,
                country_source="genbank_source_feature" if country_raw else None,
                collection_date=collection_date,
                collection_year=collection_date.year if collection_date else None,
                subtype=subtype,
                subtype_source="genbank" if subtype is not Subtype.UNKNOWN else None,
                pubmed_id=_extract_pubmed_id(record),
                ingested_at=date.today(),
                raw_record_hash=record_hash,
            )
        except ValidationError as e:
            logger.warning("dropping %s: schema validation failed: %s", record.id, e)
            return None
        except Exception as e:
            logger.exception("unexpected parse failure on %s: %s", record.id, e)
            return None


# ---- small pure helpers, testable in isolation ----


def _looks_valid_gb(path: Path) -> bool:
    """A cached GenBank flat file is valid iff its last non-blank line
    is exactly `//` (the record terminator). Complete responses always
    end that way; truncated responses don't.

    Only reads the tail of the file (last 32 bytes), so this is cheap
    even on multi-MB caches.
    """
    try:
        size = path.stat().st_size
        if size < 5:
            return False
        with path.open("rb") as fh:
            fh.seek(max(0, size - 32))
            tail = fh.read().decode("utf-8", errors="replace").rstrip()
        return tail.endswith("//")
    except OSError:
        return False


def _text_looks_valid_gb(text: str) -> bool:
    """In-memory equivalent of `_looks_valid_gb` — used before persisting."""
    return text.rstrip().endswith("//")


def _extract_source_features(record: SeqRecord) -> dict[str, Any]:
    """Flatten the `source` feature qualifiers into a plain dict."""
    for feature in record.features:
        if feature.type == "source":
            return {k: v[0] if v else None for k, v in feature.qualifiers.items()}
    return {}


def _parse_collection_date(value: str | None) -> date | None:
    """GenBank collection_date is a zoo: '2019', 'Jun-2019', '15-Jun-2019'."""
    if not value:
        return None
    from datetime import datetime

    for fmt in ("%d-%b-%Y", "%b-%Y", "%Y-%m-%d", "%Y"):
        try:
            parsed = datetime.strptime(value, fmt).date()
            # For year-only, pick July 1 as a neutral midpoint
            if fmt == "%Y":
                parsed = parsed.replace(month=7, day=1)
            elif fmt == "%b-%Y":
                parsed = parsed.replace(day=15)
            return parsed
        except ValueError:
            continue
    logger.debug("unparseable collection_date: %r", value)
    return None


def _stated_country(metadata: dict[str, Any]) -> str | None:
    """The country the record itself states, verbatim, or None.

    NCBI renamed the source qualifier `/country=` to `/geo_loc_name=` and
    `/geo_loc_name=` is the *only* form present in SAWA's 36,510-accession
    cache, so reading `country` alone found nothing at all.
    """
    for qualifier in ("geo_loc_name", "country"):
        value = metadata.get(qualifier)
        if value:
            return str(value).strip()
    return None


def _normalize_country(value: str | None) -> str | None:
    """GenBank country field is 'Kenya' or 'Kenya: Nairobi'. Map to ISO code.

    Covers every country, not just the East African ones: a record that
    states Thailand must not come out of the ingester looking Kenyan.
    """
    if not value:
        return None
    return country_name_to_code(value.split(":")[0])


# Keys depositors use to declare a subtype. `\b` on both sides so
# "intersubtype" and "Genotyping" are not keys.
_SUBTYPE_KEY = re.compile(r"\b(?:sub-?types?|genotypes?|clades?)\b[\s:=]*", re.IGNORECASE)

# Asserts a mosaic genome without a parseable designation.
_MOSAIC_HINT = re.compile(r"\b(?:recombinant|mosaic|CRF|URF)", re.IGNORECASE)

# A note packs several facts into one string, separated by semicolons; a
# designation is compact, so cap it rather than reading into the prose.
_MAX_DESIGNATION_CHARS = 48


def _normalize_subtype(note: str) -> Subtype:
    """Extract subtype from free-text 'note' qualifier. Best-effort only —
    the LANL curated annotation supersedes this once joined.

    The declaration is isolated first, then classified by
    `subtype_from_designation`. The previous substring match meant
    `"SUBTYPE: C"` matched inside `"SUBTYPE: CRF01_AE"`, so every CRF in
    the corpus was stored as subtype C.
    """
    if not note:
        return Subtype.UNKNOWN
    key = _SUBTYPE_KEY.search(note)
    if key is None:
        return Subtype.RECOMBINANT if _MOSAIC_HINT.search(note) else Subtype.UNKNOWN
    designation = note[key.end() :].split(";", 1)[0][:_MAX_DESIGNATION_CHARS]
    resolved = subtype_from_designation(designation)
    if resolved is Subtype.UNKNOWN and _MOSAIC_HINT.search(note):
        return Subtype.RECOMBINANT
    return resolved


def _infer_gene_region(length: int) -> str:
    """Label the extent the deposit covers, rather than assuming pol.

    GenBank gene annotation is inconsistent across depositors — which is
    why `_build_query` has to OR across gene-name variants — so length is
    the signal that holds for every record. Anything materially longer
    than pol's 3,012 nt is a multi-gene deposit however it is annotated;
    anything shorter is a pol subregion, which is still pol.
    """
    return "genome" if length > MAX_POL_LENGTH else "pol"


def _extract_pubmed_id(record: SeqRecord) -> str | None:
    for ref in record.annotations.get("references", []):
        if getattr(ref, "pubmed_id", None):
            return str(ref.pubmed_id)
    return None
