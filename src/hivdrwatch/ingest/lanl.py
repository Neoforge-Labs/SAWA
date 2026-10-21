"""Parse LANL HIV Sequence Database exports.

LANL has no bulk API. The workflow is:
  1. Search on https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html
  2. Download results as "Table" (TSV) — check ALL metadata columns, including
     the Sequence column
  3. Optionally download the matching alignment FASTA
  4. Save under data/raw/lanl/<cc>/ with the LANL default filenames:
         data/raw/lanl/<cc>/results.tsv    (or results.txt — we accept both)
         data/raw/lanl/<cc>/hiv-db.fasta   (optional)
     where <cc> is the ISO 3166-1 alpha-2 country code (lowercased folder name).
  5. Run `LANLIngester(...).load_directory(path)` to get a merged DataFrame,
     or `.to_sequence_records()` for canonical SequenceRecord objects.

Two things about LANL Table exports that trip up naive parsers:
  * Line 1 is a preamble ("Number of records retrieved: N") — the real header
    is on line 2. We skip the first row.
  * The header row starts with a `#` column that is just a 1-based row index.
    It's ignored via the canonical-column allowlist.

The LANL TSV column names shift between releases (they periodically rename
"Sampling Country" ↔ "Country" ↔ "Sample country"). We normalize via a
column alias map instead of hardcoding one spelling.

LANL curation is authoritative for `subtype` — it beats the free-text
GenBank `note` field every time. Downstream code should:

    merged = genbank_df.join(lanl_df.select("accession", "subtype_lanl"),
                             on="accession", how="left")
    merged = merged.with_columns(
        subtype=pl.coalesce(pl.col("subtype_lanl"), pl.col("subtype"))
    )
"""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Iterator
from datetime import date
from pathlib import Path

import polars as pl
from pydantic import ValidationError

from .schema import (
    SCHEMA_VERSION,
    SequenceRecord,
    Source,
    Subtype,
    country_name_to_code,
    subtype_from_designation,
)

logger = logging.getLogger(__name__)


# Normalize LANL's shifting column names to a canonical schema.
# Keys are lowercased/space-stripped for robust matching.
LANL_COLUMN_ALIASES: dict[str, str] = {
    "accession": "accession",
    "genbankaccession": "accession",
    "accessionnumber": "accession",
    "seid": "lanl_se_id",
    "sequenceid": "lanl_se_id",
    "name": "sequence_name",
    "sequencename": "sequence_name",
    "subtype": "subtype_lanl",
    "genotype": "subtype_lanl",
    "country": "country_name",
    "samplingcountry": "country_name",
    "samplecountry": "country_name",
    "samplingyear": "collection_year",
    "sampleyear": "collection_year",
    "yearofsampling": "collection_year",
    "sampledate": "collection_date_raw",
    "collectiondate": "collection_date_raw",
    "drugnaive": "is_drug_naive_raw",
    "patientid": "patient_id_raw",
    "patientcode": "patient_id_raw",
    "riskfactor": "risk_group",
    "pubmedid": "pubmed_id",
    "pubmed": "pubmed_id",
    "sequencelength": "length",
    "seqlength": "length",
    # sequence column variants (LANL Table exports include the full
    # nucleotide sequence when the "Sequence" metadata box is checked).
    "sequence": "sequence",
    "seq": "sequence",
}


CANONICAL_COLUMNS: list[str] = [
    "accession",
    "lanl_se_id",
    "sequence_name",
    # LANL Table export can include the full sequence directly in the TSV.
    # We keep it in the canonical frame so downstream code can use it or join
    # to a separate FASTA-derived sequence table.
    "sequence",
    "subtype_lanl",
    "country_name",
    "country_code",
    "collection_year",
    "collection_date_raw",
    "is_drug_naive",
    "patient_id_raw",
    "risk_group",
    "pubmed_id",
    "length",
    "source_file",
    "schema_version",
]


# Filename patterns we accept per country folder.
# LANL's default "Table" download is results.txt; we also accept results.tsv
# for anyone who renamed the extension after download.
RESULTS_GLOBS = ("results.tsv", "results.txt")
FASTA_GLOB = "hiv-db.fasta"


class LANLIngester:
    """Parse LANL Table exports into a normalized Polars DataFrame.

    Expected layout under raw_dir:
        <raw_dir>/
            ke/
                results.tsv        (or results.txt)
                hiv-db.fasta       (optional; alignment FASTA)
            ug/
                results.txt
                hiv-db.fasta
            tz/
                ...
    """

    def __init__(self, raw_dir: Path) -> None:
        self.raw_dir = raw_dir

    # ------------------------------------------------------------------
    # DataFrame path — for joining against GenBank
    # ------------------------------------------------------------------

    def load_directory(self, *, min_length: int = 800) -> pl.DataFrame:
        """Merge every LANL results file under raw_dir into one deduplicated frame.

        Walks country subfolders (`raw_dir/<cc>/results.{tsv,txt}`) as well as
        any legacy flat files (`raw_dir/*.tsv`). On duplicate accessions the
        *later*-modified file wins, so re-downloading a country next month
        supersedes the old export.
        """
        results_files = self._discover_results_files()
        if not results_files:
            raise FileNotFoundError(
                f"no LANL results.tsv/results.txt files under {self.raw_dir}. "
                "download from https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html "
                "and save as <raw_dir>/<cc>/results.tsv"
            )

        frames = [self._load_one(p) for p in results_files]
        merged = pl.concat(frames, how="diagonal_relaxed")

        # Deduplicate: keep the last occurrence (later file wins by mtime ordering)
        merged = merged.unique(subset=["accession"], keep="last")

        before = merged.height
        if "length" in merged.columns:
            merged = merged.filter(
                pl.col("length").is_null() | (pl.col("length") >= min_length)
            )
        logger.info(
            "LANL: %d records after length>=%d filter (dropped %d)",
            merged.height,
            min_length,
            before - merged.height,
        )
        return merged

    # ------------------------------------------------------------------
    # SequenceRecord path — for LANL-only pipelines (no GenBank)
    # ------------------------------------------------------------------

    def to_sequence_records(self, *, min_length: int = 800) -> Iterator[SequenceRecord]:
        """Yield canonical SequenceRecords built from LANL data alone.

        Uses the TSV `Sequence` column when present; falls back to the
        matching entry in hiv-db.fasta (joined on accession) when the TSV
        omits sequences. Rows without any sequence are dropped with a warning.
        """
        df = self.load_directory(min_length=min_length)

        # Build accession → sequence map from per-country FASTAs when present.
        fasta_map: dict[str, str] = {}
        for fasta_path in self.raw_dir.rglob(FASTA_GLOB):
            fasta_map.update(parse_lanl_fasta(fasta_path))
        if fasta_map:
            logger.info("LANL: loaded %d FASTA sequences for fallback", len(fasta_map))

        today = date.today()
        yielded = 0
        skipped = 0

        for row in df.iter_rows(named=True):
            accession = row.get("accession")
            if not accession:
                skipped += 1
                continue

            seq = row.get("sequence") or fasta_map.get(accession)
            if not seq:
                skipped += 1
                continue
            seq = seq.upper().replace("U", "T")

            try:
                record = SequenceRecord(
                    accession=accession,
                    source=Source.LANL,
                    sequence=seq,
                    length=len(seq),
                    gene_region="pol",
                    country=row.get("country_code"),
                    country_raw=row.get("country_name"),
                    country_source="lanl_curated" if row.get("country_code") else None,
                    collection_year=row.get("collection_year"),
                    subtype=_subtype_enum(row.get("subtype_lanl")),
                    subtype_source="lanl_curated" if row.get("subtype_lanl") else None,
                    is_drug_naive=row.get("is_drug_naive"),
                    risk_group=row.get("risk_group"),
                    pubmed_id=str(row["pubmed_id"]) if row.get("pubmed_id") else None,
                    ingested_at=today,
                    raw_record_hash=hashlib.sha256(seq.encode("utf-8")).hexdigest(),
                )
            except ValidationError as e:
                logger.warning("dropping LANL %s: %s", accession, e)
                skipped += 1
                continue

            yielded += 1
            yield record

        logger.info("LANL→SequenceRecord: %d yielded, %d skipped", yielded, skipped)

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------

    def _discover_results_files(self) -> list[Path]:
        """Find every LANL results.{tsv,txt} under raw_dir, sorted by mtime."""
        found: set[Path] = set()
        # Per-country layout: raw_dir/<cc>/results.{tsv,txt}
        for name in RESULTS_GLOBS:
            found.update(self.raw_dir.rglob(name))
        # Legacy flat layout: raw_dir/*.tsv (kept for backwards compatibility)
        found.update(self.raw_dir.glob("*.tsv"))
        return sorted(found, key=lambda p: p.stat().st_mtime)

    def _load_one(self, path: Path) -> pl.DataFrame:
        logger.info("parsing LANL results: %s", path.relative_to(self.raw_dir))

        skip = _preamble_row_count(path)

        df = pl.read_csv(
            path,
            separator="\t",
            skip_rows=skip,
            infer_schema_length=10_000,
            null_values=["", "NA", "N/A", "unknown", "Unknown", "?", "-"],
            truncate_ragged_lines=True,
            ignore_errors=True,
        )

        df = _rename_columns(df)
        df = self._normalize_types(df, country_hint=_country_from_folder(path, self.raw_dir))
        df = df.with_columns(
            source_file=pl.lit(str(path.relative_to(self.raw_dir))),
            schema_version=pl.lit(SCHEMA_VERSION),
        )

        # Keep only canonical columns that we found; missing ones stay absent.
        # This also drops the leading `#` row-index column from LANL headers.
        present = [c for c in CANONICAL_COLUMNS if c in df.columns]
        df = df.select(present)

        logger.info("  → %d rows, columns: %s", df.height, present)
        return df

    @staticmethod
    def _normalize_types(df: pl.DataFrame, *, country_hint: str | None = None) -> pl.DataFrame:
        expressions = []

        if "subtype_lanl" in df.columns:
            expressions.append(
                pl.col("subtype_lanl")
                .cast(pl.Utf8)
                .map_elements(_coerce_subtype, return_dtype=pl.Utf8)
                .alias("subtype_lanl")
            )

        if "country_name" in df.columns:
            expressions.append(
                pl.col("country_name")
                .cast(pl.Utf8)
                .map_elements(_country_name_to_code, return_dtype=pl.Utf8)
                .alias("country_code")
            )

        if "collection_year" in df.columns:
            expressions.append(pl.col("collection_year").cast(pl.Int32, strict=False))

        if "is_drug_naive_raw" in df.columns:
            expressions.append(
                pl.col("is_drug_naive_raw")
                .cast(pl.Utf8)
                .map_elements(_coerce_bool, return_dtype=pl.Boolean)
                .alias("is_drug_naive")
            )

        if "length" in df.columns:
            expressions.append(pl.col("length").cast(pl.Int32, strict=False))

        if "accession" in df.columns:
            expressions.append(
                pl.col("accession")
                .cast(pl.Utf8)
                .str.strip_chars()
                .str.to_uppercase()
                .alias("accession")
            )

        if "sequence" in df.columns:
            expressions.append(
                pl.col("sequence")
                .cast(pl.Utf8)
                .str.strip_chars()
                .str.to_uppercase()
                .alias("sequence")
            )

        if expressions:
            df = df.with_columns(expressions)

        # Fall back to the folder-derived country if TSV Country was missing/unmapped.
        if country_hint:
            if "country_code" in df.columns:
                df = df.with_columns(
                    pl.col("country_code").fill_null(country_hint).alias("country_code")
                )
            else:
                df = df.with_columns(country_code=pl.lit(country_hint))

        return df


# ---- pure helpers ----


def _rename_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Map LANL's variable column names to our canonical schema."""
    rename_map: dict[str, str] = {}
    for col in df.columns:
        key = _normalize_col_key(col)
        canonical = LANL_COLUMN_ALIASES.get(key)
        if canonical:
            rename_map[col] = canonical
    if rename_map:
        df = df.rename(rename_map)
    return df


def _normalize_col_key(col: str) -> str:
    """Lowercase and strip whitespace/underscore/hyphen/`#` for alias lookup."""
    return col.lower().replace(" ", "").replace("_", "").replace("-", "").replace("#", "")


def _preamble_row_count(path: Path) -> int:
    """Return the number of rows before the real TSV header.

    LANL Table exports begin with a line like `Number of records retrieved: N`.
    We detect that pattern and skip it; otherwise assume no preamble.
    """
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
    if first.lower().startswith("number of records"):
        return 1
    return 0


def _country_from_folder(path: Path, raw_dir: Path) -> str | None:
    """Extract an ISO-alpha-2 country code from the immediate parent folder.

    `raw_dir/ke/results.tsv` → 'KE'. Returns None if the file sits directly in
    raw_dir (legacy flat layout) or the folder name is not a 2-letter code.
    """
    try:
        rel = path.relative_to(raw_dir)
    except ValueError:
        return None
    if len(rel.parts) < 2:
        return None
    candidate = rel.parts[0]
    if len(candidate) != 2 or not candidate.isalpha():
        return None
    return candidate.upper()


def parse_lanl_fasta(path: Path) -> dict[str, str]:
    """Parse a LANL alignment FASTA into an {accession: sequence} map.

    LANL FASTA headers follow the dot-delimited convention
    `subtype.country.year.name.accession`, e.g.
        >A1.KE.1994.Q23_17.AF004885
    We extract the final dot-delimited token as the accession.
    """
    accession_to_seq: dict[str, str] = {}
    current_accession: str | None = None
    buffer: list[str] = []

    def flush() -> None:
        if current_accession and buffer:
            # LANL ships pre-aligned FASTAs (gaps as '-'), but this
            # function is used as a raw-sequence source for
            # `to_sequence_records()`; downstream we re-align to HXB2
            # ourselves. Strip alignment gaps here so the SequenceRecord
            # holds an ungapped nucleotide string.
            seq = "".join(buffer).upper().replace("U", "T").replace("-", "")
            accession_to_seq[current_accession] = seq

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                buffer = []
                header = line[1:].strip()
                # Accession is the last dot-delimited token (LANL convention).
                token = header.rsplit(".", 1)[-1].strip()
                current_accession = token.upper() if token else None
            else:
                buffer.append(line.strip())
        flush()

    return accession_to_seq


def _coerce_subtype(raw: str | None) -> str:
    """LANL Subtype column → canonical subtype value.

    Delegates to the shared classifier so a LANL "A1/D" and a GenBank
    note saying "subtype: A1/D" cannot disagree.
    """
    if raw is None:
        return Subtype.UNKNOWN.value
    token = raw.strip().removeprefix("HIV-1_").removeprefix("HIV1_")
    return subtype_from_designation(token).value


def _subtype_enum(raw: str | None) -> Subtype:
    """String → Subtype enum, defaulting to UNKNOWN."""
    if not raw:
        return Subtype.UNKNOWN
    try:
        return Subtype(raw)
    except ValueError:
        return Subtype.UNKNOWN


def _coerce_bool(raw: str | None) -> bool | None:
    if raw is None:
        return None
    v = raw.strip().lower()
    if v in {"yes", "y", "true", "1", "naive", "drug naive"}:
        return True
    if v in {"no", "n", "false", "0", "experienced", "treated"}:
        return False
    return None


def _country_name_to_code(raw: str | None) -> str | None:
    """LANL Country column → ISO code, via the shared name map."""
    return country_name_to_code(raw)
