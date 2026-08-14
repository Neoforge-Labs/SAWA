"""Canonical schema for HIV sequence records.

Every ingestion source (GenBank, LANL, KEMRI) must produce records that
validate against SequenceRecord. This is the contract that keeps downstream
pipelines source-agnostic.

Schema versioning: bump SCHEMA_VERSION on any breaking change. The version
is written into every Parquet output so future analyses can detect drift.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, Field, StringConstraints, field_validator

SCHEMA_VERSION = "0.1.0"


class Subtype(StrEnum):
    """HIV-1 subtypes and common East African CRFs.

    Not exhaustive — extend as needed. Unknown falls through to `UNKNOWN`
    rather than raising, so ingestion never blocks on annotation gaps.
    """

    A1 = "A1"
    A2 = "A2"
    B = "B"
    C = "C"
    D = "D"
    F1 = "F1"
    F2 = "F2"
    G = "G"
    H = "H"
    J = "J"
    K = "K"
    CRF01_AE = "CRF01_AE"
    CRF02_AG = "CRF02_AG"
    CRF10_CD = "CRF10_CD"
    CRF16_A2D = "CRF16_A2D"
    CRF21_A2D = "CRF21_A2D"
    CRF35_AD = "CRF35_AD"
    RECOMBINANT = "recombinant"  # catch-all for other CRFs/URFs
    UNKNOWN = "unknown"


class Source(StrEnum):
    GENBANK = "genbank"
    LANL = "lanl"
    KEMRI = "kemri"


# ISO 3166-1 alpha-2. Not validated against the full list here — the ingester
# normalizes country names to codes and warns on unknowns.
CountryCode = Annotated[str, StringConstraints(min_length=2, max_length=2, to_upper=True)]

# GenBank accession pattern: 1-2 letters + 5-8 digits, optional .version
AccessionStr = Annotated[
    str, StringConstraints(pattern=r"^[A-Z]{1,2}\d{5,8}(\.\d+)?$", to_upper=True)
]


class SequenceRecord(BaseModel):
    """One HIV sequence with the metadata phylodynamics actually needs."""

    # --- identity ---
    accession: AccessionStr
    source: Source
    schema_version: str = SCHEMA_VERSION

    # --- sequence ---
    sequence: str = Field(..., min_length=200)  # anything shorter is noise
    length: int
    gene_region: str = "pol"  # pol only for HIV-DR-Watch v1

    # --- provenance ---
    country: CountryCode | None = None
    collection_date: date | None = None
    collection_year: int | None = Field(None, ge=1980, le=2100)
    subtype: Subtype = Subtype.UNKNOWN
    subtype_source: str | None = None  # "genbank", "lanl_curated", "comet", "scueal"

    # --- clinical context (rarely available in public data) ---
    is_drug_naive: bool | None = None
    risk_group: str | None = None
    patient_id_hashed: str | None = None  # for deduping same-patient sequences

    # --- publication trail ---
    pubmed_id: str | None = None
    study_accession: str | None = None  # e.g. PANGEA cohort ID

    # --- ingestion trail ---
    ingested_at: date
    raw_record_hash: str  # SHA-256 of the source record for change detection

    @field_validator("sequence")
    @classmethod
    def _uppercase_and_check_alphabet(cls, v: str) -> str:
        v = v.upper().replace("U", "T")
        # IUPAC nucleotide codes including ambiguities and gap
        allowed = set("ACGTRYSWKMBDHVN-")
        bad = set(v) - allowed
        if bad:
            raise ValueError(f"sequence contains non-IUPAC characters: {sorted(bad)}")
        return v

    @field_validator("length")
    @classmethod
    def _length_matches_sequence(cls, v: int, info) -> int:
        seq = info.data.get("sequence")
        if seq is not None and len(seq) != v:
            raise ValueError(f"length {v} does not match sequence length {len(seq)}")
        return v
