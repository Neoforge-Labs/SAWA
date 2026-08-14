"""Canonical schema for the curated silver layer.

The curate stage produces:
    data/interim/curated/aligned.parquet   — one row per input sequence with
                                              gene-region slices + QC flags
    data/interim/curated/patients.parquet  — patient-group mapping (dedup key)

Every row carries `curate_schema_version` so downstream features/models
can pin against a specific vintage.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum

from pydantic import BaseModel, Field

from ..ingest.schema import AccessionStr

CURATE_SCHEMA_VERSION = "0.1.0"


class QCFlag(StrEnum):
    """Reasons a sequence might be flagged (multiple can co-occur)."""

    HIGH_AMBIGUITY = "high_ambiguity"       # >5% non-ACGT
    APOBEC_HYPERMUT = "apobec_hypermut"     # G→A hypermutation signature
    STOP_CODON = "stop_codon"                # premature stop in PR/RT/IN
    SHORT_COVERAGE = "short_coverage"        # region missing >20% of sites
    LOW_COMPLEXITY = "low_complexity"        # >80% single nt


class CuratedSequence(BaseModel):
    """One sequence after alignment, trim, and QC — ready for feature building."""

    # --- identity ---
    accession: AccessionStr
    sequence_sha256: str = Field(..., min_length=64, max_length=64)

    # --- alignment ---
    aligned_pol: str        # gap-preserving alignment to HXB2 pol (3012 nt)
    aligned_length: int     # should equal HXB2_POL_LENGTH after --keeplength

    # --- gene-region slices (aligned, gap-preserving; empty string if uncovered) ---
    pr_nt: str
    rt_nt: str
    in_nt: str

    # --- QC ---
    ambiguity_fraction: float = Field(..., ge=0.0, le=1.0)
    apobec_signature_ratio: float | None = None  # None when no G→A calls possible
    has_pr_stop: bool = False
    has_rt_stop: bool = False
    has_in_stop: bool = False
    qc_flags: list[QCFlag] = Field(default_factory=list)
    passes_qc: bool

    # --- dedup ---
    patient_group_id: str            # sha256(country + patient_token + first_year)[:16]
    is_earliest_per_patient: bool

    # --- provenance ---
    country: str | None = None
    collection_year: int | None = None
    subtype: str | None = None
    source_file: str | None = None
    curated_at: date
    schema_version: str = CURATE_SCHEMA_VERSION
