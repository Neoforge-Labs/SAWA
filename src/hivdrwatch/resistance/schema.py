"""Canonical schema for Stanford HIVdb Sierra drug-resistance results.

One `DRResult` per (accession, drug) pair. Kept flat and long so it slots
into Polars group-bys and joins without pivoting gymnastics.

Every row pins the Sierra `algorithm_version` at the time of the call so a
Stanford algorithm bump doesn't silently change your training set. When
re-training a model, filter to a specific algorithm_version to keep labels
consistent.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum

from pydantic import BaseModel, Field

from ..ingest.schema import AccessionStr

DR_SCHEMA_VERSION = "0.1.0"


class DrugClass(StrEnum):
    """HIV drug classes we currently score."""

    PI = "PI"       # protease inhibitor
    NRTI = "NRTI"   # nucleoside RT inhibitor
    NNRTI = "NNRTI"  # non-nucleoside RT inhibitor
    INSTI = "INSTI"  # integrase strand-transfer inhibitor


class ResistanceLevel(StrEnum):
    """Stanford's 5-level interpretation scale."""

    SUSCEPTIBLE = "susceptible"
    POTENTIAL_LOW = "potential_low"      # potential low-level resistance
    LOW = "low"                          # low-level
    INTERMEDIATE = "intermediate"
    HIGH = "high"

    @classmethod
    def from_sierra_text(cls, text: str) -> "ResistanceLevel":
        """Sierra returns strings like 'Susceptible', 'High-Level Resistance'."""
        norm = text.strip().lower()
        return _SIERRA_LEVEL_MAP.get(norm, cls.SUSCEPTIBLE)


_SIERRA_LEVEL_MAP: dict[str, ResistanceLevel] = {
    "susceptible": ResistanceLevel.SUSCEPTIBLE,
    "potential low-level resistance": ResistanceLevel.POTENTIAL_LOW,
    "low-level resistance": ResistanceLevel.LOW,
    "intermediate resistance": ResistanceLevel.INTERMEDIATE,
    "high-level resistance": ResistanceLevel.HIGH,
}


class DRResult(BaseModel):
    """One resistance call: sequence × drug × algorithm version."""

    # --- identity ---
    accession: AccessionStr
    sequence_sha256: str = Field(..., min_length=64, max_length=64)

    # --- algorithm provenance ---
    algorithm_name: str = "HIVDB"
    algorithm_version: str  # e.g. "9.6"
    algorithm_publish_date: date | None = None

    # --- drug scoring ---
    drug_class: DrugClass
    drug: str                    # short name, e.g. "3TC", "DTG"
    drug_display: str | None = None  # long name, e.g. "lamivudine"
    score: float                 # Stanford's numeric score (higher = more resistant)
    level: ResistanceLevel

    # --- mutation context ---
    gene: str                    # "PR", "RT", or "IN"
    mutations: list[str] = Field(default_factory=list)  # e.g. ["M184V", "K103N"]

    # --- ingestion trail ---
    called_at: date
    schema_version: str = DR_SCHEMA_VERSION


class SubtypeCall(BaseModel):
    """Sierra's per-sequence subtype call, kept separate from DRResult
    because it's per-sequence, not per-drug."""

    accession: AccessionStr
    sequence_sha256: str = Field(..., min_length=64, max_length=64)
    subtype: str                # Sierra's free-form subtype text
    algorithm_version: str
    called_at: date
    schema_version: str = DR_SCHEMA_VERSION
