"""Patient-level deduplication.

Phylodynamics wants one sequence per patient (or per infection episode);
otherwise longitudinal sequences from a single patient bias tree topology
and cluster size. LANL's `Name` field encodes cohort/patient/visit for
most major studies — we parse a stable `patient_token` from it and hash
to a `patient_group_id`.

Sequences that can't be tied to a patient (Name doesn't match any known
cohort pattern) get a per-accession group_id so they're treated as
independent patients. Better false-independent than false-shared.

We keep *all* sequences and add an `is_earliest_per_patient` boolean;
downstream code chooses whether to drop later timepoints.
"""

from __future__ import annotations

import hashlib
import logging
import re
from collections.abc import Iterable
from dataclasses import dataclass

logger = logging.getLogger(__name__)


# Cohort-name patterns → capture group is the patient identifier.
# Extend as new cohorts appear; each pattern must have exactly one group.
_COHORT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("RV217", re.compile(r"^RV217[_\-]?(\d+)")),
    ("MBAMARCH", re.compile(r"^MBAMARCH[_\-]?(\w+?)(?:[_\-]|$)")),
    ("KER", re.compile(r"^KER(\d+)")),
    ("UG", re.compile(r"^UG(\d+)")),
    ("Q23", re.compile(r"^(Q23)[_\-]")),
]


@dataclass(frozen=True)
class PatientKey:
    country: str | None
    cohort: str
    patient_token: str
    first_year: int | None

    def hashed(self) -> str:
        """16-char sha256 prefix — collision-safe within any realistic cohort."""
        payload = "|".join(
            [
                self.country or "?",
                self.cohort,
                self.patient_token,
                str(self.first_year) if self.first_year else "?",
            ]
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def extract_patient_token(sequence_name: str | None) -> tuple[str, str]:
    """Return (cohort, patient_token) from a LANL Name field.

    Falls back to ('unknown', first-two-underscore-tokens) when nothing
    matches — that treats each stranger-sequence as its own patient.
    """
    if not sequence_name:
        return ("unknown", "no_name")
    name = sequence_name.strip()
    for cohort, pattern in _COHORT_PATTERNS:
        m = pattern.search(name)
        if m:
            return (cohort, m.group(1))
    # Fallback: first two tokens as an opaque identifier.
    tokens = re.split(r"[_\-\s]+", name, maxsplit=2)
    return ("unknown", "_".join(tokens[:2]) or name)


def assign_patient_groups(
    rows: Iterable[dict],
) -> list[dict]:
    """Add `patient_group_id` and `is_earliest_per_patient` to each row.

    Input rows must have `accession`, `sequence_name`, `country`,
    `collection_year`. Order-independent. Ties on year within a patient
    are broken by accession sort order for determinism.
    """
    materialized = list(rows)

    # First pass: compute keys + first year per key.
    keys_by_row: list[PatientKey] = []
    first_year: dict[tuple[str, str, str | None], int] = {}
    for row in materialized:
        cohort, token = extract_patient_token(row.get("sequence_name"))
        year = row.get("collection_year")
        country = row.get("country")
        # first_year is derived from the group, so we build a preliminary key
        # without it, then refine.
        group_key = (country or "?", cohort, token)
        if year is not None and (group_key not in first_year or year < first_year[group_key]):
            first_year[group_key] = year
        keys_by_row.append(
            PatientKey(country=country, cohort=cohort, patient_token=token, first_year=None)
        )

    # Second pass: attach first-year to each key + hash, and mark earliest.
    earliest_seen: dict[str, tuple[int | None, str]] = {}
    for row, key in zip(materialized, keys_by_row, strict=True):
        gk = (key.country or "?", key.cohort, key.patient_token)
        year = first_year.get(gk)
        final_key = PatientKey(
            country=key.country,
            cohort=key.cohort,
            patient_token=key.patient_token,
            first_year=year,
        )
        group_id = final_key.hashed()
        row["patient_group_id"] = group_id

        # Rank within group by (year, accession) for a deterministic winner.
        rank_key = (row.get("collection_year") or 9999, row.get("accession") or "")
        prior = earliest_seen.get(group_id)
        if prior is None or rank_key < prior:
            earliest_seen[group_id] = rank_key

    for row in materialized:
        group_id = row["patient_group_id"]
        rank_key = (row.get("collection_year") or 9999, row.get("accession") or "")
        row["is_earliest_per_patient"] = rank_key == earliest_seen[group_id]

    n_groups = len(earliest_seen)
    logger.info("patient dedup: %d sequences → %d patient groups", len(materialized), n_groups)
    return materialized
