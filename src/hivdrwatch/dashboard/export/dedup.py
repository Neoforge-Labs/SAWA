"""Deduplication of Sierra drug-resistance results, with structural invariants.

`data/interim/dr_results.parquet` holds repeated copies of the same
resistance call (up to three each). The copies carry no information, so
dropping them is lossless, but only if two things hold: a copy is
identical in every column, and each accession maps to one sequence. Both
are checked here rather than assumed, because every count the dashboard
publishes is taken from the frame this returns.

The converse is not true and is reported, not asserted: several accessions can
carry the same sequence (identical submissions inside one study, mostly), so
`DedupReport.n_distinct_sequences` can be smaller than the accession count.
Every published count is of accessions; this number is what lets the page say
how far that differs from distinct sequences.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import polars as pl

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS: frozenset[str] = frozenset({"accession", "drug", "sequence_sha256"})

_KEY_COLUMNS: list[str] = ["accession", "drug"]


class NonRedundantDuplicateError(ValueError):
    """Rows sharing an (accession, drug) key disagree on some other column."""


class AmbiguousSequenceHashError(ValueError):
    """One accession is associated with more than one sequence hash."""


@dataclass(frozen=True)
class DedupReport:
    """Row counts before and after dropping exact duplicates, and distinct sequences kept."""

    n_rows_in: int
    n_rows_out: int
    n_redundant: int
    n_distinct_sequences: int


def dedup_dr_results(df: pl.DataFrame) -> tuple[pl.DataFrame, DedupReport]:
    """Drop exact duplicate rows and verify that doing so lost nothing.

    Raises `NonRedundantDuplicateError` if rows that share an
    (accession, drug) key survive deduplication, which means the copies
    were not identical and dropping them would have been a choice, not a
    clean-up. Raises `AmbiguousSequenceHashError` if an accession carries
    more than one `sequence_sha256`.
    """
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"dr_results missing required columns: {sorted(missing)}")

    # Stable order keeps two runs over the same input bit-identical.
    deduped = df.unique(maintain_order=True)

    n_keys = deduped.select(_KEY_COLUMNS).n_unique()
    if n_keys != deduped.height:
        n_conflicting = deduped.group_by(_KEY_COLUMNS).len().filter(pl.col("len") > 1).height
        raise NonRedundantDuplicateError(
            f"{n_conflicting} (accession, drug) keys have rows that differ in some "
            f"column ({deduped.height - n_keys} surplus rows); they are not exact "
            "duplicates and cannot be dropped losslessly"
        )

    hashes_per_accession = deduped.group_by("accession").agg(
        pl.col("sequence_sha256").n_unique().alias("n_hashes")
    )
    ambiguous = hashes_per_accession.filter(pl.col("n_hashes") > 1)
    if ambiguous.height:
        raise AmbiguousSequenceHashError(
            f"{ambiguous.height} accessions map to more than one sequence_sha256, "
            f"e.g. {ambiguous['accession'].sort()[0]}"
        )

    report = DedupReport(
        n_rows_in=df.height,
        n_rows_out=deduped.height,
        n_redundant=df.height - deduped.height,
        n_distinct_sequences=deduped["sequence_sha256"].n_unique(),
    )
    logger.info(
        "dr_results dedup: %d rows in, %d out, %d redundant, %d distinct sequences",
        report.n_rows_in,
        report.n_rows_out,
        report.n_redundant,
        report.n_distinct_sequences,
    )
    return deduped, report
