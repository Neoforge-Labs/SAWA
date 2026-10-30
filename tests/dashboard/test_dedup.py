"""Tests for the Sierra-result deduplication step of the dashboard export.

The dedup layer is the guard that stops inflated counts reaching the
published page, so these tests pin both the lossless-drop behaviour and
the two ways the input can violate the assumptions that make it lossless.
"""

from __future__ import annotations

from datetime import date
from typing import Any

import polars as pl
import pytest

from hivdrwatch.dashboard.export.dedup import (
    AmbiguousSequenceHashError,
    DedupReport,
    NonRedundantDuplicateError,
    dedup_dr_results,
)


def _row(
    accession: str,
    drug: str,
    *,
    level: str = "susceptible",
    sha: str = "a" * 64,
    mutations: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "accession": accession,
        "sequence_sha256": sha,
        "algorithm_name": "HIVDB",
        "algorithm_version": "9.6",
        "algorithm_publish_date": date(2024, 1, 1),
        "drug_class": "NRTI",
        "drug": drug,
        "drug_display": drug.lower(),
        "score": 0.0,
        "level": level,
        "gene": "RT",
        "mutations": mutations if mutations is not None else ["M184V"],
        "called_at": date(2024, 6, 1),
        "schema_version": "0.1.0",
    }


def _synthetic_dr_frame(rows: list[dict[str, Any]]) -> pl.DataFrame:
    return pl.from_dicts(
        rows,
        schema_overrides={"mutations": pl.List(pl.String)},
    )


def test_dedup_removes_exact_duplicate_rows() -> None:
    df = _synthetic_dr_frame(rows=[_row("AB1", "3TC"), _row("AB1", "3TC"), _row("AB1", "AZT")])
    out, report = dedup_dr_results(df)
    assert out.height == 2
    assert report == DedupReport(n_rows_in=3, n_rows_out=2, n_redundant=1)


def test_dedup_raises_when_duplicates_carry_different_values() -> None:
    rows = [_row("AB1", "3TC", level="high"), _row("AB1", "3TC", level="susceptible")]
    with pytest.raises(NonRedundantDuplicateError):
        dedup_dr_results(_synthetic_dr_frame(rows=rows))


def test_dedup_raises_when_accession_has_two_sequence_hashes() -> None:
    rows = [_row("AB1", "3TC", sha="a" * 64), _row("AB1", "AZT", sha="b" * 64)]
    with pytest.raises(AmbiguousSequenceHashError):
        dedup_dr_results(_synthetic_dr_frame(rows=rows))


def test_dedup_is_idempotent() -> None:
    df = _synthetic_dr_frame(rows=[_row("AB1", "3TC"), _row("AB1", "3TC")])
    once, _ = dedup_dr_results(df)
    twice, report = dedup_dr_results(once)
    assert twice.equals(once)
    assert report.n_redundant == 0


def test_dedup_keeps_a_frame_with_no_duplicates_unchanged() -> None:
    df = _synthetic_dr_frame(rows=[_row("AB1", "3TC"), _row("AB2", "3TC", mutations=[])])
    out, report = dedup_dr_results(df)
    assert out.equals(df)
    assert report == DedupReport(n_rows_in=2, n_rows_out=2, n_redundant=0)


def test_dedup_treats_rows_differing_only_in_mutations_as_conflicting() -> None:
    rows = [_row("AB1", "3TC", mutations=["M184V"]), _row("AB1", "3TC", mutations=["K65R"])]
    with pytest.raises(NonRedundantDuplicateError):
        dedup_dr_results(_synthetic_dr_frame(rows=rows))


def test_dedup_error_names_how_many_rows_conflict() -> None:
    rows = [
        _row("AB1", "3TC", level="high"),
        _row("AB1", "3TC", level="low"),
        _row("AB2", "AZT"),
    ]
    with pytest.raises(NonRedundantDuplicateError, match=r"^1 \(accession, drug\) keys"):
        dedup_dr_results(_synthetic_dr_frame(rows=rows))


def test_dedup_names_missing_columns() -> None:
    df = _synthetic_dr_frame(rows=[_row("AB1", "3TC")]).drop("sequence_sha256")
    with pytest.raises(ValueError, match="sequence_sha256"):
        dedup_dr_results(df)
