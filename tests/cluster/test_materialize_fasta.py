"""`materialize_aligned_fasta` must not hand the distance engine junk.

The FASTA this task writes is the whole input to cluster construction.
Before this was fixed it wrote every row whose sequence string was
truthy, and a row of 3,012 dashes is truthy, so all 21,358 all-gap rows
in the shipped curated vintage went into it
(`docs/evaluations/align-coverage.md` §6). `phylo/build.py` has gated
its tree input on `passes_qc` and `MIN_ALIGNED_COVERAGE` from the
start; these tests pin the cluster path to the same two gates.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from prefect.logging import disable_run_logger

from hivdrwatch.phylo.build import MIN_ALIGNED_COVERAGE
from hivdrwatch.pipeline.cluster import eligible_accessions, materialize_aligned_fasta

ALIGNMENT_WIDTH = 3012


def _row(
    accession: str,
    seq: str,
    *,
    passes_qc: bool = True,
) -> dict[str, object]:
    return {"accession": accession, "aligned_pol": seq, "passes_qc": passes_qc}


def _covered(fraction: float) -> str:
    """An aligned_pol string with `fraction` of its columns called."""
    n = int(ALIGNMENT_WIDTH * fraction)
    return ("ACGT" * ((n // 4) + 1))[:n] + "-" * (ALIGNMENT_WIDTH - n)


def _curated(rows: list[dict[str, object]]) -> pl.DataFrame:
    return pl.DataFrame(
        rows,
        schema={"accession": pl.String, "aligned_pol": pl.String, "passes_qc": pl.Boolean},
    )


def _materialize(*args: object, **kwargs: object) -> Path:
    """Call the task body outside a flow run, with Prefect's logger stubbed."""
    with disable_run_logger():
        return materialize_aligned_fasta.fn(*args, **kwargs)


def _ids(fasta: Path) -> list[str]:
    return [
        line[1:].strip() for line in fasta.read_text().splitlines() if line.startswith(">")
    ]


class TestEligibility:
    def test_all_gap_rows_are_excluded(self) -> None:
        curated = _curated(
            [
                _row("GAP_1", "-" * ALIGNMENT_WIDTH),
                _row("GAP_2", "-" * ALIGNMENT_WIDTH),
                _row("REAL_1", _covered(0.9)),
            ]
        )
        assert eligible_accessions(curated) == ["REAL_1"]

    def test_all_gap_rows_are_excluded_even_when_qc_passed(self) -> None:
        # The shipped vintage's all-gap rows all carried passes_qc=True,
        # so the coverage gate has to stand on its own.
        curated = _curated([_row("GAP_1", "-" * ALIGNMENT_WIDTH, passes_qc=True)])
        assert eligible_accessions(curated) == []

    def test_qc_failed_rows_are_excluded(self) -> None:
        curated = _curated(
            [
                _row("FAILED", _covered(0.9), passes_qc=False),
                _row("PASSED", _covered(0.9), passes_qc=True),
            ]
        )
        assert eligible_accessions(curated) == ["PASSED"]

    def test_qc_gate_can_be_lifted_for_exploration(self) -> None:
        curated = _curated([_row("FAILED", _covered(0.9), passes_qc=False)])
        assert eligible_accessions(curated, require_qc_pass=False) == ["FAILED"]

    def test_coverage_gate_matches_phylo(self) -> None:
        curated = _curated(
            [
                _row("THIN", _covered(MIN_ALIGNED_COVERAGE - 0.1)),
                _row("THICK", _covered(MIN_ALIGNED_COVERAGE + 0.1)),
            ]
        )
        assert eligible_accessions(curated) == ["THICK"]

    def test_coverage_gate_is_a_parameter(self) -> None:
        curated = _curated([_row("THIN", _covered(0.2))])
        assert eligible_accessions(curated, min_coverage=0.1) == ["THIN"]
        assert eligible_accessions(curated, min_coverage=0.3) == []

    def test_duplicate_accessions_collapse_to_one_record(self) -> None:
        # The curated table stores the same record once per country
        # query, and a FASTA with repeated ids makes the distance CSV
        # ambiguous about which record a pair refers to.
        seq = _covered(0.9)
        curated = _curated([_row("DUP", seq), _row("DUP", seq), _row("OTHER", seq)])
        assert eligible_accessions(curated) == ["DUP", "OTHER"]

    def test_missing_sequence_column_raises(self) -> None:
        curated = pl.DataFrame({"accession": ["A"]})
        with pytest.raises(ValueError, match="aligned_pol"):
            eligible_accessions(curated)


class TestMaterialize:
    def test_writes_only_eligible_records(self, tmp_path: Path) -> None:
        curated = _curated(
            [
                _row("GAP", "-" * ALIGNMENT_WIDTH),
                _row("FAILED", _covered(0.9), passes_qc=False),
                _row("THIN", _covered(0.1)),
                _row("GOOD", _covered(0.9)),
            ]
        )
        parquet = tmp_path / "curated.parquet"
        curated.write_parquet(parquet)

        out = _materialize(tmp_path / "absent.fasta", parquet, tmp_path / "work")
        assert _ids(out) == ["GOOD"]

    def test_existing_fasta_is_filtered_not_passed_through(self, tmp_path: Path) -> None:
        # A FASTA the curate stage already wrote carries no QC or
        # coverage information, so the curated Parquet still decides
        # which of its records may be used.
        curated = _curated(
            [
                _row("GAP", "-" * ALIGNMENT_WIDTH),
                _row("GOOD", _covered(0.9)),
            ]
        )
        parquet = tmp_path / "curated.parquet"
        curated.write_parquet(parquet)

        source = tmp_path / "curate" / "aligned.fasta"
        source.parent.mkdir(parents=True)
        source.write_text(
            f">GAP\n{'-' * ALIGNMENT_WIDTH}\n>GOOD\n{_covered(0.9)}\n"
            f">NOT_IN_CURATED\n{_covered(0.9)}\n"
        )

        out = _materialize(source, parquet, tmp_path / "work")
        assert _ids(out) == ["GOOD"]
        assert out != source  # the source FASTA is left alone

    def test_raises_when_nothing_is_eligible(self, tmp_path: Path) -> None:
        curated = _curated([_row("GAP", "-" * ALIGNMENT_WIDTH)])
        parquet = tmp_path / "curated.parquet"
        curated.write_parquet(parquet)
        with pytest.raises(ValueError, match="no eligible records"):
            _materialize(tmp_path / "absent.fasta", parquet, tmp_path / "work")
