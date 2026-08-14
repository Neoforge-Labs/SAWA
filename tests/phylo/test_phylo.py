"""Unit tests for phylo helpers.

Docker-invoked functions (`run_iqtree`, `run_treetime`) are not tested
here — that requires the built image + a real alignment. Those live in
`tests/phylo/test_integration.py` (fixture-based, run in CI only).

This file tests:
  * Subtype partitioning + min_taxa gate
  * Alignment FASTA writer
  * Tip-date CSV builder
  * IQ-TREE log-likelihood + version parsers
  * TreeTime clock-rate + version parsers
  * TreeMetadata schema validation
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import polars as pl
import pytest

from hivdrwatch.phylo import (
    MIN_TAXA_FOR_TREE,
    SubstitutionModel,
    TreeKind,
    TreeMetadata,
    build_tip_dates_csv,
    partition_by_subtype,
    write_subtype_alignment,
)
from hivdrwatch.phylo.iqtree import (
    _extract_log_likelihood,
)
from hivdrwatch.phylo.iqtree import (
    _extract_version as _iqtree_version,
)
from hivdrwatch.phylo.treetime import (
    _extract_clock_rate,
)
from hivdrwatch.phylo.treetime import (
    _extract_version as _treetime_version,
)

# ---------- partition_by_subtype ----------


class TestPartitionBySubtype:
    def _frame(self, n_per_subtype: dict[str, int]) -> pl.DataFrame:
        rows = []
        for st, n in n_per_subtype.items():
            for i in range(n):
                rows.append(
                    {
                        "accession": f"{st}{i:03d}",
                        "subtype": st,
                        "passes_qc": True,
                        "collection_year": 2020,
                        "aligned_pol": "A" * 3012,
                    }
                )
        return pl.from_dicts(rows)

    def test_drops_below_min_taxa(self) -> None:
        df = self._frame({"A1": 30, "C": 15, "D": 40})
        parts = partition_by_subtype(df, min_taxa=20)
        assert set(parts) == {"A1", "D"}
        assert parts["A1"].height == 30
        assert parts["D"].height == 40

    def test_default_min_matches_module_constant(self) -> None:
        df = self._frame({"A1": MIN_TAXA_FOR_TREE, "C": MIN_TAXA_FOR_TREE - 1})
        parts = partition_by_subtype(df)
        assert "A1" in parts and "C" not in parts

    def test_qc_gate_applied(self) -> None:
        df = pl.from_dicts(
            [
                {"accession": f"A{i}", "subtype": "A1", "passes_qc": i % 2 == 0}
                for i in range(60)
            ]
        )
        parts = partition_by_subtype(df, min_taxa=20)
        # Only 30 pass QC — still ≥ 20 so should survive
        assert parts["A1"].height == 30

    def test_drops_null_subtype(self) -> None:
        rows = [
            {"accession": f"A{i}", "subtype": "A1", "passes_qc": True} for i in range(30)
        ] + [{"accession": f"X{i}", "subtype": None, "passes_qc": True} for i in range(30)]
        df = pl.from_dicts(rows)
        parts = partition_by_subtype(df, min_taxa=20)
        assert set(parts) == {"A1"}


# ---------- write_subtype_alignment ----------


class TestWriteSubtypeAlignment:
    def test_writes_fasta(self, tmp_path: Path) -> None:
        df = pl.DataFrame(
            {
                "accession": ["A1", "A2", "A3"],
                "aligned_pol": ["ACGT" * 3, "GGGG" * 3, "TTTT" * 3],
            }
        )
        out = tmp_path / "A1.fasta"
        write_subtype_alignment(df, out)
        text = out.read_text()
        assert text.count(">") == 3
        assert ">A1\n" in text
        assert "ACGTACGTACGT" in text

    def test_falls_back_to_sequence_column(self, tmp_path: Path) -> None:
        df = pl.DataFrame({"accession": ["A1"], "sequence": ["ACGT" * 5]})
        out = tmp_path / "x.fasta"
        write_subtype_alignment(df, out, seq_col="aligned_pol")
        assert ">A1\n" in out.read_text()

    def test_skips_null_sequences(self, tmp_path: Path) -> None:
        df = pl.DataFrame({"accession": ["A", "B"], "aligned_pol": ["ACGT" * 3, None]})
        out = tmp_path / "x.fasta"
        write_subtype_alignment(df, out)
        assert out.read_text().count(">") == 1


# ---------- build_tip_dates_csv ----------


class TestBuildTipDatesCsv:
    def test_emits_decimal_years(self, tmp_path: Path) -> None:
        df = pl.DataFrame(
            {
                "accession": ["A1", "A2", "A3"],
                "collection_year": [2018, 2019, 2020],
            }
        )
        out = tmp_path / "dates.csv"
        build_tip_dates_csv(df, ["A1", "A2", "A3"], out)
        text = out.read_text()
        lines = text.strip().split("\n")
        assert lines[0] == "name,date"
        assert lines[1].startswith("A1,2018.5")

    def test_filters_to_requested_accessions(self, tmp_path: Path) -> None:
        df = pl.DataFrame(
            {
                "accession": ["A1", "A2", "OTHER"],
                "collection_year": [2018, 2019, 2020],
            }
        )
        out = tmp_path / "dates.csv"
        build_tip_dates_csv(df, ["A1", "A2"], out)
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 rows
        assert "OTHER" not in out.read_text()

    def test_drops_null_year_rows(self, tmp_path: Path) -> None:
        df = pl.DataFrame({"accession": ["A", "B"], "collection_year": [2020, None]})
        out = tmp_path / "dates.csv"
        build_tip_dates_csv(df, ["A", "B"], out)
        assert out.read_text().count("\n") == 2  # header + 1 row + trailing

    def test_missing_columns_raises(self, tmp_path: Path) -> None:
        df = pl.DataFrame({"accession": ["A"]})
        with pytest.raises(ValueError, match="missing columns"):
            build_tip_dates_csv(df, ["A"], tmp_path / "x.csv")


# ---------- IQ-TREE log parsers ----------


class TestIqtreeParsers:
    def test_log_likelihood_extraction(self, tmp_path: Path) -> None:
        log = tmp_path / "run.log"
        log.write_text("Some preamble\nBEST SCORE FOUND : -12345.6789\nMore text\n")
        assert _extract_log_likelihood(log) == pytest.approx(-12345.6789)

    def test_log_likelihood_missing_returns_none(self, tmp_path: Path) -> None:
        log = tmp_path / "run.log"
        log.write_text("no relevant line here\n")
        assert _extract_log_likelihood(log) is None

    def test_missing_file_returns_none(self, tmp_path: Path) -> None:
        assert _extract_log_likelihood(tmp_path / "nope.log") is None

    def test_version_extraction(self, tmp_path: Path) -> None:
        log = tmp_path / "run.log"
        log.write_text("IQ-TREE multicore version 2.3.6 for Linux x86-64\n")
        assert "2.3.6" in _iqtree_version(log)


# ---------- TreeTime parsers ----------


class TestTreetimeParsers:
    def test_clock_rate_extraction(self) -> None:
        stderr = "Some log\n   clock rate: 1.234e-03 +/- 4.56e-04\nmore\n"
        rate, sd = _extract_clock_rate(stderr)
        assert rate == pytest.approx(1.234e-3)
        assert sd == pytest.approx(4.56e-4)

    def test_clock_rate_missing(self) -> None:
        assert _extract_clock_rate("no clock here") == (None, None)

    def test_version_extraction(self) -> None:
        assert "0.11.4" in _treetime_version("TreeTime 0.11.4 running...")


# ---------- TreeMetadata schema ----------


class TestTreeMetadata:
    def test_ml_row_valid(self) -> None:
        row = TreeMetadata(
            subtype="A1",
            kind=TreeKind.ML,
            n_taxa=50,
            newick_path="data/interim/trees/A1/A1.treefile",
            substitution_model=SubstitutionModel.GTR_F_I_G4,
            log_likelihood=-12345.67,
            random_seed=42,
            tool_version="iqtree2 2.3.6",
            wall_seconds=180.0,
            built_at=date(2026, 8, 13),
        )
        assert row.kind == TreeKind.ML
        assert row.substitution_model == SubstitutionModel.GTR_F_I_G4

    def test_time_scaled_row_valid(self) -> None:
        row = TreeMetadata(
            subtype="A1",
            kind=TreeKind.TIME_SCALED,
            n_taxa=50,
            newick_path="data/interim/trees/A1/timed/timetree.nexus",
            dates_csv_path="data/interim/trees/A1/dates.csv",
            clock_rate=1.5e-3,
            clock_rate_stddev=2e-4,
            time_span_years=25.0,
            random_seed=42,
            tool_version="treetime 0.11.4",
            wall_seconds=60.0,
            built_at=date(2026, 8, 13),
        )
        assert row.clock_rate == 1.5e-3

    def test_n_taxa_positive(self) -> None:
        with pytest.raises(ValueError):
            TreeMetadata(
                subtype="A1",
                kind=TreeKind.ML,
                n_taxa=0,
                newick_path="x",
                random_seed=42,
                tool_version="x",
                built_at=date(2026, 8, 13),
            )
