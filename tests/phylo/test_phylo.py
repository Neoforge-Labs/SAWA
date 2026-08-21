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
    EXCLUDED_SUBTYPES,
    MIN_ALIGNED_COVERAGE,
    MIN_TAXA_FOR_TREE,
    SubstitutionModel,
    TreeKind,
    TreeMetadata,
    build_partition_report,
    build_tip_dates_csv,
    partition_by_subtype,
    write_partition_report,
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
    _extract_clock_rate_from_file,
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

    def test_default_excluded_subtypes_covers_unknown_and_recombinant(self) -> None:
        # Regression: 'unknown' has always been excluded; 'recombinant'
        # was added because it lumps distinct CRFs into a single label
        # and IQ-TREE fails on the mixed alignment. Both should be
        # rejected by default without an explicit `excluded=` argument.
        assert "unknown" in EXCLUDED_SUBTYPES
        assert "recombinant" in EXCLUDED_SUBTYPES

        df = self._frame({"A1": 30, "unknown": 30, "recombinant": 30})
        parts = partition_by_subtype(df, min_taxa=20)
        assert set(parts) == {"A1"}


# ---------- build_partition_report ----------


class TestPartitionReport:
    """Audit-trail parquet: one row per candidate subtype label.

    The report must be deterministic given the same curated frame and
    thresholds, and must reflect every non-null subtype label — including
    those excluded by policy and those below thresholds — so that a
    reviewer can reconstruct manifest membership without reading logs.
    """

    def _frame(self, spec: dict[str, tuple[int, str]]) -> pl.DataFrame:
        """spec: {subtype: (count, sequence_pattern)} — pattern controls coverage."""
        rows = []
        for st, (n, seq) in spec.items():
            for i in range(n):
                rows.append(
                    {
                        "accession": f"{st}{i:03d}",
                        "subtype": st,
                        "passes_qc": True,
                        "collection_year": 2020,
                        "aligned_pol": seq,
                    }
                )
        return pl.from_dicts(rows)

    def test_one_row_per_candidate_label(self) -> None:
        df = self._frame(
            {
                "A1": (30, "A" * 3012),
                "C": (5, "A" * 3012),
                "unknown": (100, "A" * 3012),
            }
        )
        report = build_partition_report(df)
        assert sorted(r.subtype for r in report) == ["A1", "C", "unknown"]

    def test_decision_built_for_healthy_subtype(self) -> None:
        df = self._frame({"A1": (30, "A" * 3012)})
        report = build_partition_report(df)
        row = next(r for r in report if r.subtype == "A1")
        assert row.decision == "built"
        assert row.n_available == 30
        assert row.n_after_qc == 30
        assert row.n_after_coverage == 30

    def test_decision_excluded_label(self) -> None:
        df = self._frame({"unknown": (100, "A" * 3012), "recombinant": (50, "A" * 3012)})
        report = build_partition_report(df)
        for name in ("unknown", "recombinant"):
            row = next(r for r in report if r.subtype == name)
            assert row.decision == "excluded_label"
            assert "EXCLUDED_SUBTYPES" in row.decision_reason

    def test_decision_below_min_taxa(self) -> None:
        df = self._frame({"A1": (5, "A" * 3012)})
        report = build_partition_report(df, min_taxa=20)
        row = next(r for r in report if r.subtype == "A1")
        assert row.decision == "below_min_taxa"
        assert row.min_taxa_used == 20

    def test_decision_below_min_coverage(self) -> None:
        # 30 sequences pass QC but every one is 90% gap → all fail coverage.
        low_cov = ("-" * 2711) + ("A" * 301)  # ~10% non-gap
        df = self._frame({"A1": (30, low_cov)})
        report = build_partition_report(df, min_taxa=20, min_coverage=0.5)
        row = next(r for r in report if r.subtype == "A1")
        assert row.decision == "below_min_coverage"
        assert row.n_after_qc == 30
        assert row.n_after_coverage == 0

    def test_deterministic_order(self) -> None:
        # Order stability matters for downstream diffs / manifest audits.
        df = self._frame(
            {"D": (30, "A" * 3012), "A1": (30, "A" * 3012), "C": (30, "A" * 3012)}
        )
        report = build_partition_report(df)
        assert [r.subtype for r in report] == ["A1", "C", "D"]

    def test_thresholds_recorded(self) -> None:
        df = self._frame({"A1": (30, "A" * 3012)})
        report = build_partition_report(df, min_taxa=25, min_coverage=0.7)
        row = report[0]
        assert row.min_taxa_used == 25
        assert row.min_coverage_used == 0.7

    def test_write_partition_report_roundtrip(self, tmp_path: Path) -> None:
        df = self._frame(
            {"A1": (30, "A" * 3012), "unknown": (10, "A" * 3012)}
        )
        report = build_partition_report(df)
        out = tmp_path / "partition_report.parquet"
        write_partition_report(report, out)
        assert out.exists()
        loaded = pl.read_parquet(out)
        assert set(loaded["subtype"].to_list()) == {"A1", "unknown"}
        assert set(loaded["decision"].to_list()) == {"built", "excluded_label"}

    def test_default_coverage_matches_module_constant(self) -> None:
        df = self._frame({"A1": (30, "A" * 3012)})
        report = build_partition_report(df)
        assert report[0].min_coverage_used == MIN_ALIGNED_COVERAGE


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

    def test_clock_rate_from_file_positive(self, tmp_path: Path) -> None:
        # Real TreeTime `molecular_clock.txt` format seen on the D subtype
        # run at v0.11.4 — tab-indented, r² on a separate line.
        (tmp_path / "molecular_clock.txt").write_text(
            "Root-Tip-Regression:\n --rate:\t2.303e-03\n --r^2:  \t0.09\n"
        )
        rate, sd = _extract_clock_rate_from_file(tmp_path)
        assert rate == pytest.approx(2.303e-3)
        assert sd is None  # File doesn't emit SD; caller must not fabricate.

    def test_clock_rate_from_file_negative(self, tmp_path: Path) -> None:
        # Real observation on the C subtype (n=20, poor temporal signal).
        # Negative rate must be preserved so downstream can flag the fit.
        (tmp_path / "molecular_clock.txt").write_text(
            "Root-Tip-Regression:\n --rate:\t-1.103e-03\n --r^2:  \t0.001\n"
        )
        rate, sd = _extract_clock_rate_from_file(tmp_path)
        assert rate == pytest.approx(-1.103e-3)
        assert sd is None

    def test_clock_rate_from_file_missing(self, tmp_path: Path) -> None:
        rate, sd = _extract_clock_rate_from_file(tmp_path)
        assert rate is None and sd is None

    def test_clock_rate_from_file_malformed(self, tmp_path: Path) -> None:
        (tmp_path / "molecular_clock.txt").write_text("gibberish, no rate line")
        rate, sd = _extract_clock_rate_from_file(tmp_path)
        assert rate is None and sd is None


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
