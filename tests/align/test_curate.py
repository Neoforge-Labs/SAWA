"""Unit tests for the curate stage — pure helpers only.

MAFFT invocation is not tested here (requires the docker image + running
daemon). That belongs in an integration test with a fixture.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hivdrwatch.align.dedup import (
    PatientKey,
    assign_patient_groups,
    extract_patient_token,
)
from hivdrwatch.align.qc import (
    ambiguity_fraction,
    apobec_signature_ratio,
    has_premature_stop,
    passes_qc,
    qc_flags,
    translate_in_frame,
)
from hivdrwatch.align.reference import (
    HXB2_POL_LENGTH,
    extract_hxb2_from_lanl_fasta,
)
from hivdrwatch.align.schema import QCFlag
from hivdrwatch.align.trim import coverage_fraction, slice_region

# ---------- reference extraction ----------


class TestExtractHXB2FromLANLFasta:
    def test_extracts_first_hxb2_record(self, tmp_path: Path) -> None:
        fasta = tmp_path / "hiv-db.fasta"
        fasta.write_text(
            ">B.FR.1983.HXB2-LAI-IIIB-BRU.K03455\n"
            "ACGT-ACGT\n"
            "AC-GT\n"
            ">A1.KE.1994.Q23.AF004885\n"
            "TTTAAA\n"
        )
        ref = extract_hxb2_from_lanl_fasta(fasta)
        assert ref.sequence == "ACGTACGTACGT"
        assert ref.length == 12
        assert ref.source_path == fasta

    def test_raises_when_no_hxb2(self, tmp_path: Path) -> None:
        fasta = tmp_path / "hiv-db.fasta"
        fasta.write_text(">A1.KE.1994.Q23.AF004885\nACGT\n")
        with pytest.raises(ValueError, match="no HXB2"):
            extract_hxb2_from_lanl_fasta(fasta)


# ---------- trim ----------


class TestSliceRegion:
    def test_slices_by_pol_relative_coordinates(self) -> None:
        seq = "N" * HXB2_POL_LENGTH
        pr = slice_region(seq, "PR")
        rt = slice_region(seq, "RT")
        integrase = slice_region(seq, "IN")
        assert len(pr) == 465 - 169 + 1
        assert len(rt) == 1785 - 466 + 1
        assert len(integrase) == 3012 - 2146 + 1

    def test_raises_on_wrong_length(self) -> None:
        with pytest.raises(ValueError, match="keeplength"):
            slice_region("ACGT", "PR")

    def test_raises_on_unknown_region(self) -> None:
        with pytest.raises(ValueError, match="unknown region"):
            slice_region("N" * HXB2_POL_LENGTH, "vif")


class TestCoverageFraction:
    def test_all_gap_is_zero(self) -> None:
        assert coverage_fraction("---") == 0.0

    def test_full_coverage_is_one(self) -> None:
        assert coverage_fraction("ACGT") == 1.0

    def test_partial(self) -> None:
        assert coverage_fraction("AC--GT") == pytest.approx(4 / 6)

    def test_empty_is_zero(self) -> None:
        assert coverage_fraction("") == 0.0


# ---------- QC ----------


class TestAmbiguityFraction:
    def test_ignores_gaps(self) -> None:
        assert ambiguity_fraction("ACGT---") == 0.0

    def test_counts_ambiguous(self) -> None:
        # 2 ambiguous (N, R) out of 6 non-gap chars
        assert ambiguity_fraction("ACGTNR") == pytest.approx(2 / 6)

    def test_all_gap_returns_zero(self) -> None:
        assert ambiguity_fraction("----") == 0.0


class TestTranslateInFrame:
    def test_translates_standard_codons(self) -> None:
        assert translate_in_frame("ATGGCA") == "MA"
        assert translate_in_frame("TAA") == "*"

    def test_gaps_and_ambiguity_become_x(self) -> None:
        assert translate_in_frame("A-G") == "X"
        assert translate_in_frame("NNN") == "X"


class TestHasPrematureStop:
    def test_true_when_internal_stop(self) -> None:
        # ATG (M) + TAA (*) + GCA (A) → stop in the middle
        assert has_premature_stop("ATGTAAGCA") is True

    def test_false_when_no_stops(self) -> None:
        assert has_premature_stop("ATGGCAGCG") is False

    def test_ignores_terminal_stop(self) -> None:
        # ATG + GCA + TAA — the final stop is fine
        assert has_premature_stop("ATGGCATAA") is False


class TestApobecSignatureRatio:
    def test_hypermutation_signal(self) -> None:
        # Reference has GG at pos 0 and GA at pos 2. Query mutates GG→AG (pos 0)
        # but leaves GA alone. Expected: gg_rate > 0, ga_rate = 0 → inf.
        ref = "GGGA"
        qry = "AGGA"
        assert apobec_signature_ratio(qry, ref) == float("inf")

    def test_no_signal_when_symmetric(self) -> None:
        # Both GG and GA sites unchanged — no APOBEC signal → 0.0
        assert apobec_signature_ratio("GGGA", "GGGA") == 0.0

    def test_returns_none_when_no_eligible_sites(self) -> None:
        # Reference with no GG or GA dinucleotides → can't compute
        assert apobec_signature_ratio("ACAC", "ACAC") is None

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="length mismatch"):
            apobec_signature_ratio("ACGT", "ACGTAA")


class TestQCFlags:
    def test_flags_combine(self) -> None:
        flags = qc_flags(
            ambiguity=0.10,  # > 0.05 threshold
            apobec_ratio=5.0,  # > 2.0 threshold
            region_coverage={"PR": 1.0, "RT": 0.5, "IN": 1.0},  # RT flags SHORT
            region_stops={"PR": True, "RT": False, "IN": False},  # PR flags STOP
        )
        assert QCFlag.HIGH_AMBIGUITY in flags
        assert QCFlag.APOBEC_HYPERMUT in flags
        assert QCFlag.SHORT_COVERAGE in flags
        assert QCFlag.STOP_CODON in flags

    def test_clean_sequence(self) -> None:
        flags = qc_flags(
            ambiguity=0.01,
            apobec_ratio=1.0,
            region_coverage={"PR": 1.0, "RT": 1.0, "IN": 1.0},
            region_stops={"PR": False, "RT": False, "IN": False},
        )
        assert flags == []
        assert passes_qc(flags) is True


# ---------- dedup ----------


class TestExtractPatientToken:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("RV217_20382_v06_GP_89_1", ("RV217", "20382")),
            ("MBAMARCH054_00", ("MBAMARCH", "054")),
            ("KER2008", ("KER", "2008")),
            ("Q23_17", ("Q23", "Q23")),
        ],
    )
    def test_matches_known_cohorts(self, name: str, expected: tuple[str, str]) -> None:
        assert extract_patient_token(name) == expected

    def test_unknown_falls_back(self) -> None:
        cohort, token = extract_patient_token("random_name_here_1")
        assert cohort == "unknown"
        assert token == "random_name"

    def test_empty_returns_no_name(self) -> None:
        assert extract_patient_token(None) == ("unknown", "no_name")
        assert extract_patient_token("") == ("unknown", "no_name")


class TestAssignPatientGroups:
    def test_same_patient_shares_group_id(self) -> None:
        rows = [
            {
                "accession": "A1",
                "sequence_name": "RV217_20382_v04",
                "country": "KE",
                "collection_year": 2015,
            },
            {
                "accession": "A2",
                "sequence_name": "RV217_20382_v06",
                "country": "KE",
                "collection_year": 2017,
            },
        ]
        out = assign_patient_groups(rows)
        assert out[0]["patient_group_id"] == out[1]["patient_group_id"]
        assert out[0]["is_earliest_per_patient"] is True
        assert out[1]["is_earliest_per_patient"] is False

    def test_different_patients_different_groups(self) -> None:
        rows = [
            {
                "accession": "A1",
                "sequence_name": "RV217_1",
                "country": "KE",
                "collection_year": 2015,
            },
            {
                "accession": "A2",
                "sequence_name": "RV217_2",
                "country": "KE",
                "collection_year": 2015,
            },
        ]
        out = assign_patient_groups(rows)
        assert out[0]["patient_group_id"] != out[1]["patient_group_id"]
        # Both are earliest of their own group.
        assert out[0]["is_earliest_per_patient"] is True
        assert out[1]["is_earliest_per_patient"] is True

    def test_same_year_ties_broken_by_accession(self) -> None:
        rows = [
            {
                "accession": "B",
                "sequence_name": "RV217_1",
                "country": "KE",
                "collection_year": 2020,
            },
            {
                "accession": "A",
                "sequence_name": "RV217_1",
                "country": "KE",
                "collection_year": 2020,
            },
        ]
        out = assign_patient_groups(rows)
        by_acc = {r["accession"]: r for r in out}
        assert by_acc["A"]["is_earliest_per_patient"] is True
        assert by_acc["B"]["is_earliest_per_patient"] is False


class TestPatientKeyHash:
    def test_stable(self) -> None:
        k1 = PatientKey("KE", "RV217", "20382", 2015)
        k2 = PatientKey("KE", "RV217", "20382", 2015)
        assert k1.hashed() == k2.hashed()
        assert len(k1.hashed()) == 16

    def test_differs_by_country(self) -> None:
        assert (
            PatientKey("KE", "RV217", "20382", 2015).hashed()
            != PatientKey("UG", "RV217", "20382", 2015).hashed()
        )
