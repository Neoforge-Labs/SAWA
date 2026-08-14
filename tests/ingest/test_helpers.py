"""Unit tests for pure helpers in genbank.py and lanl.py.

Golden-file tests for full ingestion runs live in test_integration.py
(not included in the day-1 commit — add once you have a fixture GenBank
response cached).
"""

from __future__ import annotations

from datetime import date

import pytest

from hivdrwatch.ingest.genbank import (
    _normalize_country,
    _normalize_subtype,
    _parse_collection_date,
)
from hivdrwatch.ingest.lanl import (
    _coerce_bool,
    _coerce_subtype,
    _country_from_folder,
    _country_name_to_code,
    _preamble_row_count,
    parse_lanl_fasta,
)
from hivdrwatch.ingest.schema import SequenceRecord, Source, Subtype


class TestParseCollectionDate:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("2019", date(2019, 7, 1)),
            ("Jun-2019", date(2019, 6, 15)),
            ("15-Jun-2019", date(2019, 6, 15)),
            ("2019-06-15", date(2019, 6, 15)),
        ],
    )
    def test_parses_common_formats(self, raw: str, expected: date) -> None:
        assert _parse_collection_date(raw) == expected

    def test_returns_none_on_garbage(self) -> None:
        assert _parse_collection_date("sometime in the 90s") is None

    def test_returns_none_on_empty(self) -> None:
        assert _parse_collection_date(None) is None
        assert _parse_collection_date("") is None


class TestNormalizeCountry:
    def test_maps_kenya(self) -> None:
        assert _normalize_country("Kenya") == "KE"

    def test_strips_subregion(self) -> None:
        assert _normalize_country("Kenya: Nairobi") == "KE"

    def test_returns_none_for_unknown(self) -> None:
        assert _normalize_country("Atlantis") is None


class TestNormalizeSubtype:
    def test_extracts_from_note(self) -> None:
        assert _normalize_subtype("subtype: A1") == Subtype.A1

    def test_extracts_without_colon(self) -> None:
        assert _normalize_subtype("HIV-1 subtype C from Kenya") == Subtype.C

    def test_recombinant_hint(self) -> None:
        assert _normalize_subtype("CRF10_CD recombinant") == Subtype.RECOMBINANT

    def test_unknown_default(self) -> None:
        assert _normalize_subtype("") == Subtype.UNKNOWN
        assert _normalize_subtype("no subtype info here") == Subtype.UNKNOWN


class TestLANLSubtypeCoercion:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("A1", "A1"),
            ("a1", "A1"),
            ("HIV-1_C", "C"),
            (" B ", "B"),
            ("CRF35_AD", "CRF35_AD"),
            ("A1/D", "recombinant"),
            ("CRF99_XX", "recombinant"),  # unknown CRF
            (None, "unknown"),
            ("", "unknown"),
        ],
    )
    def test_coerces(self, raw: str | None, expected: str) -> None:
        assert _coerce_subtype(raw) == expected


class TestLANLBoolCoercion:
    @pytest.mark.parametrize("raw", ["yes", "Y", "TRUE", "1", "naive"])
    def test_truthy(self, raw: str) -> None:
        assert _coerce_bool(raw) is True

    @pytest.mark.parametrize("raw", ["no", "N", "false", "0", "treated"])
    def test_falsy(self, raw: str) -> None:
        assert _coerce_bool(raw) is False

    def test_unknown_returns_none(self) -> None:
        assert _coerce_bool("maybe") is None
        assert _coerce_bool(None) is None


class TestCountryNameToCode:
    def test_common(self) -> None:
        assert _country_name_to_code("Kenya") == "KE"
        assert _country_name_to_code("united republic of tanzania") == "TZ"

    def test_unknown(self) -> None:
        assert _country_name_to_code("Wakanda") is None


class TestSequenceRecordValidation:
    def _base(self, **overrides) -> dict:
        d = {
            "accession": "MK123456",
            "source": Source.GENBANK,
            "sequence": "ACGT" * 100,
            "length": 400,
            "ingested_at": date(2026, 8, 10),
            "raw_record_hash": "a" * 64,
        }
        d.update(overrides)
        return d

    def test_valid_record(self) -> None:
        rec = SequenceRecord(**self._base())
        assert rec.accession == "MK123456"
        assert rec.subtype == Subtype.UNKNOWN

    def test_rejects_short_sequence(self) -> None:
        with pytest.raises(ValueError):
            SequenceRecord(**self._base(sequence="ACGT", length=4))

    def test_rejects_bad_alphabet(self) -> None:
        with pytest.raises(ValueError, match="non-IUPAC"):
            SequenceRecord(**self._base(sequence="ACGTZZZZ" * 50, length=400))

    def test_length_must_match_sequence(self) -> None:
        with pytest.raises(ValueError, match="does not match"):
            SequenceRecord(**self._base(sequence="ACGT" * 100, length=999))

    def test_uppercases_sequence(self) -> None:
        rec = SequenceRecord(**self._base(sequence="acgt" * 100))
        assert rec.sequence == "ACGT" * 100

    def test_converts_u_to_t(self) -> None:
        rec = SequenceRecord(**self._base(sequence="acgu" * 100))
        assert "U" not in rec.sequence
        assert rec.sequence.count("T") == 100  # ACGT × 100 → one T per repeat


class TestLANLPreambleDetection:
    def test_skips_number_of_records_preamble(self, tmp_path) -> None:
        p = tmp_path / "results.tsv"
        p.write_text(
            "Number of records retrieved: 42\n#\tAccession\tSubtype\n1\tAF004885\tA1\n"
        )
        assert _preamble_row_count(p) == 1

    def test_no_preamble(self, tmp_path) -> None:
        p = tmp_path / "results.tsv"
        p.write_text("#\tAccession\tSubtype\n1\tAF004885\tA1\n")
        assert _preamble_row_count(p) == 0


class TestCountryFromFolder:
    def test_extracts_two_letter_code(self, tmp_path) -> None:
        (tmp_path / "ke").mkdir()
        f = tmp_path / "ke" / "results.tsv"
        f.touch()
        assert _country_from_folder(f, tmp_path) == "KE"

    def test_uppercases_folder(self, tmp_path) -> None:
        (tmp_path / "ug").mkdir()
        f = tmp_path / "ug" / "results.txt"
        f.touch()
        assert _country_from_folder(f, tmp_path) == "UG"

    def test_none_for_flat_layout(self, tmp_path) -> None:
        f = tmp_path / "results.tsv"
        f.touch()
        assert _country_from_folder(f, tmp_path) is None

    def test_none_for_non_country_folder(self, tmp_path) -> None:
        (tmp_path / "kenya").mkdir()
        f = tmp_path / "kenya" / "results.tsv"
        f.touch()
        assert _country_from_folder(f, tmp_path) is None


class TestParseLANLFasta:
    def test_extracts_accession_from_dotted_header(self, tmp_path) -> None:
        p = tmp_path / "hiv-db.fasta"
        p.write_text(
            ">A1.KE.1994.Q23_17.AF004885\n"
            "ACGT-ACGT\n"
            "ACGT\n"
            ">21_A2D.KE.1999.KER2003.AF457051\n"
            "TTTAAA\n"
        )
        result = parse_lanl_fasta(p)
        assert result == {"AF004885": "ACGTACGTACGT", "AF457051": "TTTAAA"}

    def test_uppercases_and_converts_u(self, tmp_path) -> None:
        p = tmp_path / "hiv-db.fasta"
        p.write_text(">x.y.z.name.AB123456\nacgu\n")
        assert parse_lanl_fasta(p) == {"AB123456": "ACGT"}

    def test_empty_file(self, tmp_path) -> None:
        p = tmp_path / "hiv-db.fasta"
        p.write_text("")
        assert parse_lanl_fasta(p) == {}
