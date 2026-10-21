"""The LANL path end to end, from a Table export on disk to SequenceRecords.

This is the branch `--skip-genbank` and the offline pipeline run on, and it
is where LANL's `Country` column and `Subtype` column become the canonical
`country` / `subtype` fields — so it needs the same scrutiny as the GenBank
path. Fixtures are written into `tmp_path` and named `_synthetic_*` so they
can never be mistaken for real records.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import pytest

from hivdrwatch.ingest.lanl import LANLIngester
from hivdrwatch.ingest.schema import SCHEMA_VERSION, Source, Subtype

_HEADER = "#\tAccession\tName\tSubtype\tCountry\tSampling Year\tDrug Naive\tSequence"


def _table(*rows: str, preamble: bool = True) -> str:
    lines = [f"Number of records retrieved: {len(rows)}"] if preamble else []
    return "\n".join([*lines, _HEADER, *rows]) + "\n"


def _row(
    index: int,
    accession: str,
    name: str,
    subtype: str,
    country: str,
    year: str = "2010",
    naive: str = "yes",
    sequence: str | None = "ACGT",
) -> str:
    seq = (sequence * 300) if sequence else ""
    return f"{index}\t{accession}\t{name}\t{subtype}\t{country}\t{year}\t{naive}\t{seq}"


@pytest.fixture
def raw_dir(tmp_path: Path) -> Path:
    ke = tmp_path / "ke"
    ke.mkdir()
    (ke / "results.tsv").write_text(
        _table(
            _row(1, "AF004885", "_synthetic_KER2003_1", "A1", "Kenya"),
            _row(2, "AF457051", "_synthetic_KER2004_1", "CRF21_A2D", "Kenya"),
            # A mosaic shorthand the old coercer mapped to `unknown`.
            _row(3, "AF457052", "_synthetic_KER2005_1", "AD", "Kenya"),
            # Bare A: no enum member before 0.2.0.
            _row(4, "AF457053", "_synthetic_KER2006_1", "A", "Kenya"),
        )
    )
    ug = tmp_path / "ug"
    ug.mkdir()
    (ug / "results.txt").write_text(
        _table(
            _row(1, "AF457054", "_synthetic_UG1001_1", "D", "Uganda"),
            # No sequence in the TSV — must come from the FASTA fallback.
            _row(2, "AF457055", "_synthetic_UG1002_1", "D", "Uganda", sequence=None),
        )
    )
    (ug / "hiv-db.fasta").write_text(
        ">D.UG.2010._synthetic_UG1002_1.AF457055\n" + "TTTA" * 300 + "\n"
    )
    return tmp_path


class TestLoadDirectory:
    def test_reads_every_country_folder(self, raw_dir: Path) -> None:
        df = LANLIngester(raw_dir).load_directory()
        assert df.height == 6
        assert sorted(df["country_code"].unique().to_list()) == ["KE", "UG"]

    def test_skips_the_record_count_preamble(self, raw_dir: Path) -> None:
        df = LANLIngester(raw_dir).load_directory()
        assert "AF004885" in df["accession"].to_list()
        assert df["accession"].null_count() == 0

    def test_stamps_the_schema_version_and_source_file(self, raw_dir: Path) -> None:
        df = LANLIngester(raw_dir).load_directory()
        assert df["schema_version"].unique().to_list() == [SCHEMA_VERSION]
        assert sorted(df["source_file"].unique().to_list()) == [
            "ke/results.tsv",
            "ug/results.txt",
        ]

    def test_country_comes_from_the_folder_or_the_column(self, raw_dir: Path) -> None:
        df = LANLIngester(raw_dir).load_directory()
        by_accession = {r["accession"]: r for r in df.iter_rows(named=True)}
        assert by_accession["AF004885"]["country_code"] == "KE"
        assert by_accession["AF457054"]["country_code"] == "UG"

    def test_subtypes_use_the_shared_classifier(self, raw_dir: Path) -> None:
        df = LANLIngester(raw_dir).load_directory()
        by_accession = {r["accession"]: r for r in df.iter_rows(named=True)}
        assert by_accession["AF004885"]["subtype_lanl"] == Subtype.A1.value
        assert by_accession["AF457051"]["subtype_lanl"] == Subtype.CRF21_A2D.value
        # Both of these were `unknown` before 0.2.0.
        assert by_accession["AF457052"]["subtype_lanl"] == Subtype.RECOMBINANT.value
        assert by_accession["AF457053"]["subtype_lanl"] == Subtype.A.value

    def test_empty_directory_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            LANLIngester(tmp_path / "nothing").load_directory()


class TestToSequenceRecords:
    def test_builds_canonical_records(self, raw_dir: Path) -> None:
        records = list(LANLIngester(raw_dir).to_sequence_records())
        assert len(records) == 6
        assert {r.source for r in records} == {Source.LANL}
        assert {r.gene_region for r in records} == {"pol"}
        assert all(r.ingested_at == date.today() for r in records)
        assert all(len(r.raw_record_hash) == 64 for r in records)

    def test_country_provenance_is_recorded(self, raw_dir: Path) -> None:
        by_accession = {r.accession: r for r in LANLIngester(raw_dir).to_sequence_records()}
        kenyan = by_accession["AF004885"]
        assert kenyan.country == "KE"
        assert kenyan.country_raw == "Kenya"
        assert kenyan.country_source == "lanl_curated"

    def test_subtype_provenance_is_recorded(self, raw_dir: Path) -> None:
        by_accession = {r.accession: r for r in LANLIngester(raw_dir).to_sequence_records()}
        assert by_accession["AF457051"].subtype == Subtype.CRF21_A2D
        assert by_accession["AF457051"].subtype_source == "lanl_curated"

    def test_sequence_falls_back_to_the_fasta(self, raw_dir: Path) -> None:
        by_accession = {r.accession: r for r in LANLIngester(raw_dir).to_sequence_records()}
        assert by_accession["AF457055"].sequence.startswith("TTTA")
        assert by_accession["AF457055"].length == 1200

    def test_rows_with_no_sequence_anywhere_are_dropped(self, tmp_path: Path) -> None:
        ke = tmp_path / "ke"
        ke.mkdir()
        (ke / "results.tsv").write_text(
            _table(
                _row(1, "AF004885", "_synthetic_KER2003_1", "A1", "Kenya"),
                _row(2, "AF457051", "_synthetic_KER2004_1", "A1", "Kenya", sequence=None),
            )
        )
        records = list(LANLIngester(tmp_path).to_sequence_records())
        assert [r.accession for r in records] == ["AF004885"]
