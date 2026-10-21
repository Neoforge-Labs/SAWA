"""GenBank ingestion without a network, driven off the flat-file cache.

The cache-replay path is how a vintage gets rebuilt after a parser change
(`docs/evaluations/ingest-provenance.md` §8), so it is worth a test of its
own. `_esearch_all` is the only thing stubbed; everything downstream of it
runs for real against flat files written into `tmp_path`.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import pytest

from hivdrwatch.ingest.genbank import (
    GenBankIngester,
    _looks_valid_gb,
    _text_looks_valid_gb,
)
from hivdrwatch.ingest.schema import Subtype


def _flat_file(accession: str, geo: str, note: str, length: int = 1200) -> str:
    """A minimal but genuine GenBank flat file, in the cache's own shape."""
    bases = ("acgt" * ((length // 4) + 1))[:length]
    lines = [
        f"LOCUS       {accession.split('.')[0]}{length:>21} bp    RNA     linear   VRL 02-JUL-2016",
        "DEFINITION  HIV-1 pol protein (pol) gene, partial cds.",
        f"ACCESSION   {accession.split('.')[0]}",
        f"VERSION     {accession}",
        "SOURCE      Human immunodeficiency virus 1 (HIV-1)",
        "  ORGANISM  Human immunodeficiency virus 1",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{length}",
        '                     /organism="Human immunodeficiency virus 1"',
        f'                     /geo_loc_name="{geo}"',
        '                     /collection_date="2010"',
        f'                     /note="{note}"',
        "ORIGIN      ",
    ]
    for offset in range(0, length, 60):
        chunk = bases[offset : offset + 60]
        groups = " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10))
        lines.append(f"{offset + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


class TestBuildQuery:
    """The query that over-matches. Pure; exercised here so a change to it
    is visible in a diff rather than only at the next pull."""

    def test_covers_pol_and_its_subregions(self) -> None:
        query = GenBankIngester._build_query("Kenya", None, None, 800)
        assert '"Human immunodeficiency virus 1"[Organism]' in query
        assert '"Kenya"[Country]' in query
        assert "800:100000[SLEN]" in query
        for gene in ("pol[Gene]", "protease[Gene]", "integrase[Gene]"):
            assert gene in query
        assert '"reverse transcriptase"[Product]' in query
        assert "[MDAT]" not in query

    def test_date_window_uses_modification_date(self) -> None:
        query = GenBankIngester._build_query(
            "Uganda", date(2020, 1, 1), date(2024, 6, 30), 800
        )
        assert '("2020/01/01"[MDAT] : "2024/06/30"[MDAT])' in query

    @pytest.mark.parametrize(
        "since,until,expected",
        [
            (date(2020, 1, 1), None, '("2020/01/01"[MDAT] : "3000/01/01"[MDAT])'),
            (None, date(2020, 1, 1), '("1980/01/01"[MDAT] : "2020/01/01"[MDAT])'),
        ],
    )
    def test_open_ended_windows(
        self, since: date | None, until: date | None, expected: str
    ) -> None:
        assert expected in GenBankIngester._build_query("Kenya", since, until, 800)


class TestCacheValidation:
    def test_complete_file_is_valid(self, tmp_path: Path) -> None:
        path = tmp_path / "batch.gb"
        path.write_text(_flat_file("MK123456.1", "Kenya", "subtype: A1"))
        assert _looks_valid_gb(path) is True

    def test_truncated_file_is_not(self, tmp_path: Path) -> None:
        path = tmp_path / "batch.gb"
        path.write_text(_flat_file("MK123456.1", "Kenya", "subtype: A1")[:-40])
        assert _looks_valid_gb(path) is False

    def test_empty_and_missing_files_are_not(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.gb"
        empty.write_text("")
        assert _looks_valid_gb(empty) is False
        assert _looks_valid_gb(tmp_path / "absent.gb") is False

    @pytest.mark.parametrize(
        "text,valid",
        [("LOCUS x\n//\n", True), ("LOCUS x\n//  \n\n", True), ("LOCUS x\nORIGIN\n", False)],
    )
    def test_in_memory_equivalent(self, text: str, valid: bool) -> None:
        assert _text_looks_valid_gb(text) is valid


class TestFetchCountryFromCache:
    """A full `fetch_country` with no network: the batch file is already in
    `raw_dir`, so `_efetch_batch` takes its cache-hit branch."""

    @pytest.fixture
    def raw_dir(self, tmp_path: Path) -> Path:
        d = tmp_path / "genbank"
        d.mkdir()
        # `_efetch_batch` names the cache file after the first and last UID.
        (d / "batch_1001_1003.gb").write_text(
            _flat_file("MK000001.1", "Kenya: Nairobi", "subtype: CRF01_AE")
            + _flat_file("MK000002.1", "Thailand", "subtype: CRF01_AE", length=1600)
            + _flat_file("MK000003.1", "Uganda", "subtype: A", length=8800)
        )
        return d

    @pytest.fixture
    def ingester(self, raw_dir: Path, monkeypatch: pytest.MonkeyPatch) -> GenBankIngester:
        ing = GenBankIngester(email="test@example.org", raw_dir=raw_dir, batch_size=50)
        monkeypatch.setattr(
            GenBankIngester, "_esearch_all", lambda self, query: ["1001", "1002", "1003"]
        )
        return ing

    def test_yields_every_cached_record(self, ingester: GenBankIngester) -> None:
        records = list(ingester.fetch_country("KE"))
        assert [r.accession for r in records] == ["MK000001.1", "MK000002.1", "MK000003.1"]

    def test_country_comes_from_the_record_not_the_query(
        self, ingester: GenBankIngester
    ) -> None:
        by_accession = {r.accession: r for r in ingester.fetch_country("KE")}
        assert by_accession["MK000001.1"].country == "KE"
        assert by_accession["MK000001.1"].country_raw == "Kenya: Nairobi"
        # Returned by the Kenya query, but the record says Thailand.
        assert by_accession["MK000002.1"].country == "TH"
        assert by_accession["MK000003.1"].country == "UG"
        assert all(r.country_source == "genbank_source_feature" for r in by_accession.values())

    def test_crf_is_not_relabelled_subtype_c(self, ingester: GenBankIngester) -> None:
        by_accession = {r.accession: r for r in ingester.fetch_country("KE")}
        assert by_accession["MK000001.1"].subtype == Subtype.CRF01_AE
        assert by_accession["MK000003.1"].subtype == Subtype.A

    def test_whole_genome_is_not_labelled_pol(self, ingester: GenBankIngester) -> None:
        by_accession = {r.accession: r for r in ingester.fetch_country("KE")}
        assert by_accession["MK000001.1"].gene_region == "pol"
        assert by_accession["MK000003.1"].gene_region == "genome"

    def test_truncated_cache_entry_is_not_replayed(
        self, raw_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A corrupt cache file must be discarded, not parsed. With the
        network stubbed out there is nothing to refetch with, so the call
        raising rather than yielding half a record is the assertion."""
        path = raw_dir / "batch_1001_1003.gb"
        path.write_text(path.read_text()[:-40])
        ing = GenBankIngester(email="test@example.org", raw_dir=raw_dir)
        monkeypatch.setattr(
            GenBankIngester, "_esearch_all", lambda self, query: ["1001", "1002", "1003"]
        )
        monkeypatch.setattr(
            GenBankIngester,
            "_rate_limit",
            lambda self: (_ for _ in ()).throw(OSError("no network in tests")),
        )
        with pytest.raises(OSError):
            list(ing.fetch_country("KE"))
        assert not path.exists()

    def test_unknown_country_code_is_rejected(self, ingester: GenBankIngester) -> None:
        with pytest.raises(ValueError, match="unknown country code"):
            list(ingester.fetch_country("XX"))

    def test_empty_search_yields_nothing(
        self, raw_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        ing = GenBankIngester(email="test@example.org", raw_dir=raw_dir)
        monkeypatch.setattr(GenBankIngester, "_esearch_all", lambda self, query: [])
        assert list(ing.fetch_country("KE")) == []


class TestIngesterConstruction:
    def test_requires_a_real_email(self, tmp_path: Path) -> None:
        for bad in ("", "not-an-email"):
            with pytest.raises(ValueError, match="real email address"):
                GenBankIngester(email=bad, raw_dir=tmp_path / "raw")

    def test_creates_the_raw_directory(self, tmp_path: Path) -> None:
        raw = tmp_path / "nested" / "genbank"
        GenBankIngester(email="test@example.org", raw_dir=raw)
        assert raw.is_dir()

    def test_rate_limit_tightens_with_an_api_key(self, tmp_path: Path) -> None:
        with_key = GenBankIngester(
            email="test@example.org", raw_dir=tmp_path / "a", api_key="k"
        )
        assert with_key._min_interval < 0.2
