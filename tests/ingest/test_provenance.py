"""Regression tests for the three ingest provenance bugs.

Each bug silently rewrote metadata that two papers stratify on, so the
free-text cases here are the *real* `/note=` and `/geo_loc_name=` strings
that occur in `data/raw/genbank` (surveyed over all 36,510 cached
accessions), not invented ones. All offline — no network, no cache.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import polars as pl
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from hivdrwatch.ingest.genbank import (
    GenBankIngester,
    _normalize_country,
    _normalize_subtype,
    _stated_country,
)
from hivdrwatch.ingest.merge import merge_sources
from hivdrwatch.ingest.schema import Subtype

# ---------------------------------------------------------------------------
# Bug 1 — `_normalize_subtype` mapped every CRF to subtype C
# ---------------------------------------------------------------------------


class TestSubtypeEnumRoundTrip:
    """Every concrete enum member must survive `subtype: <member>`.

    This is the test the substring match failed: `"SUBTYPE: C"` is a
    substring of `"SUBTYPE: CRF01_AE"`, so every CRF came back as C.
    """

    @pytest.mark.parametrize(
        "member",
        [st for st in Subtype if st not in (Subtype.UNKNOWN, Subtype.RECOMBINANT)],
        ids=lambda st: st.value,
    )
    def test_member_round_trips(self, member: Subtype) -> None:
        assert _normalize_subtype(f"subtype: {member.value}") == member

    @pytest.mark.parametrize(
        "member",
        [st for st in Subtype if st.value.startswith("CRF")],
        ids=lambda st: st.value,
    )
    def test_crf_is_not_subtype_c(self, member: Subtype) -> None:
        assert _normalize_subtype(f"subtype: {member.value}") is not Subtype.C


class TestSubtypeFreeTextFromCache:
    """Forms that actually occur in the cached flat files."""

    @pytest.mark.parametrize(
        "note,expected",
        [
            # --- unambiguous single subtypes ---
            ("subtype: A1", Subtype.A1),
            ("subtype: C", Subtype.C),
            ("subtype: D", Subtype.D),
            ("subtype: B", Subtype.B),
            ("subtype: A2", Subtype.A2),
            ("subtype: G", Subtype.G),
            ("subtype: J", Subtype.J),
            ("subtype:A1", Subtype.A1),  # no space after the colon
            ("subtype: A1(100.0%)", Subtype.A1),  # COMET-style confidence
            ("subtype: C(100.0%)", Subtype.C),
            ("subtype: A1; group: M", Subtype.A1),
            ("IAVI-CHAVI; subtype: D", Subtype.D),
            ("LPV/r = 2nd or higher PI; subtype: C", Subtype.C),
            ("full-length genome; subtype: D", Subtype.D),
            ("genotype: D", Subtype.D),  # depositors also say "genotype"
            ("genotype: A1", Subtype.A1),
            ("genotype: HIV-1C", Subtype.C),
            ("genotype: 1C", Subtype.C),
            ("clade A1", Subtype.A1),
            # --- parent clades the depositor did not sub-classify ---
            ("subtype: A", Subtype.A),
            ("subtype: A; group: M", Subtype.A),
            ("IAVI-CHAVI; subtype: A", Subtype.A),
            ("subtype: A (A1)", Subtype.A1),  # parent + specific → specific
            ("subtype: F", Subtype.F),
            ("subtype: F (F1)", Subtype.F1),
            # --- named CRFs: the whole point of bug 1 ---
            ("subtype: CRF01_AE", Subtype.CRF01_AE),
            ("subtype: CRF02_AG", Subtype.CRF02_AG),
            ("subtype: CRF35_AD", Subtype.CRF35_AD),
            ("subtype: CRF10_CD", Subtype.CRF10_CD),
            ("subtype: CRF16_A2D", Subtype.CRF16_A2D),
            ("subtype: CRF21_A2D", Subtype.CRF21_A2D),
            ("subtype: CRF 10_CD", Subtype.CRF10_CD),  # spaced
            ("subtype: CRF_10_CD", Subtype.CRF10_CD),  # extra underscore
            ("subtype: CRF01_AE; group: M", Subtype.CRF01_AE),
            # --- CRFs with no enum member of their own ---
            ("subtype: CRF12_BF", Subtype.RECOMBINANT),
            ("subtype: CRF45_CPX", Subtype.RECOMBINANT),
            ("subtype: CRF10", Subtype.RECOMBINANT),  # truncated designation
            # --- explicit recombinants and mosaics ---
            ("subtype: A1/D recombinant", Subtype.RECOMBINANT),
            ("subtype: A1/C", Subtype.RECOMBINANT),
            ("subtype: A1/C/D recombinant", Subtype.RECOMBINANT),
            ("subtype: A1,D", Subtype.RECOMBINANT),
            ("subtype: A1 D", Subtype.RECOMBINANT),
            ("subtype: D, A1", Subtype.RECOMBINANT),
            ("subtype: recombinant C,A1", Subtype.RECOMBINANT),
            ("subtype: A1/A2/D recombinant", Subtype.RECOMBINANT),
            ("subtype: A2/C/A2", Subtype.RECOMBINANT),
            ("subtype: CRF35_AD/A/CRF35_AD", Subtype.RECOMBINANT),
            ("subtype: B/CRF01_AE recombinant", Subtype.RECOMBINANT),
            ("subtype: URF_CD", Subtype.RECOMBINANT),
            ("subtype: URF_A1C", Subtype.RECOMBINANT),
            ("subtype A and D", Subtype.RECOMBINANT),
            ("subtype: AD", Subtype.RECOMBINANT),
            ("subtype: A1D", Subtype.RECOMBINANT),
            ("subtype: ACD", Subtype.RECOMBINANT),
            ("subtype: AE", Subtype.RECOMBINANT),
            ("subtype: AG", Subtype.RECOMBINANT),
            ("genotype: AC", Subtype.RECOMBINANT),
            (
                "intersubtype recombinant virus which has been vertically transmitted",
                Subtype.RECOMBINANT,
            ),
            # --- genuinely unknown ---
            ("subtype: U", Subtype.UNKNOWN),  # LANL "unclassified"
            ("pre-ART drug resistance", Subtype.UNKNOWN),
            ("group: M", Subtype.UNKNOWN),
            ("collected during 2009-2010", Subtype.UNKNOWN),
            ("ThermoFisher HIV-1 Genotyping Kit", Subtype.UNKNOWN),
            ("", Subtype.UNKNOWN),
        ],
    )
    def test_real_cache_forms(self, note: str, expected: Subtype) -> None:
        assert _normalize_subtype(note) == expected


class TestSubtypeDoesNotReadProse:
    """The declaration ends where the designation does."""

    def test_trailing_prose_is_ignored(self) -> None:
        assert _normalize_subtype("HIV-1 subtype C from Kenya") == Subtype.C

    def test_lowercase_prose_is_not_a_mosaic(self) -> None:
        # "had" is A, D and H concatenated; only uppercase shorthand counts.
        assert _normalize_subtype("subtype had not been determined") == Subtype.UNKNOWN

    def test_substring_of_a_longer_word_is_not_a_key(self) -> None:
        assert _normalize_subtype("intersubtype variation was assessed") == Subtype.UNKNOWN


# ---------------------------------------------------------------------------
# Bug 2 — `_normalize_country` fabricated geography
# ---------------------------------------------------------------------------


class TestStatedCountry:
    """NCBI renamed `/country=` to `/geo_loc_name=`; both must be read.

    Not one record in the 36,510-accession local cache carries
    `/country=`, so reading only that qualifier made *every* country in
    the corpus the query fallback.
    """

    def test_reads_geo_loc_name(self) -> None:
        assert _stated_country({"geo_loc_name": "Kenya"}) == "Kenya"

    def test_reads_legacy_country(self) -> None:
        assert _stated_country({"country": "Kenya"}) == "Kenya"

    def test_prefers_geo_loc_name(self) -> None:
        assert _stated_country({"geo_loc_name": "Uganda", "country": "Kenya"}) == "Uganda"

    def test_none_when_absent(self) -> None:
        assert _stated_country({"organism": "Human immunodeficiency virus 1"}) is None


class TestNormalizeCountry:
    @pytest.mark.parametrize(
        "name,code",
        [
            ("Kenya", "KE"),
            ("Kenya: Nairobi", "KE"),
            ("Uganda", "UG"),
            ("Tanzania", "TZ"),
            ("Ethiopia", "ET"),
            ("Rwanda", "RW"),
            ("Burundi", "BI"),
            # Non-East-African names that occur in the cache: these used to
            # return None and inherit the query's country.
            ("Thailand", "TH"),
            ("South Africa", "ZA"),
            ("Brazil", "BR"),
            ("Malawi", "MW"),
            ("India", "IN"),
            ("Zimbabwe", "ZW"),
            ("Zambia", "ZM"),
            ("Israel", "IL"),
            ("Haiti", "HT"),
            ("Nigeria", "NG"),
            ("USA", "US"),
            ("Republic of the Congo", "CG"),
            ("Dominican Republic", "DO"),
            ("Cameroon", "CM"),
            ("Cote d'Ivoire", "CI"),
            ("Peru", "PE"),
            ("Argentina", "AR"),
            ("Spain", "ES"),
            ("Italy", "IT"),
            ("Cuba", "CU"),
            ("United Kingdom", "GB"),
            ("Finland", "FI"),
            ("Sweden", "SE"),
        ],
    )
    def test_maps_every_name_in_the_cache(self, name: str, code: str) -> None:
        assert _normalize_country(name) == code

    def test_distinguishes_the_two_congos(self) -> None:
        assert _normalize_country("Republic of the Congo") == "CG"
        assert _normalize_country("Democratic Republic of the Congo") == "CD"

    def test_unmapped_name_returns_none(self) -> None:
        assert _normalize_country("Atlantis") is None


def _record(accession: str, source_quals: dict[str, list[str]] | None) -> SeqRecord:
    """A minimal SeqRecord with (optionally) a `source` feature."""
    rec = SeqRecord(Seq("ACGT" * 300), id=accession)
    if source_quals is not None:
        rec.features.append(
            SeqFeature(SimpleLocation(0, 1200), type="source", qualifiers=source_quals)
        )
    return rec


class TestParseDoesNotInheritRequestedCountry:
    @pytest.fixture
    def ingester(self, tmp_path: Path) -> GenBankIngester:
        return GenBankIngester(email="test@example.org", raw_dir=tmp_path / "raw")

    def test_keeps_the_records_own_country(self, ingester: GenBankIngester) -> None:
        parsed = ingester._parse(
            _record("MK123456", {"geo_loc_name": ["Thailand"]}), requested_country="KE"
        )
        assert parsed is not None
        assert parsed.country == "TH"
        assert parsed.country_source == "genbank_source_feature"
        assert parsed.country_raw == "Thailand"

    def test_keeps_subdivision_verbatim(self, ingester: GenBankIngester) -> None:
        parsed = ingester._parse(
            _record("MK123457", {"geo_loc_name": ["Kenya: Nairobi"]}), requested_country="UG"
        )
        assert parsed is not None
        assert parsed.country == "KE"
        assert parsed.country_raw == "Kenya: Nairobi"

    def test_no_country_field_means_no_country(self, ingester: GenBankIngester) -> None:
        parsed = ingester._parse(
            _record("MK123458", {"organism": ["Human immunodeficiency virus 1"]}),
            requested_country="KE",
        )
        assert parsed is not None
        assert parsed.country is None
        assert parsed.country_source is None
        assert parsed.country_raw is None

    def test_stated_but_unmappable_country_is_not_dropped(
        self, ingester: GenBankIngester
    ) -> None:
        parsed = ingester._parse(
            _record("MK123459", {"geo_loc_name": ["Atlantis"]}), requested_country="KE"
        )
        assert parsed is not None
        assert parsed.country is None
        assert parsed.country_raw == "Atlantis"

    def test_subtype_source_only_claimed_when_resolved(
        self, ingester: GenBankIngester
    ) -> None:
        resolved = ingester._parse(
            _record("MK123460", {"note": ["subtype: CRF01_AE"]}), requested_country="KE"
        )
        assert resolved is not None
        assert resolved.subtype == Subtype.CRF01_AE
        assert resolved.subtype_source == "genbank"

        unresolved = ingester._parse(
            _record("MK123461", {"note": ["pre-ART drug resistance"]}),
            requested_country="KE",
        )
        assert unresolved is not None
        assert unresolved.subtype == Subtype.UNKNOWN
        assert unresolved.subtype_source is None


# ---------------------------------------------------------------------------
# Bug 3 — one row per country query, and a single patient group per country
# ---------------------------------------------------------------------------


def _merge_row(accession: str, country: str) -> dict:
    return {
        "accession": accession,
        "source": "genbank",
        "sequence": "ACGT" * 300,
        "length": 1200,
        "gene_region": "pol",
        "country": country,
        "collection_year": 2019,
        "subtype": "A1",
        "subtype_source": "genbank",
        "ingested_at": "2026-10-21",
        "raw_record_hash": "a" * 64,
    }


class TestMergeCollapsesDuplicateAccessions:
    """`fetch_country` is called once per country, so a record returned by
    several country queries was appended several times and merged as
    several rows."""

    def test_same_accession_from_three_queries_becomes_one_row(self, tmp_path: Path) -> None:
        genbank_df = pl.from_dicts(
            [
                _merge_row("MK123456", "KE"),
                _merge_row("MK123456", "UG"),
                _merge_row("MK123456", "TZ"),
                _merge_row("MK999999", "KE"),
            ]
        )
        merged = merge_sources(
            genbank_df,
            pl.DataFrame({"accession": ["MK123456"], "subtype_lanl": ["CRF01_AE"]}),
            tmp_path / "merged.parquet",
        )
        assert merged.height == 2
        assert merged["accession"].n_unique() == 2

    def test_duplicate_lanl_rows_do_not_multiply_genbank_rows(self, tmp_path: Path) -> None:
        merged = merge_sources(
            pl.from_dicts([_merge_row("MK123456", "KE")]),
            pl.DataFrame(
                {
                    "accession": ["MK123456", "MK123456"],
                    "subtype_lanl": ["A1", "A1"],
                }
            ),
            tmp_path / "merged.parquet",
        )
        assert merged.height == 1


class TestGeneRegionReflectsSequenceExtent:
    """A gag→nef deposit is 8.7 kb; calling it "pol" is a provenance claim
    the record does not support."""

    @pytest.fixture
    def ingester(self, tmp_path: Path) -> GenBankIngester:
        return GenBankIngester(email="test@example.org", raw_dir=tmp_path / "raw")

    def test_pol_sized_deposit_is_pol(self, ingester: GenBankIngester) -> None:
        parsed = ingester._parse(_record("MK123462", {}), requested_country="KE")
        assert parsed is not None
        assert parsed.gene_region == "pol"

    def test_whole_genome_deposit_is_not_labelled_pol(self, ingester: GenBankIngester) -> None:
        rec = SeqRecord(Seq("ACGT" * 2200), id="MK123463")  # 8.8 kb
        rec.features.append(SeqFeature(SimpleLocation(0, 8800), type="source", qualifiers={}))
        parsed = ingester._parse(rec, requested_country="KE")
        assert parsed is not None
        assert parsed.gene_region == "genome"


class TestIngestedAtIsStamped:
    """Guard the one thing the fixes must not disturb."""

    def test_ingested_at_is_today(self, tmp_path: Path) -> None:
        ing = GenBankIngester(email="test@example.org", raw_dir=tmp_path / "raw")
        parsed = ing._parse(_record("MK123464", {}), requested_country="KE")
        assert parsed is not None
        assert parsed.ingested_at == date.today()
