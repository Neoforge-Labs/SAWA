"""`merge_sources` — the LANL-only branch and the coalescing contract.

The LANL-only branch is what `--skip-genbank` and the offline pipeline
use, and it is the branch that writes the canonical dataset when GenBank
is unavailable.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from hivdrwatch.ingest.merge import merge_sources
from hivdrwatch.ingest.schema import SCHEMA_VERSION


def _lanl_frame(**overrides: object) -> pl.DataFrame:
    data: dict[str, object] = {
        "accession": ["AF004885", "AF457051"],
        "sequence": ["ACGT" * 300, "TTTA" * 300],
        "subtype_lanl": ["A1", "CRF21_A2D"],
        "country_code": ["KE", "KE"],
        "collection_year": [1994, 1999],
    }
    data.update(overrides)
    return pl.DataFrame(data)


class TestLANLOnlyPromotion:
    def test_promotes_lanl_to_the_canonical_layout(self, tmp_path: Path) -> None:
        merged = merge_sources(
            pl.DataFrame({"accession": []}), _lanl_frame(), tmp_path / "out.parquet"
        )
        assert merged.height == 2
        assert merged["source"].unique().to_list() == ["lanl"]
        assert merged["gene_region"].unique().to_list() == ["pol"]
        assert merged["schema_version"].unique().to_list() == [SCHEMA_VERSION]
        assert merged["merged_schema_version"].unique().to_list() == [SCHEMA_VERSION]
        assert sorted(merged["subtype"].to_list()) == ["A1", "CRF21_A2D"]
        assert merged["subtype_source"].unique().to_list() == ["lanl_curated"]
        assert merged["country"].unique().to_list() == ["KE"]
        assert "subtype_lanl" not in merged.columns
        assert "country_code" not in merged.columns

    def test_writes_the_parquet(self, tmp_path: Path) -> None:
        out = tmp_path / "nested" / "out.parquet"
        merge_sources(pl.DataFrame({"accession": []}), _lanl_frame(), out)
        assert pl.read_parquet(out).height == 2

    def test_rows_without_a_sequence_are_dropped(self, tmp_path: Path) -> None:
        merged = merge_sources(
            pl.DataFrame({"accession": []}),
            _lanl_frame(sequence=["ACGT" * 300, None]),
            tmp_path / "out.parquet",
        )
        assert merged["accession"].to_list() == ["AF004885"]

    def test_duplicate_lanl_accessions_collapse(self, tmp_path: Path) -> None:
        merged = merge_sources(
            pl.DataFrame({"accession": []}),
            _lanl_frame(
                accession=["AF004885", "AF004885"],
                subtype_lanl=["A1", "A1"],
                country_code=["KE", "UG"],
                collection_year=[1994, 1994],
            ),
            tmp_path / "out.parquet",
        )
        assert merged.height == 1

    def test_requires_a_sequence_column(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="requires a 'sequence' column"):
            merge_sources(
                pl.DataFrame({"accession": []}),
                _lanl_frame().drop("sequence"),
                tmp_path / "out.parquet",
            )

    def test_requires_an_accession_column(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="missing 'accession' column"):
            merge_sources(
                pl.DataFrame({"accession": []}),
                _lanl_frame().drop("accession"),
                tmp_path / "out.parquet",
            )


class TestGenBankLANLCoalescing:
    """LANL curation wins; GenBank fills the gaps."""

    def _genbank_frame(self) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "accession": ["AF004885", "MK999999"],
                "source": ["genbank", "genbank"],
                "sequence": ["ACGT" * 300, "ACGT" * 300],
                "length": [1200, 1200],
                "country": ["KE", "KE"],
                "country_raw": ["Kenya", "Kenya"],
                "country_source": ["genbank_source_feature"] * 2,
                "subtype": ["CRF01_AE", "D"],
                "subtype_source": ["genbank", "genbank"],
                "collection_year": [1994, 2019],
            }
        )

    def test_lanl_subtype_overrides_genbank(self, tmp_path: Path) -> None:
        merged = merge_sources(
            self._genbank_frame(),
            pl.DataFrame({"accession": ["AF004885"], "subtype_lanl": ["A1"]}),
            tmp_path / "out.parquet",
        ).sort("accession")
        assert merged["subtype"].to_list() == ["A1", "D"]
        assert merged["subtype_source"].to_list() == ["lanl_curated", "genbank"]

    def test_genbank_subtype_kept_when_lanl_has_none(self, tmp_path: Path) -> None:
        merged = merge_sources(
            self._genbank_frame(),
            pl.DataFrame({"accession": ["ZZ000000"], "subtype_lanl": ["A1"]}),
            tmp_path / "out.parquet",
        ).sort("accession")
        assert merged["subtype"].to_list() == ["CRF01_AE", "D"]

    def test_genbank_missing_accession_column_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="genbank_df missing 'accession' column"):
            merge_sources(
                pl.DataFrame({"sequence": ["ACGT" * 300]}),
                _lanl_frame(),
                tmp_path / "out.parquet",
            )
