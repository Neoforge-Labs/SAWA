"""Tests for Auspice v2 export."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest


@pytest.fixture
def toy_tree_newick(tmp_path: Path) -> Path:
    """Minimal Newick — 3 tips A/B/C, one internal node."""
    p = tmp_path / "toy.nwk"
    p.write_text("(A:0.1,B:0.2,C:0.15):0.0;\n")
    return p


@pytest.fixture
def toy_curated(tmp_path: Path) -> Path:
    p = tmp_path / "curated.parquet"
    pl.DataFrame(
        {
            "accession": ["A", "B", "C"],
            "country": ["KE", "UG", "TZ"],
            "subtype": ["A1", "A1", "D"],  # C is a different subtype
            "collection_year": [2018, 2021, 2016],
        }
    ).write_parquet(p)
    return p


@pytest.fixture
def toy_dr(tmp_path: Path) -> Path:
    p = tmp_path / "dr.parquet"
    pl.DataFrame(
        {
            "accession": ["A", "A", "B", "C"],
            "drug": ["3TC", "DTG", "3TC", "DTG"],
            "drug_class": ["NRTI", "INSTI", "NRTI", "INSTI"],
            "level": ["high", "susceptible", "susceptible", "intermediate"],
        }
    ).write_parquet(p)
    return p


class TestExportSubtypeDataset:
    def test_writes_valid_auspice_v2_json(
        self,
        toy_tree_newick: Path,
        toy_curated: Path,
        toy_dr: Path,
        tmp_path: Path,
    ) -> None:
        from hivdrwatch.phylo.auspice import export_subtype_dataset

        result = export_subtype_dataset(
            subtype="A1",
            tree_newick=toy_tree_newick,
            curated_parquet=toy_curated,
            dr_parquet=toy_dr,
            out_dir=tmp_path / "auspice",
            sawa_version="v0.1.0-test",
            git_sha="deadbeef",
        )

        assert result.n_tips_in_tree == 3
        out = Path(result.output_path)
        assert out.exists()
        payload = json.loads(out.read_text())
        assert payload["version"] == 2
        assert "tree" in payload
        assert "meta" in payload
        assert payload["meta"]["title"].startswith("SAWA HIV-1 subtype A1")
        assert payload["meta"]["extensions"]["sawa"]["git_sha"] == "deadbeef"

    def test_resistant_flag_from_any_class(
        self,
        toy_tree_newick: Path,
        toy_curated: Path,
        toy_dr: Path,
        tmp_path: Path,
    ) -> None:
        from hivdrwatch.phylo.auspice import export_subtype_dataset

        result = export_subtype_dataset(
            subtype="A1",
            tree_newick=toy_tree_newick,
            curated_parquet=toy_curated,
            dr_parquet=toy_dr,
            out_dir=tmp_path / "auspice",
            sawa_version="v0.1.0-test",
            git_sha="deadbeef",
        )
        # A has high NRTI resistance; C has intermediate INSTI
        assert result.n_resistant == 2

    def test_resistant_flag_scoped_to_drug_class(
        self,
        toy_tree_newick: Path,
        toy_curated: Path,
        toy_dr: Path,
        tmp_path: Path,
    ) -> None:
        from hivdrwatch.phylo.auspice import export_subtype_dataset

        result = export_subtype_dataset(
            subtype="A1",
            tree_newick=toy_tree_newick,
            curated_parquet=toy_curated,
            dr_parquet=toy_dr,
            out_dir=tmp_path / "auspice",
            sawa_version="v0.1.0-test",
            git_sha="deadbeef",
            resistance_drug_class="INSTI",
        )
        # Only C carries INSTI resistance in the fixture
        assert result.n_resistant == 1

    def test_missing_subtype_metadata_still_annotates_tip(
        self,
        toy_tree_newick: Path,
        toy_curated: Path,
        toy_dr: Path,
        tmp_path: Path,
    ) -> None:
        """Tip C has subtype D in curated but we export as A1 — shouldn't drop it."""
        from hivdrwatch.phylo.auspice import export_subtype_dataset

        result = export_subtype_dataset(
            subtype="A1",
            tree_newick=toy_tree_newick,
            curated_parquet=toy_curated,
            dr_parquet=toy_dr,
            out_dir=tmp_path / "auspice",
            sawa_version="v0.1.0-test",
            git_sha="deadbeef",
        )
        assert result.n_tips_in_tree == 3
        # All three tips get metadata annotation (via all_metadata fallback)
        assert result.n_tips_with_metadata >= 2


class TestSlug:
    def test_alphanumeric_preserved(self) -> None:
        from hivdrwatch.phylo.auspice import _slug

        assert _slug("A1") == "A1"
        assert _slug("CRF02_AG") == "CRF02_AG"

    def test_special_chars_replaced(self) -> None:
        from hivdrwatch.phylo.auspice import _slug

        # Consecutive non-alnum runs collapse into consecutive underscores;
        # leading/trailing underscores are stripped. "A1 (4.17%)" → "A1__4_17"
        # (space and paren both become underscores, . and % become underscores).
        assert _slug("A1 (4.17%)") == "A1__4_17"
        # Empty / all-junk should still produce a stable slug
        assert _slug("") == "unknown"


class TestStripAccession:
    def test_lanl_dot_delimited_header_returns_position_7(self) -> None:
        """Real LANL FASTA header format:
        Subtype.Country.Year.Name.PatientID.CDID.SEID.Accession.Length.Start.End.HIV-1
        """
        from hivdrwatch.phylo.auspice import _strip_accession

        assert (
            _strip_accession("A1.RW.1992.92RW008.RW008.102981.606.AB253421.9632.1.9643.HIV-1")
            == "AB253421"
        )
        assert (
            _strip_accession("D.BI.2002.U0106.U0106.84925.45071.AM260218.3843.2085.1783.HIV-1")
            == "AM260218"
        )
        assert (
            _strip_accession("C.ET.2003.ETH_G_001.-.68590.-.AB285752.3347.2253.1095.HIV-1")
            == "AB285752"
        )

    def test_pipe_delimited_fallback(self) -> None:
        """Older SAWA / test-fixture format kept alive as fallback."""
        from hivdrwatch.phylo.auspice import _strip_accession

        assert _strip_accession("AB253421|RW|1992") == "AB253421"

    def test_bare_accession_returned_unchanged(self) -> None:
        from hivdrwatch.phylo.auspice import _strip_accession

        assert _strip_accession("AB253421") == "AB253421"

    def test_empty_returns_empty(self) -> None:
        from hivdrwatch.phylo.auspice import _strip_accession

        assert _strip_accession("") == ""
