"""Row-order regression tests for the ingest merge path.

Polars' `unique()` returns rows in an unspecified order unless told
otherwise, and two of them sit between the LANL export and
`east_africa_pol.parquet`. The Parquet's md5 is the DVC hash cited in the
manuscript, so an unstable row order churns that hash on every re-run even
though the data has not changed.

See docs/evaluations/unique-row-ordering.md.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import polars as pl

from hivdrwatch.ingest.lanl import LANLIngester
from hivdrwatch.ingest.merge import merge_sources

# Enough distinct keys that Polars' parallel hash path actually reorders.
# A frame with two or three keys passes with the flag removed, purely by
# luck — measured in the evaluation doc's §2.
N_ACCESSIONS = 40
N_RUNS = 8


def _write_lanl_export(root: Path, country: str, accessions: list[str]) -> None:
    """Minimal LANL Table TSV, one row per accession."""
    folder = root / country
    folder.mkdir(parents=True, exist_ok=True)
    header = "Accession\tSE id(SA)\tSubtype\tCountry\tSampling Year\tSequence Length\tSequence"
    lines = [header]
    for i, acc in enumerate(accessions):
        lines.append(f"{acc}\tSE{i:05d}\tA1\tKenya\t2019\t1200\t{'ACGT' * 300}")
    (folder / "results.tsv").write_text("\n".join(lines) + "\n")


def _accessions() -> list[str]:
    return [f"AF{i:05d}" for i in range(N_ACCESSIONS)]


class TestLanlLoadDirectoryOrder:
    def test_duplicate_accessions_across_files_resolve_in_a_stable_order(
        self, tmp_path: Path
    ) -> None:
        """`load_directory` must return the same row order on every call.

        Two country folders carry every accession, so the dedup in
        `load_directory` has a duplicate to resolve for each one — which is
        exactly the call that reorders.
        """
        accs = _accessions()
        _write_lanl_export(tmp_path, "ke", accs)
        _write_lanl_export(tmp_path, "ug", accs)

        ingester = LANLIngester(tmp_path)
        orders = {
            tuple(ingester.load_directory()["accession"].to_list()) for _ in range(N_RUNS)
        }
        assert len(orders) == 1, (
            f"load_directory returned {len(orders)} distinct row orders in "
            f"{N_RUNS} calls; the dedup needs maintain_order=True"
        )

    def test_keeps_one_row_per_accession(self, tmp_path: Path) -> None:
        accs = _accessions()
        _write_lanl_export(tmp_path, "ke", accs)
        _write_lanl_export(tmp_path, "ug", accs)

        df = LANLIngester(tmp_path).load_directory()
        assert df.height == N_ACCESSIONS
        assert sorted(df["accession"].to_list()) == sorted(accs)


class TestMergeSourcesParquetBytes:
    def _genbank(self, accs: list[str]) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "accession": accs,
                "sequence": ["ACGT" * 300] * len(accs),
                "subtype": ["unknown"] * len(accs),
                "subtype_source": ["genbank_note"] * len(accs),
                "country": ["KE"] * len(accs),
                "source": ["genbank"] * len(accs),
            }
        )

    def _lanl(self, accs: list[str]) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "accession": accs,
                "subtype_lanl": ["A1"] * len(accs),
                "country_code": ["KE"] * len(accs),
                "collection_year": [2019] * len(accs),
            }
        )

    def test_merged_parquet_is_byte_identical_across_runs(self, tmp_path: Path) -> None:
        """The DVC hash must not move when only the row order could.

        This is the claim the manuscript rests on: same input, same bytes.
        Guards the `maintain_order=True` in `_collapse_duplicate_accessions`,
        which is the last operation before the write.
        """
        accs = _accessions()
        gb, lanl = self._genbank(accs), self._lanl(accs)

        digests = set()
        for run in range(N_RUNS):
            out = tmp_path / f"merged_{run}.parquet"
            merge_sources(gb, lanl, out)
            digests.add(hashlib.md5(out.read_bytes()).hexdigest())

        assert len(digests) == 1, (
            f"merge_sources produced {len(digests)} distinct Parquet digests in "
            f"{N_RUNS} runs on identical input; the DVC hash cited in the "
            "manuscript would churn on every re-run"
        )

    def test_merge_keeps_every_accession_once(self, tmp_path: Path) -> None:
        accs = _accessions()
        merged = merge_sources(
            self._genbank(accs), self._lanl(accs), tmp_path / "merged.parquet"
        )
        assert merged.height == N_ACCESSIONS
        assert set(merged["accession"].to_list()) == set(accs)
