"""Row-order regression test for the subtype partition.

`build_subtype_trees` extends one manifest list in the iteration order of
`partition_by_subtype`'s dict, so an unordered `unique()` upstream gives
`tree_manifest.parquet` a different row order on every run.

See docs/evaluations/unique-row-ordering.md.
"""

from __future__ import annotations

import polars as pl

from hivdrwatch.phylo.build import EXCLUDED_SUBTYPES, partition_by_subtype

# Enough distinct labels that Polars' parallel hash path actually reorders;
# a two- or three-label frame passes with the fix removed, by luck.
SUBTYPES = (
    "A1",
    "C",
    "D",
    "B",
    "G",
    "A2",
    "F1",
    "J",
    "H",
    "CRF01_AE",
    "CRF02_AG",
    "recombinant",
)
N_PER_SUBTYPE = 4
N_RUNS = 8


def _curated() -> pl.DataFrame:
    rows = [
        {"accession": f"{s}_{i:02d}", "subtype": s, "passes_qc": True}
        for s in SUBTYPES
        for i in range(N_PER_SUBTYPE)
    ]
    return pl.from_dicts(rows)


class TestPartitionBySubtypeOrder:
    def test_partition_keys_are_in_a_stable_order(self) -> None:
        df = _curated()
        orders = {tuple(partition_by_subtype(df, min_taxa=1).keys()) for _ in range(N_RUNS)}
        assert len(orders) == 1, (
            f"partition_by_subtype returned {len(orders)} distinct key orders in "
            f"{N_RUNS} calls; the tree manifest's row order follows this"
        )

    def test_partition_keys_are_sorted(self) -> None:
        keys = list(partition_by_subtype(_curated(), min_taxa=1).keys())
        assert keys == sorted(keys)

    def test_every_eligible_subtype_is_still_present(self) -> None:
        parts = partition_by_subtype(_curated(), min_taxa=1)
        assert set(parts) == set(SUBTYPES) - EXCLUDED_SUBTYPES
        assert all(p.height == N_PER_SUBTYPE for p in parts.values())
