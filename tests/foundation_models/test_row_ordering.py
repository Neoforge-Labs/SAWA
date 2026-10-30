"""Row-order regression tests for Paper 3's train/val/test splits.

`_subtype_holdout_split` gathers rows by position from a seeded
permutation. A fixed seed only reproduces a split if the positions mean
the same thing on every run — and two `unique()` calls plus an inner join
sit between the Parquet and the pool, none of which promises an order.
Without the sort in `build_fm_dataset`, a run puts *different sequences*
in train and val while reporting the same seed.

See docs/evaluations/unique-row-ordering.md.
"""

from __future__ import annotations

import polars as pl

from foundation_models.data import build_fm_dataset

N_ACCESSIONS = 60
N_RUNS = 8
SUBTYPES = ("B", "A1", "C", "D")


def _toy() -> tuple[pl.DataFrame, pl.DataFrame]:
    """A cohort big enough that Polars' parallel hash path actually reorders.

    Each accession carries two drugs, so the `(accession, drug)` dedup has a
    duplicate to resolve for every row.
    """
    dr_rows: list[dict] = []
    cov_rows: list[dict] = []
    for i in range(N_ACCESSIONS):
        acc = f"A{i:05d}"
        for drug in ("3TC", "AZT"):
            dr_rows.append(
                {
                    "accession": acc,
                    "drug": drug,
                    "drug_class": "NRTI",
                    "level": "high" if i % 3 == 0 else "susceptible",
                }
            )
            # A second, redundant Sierra row per (accession, drug) — the
            # shape the keep="last" dedup exists for.
            dr_rows.append(
                {
                    "accession": acc,
                    "drug": drug,
                    "drug_class": "NRTI",
                    "level": "high" if i % 3 == 0 else "susceptible",
                }
            )
        cov_rows.append(
            {
                "accession": acc,
                "sequence": "ACGT" * 100,
                "subtype": SUBTYPES[i % len(SUBTYPES)],
                "country": "KE",
                "collection_year": 2019,
            }
        )
    return pl.from_dicts(dr_rows), pl.from_dicts(cov_rows)


def _membership(df: pl.DataFrame) -> frozenset[tuple[str, str]]:
    return frozenset(zip(df["accession"].to_list(), df["drug"].to_list(), strict=True))


def _row_order(df: pl.DataFrame) -> tuple[tuple[str, str], ...]:
    return tuple(zip(df["accession"].to_list(), df["drug"].to_list(), strict=True))


class TestSubtypeHoldoutSplitDeterminism:
    """Paper 3's thesis split — train on B, test on non-B."""

    def _build(self):
        dr, curated = _toy()
        return build_fm_dataset(
            dr,
            curated,
            strategy="subtype_holdout",
            holdout_train_subtypes=("B",),
            seed=42,
        )

    def test_same_seed_selects_the_same_train_rows(self) -> None:
        """The headline claim: a seed is only a seed if the order holds."""
        memberships = {_membership(self._build().train) for _ in range(N_RUNS)}
        assert len(memberships) == 1, (
            f"seed=42 produced {len(memberships)} distinct train sets in "
            f"{N_RUNS} runs; the split is drawing by position from a pool "
            "whose row order moves"
        )

    def test_same_seed_selects_the_same_val_rows(self) -> None:
        memberships = {_membership(self._build().val) for _ in range(N_RUNS)}
        assert len(memberships) == 1

    def test_train_and_val_row_order_is_stable(self) -> None:
        """Row order decides the Parquet's bytes, so DVC sees it too."""
        orders = {
            (_row_order(d.train), _row_order(d.val))
            for d in (self._build() for _ in range(N_RUNS))
        }
        assert len(orders) == 1

    def test_test_arm_row_order_is_stable(self) -> None:
        orders = {_row_order(self._build().test) for _ in range(N_RUNS)}
        assert len(orders) == 1

    def test_train_and_test_subtypes_are_still_disjoint(self) -> None:
        ds = self._build()
        assert set(ds.train["subtype"].to_list()) == {"B"}
        assert "B" not in set(ds.test["subtype"].to_list())


class TestStratifiedSplitDeterminism:
    def _build(self):
        dr, curated = _toy()
        return build_fm_dataset(dr, curated, strategy="stratified", seed=42)

    def test_same_seed_selects_the_same_train_rows(self) -> None:
        memberships = {_membership(self._build().train) for _ in range(N_RUNS)}
        assert len(memberships) == 1

    def test_split_row_order_is_stable(self) -> None:
        orders = {
            (_row_order(d.train), _row_order(d.val), _row_order(d.test))
            for d in (self._build() for _ in range(N_RUNS))
        }
        assert len(orders) == 1

    def test_splits_partition_the_cohort(self) -> None:
        ds = self._build()
        total = ds.train.height + ds.val.height + ds.test.height
        assert total == N_ACCESSIONS * 2  # one row per (accession, drug)
