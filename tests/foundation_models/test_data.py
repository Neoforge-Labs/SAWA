"""Unit tests for foundation_models.data — Paper 1 → benchmark dataset."""

from __future__ import annotations

import polars as pl
import pytest

from foundation_models.data import (
    RESISTANT_LEVELS,
    build_fm_dataset,
    stratified_train_val_test_split,
)


def _dr(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


def _curated(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


DEFAULT_DR = {"accession": "A1", "drug": "3TC", "drug_class": "NRTI", "level": "high"}
DEFAULT_COV = {
    "accession": "A1",
    "sequence": "ACGT" * 100,
    "subtype": "B",
    "country": "US",
    "collection_year": 2020,
}


def _make_toy(
    n: int = 40, subtypes: tuple[str, ...] = ("B", "A1", "C")
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Balanced toy set with `n` accessions across the given subtypes."""
    dr_rows: list[dict] = []
    cov_rows: list[dict] = []
    for i in range(n):
        acc = f"A{i:05d}"
        subtype = subtypes[i % len(subtypes)]
        level = "high" if i % 3 == 0 else "susceptible"
        dr_rows.append({**DEFAULT_DR, "accession": acc, "level": level})
        cov_rows.append(
            {**DEFAULT_COV, "accession": acc, "subtype": subtype, "sequence": "ACGT" * 100}
        )
    return _dr(dr_rows), _curated(cov_rows)


class TestConstants:
    def test_resistant_levels(self) -> None:
        assert frozenset({"high", "intermediate"}) == RESISTANT_LEVELS


class TestBuildFmDataset:
    def test_stratified_split_sums_to_input(self) -> None:
        dr, curated = _make_toy(n=60)
        ds = build_fm_dataset(dr, curated, strategy="stratified", seed=42)
        assert ds.train.height + ds.val.height + ds.test.height == 60
        assert ds.counts.n_train + ds.counts.n_val + ds.counts.n_test == 60
        assert ds.counts.strategy == "stratified"

    def test_stratified_is_deterministic(self) -> None:
        dr, curated = _make_toy(n=60)
        ds1 = build_fm_dataset(dr, curated, strategy="stratified", seed=7)
        ds2 = build_fm_dataset(dr, curated, strategy="stratified", seed=7)
        assert ds1.train["accession"].to_list() == ds2.train["accession"].to_list()
        assert ds1.test["accession"].to_list() == ds2.test["accession"].to_list()

    def test_stratified_different_seeds_differ(self) -> None:
        dr, curated = _make_toy(n=60)
        ds1 = build_fm_dataset(dr, curated, strategy="stratified", seed=1)
        ds2 = build_fm_dataset(dr, curated, strategy="stratified", seed=2)
        assert ds1.train["accession"].to_list() != ds2.train["accession"].to_list()

    def test_subtype_holdout_puts_all_non_train_subtypes_in_test(self) -> None:
        dr, curated = _make_toy(n=60, subtypes=("B", "A1", "C"))
        ds = build_fm_dataset(
            dr,
            curated,
            strategy="subtype_holdout",
            holdout_train_subtypes=("B",),
            seed=42,
        )
        train_subtypes = set(ds.train["subtype"].to_list())
        val_subtypes = set(ds.val["subtype"].to_list())
        test_subtypes = set(ds.test["subtype"].to_list())
        assert train_subtypes == {"B"}
        assert val_subtypes == {"B"}
        assert test_subtypes == {"A1", "C"}
        assert ds.counts.strategy == "subtype_holdout"

    def test_subtype_holdout_train_pool_empty_raises(self) -> None:
        dr, curated = _make_toy(n=30, subtypes=("A1", "C"))  # no B
        with pytest.raises(ValueError, match="no rows in train subtypes"):
            build_fm_dataset(
                dr,
                curated,
                strategy="subtype_holdout",
                holdout_train_subtypes=("B",),
            )

    def test_missing_columns_raise(self) -> None:
        with pytest.raises(ValueError, match="dr_results missing"):
            build_fm_dataset(pl.DataFrame({"x": [1]}), _curated([DEFAULT_COV]))
        with pytest.raises(ValueError, match="curated missing"):
            build_fm_dataset(_dr([DEFAULT_DR]), pl.DataFrame({"x": [1]}))

    def test_unknown_strategy_raises(self) -> None:
        dr, curated = _make_toy(n=10)
        with pytest.raises(ValueError, match="unknown split strategy"):
            build_fm_dataset(dr, curated, strategy="magic")  # type: ignore[arg-type]

    def test_label_binary_computed_from_level(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR, "accession": "A1", "level": "high"},
                {**DEFAULT_DR, "accession": "A2", "level": "intermediate"},
                {**DEFAULT_DR, "accession": "A3", "level": "low"},
                {**DEFAULT_DR, "accession": "A4", "level": "susceptible"},
            ]
        )
        curated = _curated(
            [{**DEFAULT_COV, "accession": f"A{i}", "subtype": "B"} for i in range(1, 5)]
        )
        ds = build_fm_dataset(dr, curated, strategy="stratified", seed=0)
        all_rows = pl.concat([ds.train, ds.val, ds.test]).sort("accession")
        labels = dict(
            zip(
                all_rows["accession"].to_list(),
                all_rows["label_binary"].to_list(),
                strict=True,
            )
        )
        assert labels == {"A1": 1, "A2": 1, "A3": 0, "A4": 0}

    def test_missing_covariates_dropped(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR, "accession": "A1"},
                {**DEFAULT_DR, "accession": "A2"},
            ]
        )
        curated = _curated(
            [
                {**DEFAULT_COV, "accession": "A1", "subtype": "B"},
                {**DEFAULT_COV, "accession": "A2", "subtype": None},
            ]
        )
        ds = build_fm_dataset(dr, curated, strategy="stratified", seed=0)
        all_accessions = set(
            ds.train["accession"].to_list()
            + ds.val["accession"].to_list()
            + ds.test["accession"].to_list()
        )
        assert all_accessions == {"A1"}


class TestSplitHelper:
    def test_bad_fractions_raise(self) -> None:
        _make_toy(n=10)
        with pytest.raises(ValueError, match="bad fractions"):
            stratified_train_val_test_split(
                pl.DataFrame({"subtype": ["A"], "drug": ["x"], "label_binary": [0]}),
                train_frac=1.5,
            )

    def test_no_room_for_test_raises(self) -> None:
        with pytest.raises(ValueError, match="leave room for test"):
            stratified_train_val_test_split(
                pl.DataFrame({"subtype": ["A"], "drug": ["x"], "label_binary": [0]}),
                train_frac=0.6,
                val_frac=0.5,
            )
