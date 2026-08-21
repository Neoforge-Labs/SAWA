"""Data preparation for the Paper 3 foundation-model benchmark.

Turns SAWA's Paper 1 Parquet outputs into an `FMDataset` — one
sequence per row, one label per drug, plus the covariates the
benchmark stratifies on.

## Split contract

Two supported split strategies:

  * **`stratified`** (default): 70/15/15 train/val/test, stratified
    on `(subtype, label_per_drug)` so every subtype × drug cell is
    represented in every split. This is the fair baseline that
    classical models are expected to do well on.

  * **`subtype_holdout`**: train on subtype B (or specified pool),
    validate on remaining B, TEST on non-B subtypes only. This is
    Paper 3's thesis-testing split — the whole point is to measure
    the generalization gap. Non-B subtypes are entirely unseen
    during training.

Splits are deterministic given the seed. Written to Parquet so the
sklearn baselines, PyTorch fine-tuners, and evaluators all read the
exact same rows.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

# Same convention as causal.data_prep and paper1.fig02_aggregate.
RESISTANT_LEVELS: frozenset[str] = frozenset({"high", "intermediate"})

SplitStrategy = Literal["stratified", "subtype_holdout"]


@dataclass(frozen=True)
class SplitCounts:
    """Row counts per split for the summary log."""

    n_train: int
    n_val: int
    n_test: int
    n_unique_accessions: int
    n_drugs: int
    n_subtypes: int
    strategy: str


@dataclass(frozen=True)
class FMDataset:
    """Split dataset with everything the benchmark needs.

    * `train`, `val`, `test` — one row per (accession, drug) with:
      `accession, sequence, subtype, country, collection_year,
       drug, drug_class, label_binary`.
    * `counts` — reproducibility summary.

    Every split is a Polars DataFrame; use `.to_pandas()` or
    `.to_numpy()` at the model boundary as needed.
    """

    train: pl.DataFrame
    val: pl.DataFrame
    test: pl.DataFrame
    counts: SplitCounts


def build_fm_dataset(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    strategy: SplitStrategy = "stratified",
    train_frac: float = 0.70,
    val_frac: float = 0.15,
    holdout_train_subtypes: tuple[str, ...] = ("B",),
    seed: int = 42,
    sequence_column: str | None = None,
    strip_gaps: bool = True,
    min_sequence_length: int = 300,
) -> FMDataset:
    """Build a foundation-model dataset from SAWA Paper 1 Parquets.

    `strategy="stratified"` — proportional split stratified on
    `(subtype, drug, label)` for balanced baselines.

    `strategy="subtype_holdout"` — train + val on
    `holdout_train_subtypes` only; test on ALL other subtypes.
    Use this to reproduce Paper 3's headline generalization gap.

    Parameters
    ----------
    sequence_column
        Column in `curated` holding the sequence. Default `None` = auto-
        detect: prefer `"sequence"` (raw, from the ingest layer) then
        `"aligned_pol"` (from the curate layer). The output frame always
        exposes the sequence under the canonical name `sequence`.
    strip_gaps
        When the source column is aligned (e.g. `aligned_pol`), remove
        gap characters (`-`) before emitting. Default `True` because
        every FM tokenizer we support treats `-` as an unknown token,
        which is a no-op that just wastes context length. Set to
        `False` if a downstream evaluator needs alignment-preserving
        positions.
    min_sequence_length
        Drop rows whose gap-stripped sequence is shorter than this
        many characters. Default 300 — below that, HIV pol is missing
        so much of RT that resistance calls in the drug-class labels
        aren't well-defined from sequence alone. This is a coarse
        guard; per-gene coverage filtering happens upstream in
        `curate/qc.py`.
    """
    _require_cols(
        "dr_results",
        dr_results,
        {"accession", "drug", "drug_class", "level"},
    )
    _require_cols(
        "curated",
        curated,
        {"accession", "subtype", "country", "collection_year"},
    )

    if sequence_column is None:
        if "sequence" in curated.columns:
            sequence_column = "sequence"
        elif "aligned_pol" in curated.columns:
            sequence_column = "aligned_pol"
        else:
            raise ValueError(
                "curated has neither 'sequence' nor 'aligned_pol'; supply "
                "sequence_column= explicitly with a column name that holds "
                "the pol sequences"
            )
    elif sequence_column not in curated.columns:
        raise ValueError(
            f"sequence_column={sequence_column!r} not in curated columns: "
            f"{curated.columns}"
        )

    logger.info("using sequence column %r (strip_gaps=%s)", sequence_column, strip_gaps)

    # Per-(accession, drug) binary label from Sierra levels.
    labeled = (
        dr_results.with_columns(
            pl.col("level").is_in(list(RESISTANT_LEVELS)).cast(pl.Int8).alias("label_binary"),
        )
        .select(["accession", "drug", "drug_class", "label_binary"])
        .unique(subset=["accession", "drug"], keep="last")
    )

    seq_expr = pl.col(sequence_column)
    if strip_gaps and sequence_column != "sequence":
        # Gaps in aligned_pol shouldn't survive into the FM tokenizer.
        seq_expr = seq_expr.str.replace_all("-", "")
    seq_expr = seq_expr.alias("sequence")

    covs = (
        curated.with_columns(seq_expr)
        .select(["accession", "sequence", "subtype", "country", "collection_year"])
        .unique(subset=["accession"], keep="first")
        .filter(
            pl.col("subtype").is_not_null()
            & pl.col("sequence").is_not_null()
            & (pl.col("sequence").str.len_chars() >= min_sequence_length)
        )
    )

    joined = labeled.join(covs, on="accession", how="inner")

    if joined.is_empty():
        raise ValueError(
            "no rows survived the (accession, drug) join — check that "
            "dr_results.accession overlaps with curated.accession"
        )

    if strategy == "stratified":
        train, val, test = _stratified_split(
            joined,
            train_frac=train_frac,
            val_frac=val_frac,
            seed=seed,
        )
    elif strategy == "subtype_holdout":
        train, val, test = _subtype_holdout_split(
            joined,
            holdout_train_subtypes=holdout_train_subtypes,
            train_frac_of_pool=train_frac / (train_frac + val_frac),
            seed=seed,
        )
    else:
        raise ValueError(f"unknown split strategy: {strategy!r}")

    counts = SplitCounts(
        n_train=train.height,
        n_val=val.height,
        n_test=test.height,
        n_unique_accessions=joined["accession"].n_unique(),
        n_drugs=joined["drug"].n_unique(),
        n_subtypes=joined["subtype"].n_unique(),
        strategy=strategy,
    )
    logger.info(
        "FM dataset (%s): train=%d val=%d test=%d — %d accessions × %d drugs × %d subtypes",
        counts.strategy,
        counts.n_train,
        counts.n_val,
        counts.n_test,
        counts.n_unique_accessions,
        counts.n_drugs,
        counts.n_subtypes,
    )
    return FMDataset(train=train, val=val, test=test, counts=counts)


def stratified_train_val_test_split(
    df: pl.DataFrame,
    *,
    stratify_cols: tuple[str, ...] = ("subtype", "drug", "label_binary"),
    train_frac: float = 0.70,
    val_frac: float = 0.15,
    seed: int = 42,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Reusable helper — stratified split on arbitrary column tuple.

    Deterministic given seed. Every stratum is split independently so
    small strata still get at least one row per split when possible.
    """
    return _stratified_split(
        df, train_frac=train_frac, val_frac=val_frac, seed=seed, stratify_cols=stratify_cols
    )


# ------------------------------------------------------------------
# internals
# ------------------------------------------------------------------


def _stratified_split(
    df: pl.DataFrame,
    *,
    train_frac: float,
    val_frac: float,
    seed: int,
    stratify_cols: tuple[str, ...] = ("subtype", "drug", "label_binary"),
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    if not 0 < train_frac < 1 or not 0 <= val_frac < 1:
        raise ValueError(f"bad fractions: train={train_frac}, val={val_frac}")
    if train_frac + val_frac >= 1.0:
        raise ValueError("train_frac + val_frac must leave room for test")

    rng = np.random.default_rng(seed)

    train_rows: list[pl.DataFrame] = []
    val_rows: list[pl.DataFrame] = []
    test_rows: list[pl.DataFrame] = []

    # Iterate strata and rows in a fixed order so the RNG produces
    # identical splits across runs — Polars' group_by iteration order
    # is not guaranteed stable.
    stratify_list = list(stratify_cols)
    sort_cols = stratify_list + (["accession"] if "accession" in df.columns else [])
    df_sorted = df.sort(sort_cols)
    for key, group in df_sorted.group_by(stratify_list, maintain_order=True):
        del key
        n = group.height
        if n == 0:
            continue
        perm = rng.permutation(n)
        n_train = int(n * train_frac)
        n_val = int(n * val_frac)
        # Guarantee test set gets whatever's left, at least 0 rows.
        train_ix = perm[:n_train]
        val_ix = perm[n_train : n_train + n_val]
        test_ix = perm[n_train + n_val :]
        # Polars gather via Series index.
        train_rows.append(group[train_ix.tolist()])
        val_rows.append(group[val_ix.tolist()])
        test_rows.append(group[test_ix.tolist()])

    train = pl.concat(train_rows) if train_rows else df.head(0)
    val = pl.concat(val_rows) if val_rows else df.head(0)
    test = pl.concat(test_rows) if test_rows else df.head(0)
    return train, val, test


def _subtype_holdout_split(
    df: pl.DataFrame,
    *,
    holdout_train_subtypes: tuple[str, ...],
    train_frac_of_pool: float,
    seed: int,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Train + val drawn from holdout_train_subtypes; test = all others."""
    train_pool = df.filter(pl.col("subtype").is_in(list(holdout_train_subtypes)))
    test = df.filter(~pl.col("subtype").is_in(list(holdout_train_subtypes)))

    if train_pool.is_empty():
        raise ValueError(
            f"subtype_holdout: no rows in train subtypes {holdout_train_subtypes}. "
            f"Available: {sorted(df['subtype'].unique().to_list())}"
        )

    rng = np.random.default_rng(seed)
    perm = rng.permutation(train_pool.height)
    n_train = int(train_pool.height * train_frac_of_pool)
    train_ix = perm[:n_train]
    val_ix = perm[n_train:]
    train = train_pool[train_ix.tolist()]
    val = train_pool[val_ix.tolist()]

    logger.info(
        "subtype_holdout: train pool=%s (train=%d, val=%d), test subtypes=%s (n=%d)",
        list(holdout_train_subtypes),
        train.height,
        val.height,
        sorted(test["subtype"].unique().to_list()) if test.height else [],
        test.height,
    )
    return train, val, test


def _require_cols(name: str, df: pl.DataFrame, cols: set[str]) -> None:
    missing = cols - set(df.columns)
    if missing:
        raise ValueError(f"{name} missing required columns: {sorted(missing)}")


__all__ = [
    "RESISTANT_LEVELS",
    "FMDataset",
    "SplitCounts",
    "SplitStrategy",
    "build_fm_dataset",
    "stratified_train_val_test_split",
]
