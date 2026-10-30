"""Measure how much Polars' `unique()` moves row order, and where that escapes.

Exploratory, and deliberately not part of the package. Its numbers are the
ones quoted in `docs/evaluations/unique-row-ordering.md`.

Two sections. `polars` measures the primitive on synthetic frames — how
many distinct row orders a repeated call produces, whether the *retained*
row changes with it, and how a join's output order depends on its
operands. `sites` measures the repo's own functions, running each a number
of times and counting distinct outputs.

    python scripts/unique_ordering_probe.py polars --runs 12
    python scripts/unique_ordering_probe.py sites --runs 8
"""

from __future__ import annotations

import argparse
import hashlib
import io
import logging
import tempfile
from collections.abc import Callable
from functools import partial
from pathlib import Path
from typing import Any, Literal

import numpy as np
import polars as pl

logger = logging.getLogger("unique_ordering_probe")

DEFAULT_RUNS = 10


def _frame(n_keys: int, reps: int) -> pl.DataFrame:
    """`n_keys` distinct accessions, each repeated `reps` times.

    Every duplicate carries a distinguishable payload, so a caller can tell
    *which* row survived a `keep=` as well as where it landed.
    """
    return pl.DataFrame(
        {
            "accession": [f"ACC{i % n_keys:04d}" for i in range(n_keys * reps)],
            "payload": [f"p{i:05d}" for i in range(n_keys * reps)],
        }
    )


def _distinct(fn: Callable[[], Any], runs: int) -> int:
    return len({fn() for _ in range(runs)})


def _distinct_orders(fn: Callable[[], pl.DataFrame], runs: int) -> int:
    """`_distinct` over `_order`, as its own function so the loops below
    can bind the frame-producing callable without an inline lambda."""
    return len({_order(fn()) for _ in range(runs)})


def _distinct_retained(fn: Callable[[], pl.DataFrame], runs: int) -> int:
    return len({_retained(fn()) for _ in range(runs)})


def _order(df: pl.DataFrame) -> tuple[str, ...]:
    return tuple(df["accession"].to_list())


def _retained(df: pl.DataFrame) -> frozenset[tuple[str, str]]:
    """Which rows survived, ignoring where they landed."""
    return frozenset(zip(df["accession"].to_list(), df["payload"].to_list(), strict=True))


# ------------------------------------------------------------------
# section 1 — the primitive
# ------------------------------------------------------------------


def probe_polars(runs: int) -> None:
    small, big, tiny = _frame(6, 2), _frame(40, 3), _frame(3, 2)

    print(f"polars {pl.__version__}, {runs} runs per row\n")
    print("-- row order --")
    rows: list[tuple[str, Callable[[], pl.DataFrame]]] = [
        ("12 rows / 6 keys   df.unique()", lambda: small.unique()),
        ("12 rows / 6 keys   + maintain_order", lambda: small.unique(maintain_order=True)),
        (
            "120 rows / 40 keys unique(subset, keep=last)",
            lambda: big.unique(subset=["accession"], keep="last"),
        ),
        (
            "120 rows / 40 keys unique(subset, keep=first)",
            lambda: big.unique(subset=["accession"], keep="first"),
        ),
        (
            "120 rows / 40 keys + maintain_order",
            lambda: big.unique(subset=["accession"], keep="last", maintain_order=True),
        ),
        (
            "6 rows / 3 keys    unique(subset, keep=last)",
            lambda: tiny.unique(subset=["accession"], keep="last"),
        ),
    ]
    for label, fn in rows:
        print(f"  {label:<48} {_distinct_orders(fn, runs):>3} distinct order(s)")

    print("\n-- which row is retained (the negative result) --")
    keeps: tuple[tuple[str, Literal["first", "last"]], ...] = (
        ("keep=last", "last"),
        ("keep=first", "first"),
    )
    for label, keep in keeps:
        n = _distinct_retained(partial(big.unique, subset=["accession"], keep=keep), runs)
        print(f"  120 rows / 40 keys {label:<30} {n:>3} distinct retained-value set(s)")

    print("\n-- Series.unique() --")
    print(
        f"  {'no flag':<48} "
        f"{_distinct(lambda: tuple(big['accession'].unique().to_list()), runs):>3} distinct order(s)"
    )
    print(
        f"  {'maintain_order=True':<48} "
        f"{_distinct(lambda: tuple(big['accession'].unique(maintain_order=True).to_list()), runs):>3} distinct order(s)"
    )

    print("\n-- group_by, the sibling defect --")
    print(
        f"  {'group_by(..).agg(..)':<48} "
        f"{_distinct(lambda: _order(big.group_by('accession').agg(pl.col('payload').max())), runs):>3} distinct order(s)"
    )
    print(
        f"  {'group_by(maintain_order=True)':<48} "
        f"{_distinct(lambda: _order(big.group_by('accession', maintain_order=True).agg(pl.col('payload').max())), runs):>3} distinct order(s)"
    )

    print("\n-- what a join's output order follows --")
    n = 40
    ordered_l = pl.DataFrame(
        {"accession": [f"ACC{i:04d}" for i in range(n)], "x": list(range(n))}
    )
    ordered_r = pl.DataFrame(
        {"accession": [f"ACC{i:04d}" for i in range(n)], "y": list(range(n))}
    )
    scrambled = lambda: big.unique(subset=["accession"])  # noqa: E731
    joins: list[tuple[str, Callable[[], pl.DataFrame]]] = [
        (
            "ordered L x ordered R, inner",
            lambda: ordered_l.join(ordered_r, on="accession", how="inner"),
        ),
        (
            "ordered L x ordered R, left",
            lambda: ordered_l.join(ordered_r, on="accession", how="left"),
        ),
        (
            "ordered L x scrambled R, inner",
            lambda: ordered_l.join(scrambled(), on="accession", how="inner"),
        ),
        (
            "ordered L x scrambled R, left",
            lambda: ordered_l.join(scrambled(), on="accession", how="left"),
        ),
        (
            "scrambled L x ordered R, inner",
            lambda: scrambled().join(ordered_r, on="accession", how="inner"),
        ),
    ]
    for label, fn in joins:
        print(f"  {label:<48} {_distinct_orders(fn, runs):>3} distinct order(s)")

    print("\n-- Parquet bytes, which is what DVC hashes --")
    for label, flag in (("no flag", False), ("maintain_order=True", True)):

        def _md5(flag: bool = flag) -> str:
            buf = io.BytesIO()
            big.unique(subset=["accession"], keep="last", maintain_order=flag).write_parquet(
                buf, compression="zstd", compression_level=6
            )
            return hashlib.md5(buf.getvalue()).hexdigest()

        print(f"  {label:<48} {_distinct(_md5, runs):>3} distinct md5(s)")

    print("\n-- a seeded split over a pool whose order moves --")
    for label, flag in (("no flag", False), ("maintain_order=True", True)):

        def _train(flag: bool = flag) -> frozenset[str]:
            pool = big.unique(subset=["accession"], keep="last", maintain_order=flag)
            perm = np.random.default_rng(42).permutation(pool.height)
            return frozenset(pool[perm[: pool.height // 2].tolist()]["accession"].to_list())

        print(f"  seed=42, {label:<39} {_distinct(_train, runs):>3} distinct train set(s)")


# ------------------------------------------------------------------
# section 2 — the repo's own call sites
# ------------------------------------------------------------------


def probe_sites(runs: int) -> None:
    from causal.data_prep import build_dr_causal_frames
    from foundation_models.data import build_fm_dataset
    from hivdrwatch.ingest.lanl import LANLIngester
    from hivdrwatch.ingest.merge import merge_sources
    from hivdrwatch.phylo.build import partition_by_subtype

    print(f"polars {pl.__version__}, {runs} runs per site\n")
    subtypes = ("A1", "C", "D", "B", "G", "A2", "F1", "J", "H", "CRF01_AE", "CRF02_AG")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)

        # ingest/lanl.py — two country exports carrying the same accessions
        accs = [f"AF{i:05d}" for i in range(40)]
        header = (
            "Accession\tSE id(SA)\tSubtype\tCountry\tSampling Year\tSequence Length\tSequence"
        )
        for cc in ("ke", "ug"):
            (root / cc).mkdir(parents=True, exist_ok=True)
            body = [header] + [
                f"{a}\tSE{i:05d}\tA1\tKenya\t2019\t1200\t{'ACGT' * 300}"
                for i, a in enumerate(accs)
            ]
            (root / cc / "results.tsv").write_text("\n".join(body) + "\n")
        ingester = LANLIngester(root)
        n = _distinct(lambda: _order(ingester.load_directory()), runs)
        print(f"  {'ingest/lanl.py  load_directory row order':<52} {n:>3} distinct order(s)")

        # ingest/merge.py — the Parquet md5 DVC would hash
        gb = pl.DataFrame(
            {
                "accession": accs,
                "sequence": ["ACGT" * 300] * len(accs),
                "subtype": ["unknown"] * len(accs),
                "subtype_source": ["genbank_note"] * len(accs),
                "country": ["KE"] * len(accs),
                "source": ["genbank"] * len(accs),
            }
        )
        lanl = pl.DataFrame(
            {
                "accession": accs,
                "subtype_lanl": ["A1"] * len(accs),
                "country_code": ["KE"] * len(accs),
                "collection_year": [2019] * len(accs),
            }
        )

        def _merged_md5() -> str:
            out = root / "merged.parquet"
            merge_sources(gb, lanl, out)
            return hashlib.md5(out.read_bytes()).hexdigest()

        print(
            f"  {'ingest/merge.py merged Parquet md5':<52} "
            f"{_distinct(_merged_md5, runs):>3} distinct md5(s)"
        )

    # phylo/build.py — the order the tree manifest is appended in
    curated = pl.from_dicts(
        [
            {"accession": f"{s}_{i:02d}", "subtype": s, "passes_qc": True}
            for s in subtypes
            for i in range(4)
        ]
    )
    print(
        f"  {'phylo/build.py  partition key order':<52} "
        f"{_distinct(lambda: tuple(partition_by_subtype(curated, min_taxa=1)), runs):>3} distinct order(s)"
    )

    # causal/data_prep.py — the frame _bootstrap_ci resamples by position
    dr = pl.from_dicts(
        [
            {
                "accession": f"A{i:05d}",
                "gene": "RT",
                "mutations": ["M184V"],
                "drug_class": "NRTI",
                "drug": "3TC",
                "level": "high" if i % 3 == 0 else "susceptible",
            }
            for i in range(60)
        ]
    )
    cov = pl.from_dicts(
        [
            {
                "accession": f"A{i:05d}",
                "subtype": subtypes[i % len(subtypes)],
                "country": "KE",
                "collection_year": 2019,
            }
            for i in range(60)
            for _ in range(2)
        ]
    )
    print(
        f"  {'causal/data_prep.py  frame row order':<52} "
        f"{_distinct(lambda: _order(build_dr_causal_frames(dr, cov)[0]['NRTI']), runs):>3} distinct order(s)"
    )

    # foundation_models/data.py — Paper 3's holdout split membership
    fm_dr = pl.from_dicts(
        [
            {
                "accession": f"A{i:05d}",
                "drug": d,
                "drug_class": "NRTI",
                "level": "high" if i % 3 == 0 else "susceptible",
            }
            for i in range(60)
            for d in ("3TC", "AZT")
            for _ in range(2)
        ]
    )
    fm_cov = pl.from_dicts(
        [
            {
                "accession": f"A{i:05d}",
                "sequence": "ACGT" * 100,
                "subtype": ("B", "A1", "C", "D")[i % 4],
                "country": "KE",
                "collection_year": 2019,
            }
            for i in range(60)
        ]
    )

    def _train_set() -> frozenset[tuple[str, str]]:
        ds = build_fm_dataset(
            fm_dr,
            fm_cov,
            strategy="subtype_holdout",
            holdout_train_subtypes=("B",),
            seed=42,
        )
        return frozenset(
            zip(ds.train["accession"].to_list(), ds.train["drug"].to_list(), strict=True)
        )

    print(
        f"  {'foundation_models/data.py  seed=42 train set':<52} "
        f"{_distinct(_train_set, runs):>3} distinct membership(s)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("section", choices=["polars", "sites"])
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    parser.add_argument("--log-level", default="WARNING")
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s %(name)s: %(message)s")

    if args.section == "polars":
        probe_polars(args.runs)
    else:
        probe_sites(args.runs)


if __name__ == "__main__":
    main()
