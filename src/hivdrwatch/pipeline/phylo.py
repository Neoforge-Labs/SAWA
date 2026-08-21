"""Prefect flow: curated alignment parquet → per-subtype phylogenies.

Inputs:
    data/interim/curated/aligned.parquet  (from pipeline/curate.py)

Outputs:
    data/interim/trees/<subtype>/         (raw IQ-TREE + TreeTime artifacts)
    data/gold/trees/manifest.parquet      (TreeMetadata index)

Run locally:
    make restart-worker
    python -m hivdrwatch.pipeline.phylo \\
        --curated data/interim/curated/aligned.parquet \\
        --trees-root data/interim/trees \\
        --manifest data/gold/trees/manifest.parquet
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl
from prefect import flow, get_run_logger, task
from prefect.tasks import task_input_hash

from ..phylo import (
    MIN_TAXA_FOR_TREE,
    build_all_trees,
    write_manifest,
)
from ..phylo.build import MIN_ALIGNED_COVERAGE
from ..tracking import tracked_run


@task(cache_key_fn=task_input_hash)
def load_curated(path: Path) -> pl.DataFrame:
    log = get_run_logger()
    df = pl.read_parquet(path)
    log.info("loaded %d curated rows from %s", df.height, path)
    return df


@task(retries=1)
def build_trees(
    curated: pl.DataFrame,
    trees_root: Path,
    seed: int,
    bootstrap: int,
    fast: bool,
    skip_treetime: bool,
    subtypes: list[str] | None,
    min_taxa: int,
    min_coverage: float,
    partition_report_path: Path,
) -> list:
    log = get_run_logger()
    # Pass the Prefect run logger down so IQ-TREE/TreeTime subprocess
    # output streams into the Prefect UI in real time (instead of only
    # showing up via stdlib logging, which Prefect doesn't forward).
    meta = build_all_trees(
        curated,
        trees_root,
        subtypes=subtypes,
        seed=seed,
        bootstrap=bootstrap,
        fast=fast,
        skip_treetime=skip_treetime,
        min_taxa=min_taxa,
        min_coverage=min_coverage,
        stream_logger=log,
        partition_report_path=partition_report_path,
    )
    log.info(
        "built %d tree artifacts; partition audit trail → %s",
        len(meta),
        partition_report_path,
    )
    return meta


@task
def write_and_track(rows: list, manifest_path: Path) -> Path:
    log = get_run_logger()
    out = write_manifest(rows, manifest_path)
    log.info("manifest → %s", out)

    with tracked_run(
        experiment="phylo",
        run_name=f"phylo-{len(rows)}-trees",
        tags={"flow": "phylo"},
    ) as run:
        run.log_metrics(
            {
                "n_tree_rows": float(len(rows)),
                "n_subtypes": float(len({r.subtype for r in rows})),
            }
        )
        for r in rows:
            run.log_metrics(
                {
                    f"n_taxa_{r.subtype}_{r.kind.value}": float(r.n_taxa),
                    f"wall_seconds_{r.subtype}_{r.kind.value}": float(r.wall_seconds or 0),
                }
            )
        if out.exists():
            run.log_artifact(out)
    return out


@flow(name="phylo")
def phylo_flow(
    curated_parquet: Path,
    trees_root: Path,
    manifest_path: Path,
    *,
    seed: int = 42,
    bootstrap: int = 0,
    fast: bool = True,
    skip_treetime: bool = False,
    subtypes: list[str] | None = None,
    min_taxa: int = MIN_TAXA_FOR_TREE,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
    partition_report_path: Path | None = None,
) -> Path:
    curated = load_curated(curated_parquet)
    # Default the partition report next to the manifest — one gold-tier
    # artifact directory, one audit trail.
    report_path = partition_report_path or (manifest_path.parent / "partition_report.parquet")
    rows = build_trees(
        curated,
        trees_root,
        seed=seed,
        bootstrap=bootstrap,
        fast=fast,
        skip_treetime=skip_treetime,
        subtypes=subtypes,
        min_taxa=min_taxa,
        min_coverage=min_coverage,
        partition_report_path=report_path,
    )
    return write_and_track(rows, manifest_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.pipeline.phylo")
    parser.add_argument(
        "--curated",
        type=Path,
        default=Path("data/interim/curated/aligned.parquet"),
    )
    parser.add_argument(
        "--trees-root",
        type=Path,
        default=Path("data/interim/trees"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/gold/trees/manifest.parquet"),
    )
    parser.add_argument(
        "--partition-report",
        type=Path,
        default=None,
        help=(
            "audit trail: per-subtype accounting parquet "
            "(default: <manifest-dir>/partition_report.parquet)"
        ),
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help="IQ-TREE ultrafast-bootstrap replicates (0=off, 1000 recommended for paper)",
    )
    parser.add_argument(
        "--skip-treetime",
        action="store_true",
        help="only build ML trees; skip the time-scaled TreeTime step",
    )
    parser.add_argument(
        "--subtypes",
        nargs="+",
        default=None,
        help="restrict to a subset of subtypes (default: all with >=min_taxa)",
    )
    parser.add_argument(
        "--min-taxa",
        type=int,
        default=MIN_TAXA_FOR_TREE,
        help=f"minimum taxa per subtype to build a tree (default {MIN_TAXA_FOR_TREE})",
    )
    parser.add_argument(
        "--no-fast",
        action="store_true",
        help=(
            "disable IQ-TREE 2's -fast mode. Default is fast (approximate "
            "ML search, ~5-10x speedup, quality suitable for interim / "
            "phylodynamics work). Turn off for final paper-figure runs."
        ),
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=MIN_ALIGNED_COVERAGE,
        help=(
            f"drop sequences with less than this fraction of non-gap "
            f"positions in aligned_pol before tree inference "
            f"(default {MIN_ALIGNED_COVERAGE:.0%}). Lower this if you want "
            f"integrase-only sequences included; use per-gene trees for "
            f"the scientifically cleaner alternative."
        ),
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    phylo_flow(
        curated_parquet=args.curated,
        trees_root=args.trees_root,
        manifest_path=args.manifest,
        seed=args.seed,
        bootstrap=args.bootstrap,
        fast=not args.no_fast,
        skip_treetime=args.skip_treetime,
        subtypes=args.subtypes,
        min_taxa=args.min_taxa,
        min_coverage=args.min_coverage,
        partition_report_path=args.partition_report,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
