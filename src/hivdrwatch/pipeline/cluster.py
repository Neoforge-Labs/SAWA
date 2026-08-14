"""Prefect flow: aligned FASTA → TN93 distances → cluster assignments +
outbreak features parquet.

Inputs:
    data/interim/curate/aligned.fasta         (from flows/curate.py)
    data/interim/curated/aligned.parquet      (metadata)

Outputs:
    data/interim/cluster/distances.csv        (below-threshold TN93 pairs)
    data/gold/cluster_assignments.parquet     (accession → cluster_id)
    data/gold/outbreak_features.parquet       (country × year rollups)

Run locally:
    prefect server start &
    python -m hivdrwatch.pipeline.cluster \\
        --aligned-fasta data/interim/curate/aligned.fasta \\
        --curated data/interim/curated/aligned.parquet \\
        --threshold 0.015 \\
        --work-dir data/interim/cluster \\
        --gold-dir data/gold
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl
from prefect import flow, get_run_logger, task
from prefect.tasks import task_input_hash

from ..cluster import (
    DEFAULT_THRESHOLD,
    assign_clusters,
    cluster_geographic_spread,
    cluster_growth_rate,
    cluster_size_distribution,
    cluster_subtype_mix,
    cluster_time_series,
    read_distances,
    run_tn93,
)
from ..features import build_outbreak_features


@task(retries=1)
def run_distances(aligned_fasta: Path, work_dir: Path, threshold: float) -> Path:
    log = get_run_logger()
    work_dir.mkdir(parents=True, exist_ok=True)
    output_csv = work_dir / "distances.csv"
    run_tn93(aligned_fasta, output_csv, threshold=threshold)
    log.info("tn93 distances → %s", output_csv)
    return output_csv


@task(cache_key_fn=task_input_hash)
def load_curated(path: Path) -> pl.DataFrame:
    log = get_run_logger()
    df = pl.read_parquet(path)
    log.info("loaded %d curated rows", df.height)
    return df


@task
def cluster(distances_csv: Path, curated: pl.DataFrame, threshold: float) -> pl.DataFrame:
    log = get_run_logger()
    dists = read_distances(distances_csv)
    all_accessions = curated["accession"].to_list()
    result = assign_clusters(dists, all_accessions, threshold=threshold)
    sizes = cluster_size_distribution(result.df)
    log.info("cluster size distribution:\n%s", sizes)
    return result.df


@task
def outbreak_rollup(curated: pl.DataFrame, cluster_df: pl.DataFrame) -> pl.DataFrame:
    metadata_cols = [
        c
        for c in ("accession", "country", "collection_year", "subtype")
        if c in curated.columns
    ]
    metadata = curated.select(metadata_cols)

    time_series = cluster_time_series(cluster_df, metadata)
    growth = cluster_growth_rate(time_series)
    subtype_mix = cluster_subtype_mix(cluster_df, metadata)
    _ = cluster_geographic_spread(cluster_df, metadata)  # kept for future

    features = build_outbreak_features(
        sequence_metadata=metadata,
        cluster_assignments=cluster_df,
        cluster_growth=growth,
        cluster_subtypes=subtype_mix,
    )
    return features


@task
def write_gold(
    cluster_df: pl.DataFrame,
    features_df: pl.DataFrame,
    gold_dir: Path,
) -> tuple[Path, Path]:
    log = get_run_logger()
    gold_dir.mkdir(parents=True, exist_ok=True)
    cluster_path = gold_dir / "cluster_assignments.parquet"
    features_path = gold_dir / "outbreak_features.parquet"
    cluster_df.write_parquet(cluster_path, compression="zstd", compression_level=6)
    features_df.write_parquet(features_path, compression="zstd", compression_level=6)
    log.info("wrote %d cluster rows → %s", cluster_df.height, cluster_path)
    log.info("wrote %d outbreak-feature rows → %s", features_df.height, features_path)
    return cluster_path, features_path


@flow(name="cluster")
def cluster_flow(
    aligned_fasta: Path,
    curated_parquet: Path,
    work_dir: Path,
    gold_dir: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    to_postgres: bool | None = None,
) -> tuple[Path, Path]:
    curated = load_curated(curated_parquet)
    distances_csv = run_distances(aligned_fasta, work_dir, threshold)
    cluster_df = cluster(distances_csv, curated, threshold)
    features_df = outbreak_rollup(curated, cluster_df)
    cluster_path, features_path = write_gold(cluster_df, features_df, gold_dir)

    # Storage tail — never crashes the flow.
    from ._storage_tail import storage_tail

    storage_tail("cluster", cluster_path, enabled=to_postgres)

    return cluster_path, features_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.pipeline.cluster")
    parser.add_argument(
        "--aligned-fasta",
        type=Path,
        default=Path("data/interim/curate/aligned.fasta"),
    )
    parser.add_argument(
        "--curated",
        type=Path,
        default=Path("data/interim/curated/aligned.parquet"),
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("data/interim/cluster"),
    )
    parser.add_argument("--gold-dir", type=Path, default=Path("data/gold"))
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="TN93 clustering threshold; CDC HIV-TRACE convention is 0.015 (1.5%)",
    )
    parser.add_argument(
        "--to-postgres",
        action="store_true",
        help="After writing the Parquet, also load cluster memberships "
        "into Postgres (also honors HIVDRWATCH_LOAD_TO_POSTGRES=1 env).",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cluster_flow(
        aligned_fasta=args.aligned_fasta,
        curated_parquet=args.curated,
        work_dir=args.work_dir,
        gold_dir=args.gold_dir,
        threshold=args.threshold,
        to_postgres=args.to_postgres or None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
