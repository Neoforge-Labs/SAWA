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
    build_outbreak_features,
    cluster_geographic_spread,
    cluster_growth_rate,
    cluster_size_distribution,
    cluster_subtype_mix,
    cluster_time_series,
    read_distances,
    run_tn93,
)
from ..cluster.distance import DEFAULT_MIN_OVERLAP, compute_pairwise_distances_python
from ..phylo.build import MIN_ALIGNED_COVERAGE, filter_by_coverage

logger = logging.getLogger(__name__)


def eligible_accessions(
    curated: pl.DataFrame,
    *,
    seq_col: str = "aligned_pol",
    require_qc_pass: bool = True,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
) -> list[str]:
    """Accessions whose alignment is fit to compute distances from.

    `phylo/build.py` has filtered its tree input on `passes_qc` and
    `MIN_ALIGNED_COVERAGE` since it was written; the cluster path never
    did, which is how 21,358 rows of 3,012 dashes reached the distance
    engine (`docs/evaluations/align-coverage.md` §6). The same two
    gates, applied with the same helper, so the tree and the network
    are built from the same cohort.

    Returns one accession per record, first row wins. The curated table
    carries the same accession up to three times (one row per country
    query that returned it — `align-coverage.md` §7 item 3), and a
    FASTA with repeated identifiers makes the distance CSV ambiguous
    about which record a pair refers to.
    """
    if seq_col not in curated.columns:
        raise ValueError(f"curated frame lacks '{seq_col}' column; cannot materialize FASTA")

    df = curated.filter(pl.col(seq_col).is_not_null() & (pl.col(seq_col).str.len_chars() > 0))
    n_present = df.height

    if require_qc_pass and "passes_qc" in df.columns:
        df = df.filter(pl.col("passes_qc"))
    n_qc = df.height

    df = filter_by_coverage(df, seq_col=seq_col, min_coverage=min_coverage)
    n_covered = df.height

    accessions = df["accession"].to_list()
    deduped = list(dict.fromkeys(accessions))
    logger.info(
        "distance-engine eligibility: %d rows with a sequence → %d passing QC "
        "→ %d with coverage ≥ %.0f%% → %d distinct accessions",
        n_present,
        n_qc,
        n_covered,
        min_coverage * 100,
        len(deduped),
    )
    if len(deduped) < n_covered:
        logger.warning(
            "dropped %d duplicate accession rows; the curated table stores "
            "the same record once per country query",
            n_covered - len(deduped),
        )
    return deduped


@task
def materialize_aligned_fasta(
    aligned_fasta: Path,
    curated_parquet: Path,
    work_dir: Path,
    *,
    require_qc_pass: bool = True,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
) -> Path:
    """Write the aligned FASTA the distance engine consumes.

    Sequences come from the curated Parquet's `aligned_pol` column, or
    from `aligned_fasta` when the curate stage already wrote one. Either
    way the record set is restricted to `eligible_accessions` — the
    Parquet is the only place `passes_qc` and the coverage of each
    alignment are recorded, so it decides what is eligible even when the
    sequences themselves come from a file.
    """
    log = get_run_logger()
    work_dir.mkdir(parents=True, exist_ok=True)
    curated = pl.read_parquet(curated_parquet)
    eligible = eligible_accessions(
        curated,
        require_qc_pass=require_qc_pass,
        min_coverage=min_coverage,
    )
    if not eligible:
        raise ValueError(
            f"no eligible records in {curated_parquet}: every row is empty, "
            f"QC-failed, or below {min_coverage:.0%} aligned coverage"
        )
    keep = set(eligible)

    out = work_dir / "aligned.fasta"
    if aligned_fasta.exists() and aligned_fasta.stat().st_size > 0:
        n_written, n_dropped = _filter_fasta(aligned_fasta, out, keep)
        log.info(
            "filtered %s → %s (%d records kept, %d dropped as ineligible)",
            aligned_fasta,
            out,
            n_written,
            n_dropped,
        )
    else:
        n_written = _write_fasta_from_curated(curated, out, eligible)
        log.info(
            "materialized aligned FASTA from %s → %s (%d records)",
            curated_parquet,
            out,
            n_written,
        )

    if n_written == 0:
        raise ValueError(f"wrote 0 records to {out}: no eligible accession carried a sequence")
    return out


def _write_fasta_from_curated(
    curated: pl.DataFrame,
    out: Path,
    accessions: list[str],
) -> int:
    """Write one record per accession in `accessions`, in that order."""
    by_accession = {
        row["accession"]: row["aligned_pol"]
        for row in curated.select(["accession", "aligned_pol"]).iter_rows(named=True)
    }
    n_written = 0
    with out.open("w", encoding="utf-8") as fh:
        for accession in accessions:
            seq = by_accession.get(accession)
            if not seq:
                continue
            fh.write(f">{accession}\n{seq}\n")
            n_written += 1
    return n_written


def _filter_fasta(source: Path, out: Path, keep: set[str]) -> tuple[int, int]:
    """Copy `source` to `out` keeping only records whose id is in `keep`."""
    n_written = 0
    n_dropped = 0
    seen: set[str] = set()
    with source.open("r", encoding="utf-8") as fh, out.open("w", encoding="utf-8") as sink:
        keeping = False
        for raw in fh:
            if raw.startswith(">"):
                accession = raw[1:].strip().split()[0] if raw[1:].strip() else ""
                keeping = accession in keep and accession not in seen
                if keeping:
                    seen.add(accession)
                    n_written += 1
                else:
                    n_dropped += 1
            if keeping:
                sink.write(raw)
    return n_written, n_dropped


@task(retries=1)
def run_distances(
    aligned_fasta: Path,
    work_dir: Path,
    threshold: float,
    engine: str = "python",
    min_overlap: int = DEFAULT_MIN_OVERLAP,
) -> Path:
    """Compute pairwise distances → CSV of below-threshold pairs.

    Two engines:

    * ``python`` (default) — the in-process vectorized p-distance over
      the pairwise unambiguous overlap, monotone with TN93 for small
      divergences. See
      ``cluster/distance.py::compute_pairwise_distances_python``.
    * ``docker`` — the pinned `tn93` container, the reference
      implementation. At v1.0.13 it emitted zero pairs on both the
      synthetic and the real cohort; its own 100-site minimum overlap
      is one reason why, since the FASTA it was handed was half
      all-gap rows (``docs/evaluations/align-coverage.md``). Retained
      for A/B comparison.
    """
    log = get_run_logger()
    work_dir.mkdir(parents=True, exist_ok=True)
    output_csv = work_dir / "distances.csv"
    if engine == "docker":
        run_tn93(aligned_fasta, output_csv, threshold=threshold)
        log.info("tn93 (docker) distances → %s", output_csv)
    else:
        compute_pairwise_distances_python(
            aligned_fasta, output_csv, threshold=threshold, min_overlap=min_overlap
        )
        log.info("python distances → %s (min overlap %d sites)", output_csv, min_overlap)
    n_pairs = max(0, sum(1 for _ in output_csv.open()) - 1)
    log.info("distances CSV: %d below-threshold pairs at t=%.4f", n_pairs, threshold)
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
    distance_engine: str = "python",
    min_overlap: int = DEFAULT_MIN_OVERLAP,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
) -> tuple[Path, Path]:
    curated = load_curated(curated_parquet)
    aligned_fasta = materialize_aligned_fasta(
        aligned_fasta, curated_parquet, work_dir, min_coverage=min_coverage
    )
    distances_csv = run_distances(
        aligned_fasta,
        work_dir,
        threshold,
        engine=distance_engine,
        min_overlap=min_overlap,
    )
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
    parser.add_argument(
        "--distance-engine",
        choices=("python", "docker"),
        default="python",
        help=(
            "Distance backend. Default 'python' — the in-process "
            "vectorized p-distance over the pairwise unambiguous overlap. "
            "'docker' runs the pinned tn93 container as the reference "
            "implementation; at v1.0.13 it emitted zero pairs on this "
            "cohort, see docs/evaluations/align-coverage.md."
        ),
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=DEFAULT_MIN_OVERLAP,
        help=(
            "Fewest shared unambiguous sites a pair needs before a distance "
            f"is reported for it (default {DEFAULT_MIN_OVERLAP}, matching "
            "tn93's own -l default). Pairs below it are not emitted."
        ),
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=MIN_ALIGNED_COVERAGE,
        help=(
            "Least non-gap fraction of the alignment a record needs to enter "
            f"the distance engine (default {MIN_ALIGNED_COVERAGE:.0%}, the "
            "same gate phylo/ applies to tree input)."
        ),
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
        distance_engine=args.distance_engine,
        min_overlap=args.min_overlap,
        min_coverage=args.min_coverage,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
