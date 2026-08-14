"""Prefect flow: merged sequence parquet → Sierra → dr_results parquet.

This is the Phase-1 walking skeleton for the DR side of the pipeline.
It reads the merged silver-layer parquet produced by `hivdrwatch.ingest`,
sends every sequence through Stanford HIVdb Sierra (cache-first), and
writes canonical DRResult + SubtypeCall parquets to `data/interim/`.

Run locally:
    prefect worker start --pool default-agent-pool   # in another shell
    python -m hivdrwatch.pipeline.dr_ingest \\
        --input data/interim/east_africa_pol.parquet \\
        --cache-dir data/raw/hivdb \\
        --out-dir data/interim

Cache is the primary interface — re-runs of the same input do zero network
work. Use --offline to hard-fail on cache misses (useful for CI).
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl
from prefect import flow, get_run_logger, task
from prefect.tasks import task_input_hash

from ..resistance import (
    AlgorithmVersion,
    DRResult,
    SierraCache,
    SierraClient,
    SubtypeCall,
    pin_algorithm,
)

# ------------------------------------------------------------------
# Tasks
# ------------------------------------------------------------------


@task(retries=2, retry_delay_seconds=10, cache_key_fn=task_input_hash)
def load_merged_sequences(parquet_path: Path) -> pl.DataFrame:
    log = get_run_logger()
    df = pl.read_parquet(parquet_path)
    if "accession" not in df.columns or "sequence" not in df.columns:
        raise ValueError(f"{parquet_path} missing required columns (accession, sequence)")
    # Sierra needs non-gapped sequences; strip alignment gaps here.
    df = df.with_columns(pl.col("sequence").str.replace_all("-", "").alias("sequence")).filter(
        pl.col("sequence").str.len_chars() >= 200
    )
    log.info("loaded %d sequences from %s", df.height, parquet_path)
    return df


@task(retries=3, retry_delay_seconds=30)
def call_sierra(
    df: pl.DataFrame,
    cache_dir: Path,
    *,
    batch_size: int = 100,
    offline: bool = False,
) -> tuple[list[DRResult], list[SubtypeCall], AlgorithmVersion]:
    """Run every sequence through Sierra (cache-first). Returns flat lists."""
    log = get_run_logger()
    cache = SierraCache(cache_dir)
    with SierraClient(cache=cache, batch_size=batch_size, offline=offline) as client:
        algo = pin_algorithm(client)
        pairs = list(zip(df["accession"].to_list(), df["sequence"].to_list(), strict=True))
        dr_results, subtypes = client.analyze_sequences(pairs)
    log.info(
        "Sierra: %d DR rows, %d subtype calls, algorithm=%s",
        len(dr_results),
        len(subtypes),
        algo.version,
    )
    return dr_results, subtypes, algo


@task
def write_dr_parquet(
    dr_results: list[DRResult], subtypes: list[SubtypeCall], out_dir: Path
) -> tuple[Path, Path]:
    log = get_run_logger()
    out_dir.mkdir(parents=True, exist_ok=True)

    dr_path = out_dir / "dr_results.parquet"
    if dr_results:
        pl.from_dicts([r.model_dump(mode="json") for r in dr_results]).write_parquet(
            dr_path, compression="zstd", compression_level=6
        )
    else:
        pl.DataFrame().write_parquet(dr_path)
    log.info("wrote %d DR rows → %s", len(dr_results), dr_path)

    st_path = out_dir / "sierra_subtypes.parquet"
    if subtypes:
        pl.from_dicts([s.model_dump(mode="json") for s in subtypes]).write_parquet(
            st_path, compression="zstd", compression_level=6
        )
    else:
        pl.DataFrame().write_parquet(st_path)
    log.info("wrote %d subtype rows → %s", len(subtypes), st_path)

    return dr_path, st_path


# ------------------------------------------------------------------
# Flow
# ------------------------------------------------------------------


@flow(name="dr-ingest")
def dr_ingest_flow(
    input_parquet: Path,
    cache_dir: Path,
    out_dir: Path,
    *,
    batch_size: int = 100,
    offline: bool = False,
    to_postgres: bool | None = None,
) -> dict[str, Path | str]:
    df = load_merged_sequences(input_parquet)
    dr_results, subtypes, algo = call_sierra(
        df, cache_dir, batch_size=batch_size, offline=offline
    )
    dr_path, st_path = write_dr_parquet(dr_results, subtypes, out_dir)

    # Storage tail — never crashes the flow.
    from ._storage_tail import storage_tail

    storage_tail("sierra", dr_path, enabled=to_postgres)

    return {
        "dr_results": dr_path,
        "subtypes": st_path,
        "algorithm_version": algo.version,
    }


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.pipeline.dr_ingest")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/interim/east_africa_pol.parquet"),
        help="merged sequence parquet from hivdrwatch.ingest",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("data/raw/hivdb"),
        help="Sierra response cache",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/interim"),
        help="where dr_results.parquet + sierra_subtypes.parquet are written",
    )
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument(
        "--offline",
        action="store_true",
        help="hard-fail on cache misses (CI / reproducibility mode)",
    )
    parser.add_argument(
        "--to-postgres",
        action="store_true",
        help="After writing the Parquet, also load into Postgres "
        "(also honors HIVDRWATCH_LOAD_TO_POSTGRES=1 env var).",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    result = dr_ingest_flow(
        input_parquet=args.input,
        cache_dir=args.cache_dir,
        out_dir=args.out_dir,
        batch_size=args.batch_size,
        offline=args.offline,
        to_postgres=args.to_postgres or None,
    )
    print(f"algorithm: {result['algorithm_version']}")
    print(f"dr_results: {result['dr_results']}")
    print(f"subtypes:   {result['subtypes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
