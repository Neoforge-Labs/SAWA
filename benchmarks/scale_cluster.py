"""Scalability benchmark for the clustering stage.

Simulates a mixed set of true clusters + background singletons, writes
a FASTA, runs tn93 (Docker) + connected-components assembly, and
records wall time + peak RSS at doubling input sizes.

Requires the tn93 image (`make images`). If Docker isn't available,
the script logs a WARNING and skips — CI runs it on Linux with
Docker; local dev may not.

CLI:
    python -m benchmarks.scale_cluster \\
        --sizes 100 200 500 1000 2000 5000 \\
        --out benchmarks/results/scale_cluster.json \\
        --seed 42
"""

from __future__ import annotations

import argparse
import logging
import shutil
import tempfile
from pathlib import Path

from hivdrwatch.cluster import (
    DEFAULT_THRESHOLD,
    assign_clusters,
    read_distances,
    run_tn93,
    simulate_clusters,
)

from .common import (
    BenchmarkMeasurement,
    BenchmarkResult,
    HostInfo,
    measure_block,
    parse_sizes,
    today_iso,
    write_benchmark_result,
)

logger = logging.getLogger(__name__)


def _has_docker() -> bool:
    return shutil.which("docker") is not None


def run_cluster_benchmark(
    sizes: tuple[int, ...],
    seed: int,
    threshold: float = DEFAULT_THRESHOLD,
) -> list[BenchmarkMeasurement]:
    if not _has_docker():
        logger.warning("docker not on PATH; cluster benchmark cannot run")
        return []

    measurements: list[BenchmarkMeasurement] = []
    with tempfile.TemporaryDirectory(prefix="hivdrwatch-bench-") as tmp:
        work = Path(tmp)
        for n in sizes:
            n_clusters = max(1, n // 20)
            per = 5
            n_singletons = max(0, n - n_clusters * per)
            sim = simulate_clusters(
                n_clusters=n_clusters,
                cluster_size=per,
                n_singletons=n_singletons,
                seed=seed + n,
            )
            fasta = sim.write_fasta(work / f"n{n}.fasta")
            csv = work / f"n{n}.csv"
            with measure_block(n_sequences=n) as m:
                run_tn93(fasta, csv, threshold=threshold)
                dists = read_distances(csv)
                result = assign_clusters(dists, list(sim.accession_to_seq), threshold=threshold)
                m["extras"] = {
                    "n_predicted_clusters": float(result.n_clusters),
                    "n_pairs_below_threshold": float(len(dists)),
                }
            measurements.append(BenchmarkMeasurement(**m))
            logger.info(
                "cluster n=%d wall=%.2fs throughput=%.1f/s peak=%.1fMB",
                n, m["wall_seconds"], m["throughput_per_sec"], m["peak_mb"],
            )
    return measurements


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="benchmarks.scale_cluster")
    ap.add_argument("--sizes", nargs="+", default=["100", "200", "500", "1000"])
    ap.add_argument("--out", type=Path, default=Path("benchmarks/results/scale_cluster.json"))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    sizes = parse_sizes(args.sizes)
    measurements = run_cluster_benchmark(sizes, seed=args.seed, threshold=args.threshold)

    result = BenchmarkResult(
        stage="cluster",
        seed=args.seed,
        generated_at=today_iso(),
        host=HostInfo.collect(),
        measurements=tuple(measurements),
        notes=f"tn93 Docker + connected-components at threshold={args.threshold}",
    )
    write_benchmark_result(result, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
