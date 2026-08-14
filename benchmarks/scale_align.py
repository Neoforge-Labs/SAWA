"""Scalability benchmark for the alignment stage.

Generates synthetic pol-length sequences via the JC69 simulator we
use for cluster ground truth, aligns them to a synthetic reference
via `mappy`, and records wall time + peak RSS at doubling input
sizes. Deterministic given seed.

The reference here is a random 3012-nt string — not real HXB2. We're
measuring **throughput**, not resistance-call accuracy. Alignment
correctness has its own test suite.

CLI:
    python -m benchmarks.scale_align \\
        --sizes 100 200 500 1000 2000 5000 \\
        --out benchmarks/results/scale_align.json \\
        --seed 42
"""

from __future__ import annotations

import argparse
import logging
import tempfile
from pathlib import Path

from hivdrwatch.align.align import DEFAULT_ALIGNER, align_sequences
from hivdrwatch.cluster.simulate import simulate_clusters

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


def _write_synthetic_reference(path: Path, seed: int, length: int = 3012) -> Path:
    """Emit a single-record FASTA of a random 3012-nt sequence.

    Real HXB2 would be more meaningful but requires the LANL alignment
    on disk. For throughput measurement, any 3012-nt reference works —
    mappy's build time and per-query alignment cost don't depend on
    reference *identity*, only on its length.
    """
    sim = simulate_clusters(
        n_clusters=1, cluster_size=2, n_singletons=0, seq_length=length, seed=seed
    )
    ref_seq = next(iter(sim.accession_to_seq.values()))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">HXB2_synthetic_pol\n{ref_seq}\n")
    return path


def _synthetic_records(n: int, seed: int) -> list[tuple[str, str]]:
    """Make `n` synthetic pol-length (accession, sequence) tuples."""
    sim = simulate_clusters(
        n_clusters=max(1, n // 20),
        cluster_size=(2, 5),
        n_singletons=n,  # ensures at least n across clusters + singletons
        seq_length=3012,
        seed=seed,
    )
    items = list(sim.accession_to_seq.items())[:n]
    return [(f"BENCH{i:07d}", seq) for i, (_, seq) in enumerate(items)]


def run_align_benchmark(
    sizes: tuple[int, ...],
    seed: int,
    aligner: str = DEFAULT_ALIGNER,
) -> list[BenchmarkMeasurement]:
    measurements: list[BenchmarkMeasurement] = []
    with tempfile.TemporaryDirectory(prefix="hivdrwatch-align-bench-") as tmp:
        work = Path(tmp)
        ref = _write_synthetic_reference(work / "ref.fasta", seed=seed)
        for n in sizes:
            records = _synthetic_records(n, seed=seed + n)
            with measure_block(n_sequences=n) as m:
                aligned = align_sequences(
                    records,
                    reference_fasta=ref,
                    work_dir=work,
                    aligner=aligner,
                )
                m["extras"] = {"n_aligned": float(len(aligned))}
            measurements.append(BenchmarkMeasurement(**m))
            logger.info(
                "align n=%d wall=%.2fs throughput=%.1f/s peak=%.1fMB",
                n, m["wall_seconds"], m["throughput_per_sec"], m["peak_mb"],
            )
    return measurements


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="benchmarks.scale_align")
    ap.add_argument("--sizes", nargs="+", default=["100", "200", "500", "1000", "2000"])
    ap.add_argument("--out", type=Path, default=Path("benchmarks/results/scale_align.json"))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--aligner", default=DEFAULT_ALIGNER)
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    sizes = parse_sizes(args.sizes)
    measurements = run_align_benchmark(sizes, seed=args.seed, aligner=args.aligner)

    result = BenchmarkResult(
        stage="align",
        seed=args.seed,
        generated_at=today_iso(),
        host=HostInfo.collect(),
        measurements=tuple(measurements),
        notes=f"aligner={args.aligner}; synthetic JC69 sequences at pol length",
    )
    write_benchmark_result(result, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
