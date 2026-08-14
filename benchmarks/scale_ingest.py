"""Scalability benchmark for the ingest stage.

We can't hammer LANL or GenBank in CI (rate limits + reproducibility).
Instead we time the *parser* against a synthetic LANL-shaped TSV of
increasing size. That's the honest scale story: given the file, how
fast does the parser hand off canonical `SequenceRecord`s? Network
throughput is a separate story documented in the ingest README.

CLI:
    python -m benchmarks.scale_ingest \\
        --sizes 100 500 1000 5000 10000 \\
        --out benchmarks/results/scale_ingest.json \\
        --seed 42
"""

from __future__ import annotations

import argparse
import logging
import tempfile
from pathlib import Path

from hivdrwatch.cluster.simulate import simulate_clusters
from hivdrwatch.ingest.lanl import LANLIngester

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


# Minimal LANL Table column set the parser recognizes. Extra cells are
# tolerated; missing ones default to null. The parser uses `Accession`,
# `Subtype`, `Country`, `Sampling Year`, `Sequence Length`, `Sequence`.
_HEADER = [
    "Accession",
    "Name",
    "Subtype",
    "Country",
    "Sampling Year",
    "Genomic Region",
    "Sequence Length",
    "Sequence",
]


def _write_synthetic_lanl_dir(n: int, root: Path, seed: int) -> Path:
    """Emit `root/ke/results.tsv` with `n` synthetic LANL rows."""
    country_dir = root / "ke"
    country_dir.mkdir(parents=True, exist_ok=True)
    tsv = country_dir / "results.tsv"

    sim = simulate_clusters(
        n_clusters=1,
        cluster_size=max(2, min(n, 5)),
        n_singletons=max(0, n - 5),
        seq_length=1000,
        seed=seed,
    )
    seqs = list(sim.accession_to_seq.values())[:n]
    with tsv.open("w", encoding="utf-8") as fh:
        # LANL prepends a preamble line before the header — the parser
        # detects and skips lines starting with "Number of records".
        # Anything else would be misread as the header row.
        fh.write(f"Number of records retrieved: {n}\n")
        fh.write("\t".join(_HEADER) + "\n")
        # Accessions must match SequenceRecord's regex: [A-Z]{1,2}\d{5,8}.
        # Use "BE1234567" style — 2 letters + 7 digits.
        for i, seq in enumerate(seqs):
            row = [
                f"BE{i:07d}",        # Accession
                f"KE.BENCH.{i}",     # Name
                "A1",                # Subtype
                "Kenya",             # Country (name; the parser maps to KE)
                "2020",              # Sampling Year
                "pol",               # Genomic Region
                str(len(seq)),       # Sequence Length
                seq,                 # Sequence
            ]
            fh.write("\t".join(row) + "\n")
    return root


def run_ingest_benchmark(
    sizes: tuple[int, ...],
    seed: int,
) -> list[BenchmarkMeasurement]:
    measurements: list[BenchmarkMeasurement] = []
    for n in sizes:
        with tempfile.TemporaryDirectory(prefix=f"hivdrwatch-ingest-{n}-") as tmp:
            root = _write_synthetic_lanl_dir(n, Path(tmp), seed=seed + n)
            with measure_block(n_sequences=n) as m:
                ingester = LANLIngester(root)
                records = list(ingester.to_sequence_records(min_length=200))
                m["extras"] = {"n_parsed": float(len(records))}
            measurements.append(BenchmarkMeasurement(**m))
            logger.info(
                "ingest n=%d wall=%.3fs throughput=%.1f/s peak=%.1fMB",
                n, m["wall_seconds"], m["throughput_per_sec"], m["peak_mb"],
            )
    return measurements


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="benchmarks.scale_ingest")
    ap.add_argument("--sizes", nargs="+", default=["100", "500", "1000", "5000"])
    ap.add_argument("--out", type=Path, default=Path("benchmarks/results/scale_ingest.json"))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    sizes = parse_sizes(args.sizes)
    measurements = run_ingest_benchmark(sizes, seed=args.seed)

    result = BenchmarkResult(
        stage="ingest",
        seed=args.seed,
        generated_at=today_iso(),
        host=HostInfo.collect(),
        measurements=tuple(measurements),
        notes="synthetic LANL-shaped TSV parser throughput (no network)",
    )
    write_benchmark_result(result, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
