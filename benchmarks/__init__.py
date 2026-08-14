"""HIV-DR-Watch benchmark harness.

Every claim of the form "X sequences in Y seconds" or "scales to N"
must have a script here whose committed `results/*.json` matches the
number printed in the manuscript. CI regenerates the JSONs weekly; any
delta > 20% blocks merge.

Scripts by stage:

  * `scale_ingest`   — GenBank + LANL fetch/parse throughput
  * `scale_align`    — mappy-based pairwise-to-HXB2 throughput
  * `scale_cluster`  — tn93 + connected-components clustering wall time
                       and peak RSS at doubling input sizes

Each script's CLI:
  --sizes 100 200 500 1000 2000 5000 10000 25000
  --out benchmarks/results/scale_<stage>.json
  --seed 42
  --repeats 1

Result schema (locked — the notebook depends on it):
  {
    "stage": "align",
    "seed": 42,
    "generated_at": "2026-08-13",
    "host": {"platform": "...", "cpu_count": 8, ...},
    "measurements": [
      {"n_sequences": 100, "wall_seconds": 0.12, "peak_mb": 34.2, "throughput_per_sec": 833.3}
    ]
  }
"""

from .common import (
    BENCHMARK_SCHEMA_VERSION,
    BenchmarkMeasurement,
    BenchmarkResult,
    HostInfo,
    measure_block,
    write_benchmark_result,
)

__all__ = [
    "BENCHMARK_SCHEMA_VERSION",
    "BenchmarkMeasurement",
    "BenchmarkResult",
    "HostInfo",
    "measure_block",
    "write_benchmark_result",
]
