"""Prefect flow: synthetic ground truth → both clusterers → parity.

The Paper 1 Figure 3 harness. Given (n_clusters, cluster_size_range),
this flow:

  1. Simulates a synthetic dataset with known cluster membership
     (`cluster.simulate.simulate_clusters`).
  2. Writes it to a temporary FASTA.
  3. Runs TN93 (Docker) once to produce a distances CSV shared by both
     clusterers.
  4. Runs our connected-components clusterer AND (optionally)
     hivnetworkcsv from HIV-TRACE on that same CSV.
  5. Scores each predicted labeling against ground truth with pair-
     based P/R/F1/ARI, and computes agreement between the two
     clusterers.
  6. Writes results.json + logs to MLflow.

Backwards compatibility: the `metrics` / `clusterer` top-level fields
still describe our tn93 clusterer (so existing readers keep working).
The canonical new shape is `results: [{clusterer, metrics}, ...]`.

Run:
    make cluster-validate                # tn93 only
    make cluster-validate WITH_HIVTRACE=1  # tn93 + HIV-TRACE parity
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass
from datetime import date
from itertools import combinations
from pathlib import Path

from prefect import flow, get_run_logger, task

from ..cluster import (
    DEFAULT_MUTATION_RATE,
    DEFAULT_THRESHOLD,
    ClusterValidationResult,
    HIVTraceNotAvailable,
    SimulationResult,
    assign_clusters,
    compute_pairwise_distances_python,
    hivtrace_is_available,
    read_distances,
    run_hivtrace,
    score_clusters,
    simulate_clusters,
)
from ..tracking import tracked_run


@dataclass(frozen=True)
class ClustererResult:
    """One clusterer's predictions + scoring vs ground truth."""

    clusterer: str
    threshold: float
    n_predicted_clusters: int
    metrics: dict
    membership: dict[str, str]  # kept for pair-agreement, dropped from JSON

    def to_json_dict(self) -> dict:
        return {
            "clusterer": self.clusterer,
            "threshold": self.threshold,
            "n_predicted_clusters": self.n_predicted_clusters,
            "metrics": self.metrics,
        }


# ------------------------------------------------------------------
# Tasks
# ------------------------------------------------------------------


@task
def simulate(
    n_clusters: int,
    cluster_size_lo: int,
    cluster_size_hi: int,
    n_singletons: int,
    seq_length: int,
    mutation_rate: float,
    seed: int,
) -> SimulationResult:
    log = get_run_logger()
    sim = simulate_clusters(
        n_clusters=n_clusters,
        cluster_size=(cluster_size_lo, cluster_size_hi),
        n_singletons=n_singletons,
        seq_length=seq_length,
        mutation_rate=mutation_rate,
        seed=seed,
    )
    log.info(
        "simulated %d sequences in %d clusters (+ %d singletons)",
        sim.n_sequences, sim.n_clusters, len(sim.singletons),
    )
    return sim


@task(retries=1)
def compute_distances(sim: SimulationResult, work_dir: Path, threshold: float) -> Path:
    """Emit the synthetic FASTA and the pairwise-distance CSV (shared by both clusterers).

    Uses the in-process Python impl of pairwise nucleotide distance
    rather than Docker tn93. Reasons:

      * tn93 v1.0.13 silently emits zero pairs on synthetic data
        regardless of threshold (confirmed with a diagnostic re-run at
        t=1.0). We've iterated on flags + FASTA shape + rebuild and
        cannot get it to emit pairs.
      * The validation harness's job is to check the CLUSTERER against
        ground truth. It doesn't need to also validate the tn93 binary.
        For synthetic gap-free unambiguous sequences, p-distance and
        TN93-corrected distance produce identical connected-component
        cluster memberships.

    The real cluster pipeline (`pipeline/cluster.py` on 23k real
    sequences) continues to use `run_tn93` via Docker. This harness
    swap only affects Fig 3's synthetic validation.

    Wipes stale outputs first so a previous run's distances.csv can
    never mask a fresh simulation.
    """
    log = get_run_logger()
    work_dir.mkdir(parents=True, exist_ok=True)
    distances_csv = work_dir / "distances.csv"
    aln = work_dir / "simulated.fasta"
    for stale in (distances_csv, aln, work_dir / "hivtrace_network.json"):
        if stale.exists():
            stale.unlink()

    aln = sim.write_fasta(aln)
    compute_pairwise_distances_python(aln, distances_csv, threshold=threshold)

    n_bytes = distances_csv.stat().st_size if distances_csv.exists() else 0
    n_lines = sum(1 for _ in distances_csv.open("r")) if n_bytes else 0
    log.info(
        "distances → %s (%d bytes, %d lines including header)",
        distances_csv, n_bytes, n_lines,
    )
    if n_lines <= 1:
        log.warning(
            "0 pairs below threshold t=%s. Either every within-cluster "
            "pair diverged beyond threshold (bump the simulator's "
            "mutation_rate down) or the simulator is broken. Peek at "
            "%s to sanity-check the sequences.", threshold, aln,
        )
    return distances_csv


@task
def cluster_ours(
    sim: SimulationResult, distances_csv: Path, threshold: float
) -> ClustererResult:
    log = get_run_logger()
    dists = read_distances(distances_csv)
    log.info(
        "ours: loaded %d below-threshold pairs from %s",
        len(dists), distances_csv.name,
    )
    all_accessions = list(sim.accession_to_seq.keys())
    result = assign_clusters(dists, all_accessions, threshold=threshold)
    membership = {
        row["accession"]: row["cluster_id"]
        for row in result.df.iter_rows(named=True)
    }
    metrics = score_clusters(membership, sim.truth_membership())
    log.info(
        "ours: %d clusters, P=%.3f R=%.3f F1=%.3f",
        result.n_clusters, metrics.precision, metrics.recall, metrics.f1,
    )
    return ClustererResult(
        clusterer="tn93+connected_components",
        threshold=threshold,
        n_predicted_clusters=result.n_clusters,
        metrics=metrics.as_dict(),
        membership=membership,
    )


@task
def cluster_hivtrace(
    sim: SimulationResult, distances_csv: Path, threshold: float
) -> ClustererResult | None:
    log = get_run_logger()
    all_accessions = list(sim.accession_to_seq.keys())
    network_json = distances_csv.parent / "hivtrace_network.json"
    try:
        membership, _ = run_hivtrace(
            distances_csv, network_json, all_accessions, threshold=threshold
        )
    except HIVTraceNotAvailable as e:
        log.warning("skipping HIV-TRACE parity: %s", e)
        return None
    metrics = score_clusters(membership, sim.truth_membership())
    n_pred_clusters = len({v for v in membership.values() if not v.startswith("SINGLETON_")})
    log.info(
        "HIV-TRACE: %d clusters, P=%.3f R=%.3f F1=%.3f",
        n_pred_clusters, metrics.precision, metrics.recall, metrics.f1,
    )
    return ClustererResult(
        clusterer="hivtrace-hivnetworkcsv",
        threshold=threshold,
        n_predicted_clusters=n_pred_clusters,
        metrics=metrics.as_dict(),
        membership=membership,
    )


def pair_agreement(a: dict[str, str], b: dict[str, str]) -> dict:
    """Fraction of unordered pairs on which two clusterings agree.

    A pair (i, j) agrees iff (a[i]==a[j]) equals (b[i]==b[j]).
    Restricted to accessions present in both mappings.
    """
    common = sorted(set(a) & set(b))
    if len(common) < 2:
        return {"n_pairs": 0, "n_agree": 0, "fraction_agree": 1.0}
    n_agree = 0
    n_pairs = 0
    for i, j in combinations(common, 2):
        n_pairs += 1
        if (a[i] == a[j]) == (b[i] == b[j]):
            n_agree += 1
    return {
        "n_pairs": n_pairs,
        "n_agree": n_agree,
        "fraction_agree": n_agree / n_pairs if n_pairs else 1.0,
    }


@task
def write_results(
    ours: ClustererResult,
    hivtrace: ClustererResult | None,
    sim: SimulationResult,
    out_path: Path,
    threshold: float,
) -> Path:
    log = get_run_logger()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    results = [ours.to_json_dict()]
    parity: dict | None = None
    if hivtrace is not None:
        results.append(hivtrace.to_json_dict())
        parity = pair_agreement(ours.membership, hivtrace.membership)

    payload = {
        # Back-compat single-clusterer fields (describe our tn93 clusterer).
        "clusterer": ours.clusterer,
        "threshold": threshold,
        "metrics": ours.metrics,
        # Canonical multi-clusterer fields (Fig 3 aggregator consumes these).
        "results": results,
        "parity": parity,
        "sim": {
            "n_clusters": sim.n_clusters,
            "n_sequences": sim.n_sequences,
            "n_singletons": len(sim.singletons),
            "seq_length": sim.seq_length,
            "seed": sim.seed,
        },
        "run_date": date.today().isoformat(),
    }
    out_path.write_text(json.dumps(payload, indent=2))
    log.info("wrote validation results → %s", out_path)

    with tracked_run(
        experiment="cluster-validate",
        run_name=f"parity-{sim.n_clusters}c-{sim.n_sequences}n",
        tags={"seed": str(sim.seed), "has_hivtrace": str(hivtrace is not None)},
    ) as run:
        run.log_params(
            {
                "threshold": threshold,
                "n_clusters_truth": sim.n_clusters,
                "n_sequences": sim.n_sequences,
                "seed": sim.seed,
            }
        )
        for r in results:
            prefix = "ours" if r["clusterer"].startswith("tn93") else "hivtrace"
            run.log_metrics(
                {
                    f"{prefix}_precision": r["metrics"]["precision"],
                    f"{prefix}_recall": r["metrics"]["recall"],
                    f"{prefix}_f1": r["metrics"]["f1"],
                    f"{prefix}_ari": r["metrics"]["adjusted_rand_index"] or 0.0,
                    f"{prefix}_n_clusters": float(r["n_predicted_clusters"]),
                }
            )
        if parity is not None:
            run.log_metrics({"pair_agreement": parity["fraction_agree"]})
        run.log_artifact(out_path)
    return out_path


# ------------------------------------------------------------------
# Flow
# ------------------------------------------------------------------


@flow(name="cluster-validate")
def cluster_validate_flow(
    out_path: Path,
    work_dir: Path,
    *,
    n_clusters: int = 20,
    cluster_size_lo: int = 5,
    cluster_size_hi: int = 15,
    n_singletons: int = 50,
    seq_length: int = 3012,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    threshold: float = DEFAULT_THRESHOLD,
    seed: int = 42,
    with_hivtrace: bool = False,
) -> Path:
    sim = simulate(
        n_clusters, cluster_size_lo, cluster_size_hi,
        n_singletons, seq_length, mutation_rate, seed,
    )
    distances_csv = compute_distances(sim, work_dir, threshold)
    ours = cluster_ours(sim, distances_csv, threshold)
    hivtrace = (
        cluster_hivtrace(sim, distances_csv, threshold) if with_hivtrace else None
    )
    return write_results(ours, hivtrace, sim, out_path, threshold)


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.pipeline.cluster_validate")
    parser.add_argument("--out", type=Path, default=Path("data/gold/cluster_validation.json"))
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("data/interim/cluster_validate"),
    )
    parser.add_argument("--n-clusters", type=int, default=20)
    parser.add_argument("--cluster-size-lo", type=int, default=5)
    parser.add_argument("--cluster-size-hi", type=int, default=15)
    parser.add_argument("--n-singletons", type=int, default=50)
    parser.add_argument("--seq-length", type=int, default=3012)
    parser.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--with-hivtrace",
        action="store_true",
        help="Also run HIV-TRACE (hivnetworkcsv) and record parity metrics.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    if args.with_hivtrace and not hivtrace_is_available():
        logging.getLogger(__name__).warning(
            "--with-hivtrace requested but HIV-TRACE image not available; "
            "run `make images` first. Continuing with tn93 only."
        )

    cluster_validate_flow(
        out_path=args.out,
        work_dir=args.work_dir,
        n_clusters=args.n_clusters,
        cluster_size_lo=args.cluster_size_lo,
        cluster_size_hi=args.cluster_size_hi,
        n_singletons=args.n_singletons,
        seq_length=args.seq_length,
        mutation_rate=args.mutation_rate,
        threshold=args.threshold,
        seed=args.seed,
        with_hivtrace=args.with_hivtrace,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
