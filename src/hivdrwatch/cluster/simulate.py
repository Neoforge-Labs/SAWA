"""Synthetic transmission-cluster + sequence simulator.

Produces alignments with **known ground-truth cluster membership** so
we can score any clusterer's precision/recall. Pure-Python, no external
deps beyond numpy + biopython (both already in the env). If we later
want richer epidemic dynamics, `favites` on bioconda plugs into the
same interface — this module's job is to define the interface, not to
be the most sophisticated simulator.

## Model

For each of `n_clusters` clusters:
  1. Sample a **founder** sequence uniformly at random (nucleotides).
  2. Evolve **within-cluster** members from the founder under a
     Jukes-Cantor process at a low substitution rate → member
     sequences are all within ~cluster_divergence of the founder.

Then generate `n_singletons` **background** sequences — random
nucleotides at full divergence; these are the "outside any cluster"
noise a real dataset always has.

Because we generated the sequences ourselves, ground-truth cluster
membership is trivially known: `{cluster_id: [accession, ...]}`.

## Determinism

Every stochastic step draws from a `numpy.random.Generator` seeded by
the `seed` param. Two runs with the same seed produce bit-identical
output.

## Realism caveats (call out in the paper's methods)

  * JC69 is the simplest substitution model — no gamma rate variation,
    no codon-position bias, no APOBEC signature. Fine for validating
    that a *distance-based* clusterer recovers true clusters, since
    tn93 is also nucleotide-only.
  * No molecular clock — divergence controlled by rate × time product,
    not two separate parameters.
  * No transmission-network topology within a cluster — all members
    branch from the founder, not from each other. For fig-3-level
    ground truth this is sufficient; a proper coalescent within-cluster
    is a v0.2 improvement.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# HXB2 pol is ~3012 nt; we default to the same length so simulated
# sequences slot into the same tn93 threshold as real data (0.015).
DEFAULT_SEQ_LENGTH = 3012

# Substitution rate for the within-cluster evolution step. Chosen so
# that at the default `cluster_divergence=0.005` (0.5% pairwise), the
# resulting sequences are well below the 1.5% cluster threshold.
DEFAULT_MUTATION_RATE = 0.001  # per-site substitution probability


@dataclass(frozen=True)
class SimulationResult:
    """A generated synthetic dataset with ground-truth cluster labels."""

    accession_to_seq: dict[str, str]
    ground_truth: dict[str, list[str]]     # cluster_id → [accessions]
    singletons: list[str]                   # accessions with no cluster
    n_clusters: int
    seq_length: int
    seed: int

    @property
    def n_sequences(self) -> int:
        return len(self.accession_to_seq)

    def write_fasta(self, out_path: Path, *, wrap: int = 0) -> Path:
        """Emit the alignment as FASTA — ordered clusters then singletons.

        Default is single-line-per-sequence (wrap=0). This matches what
        the previous successful cluster-validate run used and what tn93
        v1.0.13 parses correctly. Pass `wrap=60` for NCBI-style wrapped
        output if you need to feed the FASTA to a tool with a fixed
        line buffer.
        """
        out_path.parent.mkdir(parents=True, exist_ok=True)

        def _write(fh, accession: str) -> None:
            seq = self.accession_to_seq[accession]
            fh.write(f">{accession}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    fh.write(seq[i : i + wrap] + "\n")
            else:
                fh.write(seq + "\n")

        with out_path.open("w", encoding="utf-8") as fh:
            for cluster_id in sorted(self.ground_truth):
                for accession in self.ground_truth[cluster_id]:
                    _write(fh, accession)
            for accession in self.singletons:
                _write(fh, accession)
        return out_path

    def truth_membership(self) -> dict[str, str]:
        """{accession: cluster_id} — singletons get their own IDs."""
        out: dict[str, str] = {}
        for cluster_id, accessions in self.ground_truth.items():
            for accession in accessions:
                out[accession] = cluster_id
        for accession in self.singletons:
            out[accession] = f"SINGLETON_{accession}"
        return out


# ------------------------------------------------------------------
# public API
# ------------------------------------------------------------------


def simulate_clusters(
    *,
    n_clusters: int = 10,
    cluster_size: int | tuple[int, int] = (5, 15),
    n_singletons: int = 20,
    seq_length: int = DEFAULT_SEQ_LENGTH,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    seed: int = 42,
) -> SimulationResult:
    """Generate a synthetic dataset with `n_clusters` known clusters.

    `cluster_size` is either a fixed int or a `(lo, hi)` range; sizes
    are drawn uniformly in `[lo, hi]` when a tuple is passed.

    Returns a `SimulationResult` — call `.write_fasta(path)` to emit
    the alignment, or use `.truth_membership()` to score any clusterer.
    """
    rng = np.random.default_rng(seed)

    if isinstance(cluster_size, tuple):
        lo, hi = cluster_size
        if lo < 2 or hi < lo:
            raise ValueError(f"cluster_size range must be lo>=2 and hi>=lo; got {cluster_size}")
        sizes = rng.integers(lo, hi + 1, size=n_clusters).tolist()
    else:
        if cluster_size < 2:
            raise ValueError(f"cluster_size must be >= 2; got {cluster_size}")
        sizes = [cluster_size] * n_clusters

    logger.info(
        "simulating: %d clusters (sizes %s), %d singletons, seq_len=%d, "
        "mutation_rate=%.4f, seed=%d",
        n_clusters,
        _summarize_sizes(sizes),
        n_singletons,
        seq_length,
        mutation_rate,
        seed,
    )

    accession_to_seq: dict[str, str] = {}
    ground_truth: dict[str, list[str]] = {}
    accession_ix = 0

    for cluster_ix, size in enumerate(sizes):
        cluster_id = f"C-{cluster_ix:03d}"
        ground_truth[cluster_id] = []
        founder = _random_sequence(rng, seq_length)

        for _ in range(size):
            accession = f"S{accession_ix:06d}"
            accession_ix += 1
            member = _mutate(rng, founder, mutation_rate)
            accession_to_seq[accession] = member
            ground_truth[cluster_id].append(accession)

    singletons: list[str] = []
    for _ in range(n_singletons):
        accession = f"S{accession_ix:06d}"
        accession_ix += 1
        # Full-random singleton — should be far from every cluster's founder.
        accession_to_seq[accession] = _random_sequence(rng, seq_length)
        singletons.append(accession)

    return SimulationResult(
        accession_to_seq=accession_to_seq,
        ground_truth=ground_truth,
        singletons=singletons,
        n_clusters=n_clusters,
        seq_length=seq_length,
        seed=seed,
    )


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------

_ALPHABET = np.array(["A", "C", "G", "T"])


def _random_sequence(rng: np.random.Generator, length: int) -> str:
    """Uniformly random ACGT string of the given length."""
    return "".join(_ALPHABET[rng.integers(0, 4, size=length)])


def _mutate(rng: np.random.Generator, template: str, rate: float) -> str:
    """Per-site independent JC69 substitution."""
    if not 0 <= rate <= 1:
        raise ValueError(f"mutation rate must be in [0, 1]; got {rate}")
    arr = np.array(list(template))
    to_mutate = rng.random(len(arr)) < rate
    if not to_mutate.any():
        return template
    # For each mutated site, replace with one of the other 3 nucleotides.
    for i in np.where(to_mutate)[0]:
        current = arr[i]
        choices = [b for b in "ACGT" if b != current]
        arr[i] = choices[int(rng.integers(0, 3))]
    return "".join(arr)


def _summarize_sizes(sizes: list[int]) -> str:
    if not sizes:
        return "[]"
    if min(sizes) == max(sizes):
        return f"[all {sizes[0]}]"
    return f"[{min(sizes)}-{max(sizes)}, mean {sum(sizes)/len(sizes):.1f}]"
