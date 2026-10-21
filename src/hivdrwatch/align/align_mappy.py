"""In-process pairwise alignment via `mappy` (minimap2 Python bindings).

Why mappy?
----------
Benchmark on 1000 synthetic HIV pol sequences to HXB2:
  * MAFFT `--add --keeplength` in Docker: ~1 hour per 1000 seqs (thrashing)
  * Biopython `PairwiseAligner`: ~200 seconds per 1000 seqs
  * **mappy**: ~1.5 seconds per 1000 seqs (~700 seq/s single-thread)

That's 100–1000x faster than MAFFT and 100x faster than Biopython.
For 23k East Africa sequences, wall time is under a minute.

mappy uses minimap2's C implementation under the hood — k-mer seeding
+ chained SIMD Smith-Waterman for divergent regions. It's the
purpose-built tool for aligning short-ish sequences to a fixed
reference, which is exactly what we're doing.

Preset choice
-------------
HXB2 is subtype B; SAWA's cohort is A1/C/D and the CRFs. Measured over
the 36,510-sequence East Africa corpus, the pol mappy recovers sits
**~13% diverged** from HXB2 — past what `asm10` tolerates. `asm10`
discards every block it cannot chain inside that budget, which
truncated 68% of sequences by more than 150 nt and returned no
alignment at all for 52% of them. `asm20` recovers those at the same
base calls (99.98% agreement on columns both presets place) and costs
coverage on none. See `docs/evaluations/align-coverage.md`.

Lower the preset only for a deliberately same-subtype cohort.

Output semantics match MAFFT `--keeplength`: every returned sequence
is exactly `HXB2_POL_LENGTH` (3012) nucleotides, insertions vs HXB2
dropped, deletions become gaps.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any

import mappy as mp

from .reference import HXB2_POL_LENGTH

logger = logging.getLogger(__name__)

# `mappy` ships no type stubs, so one alignment record is opaque here.
# The fields we read off it are r_st, q_st, mlen, strand and cigar.
type MappyHit = Any

# Minimap2 preset — "asm20" allows up to ~20% divergence, which is what
# non-B pol against a subtype-B reference needs. "asm10" sits below the
# corpus's own ~13% divergence and silently drops sequence; "asm5" finds
# nothing outside subtype B at all.
DEFAULT_PRESET = "asm20"


def align_sequences_mappy(
    records: list[tuple[str, str]],
    reference_fasta: Path,
    *,
    preset: str = DEFAULT_PRESET,
    n_workers: int | None = None,  # kept for API parity; unused
) -> dict[str, str]:
    """Align every (accession, sequence) pair to the HXB2 reference.

    Returns `{accession: aligned_sequence}` where each aligned sequence
    is exactly `HXB2_POL_LENGTH` (3012) characters — insertions vs
    HXB2 are dropped (`--keeplength` semantics), deletions become gaps.

    Failed alignments (no significant match) get an all-gap output so
    downstream QC flags SHORT_COVERAGE without crashing.

    `n_workers` is accepted for API parity with the Biopython path but
    unused here — mappy's C core is fast enough single-thread that
    process-pool overhead would only slow things down.
    """
    n = len(records)
    if n == 0:
        return {}

    aligner = mp.Aligner(str(reference_fasta), preset=preset)
    if aligner is None or not aligner.seq_names:
        raise RuntimeError(
            f"mappy failed to index reference {reference_fasta}. Is it a valid FASTA?"
        )

    logger.info(
        "aligning %d sequences to HXB2 via mappy (preset=%s, ref_len=%d)",
        n,
        preset,
        HXB2_POL_LENGTH,
    )

    t0 = time.monotonic()
    aligned: dict[str, str] = {}
    n_failed = 0

    for accession, seq in records:
        query = seq.upper().replace("U", "T").replace("-", "").strip()
        if len(query) < 100:
            n_failed += 1
            aligned[accession] = "-" * HXB2_POL_LENGTH
            continue

        hits = _usable_hits(aligner, query)
        if not hits:
            n_failed += 1
            aligned[accession] = "-" * HXB2_POL_LENGTH
            continue

        aligned[accession] = _build_keeplength(query, hits, HXB2_POL_LENGTH)

    elapsed = time.monotonic() - t0
    logger.info(
        "mappy alignment complete: %d aligned, %d failed, in %.1fs (%.0f seq/s)",
        n - n_failed,
        n_failed,
        elapsed,
        n / max(elapsed, 1e-6),
    )
    return aligned


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------


def _usable_hits(aligner: mp.Aligner, query: str) -> list[MappyHit]:
    """Every forward-strand hit for this query, best first.

    A query that minimap2 splits into more than one block — which pol
    fragments deposited with an internal deletion, and whole-genome
    records, both do — yields several hits. Keeping only the largest
    throws the rest of the sequence away, so we keep them all and let
    `_build_keeplength` lay them down best-first.

    Reverse-strand hits are dropped rather than reverse-complemented:
    the CIGAR for those walks the complement of the query, so the
    forward-coordinate projection below would write garbage. A pol
    sequence that only matches HXB2 backwards is a deposit-orientation
    problem to surface, not to silently repair.
    """
    forward = [h for h in aligner.map(query) if h.strand == 1]
    return sorted(forward, key=lambda h: h.mlen, reverse=True)


def _build_keeplength(query_seq: str, hits: list[MappyHit], ref_length: int) -> str:
    """Project one or more alignment blocks onto reference coordinates.

    Walks each CIGAR:
      * M/=/X (0/7/8): consume both — copy query base to ref position
      * I (1): insertion in query — skip in ref (this is --keeplength)
      * D/N (2/3): deletion in query — gap in ref position
      * S (4): soft clip on query — skip
      * H (5): hard clip — nothing to consume

    `hits` is applied in order, and a block only writes where the
    reference column is still a gap, so the best hit wins any overlap
    and a weaker block can fill around it but never overwrite it.

    Returns exactly `ref_length` characters — positions no block
    reaches remain as gaps.
    """
    aligned = ["-"] * ref_length
    for hit in hits:
        _project_hit(query_seq, hit, aligned, ref_length)
    return "".join(aligned)


def _project_hit(query_seq: str, hit: MappyHit, aligned: list[str], ref_length: int) -> None:
    """Lay one block into `aligned`, writing only into gap columns."""
    q_pos = hit.q_st
    r_pos = hit.r_st

    for length, op in hit.cigar:
        if op == 0 or op == 7 or op == 8:  # M / = / X
            for i in range(length):
                r, q = r_pos + i, q_pos + i
                if r < ref_length and q < len(query_seq) and aligned[r] == "-":
                    aligned[r] = query_seq[q]
            r_pos += length
            q_pos += length
        elif op == 1:  # I — insertion in query, skip in ref (--keeplength)
            q_pos += length
        elif op == 2 or op == 3:  # D / N — gap in query
            r_pos += length
        elif op == 4:  # S — soft clip
            q_pos += length
        # H (5) — hard clip, nothing to advance
