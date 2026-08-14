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

For HIV specifically we use `preset="asm10"` which tolerates the
~10–15% intersubtype divergence typical between HIV-1 subtypes.

Output semantics match MAFFT `--keeplength`: every returned sequence
is exactly `HXB2_POL_LENGTH` (3012) nucleotides, insertions vs HXB2
dropped, deletions become gaps.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import mappy as mp

from .reference import HXB2_POL_LENGTH

logger = logging.getLogger(__name__)

# Minimap2 preset — "asm10" allows up to ~10% divergence, sufficient
# for intra-HIV alignment. Bump to "asm20" if you're aligning very
# divergent subtypes (e.g. SIV outliers). "asm5" for same-subtype only.
DEFAULT_PRESET = "asm10"


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

        best = _best_hit(aligner, query)
        if best is None:
            n_failed += 1
            aligned[accession] = "-" * HXB2_POL_LENGTH
            continue

        aligned[accession] = _build_keeplength(query, best, HXB2_POL_LENGTH)

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


def _best_hit(aligner: mp.Aligner, query: str):
    """Return the mapping hit with the most matches, or None."""
    best = None
    for hit in aligner.map(query):
        if best is None or hit.mlen > best.mlen:
            best = hit
    return best


def _build_keeplength(query_seq: str, hit, ref_length: int) -> str:
    """Project the query alignment onto reference coordinates.

    Walks the CIGAR:
      * M/=/X (0/7/8): consume both — copy query base to ref position
      * I (1): insertion in query — skip in ref (this is --keeplength)
      * D/N (2/3): deletion in query — gap in ref position
      * S (4): soft clip on query — skip
      * H (5): hard clip — nothing to consume

    Returns exactly `ref_length` characters — positions outside the
    aligned interval remain as gaps.
    """
    aligned = ["-"] * ref_length
    q_pos = hit.q_st
    r_pos = hit.r_st

    for length, op in hit.cigar:
        if op == 0 or op == 7 or op == 8:  # M / = / X
            for i in range(length):
                if r_pos + i < ref_length and q_pos + i < len(query_seq):
                    aligned[r_pos + i] = query_seq[q_pos + i]
            r_pos += length
            q_pos += length
        elif op == 1:  # I — insertion in query, skip in ref (--keeplength)
            q_pos += length
        elif op == 2 or op == 3:  # D / N — gap in query
            r_pos += length
        elif op == 4:  # S — soft clip
            q_pos += length
        # H (5) — hard clip, nothing to advance

    return "".join(aligned)
