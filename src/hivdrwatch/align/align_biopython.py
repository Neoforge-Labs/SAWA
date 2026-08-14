"""In-process pairwise alignment via Biopython.

Why not MAFFT?
--------------
MAFFT `--add --keeplength` builds a joint guide-tree across all added
sequences. That's mathematically correct for MSA but wasteful for our
use case: we're aligning each sequence *independently* to HXB2 pol,
which is a highly-conserved ~3kb reference. Every added sequence has
the same reference; there's no mutual-context benefit.

Pairwise global alignment against HXB2 per sequence:
  * Trivially parallelizable — no shared state between sequences
  * Runs in ~5–20 ms per sequence via Biopython's C-backed aligner
  * No Docker container startup overhead
  * No cross-arch install pain
  * Same `--keeplength` semantics (extract reference-aligned positions)

Estimated wall time for the East Africa footprint (~25k sequences):
  * Single-threaded: ~5 min
  * With ProcessPoolExecutor and 8 cores: ~45s

Compare to Docker-MAFFT at chunk_size=1000, concurrency=2: hours.

Fallback: if a query has too much divergence for pairwise (rare for
HIV pol), the sequence gets a warning + all-gap alignment so downstream
QC drops it as SHORT_COVERAGE.
"""

from __future__ import annotations

import logging
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from Bio import Align

from .reference import HXB2_POL_LENGTH

logger = logging.getLogger(__name__)

# Batch size for ProcessPool submissions. Small enough to give good
# load-balancing across workers, large enough to amortize per-job
# process communication cost.
DEFAULT_BATCH_SIZE = 200

# Minimum alignment score fraction (score / max_possible_score) below
# which we consider the query too diverged from HXB2 to trust.
MIN_ALIGNMENT_SCORE_FRAC = 0.5


@dataclass(frozen=True)
class AlignmentStats:
    n_aligned: int
    n_failed: int
    wall_seconds: float


# ------------------------------------------------------------------
# public API
# ------------------------------------------------------------------


def align_sequences_biopython(
    records: list[tuple[str, str]],
    reference_fasta: Path,
    *,
    n_workers: int | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> dict[str, str]:
    """Align every (accession, sequence) pair to the HXB2 reference.

    Returns `{accession: aligned_sequence}` where each aligned sequence
    is exactly `HXB2_POL_LENGTH` (3012) characters — insertions vs
    HXB2 are dropped (`--keeplength` semantics), deletions become gaps.

    Failed alignments (score too low, sequence too short) are logged
    and get an all-gap output so downstream QC flags SHORT_COVERAGE
    without crashing.
    """
    if n_workers is None:
        n_workers = max(1, (os.cpu_count() or 4) - 2)

    ref_seq = _read_reference(reference_fasta)
    n = len(records)
    if n == 0:
        return {}

    logger.info(
        "aligning %d sequences pairwise to HXB2 (workers=%d, batch=%d, ref_len=%d)",
        n,
        n_workers,
        batch_size,
        len(ref_seq),
    )

    t0 = time.monotonic()
    aligned: dict[str, str] = {}
    n_failed = 0

    # Split into batches for ProcessPoolExecutor.
    batches = [
        records[i : i + batch_size] for i in range(0, n, batch_size)
    ]

    # If we only have one core or one batch, skip the process pool
    # overhead and run in the calling process.
    if n_workers <= 1 or len(batches) == 1:
        for batch in batches:
            for acc, aln in _align_batch(batch, ref_seq):
                if aln is None:
                    n_failed += 1
                    aligned[acc] = "-" * HXB2_POL_LENGTH
                else:
                    aligned[acc] = aln
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(_align_batch, batch, ref_seq): i
                for i, batch in enumerate(batches)
            }
            completed = 0
            for future in as_completed(futures):
                completed += 1
                for acc, aln in future.result():
                    if aln is None:
                        n_failed += 1
                        aligned[acc] = "-" * HXB2_POL_LENGTH
                    else:
                        aligned[acc] = aln
                if completed % 10 == 0 or completed == len(batches):
                    logger.info(
                        "batches complete: %d/%d (%d sequences aligned so far)",
                        completed,
                        len(batches),
                        len(aligned),
                    )

    elapsed = time.monotonic() - t0
    logger.info(
        "alignment complete: %d aligned, %d failed, in %.1fs (%.1f seq/s)",
        n - n_failed,
        n_failed,
        elapsed,
        n / max(elapsed, 1e-6),
    )
    return aligned


# ------------------------------------------------------------------
# worker (module-level so ProcessPoolExecutor can pickle it)
# ------------------------------------------------------------------


def _align_batch(
    records: list[tuple[str, str]], reference: str
) -> list[tuple[str, str | None]]:
    """Align one batch of sequences. Runs in a worker process.

    Reference is passed by value (each worker gets its own copy) — the
    reference is ~3 KB so this is cheap.
    """
    aligner = _make_aligner()
    ref_upper = reference.upper()
    max_score = len(ref_upper) * aligner.match_score
    threshold = max_score * MIN_ALIGNMENT_SCORE_FRAC

    out: list[tuple[str, str | None]] = []
    for accession, seq in records:
        query = _clean_sequence(seq)
        if len(query) < 100:
            out.append((accession, None))
            continue
        try:
            alignments = aligner.align(ref_upper, query)
            best = alignments[0]
        except (ValueError, Exception):  # noqa: BLE001
            out.append((accession, None))
            continue
        if best.score < threshold:
            out.append((accession, None))
            continue
        aligned = _extract_keeplength(best, len(ref_upper))
        out.append((accession, aligned))
    return out


def _make_aligner() -> Align.PairwiseAligner:
    """Configure a global aligner tuned for HIV pol vs HXB2.

    Match / mismatch scores loosely modeled on EDNAFULL nucleotide
    matrix used by MAFFT for DNA — they should give qualitatively
    similar alignment to `mafft --keeplength`.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    # Penalize long gaps more than short ones (affine gap).
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    # End gaps are cheap so we don't penalize truncated queries.
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    return aligner


def _extract_keeplength(alignment, ref_length: int) -> str:
    """Extract the query aligned to reference coordinates.

    Iterates the alignment; keeps a query character for every reference
    position that isn't a gap. Reference-gap columns (query insertions)
    are dropped — same semantic as MAFFT `--keeplength`.

    Biopython's Alignment object exposes `.aligned` as two arrays of
    coordinate ranges, one per sequence. We rebuild the projected
    query using those coordinates.
    """
    # Fast path via .aligned attribute (available on Biopython >= 1.80)
    try:
        ref_blocks, qry_blocks = alignment.aligned
    except AttributeError:
        # Very old Biopython — fall back to string parsing.
        return _extract_keeplength_from_str(str(alignment), ref_length)

    aligned = ["-"] * ref_length
    # `alignment[1]` is the query aligned string with gaps; we can
    # index by (ref_position - gap_offset) to fetch chars.
    # But blocks give us the exact intervals; simpler to walk them.
    query_seq = alignment.query
    for (r_start, r_end), (q_start, q_end) in zip(ref_blocks, qry_blocks, strict=True):
        block_len = r_end - r_start
        if block_len != (q_end - q_start):
            # Shouldn't happen from .aligned but guard against it.
            continue
        query_slice = str(query_seq)[q_start:q_end]
        for i, ch in enumerate(query_slice):
            aligned[r_start + i] = ch
    return "".join(aligned)


def _extract_keeplength_from_str(alignment_str: str, ref_length: int) -> str:
    """Fallback that parses the two-line alignment string representation."""
    lines = alignment_str.strip().split("\n")
    if len(lines) < 2:
        return "-" * ref_length
    ref_line = lines[0]
    qry_line = lines[-1]
    aligned = []
    for r_char, q_char in zip(ref_line, qry_line, strict=False):
        if r_char == "-":
            continue  # reference gap → insertion in query → drop
        aligned.append(q_char)
    result = "".join(aligned)
    if len(result) < ref_length:
        result += "-" * (ref_length - len(result))
    return result[:ref_length]


def _clean_sequence(seq: str) -> str:
    """Uppercase, U→T, strip gaps and whitespace."""
    return seq.upper().replace("U", "T").replace("-", "").strip()


def _read_reference(reference_fasta: Path) -> str:
    """Load the HXB2 pol reference as an ungapped uppercase string."""
    buffer: list[str] = []
    with reference_fasta.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            buffer.append(line.strip())
    ref = "".join(buffer).replace("-", "").upper()
    if len(ref) != HXB2_POL_LENGTH:
        logger.warning(
            "reference length %d != expected HXB2_POL_LENGTH %d — "
            "using as-is; slice coordinates may be off",
            len(ref),
            HXB2_POL_LENGTH,
        )
    return ref
