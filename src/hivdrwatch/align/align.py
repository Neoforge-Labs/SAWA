"""MAFFT alignment stage — run inside the pinned Docker image.

We use `mafft --add newseqs.fasta --keeplength --reorder hxb2_pol.fasta`.
The `--keeplength` flag guarantees the output alignment has exactly the
same number of columns as the reference (HXB2 pol, 3012 nt), so every
downstream stage can address positions in HXB2 coordinates without
re-mapping. New sequences that don't fit lose their inserts vs HXB2 —
that's intentional; we don't want free-floating inserts in the feature
space.

Docker is invoked via `subprocess.run(["docker", "run", ...])`. Simple,
works, no extra dep. If we ever need Prefect-native container tasks with
resource limits and the UI showing container logs, we swap the executor
for `prefect_docker.DockerContainer` — same interface.
"""

from __future__ import annotations

import logging
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


def _default_thread_count() -> int:
    """Auto-detect a reasonable MAFFT thread count.

    On M1/M2 Macs `os.cpu_count()` returns the total core count
    (e.g. 8, 10, 12). We use N-2 to leave headroom for the OS +
    Prefect worker + other services, floored at 4 so small machines
    still work.
    """
    n = os.cpu_count() or 4
    return max(4, n - 2)

# Tag is `latest` (not a version pin) because the arm64 build path
# uses Debian's apt package (7.505) while amd64 uses the pinned 7.526 .deb.
# The exact version baked in is recorded at /mafft.version inside the image
# (see docker/mafft/Dockerfile).
DEFAULT_IMAGE = "hivdrwatch/mafft:latest"

# `mafft --add` builds an O(n²) distance matrix + per-sequence
# alloclen growth. Empirically 1000 seqs / chunk on ~40k East Africa
# sequences fits in 4 GB comfortably; 2000 is borderline once inserts
# force reallocation ("Reallocating..done. *alloclen = ...") past 30k.
DEFAULT_CHUNK_SIZE = 1000

# Hard memory ceiling passed to `docker run --memory`. If MAFFT hits
# this, the OOM killer stops the container immediately with exit 137
# instead of dragging the host into swap.
DEFAULT_MEMORY_LIMIT = "6g"

# MAFFT default is L-INS-i (2 tree iterations). For HIV pol vs a
# highly-conserved HXB2 reference, 1 iteration is indistinguishable
# in alignment quality but 2–4x faster. Empirically this is what
# HIV-TRACE and other HIV toolchains use.
DEFAULT_MAFFT_RETREE = 1

# `--maxiterate 0` skips iterative refinement. Per MAFFT's own docs,
# recommended for >10k-scale alignments. For our --add + --keeplength
# use case (reference-based), refinement is redundant.
DEFAULT_MAFFT_MAXITERATE = 0

# Run this many MAFFT invocations concurrently. Each uses ~memory_limit
# memory and cpu_count/concurrency threads. Default 2 assumes 16 GB
# Docker Desktop; bump to 4 if you have 32 GB allocated.
DEFAULT_CHUNK_CONCURRENCY = 2


@dataclass(frozen=True)
class MafftInvocation:
    """Records the exact docker + mafft invocation for provenance."""

    image: str
    args: tuple[str, ...]
    input_fasta: Path
    reference_fasta: Path
    output_fasta: Path


def write_fasta(records: list[tuple[str, str]], path: Path) -> Path:
    """Write (accession, sequence) pairs to a FASTA file.

    Accessions are used as headers unchanged — they must be unique and
    contain no whitespace. Sequences have gaps stripped since we're
    aligning ungapped inputs against the reference alignment.
    """
    seen: set[str] = set()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for accession, seq in records:
            if accession in seen:
                raise ValueError(f"duplicate accession in FASTA input: {accession}")
            seen.add(accession)
            ungapped = seq.replace("-", "").upper()
            fh.write(f">{accession}\n{ungapped}\n")
    logger.info("wrote %d sequences → %s", len(seen), path)
    return path


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA into {accession: sequence} (gap-preserving)."""
    out: dict[str, str] = {}
    current: str | None = None
    buffer: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    out[current] = "".join(buffer).upper()
                current = line[1:].strip().split()[0]  # first token as accession
                buffer = []
            else:
                buffer.append(line.strip())
    if current is not None:
        out[current] = "".join(buffer).upper()
    return out


def run_mafft(
    input_fasta: Path,
    reference_fasta: Path,
    output_fasta: Path,
    *,
    image: str = DEFAULT_IMAGE,
    thread_count: int | None = None,
    memory_limit: str = DEFAULT_MEMORY_LIMIT,
) -> MafftInvocation:
    """Run MAFFT --add via docker; return the invocation record for provenance.

    All paths are made absolute before mounting into the container so
    relative-path resolution can't bite us. We mount two host directories
    read-only for reference and writable for the working area.

    `memory_limit` is a hard ceiling passed to `docker run --memory` —
    prevents runaway MAFFT invocations from OOM'ing the host. Container
    exit code will be 137 (SIGKILL) if hit.
    """
    input_fasta = input_fasta.resolve()
    reference_fasta = reference_fasta.resolve()
    output_fasta = output_fasta.resolve()

    if thread_count is None:
        thread_count = _default_thread_count()

    if not input_fasta.exists():
        raise FileNotFoundError(f"input FASTA missing: {input_fasta}")
    if not reference_fasta.exists():
        raise FileNotFoundError(f"reference FASTA missing: {reference_fasta}")

    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # Working directory the container sees at /work; reference at /ref
    work_dir = input_fasta.parent
    ref_dir = reference_fasta.parent
    if output_fasta.parent != work_dir:
        raise ValueError(
            f"output must live alongside input for docker volume mount "
            f"(input parent {work_dir} != output parent {output_fasta.parent})"
        )

    mafft_args = (
        "--add",
        f"/work/{input_fasta.name}",
        "--keeplength",
        "--reorder",
        # --retree 1 halves the guide-tree cost with no meaningful
        # accuracy hit on reference-based HIV pol alignment.
        "--retree",
        str(DEFAULT_MAFFT_RETREE),
        # --maxiterate 0 skips iterative refinement. Redundant for
        # --add + --keeplength; recommended by MAFFT for large sets.
        "--maxiterate",
        str(DEFAULT_MAFFT_MAXITERATE),
        "--thread",
        str(thread_count),
        f"/ref/{reference_fasta.name}",
    )
    docker_cmd = [
        "docker", "run", "--rm",
        "--memory", memory_limit,
        # Match CPU quota to thread_count so MAFFT can't oversubscribe.
        "--cpus", str(thread_count),
        "-v", f"{work_dir}:/work",
        "-v", f"{ref_dir}:/ref:ro",
        image,
        *mafft_args,
    ]
    logger.info("running MAFFT: %s", " ".join(docker_cmd))

    with output_fasta.open("wb") as fh:
        result = subprocess.run(
            docker_cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            check=False,
        )
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")[-2000:]
        hint = ""
        if result.returncode == 137:
            hint = (
                f"\n\nExit 137 = SIGKILL, likely OOM against --memory {memory_limit}. "
                f"Try a smaller --chunk-size or raise DEFAULT_MEMORY_LIMIT in "
                f"curate/align.py."
            )
        raise RuntimeError(
            f"MAFFT failed (exit {result.returncode}):\n{stderr}{hint}"
        )

    logger.info("MAFFT complete → %s (%d bytes)", output_fasta, output_fasta.stat().st_size)
    return MafftInvocation(
        image=image,
        args=mafft_args,
        input_fasta=input_fasta,
        reference_fasta=reference_fasta,
        output_fasta=output_fasta,
    )


# Aligner selection —
#   "mappy"      → in-process minimap2 via Python bindings (fastest, default)
#   "biopython"  → Bio.Align.PairwiseAligner (portable fallback, ~100x slower)
#   "mafft"      → Docker MAFFT --add --keeplength (only for MSA needs)
DEFAULT_ALIGNER = "mappy"


def align_sequences(
    records: list[tuple[str, str]],
    reference_fasta: Path,
    work_dir: Path,
    *,
    aligner: str = DEFAULT_ALIGNER,
    image: str = DEFAULT_IMAGE,
    thread_count: int | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    chunk_concurrency: int = DEFAULT_CHUNK_CONCURRENCY,
) -> dict[str, str]:
    """Dispatch to the requested aligner backend.

    * `mappy` (default): in-process minimap2 bindings. ~700 seq/s
      single-thread — the whole East Africa footprint (~25k
      sequences) aligns in under a minute.
    * `biopython`: portable Python fallback. ~5 seq/s. Only use if
      mappy fails to install for some reason.
    * `mafft`: Docker MAFFT --add --keeplength. Only use when you
      specifically need joint-guide-tree MSA across the added
      sequences (rare for reference-based work).
    """
    if aligner == "mappy":
        from .align_mappy import align_sequences_mappy

        return align_sequences_mappy(records, reference_fasta)
    if aligner == "biopython":
        from .align_biopython import align_sequences_biopython

        return align_sequences_biopython(
            records,
            reference_fasta,
            n_workers=thread_count,
        )
    if aligner != "mafft":
        raise ValueError(
            f"unknown aligner {aligner!r}; expected 'mappy', 'biopython' or 'mafft'"
        )
    return _align_sequences_mafft(
        records,
        reference_fasta,
        work_dir,
        image=image,
        thread_count=thread_count,
        chunk_size=chunk_size,
        chunk_concurrency=chunk_concurrency,
    )


def _align_sequences_mafft(
    records: list[tuple[str, str]],
    reference_fasta: Path,
    work_dir: Path,
    *,
    image: str = DEFAULT_IMAGE,
    thread_count: int | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    chunk_concurrency: int = DEFAULT_CHUNK_CONCURRENCY,
) -> dict[str, str]:
    """End-to-end MAFFT path: chunked --add via Docker.

    Kept for cases that specifically need MAFFT-quality MSA (rare for
    the DR-detection pipeline). Chunks avoid the O(n²) distance-matrix
    OOM at large n. Because `--keeplength` anchors each chunk's output
    to the same HXB2 column coordinates, chunk outputs are drop-in
    concatenable.

    Returns `{accession: aligned_sequence}`; the HXB2 reference row is
    dropped from the return dict since it's already the coordinate anchor.
    """
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")

    # Dedupe input on accession. The upstream ingest merge occasionally
    # emits the same accession twice (e.g. LANL rows whose `Accession`
    # column has full FASTA-header strings vs GenBank rows with the same
    # underlying accession). We can't have duplicates because MAFFT
    # would then emit repeated headers in chunk outputs — the assertion
    # below would fire spuriously. Keep first occurrence; log the count
    # so the root cause is investigable.
    original_n = len(records)
    seen_input: dict[str, str] = {}
    for acc, seq in records:
        if acc not in seen_input:
            seen_input[acc] = seq
    records = list(seen_input.items())
    dropped_dupes = original_n - len(records)
    if dropped_dupes:
        logger.warning(
            "input had %d duplicate accessions; dedup kept first occurrence each",
            dropped_dupes,
        )

    if thread_count is None:
        # Divide available threads across concurrent chunks so parallel
        # MAFFT processes don't oversubscribe cores.
        total = _default_thread_count()
        thread_count = max(2, total // max(1, chunk_concurrency))

    n = len(records)
    n_chunks = (n + chunk_size - 1) // chunk_size
    logger.info(
        "aligning %d unique sequences in %d chunk(s) of ≤ %d "
        "(chunk_concurrency=%d, threads/chunk=%d, memory limit=%s, --retree %d)",
        n,
        n_chunks,
        chunk_size,
        chunk_concurrency,
        thread_count,
        DEFAULT_MEMORY_LIMIT,
        DEFAULT_MAFFT_RETREE,
    )

    # Pre-write every chunk's input FASTA so the parallel workers only
    # need to invoke MAFFT + read output. Cheap sequential I/O keeps the
    # parallel section clean.
    chunk_plan: list[tuple[int, Path, Path, int]] = []
    for chunk_ix in range(n_chunks):
        start = chunk_ix * chunk_size
        stop = min(start + chunk_size, n)
        chunk_records = records[start:stop]
        input_fasta = work_dir / f"input_{chunk_ix:04d}.fasta"
        output_fasta = work_dir / f"aligned_{chunk_ix:04d}.fasta"
        write_fasta(chunk_records, input_fasta)
        chunk_plan.append((chunk_ix, input_fasta, output_fasta, len(chunk_records)))

    aligned: dict[str, str] = {}
    completed = 0
    total_start = time.monotonic()

    def _run_one(job: tuple[int, Path, Path, int]) -> tuple[int, dict[str, str], float]:
        chunk_ix, input_fasta, output_fasta, n_seqs = job
        t0 = time.monotonic()
        run_mafft(
            input_fasta=input_fasta,
            reference_fasta=reference_fasta,
            output_fasta=output_fasta,
            image=image,
            thread_count=thread_count,
        )
        elapsed = time.monotonic() - t0
        chunk_aligned = read_fasta(output_fasta)
        for h in [k for k in chunk_aligned if _is_reference_header(k)]:
            chunk_aligned.pop(h, None)
        return chunk_ix, chunk_aligned, elapsed

    with ThreadPoolExecutor(max_workers=chunk_concurrency) as pool:
        futures = {pool.submit(_run_one, job): job for job in chunk_plan}
        for future in as_completed(futures):
            chunk_ix, chunk_aligned, elapsed = future.result()
            completed += 1
            n_seqs = len(chunk_aligned)
            logger.info(
                "chunk %d/%d done in %.1fs (%.2f seq/s) — %d/%d chunks complete",
                chunk_ix + 1,
                n_chunks,
                elapsed,
                n_seqs / max(elapsed, 1e-6),
                completed,
                n_chunks,
            )
            if elapsed > 300:
                logger.warning(
                    "chunk %d took %.1fs — MAFFT is thrashing. Try "
                    "`CHUNK_SIZE=500 make curate` or lower CHUNK_CONCURRENCY.",
                    chunk_ix + 1,
                    elapsed,
                )
            overlap = set(aligned) & set(chunk_aligned)
            if overlap:
                raise RuntimeError(
                    f"chunk {chunk_ix} produced duplicate accessions: "
                    f"{sorted(overlap)[:5]}"
                )
            aligned.update(chunk_aligned)

    total_elapsed = time.monotonic() - total_start
    logger.info(
        "alignment complete: %d sequences in %.1fs wall (%.2f seq/s)",
        len(aligned),
        total_elapsed,
        len(aligned) / max(total_elapsed, 1e-6),
    )

    # Concatenate all chunk FASTAs into a single `aligned.fasta` for
    # downstream tools (tn93) that expect one file per dataset.
    combined = work_dir / "aligned.fasta"
    _concatenate_chunk_fastas(work_dir, combined, n_chunks)
    logger.info("wrote combined alignment (%d sequences) → %s", len(aligned), combined)

    return aligned


def _is_reference_header(header: str) -> bool:
    u = header.upper()
    return "HXB2" in u or u.startswith("K03455")


def _concatenate_chunk_fastas(
    work_dir: Path, combined_path: Path, n_chunks: int
) -> None:
    """Merge per-chunk aligned FASTAs into one file, keeping the reference
    from chunk 0 only. Chunk files are left in place for debugging."""
    with combined_path.open("wb") as out_fh:
        for chunk_ix in range(n_chunks):
            chunk_path = work_dir / f"aligned_{chunk_ix:04d}.fasta"
            if not chunk_path.exists():
                continue
            skipping_reference = chunk_ix > 0  # reference appears once
            in_reference_record = False
            with chunk_path.open("rb") as in_fh:
                for line in in_fh:
                    if line.startswith(b">"):
                        header = line.decode("utf-8", errors="replace")[1:].strip()
                        in_reference_record = _is_reference_header(header)
                        if skipping_reference and in_reference_record:
                            continue
                        out_fh.write(line)
                    else:
                        if skipping_reference and in_reference_record:
                            continue
                        out_fh.write(line)
