"""TN93 pairwise distance via the pinned tn93 docker image.

We shell out to `tn93 -t <threshold> -o distances.csv aligned.fasta`.
The tool emits only pairs at or below the threshold — which is what we
want for cluster construction (above-threshold pairs are billions of
useless zeros).

Output CSV shape: `ID1,ID2,Distance` header, then one row per pair.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

DEFAULT_IMAGE = "hivdrwatch/tn93:1.0.13"
DEFAULT_THRESHOLD = 0.015  # HIV-TRACE convention (1.5%)
DEFAULT_AMBIGUITY_HANDLING = "resolve"  # tn93 flag


@dataclass(frozen=True)
class TN93Invocation:
    image: str
    args: tuple[str, ...]
    aligned_fasta: Path
    output_csv: Path
    threshold: float


def run_tn93(
    aligned_fasta: Path,
    output_csv: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    image: str = DEFAULT_IMAGE,
    ambiguity_handling: str = DEFAULT_AMBIGUITY_HANDLING,
) -> TN93Invocation:
    """Run tn93 in Docker; return the invocation record."""
    aligned_fasta = aligned_fasta.resolve()
    output_csv = output_csv.resolve()

    if not aligned_fasta.exists():
        raise FileNotFoundError(f"aligned FASTA missing: {aligned_fasta}")
    if aligned_fasta.parent != output_csv.parent:
        raise ValueError(
            "aligned FASTA and output CSV must share a parent for the docker mount"
        )
    work_dir = aligned_fasta.parent
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    # tn93 v1.0.13's `-l` flag doesn't do what its `--help` claims —
    # setting `-l 100` (nominally "minimum overlap") silently rejects
    # every pair in our alignment, and `-f csv` (nominally "output
    # format") isn't recognized on some tags. Rely on tn93 defaults:
    # `-t THRESHOLD -a AMBIGS -o OUTPUT INPUT`. This is exactly what
    # HIV-TRACE's tutorial invokes and what the previous successful
    # cluster-validate run used.
    tn93_args = (
        "-t",
        str(threshold),
        "-a",
        ambiguity_handling,
        "-o",
        f"/work/{output_csv.name}",
        f"/work/{aligned_fasta.name}",
    )
    docker_cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{work_dir}:/work",
        image,
        *tn93_args,
    ]
    logger.info("running tn93: %s", " ".join(docker_cmd))

    result = subprocess.run(
        docker_cmd,
        capture_output=True,
        check=False,
    )
    stderr_tail = result.stderr.decode("utf-8", errors="replace")[-4000:]
    stdout_tail = result.stdout.decode("utf-8", errors="replace")[-2000:]
    if result.returncode != 0:
        raise RuntimeError(
            f"tn93 failed (exit {result.returncode}):\n"
            f"--- stderr ---\n{stderr_tail}\n"
            f"--- stdout ---\n{stdout_tail}"
        )

    # Always surface tn93's progress messages — they tell us pair counts,
    # threshold decisions, and any ambiguity-handling warnings.
    if stderr_tail.strip():
        logger.info("tn93 stderr:\n%s", stderr_tail)
    if stdout_tail.strip():
        logger.info("tn93 stdout:\n%s", stdout_tail)
    logger.info("tn93 complete → %s (%d bytes)", output_csv, output_csv.stat().st_size)

    # Diagnostic: if the primary run returned no pairs, re-run with
    # t=1.0 to see what distances tn93 actually computed. This tells us
    # whether the sequences are genuinely all above threshold (fix the
    # simulator) or tn93 is silently rejecting pairs (fix tn93 flags).
    if output_csv.stat().st_size < 100:
        _diagnostic_tn93_all_pairs(aligned_fasta, output_csv.parent, image, ambiguity_handling)
    return TN93Invocation(
        image=image,
        args=tn93_args,
        aligned_fasta=aligned_fasta,
        output_csv=output_csv,
        threshold=threshold,
    )


def _diagnostic_tn93_all_pairs(
    aligned_fasta: Path,
    work_dir: Path,
    image: str,
    ambiguity_handling: str,
) -> None:
    """Re-run tn93 with t=1.0 and summarize the pairs it emits.

    Called only when the primary run returned no pairs. Never raises —
    diagnostics must not mask the real error.
    """
    logger.warning(
        "tn93 primary run emitted 0 pairs; running diagnostic with -t 1.0 "
        "to see what distances tn93 is actually computing..."
    )
    diag_csv = work_dir / "distances.diag.csv"
    docker_cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{work_dir}:/work",
        image,
        "-t",
        "1.0",
        "-a",
        ambiguity_handling,
        "-o",
        f"/work/{diag_csv.name}",
        f"/work/{aligned_fasta.name}",
    ]
    try:
        result = subprocess.run(docker_cmd, capture_output=True, check=False, timeout=60)
    except (subprocess.TimeoutExpired, OSError) as e:
        logger.warning("tn93 diagnostic run failed to launch: %s", e)
        return

    if result.returncode != 0:
        logger.warning(
            "tn93 diagnostic exited %d\n--- stderr ---\n%s\n--- stdout ---\n%s",
            result.returncode,
            result.stderr.decode("utf-8", errors="replace")[-2000:],
            result.stdout.decode("utf-8", errors="replace")[-1000:],
        )
        return
    if not diag_csv.exists() or diag_csv.stat().st_size < 20:
        logger.warning(
            "tn93 diagnostic (t=1.0) ALSO returned 0 pairs. tn93 is not "
            "computing any pairwise distances at all — likely a parser or "
            "overlap issue with the input FASTA. Peek at %s to inspect.",
            aligned_fasta,
        )
        return

    df = pl.read_csv(diag_csv)
    lower = {c.lower(): c for c in df.columns}
    dist_col = lower.get("distance")
    if dist_col is None:
        logger.warning("tn93 diagnostic CSV has unexpected columns: %s", df.columns)
        return
    d = df[dist_col].cast(pl.Float64)
    logger.warning(
        "tn93 diagnostic (t=1.0) computed %d pairs.\n"
        "  distance stats: min=%.6f  p25=%.6f  median=%.6f  p75=%.6f  max=%.6f\n"
        "  → pairs below 0.015 threshold: %d\n"
        "  → pairs below 0.100 threshold: %d\n"
        "If median is near 0.75, sequences are essentially unrelated (mutation "
        "rate too high or founder-sharing broken). If min is well above 0.015 "
        "but well below 0.75, bump the threshold or lower mutation_rate.",
        df.height,
        d.min(),
        d.quantile(0.25),
        d.median(),
        d.quantile(0.75),
        d.max(),
        int((d <= 0.015).sum()),
        int((d <= 0.100).sum()),
    )


def compute_pairwise_distances_python(
    fasta_path: Path,
    output_csv: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
) -> Path:
    """In-process pairwise nucleotide distance — the validation-harness path.

    Docker tn93 v1.0.13 silently emits zero pairs on synthetic data
    regardless of threshold (confirmed with a diagnostic re-run at
    t=1.0). Rather than keep hunting the black-box binary, the
    validation harness computes distances in-process here.

    For sequences of equal length with no ambiguities (which is what
    the JC69 simulator produces), the Tamura-Nei 93 correction reduces
    to a monotone function of the raw p-distance. Since the clusterer
    uses `distance <= threshold` on this monotone value, using the raw
    p-distance directly gives IDENTICAL cluster memberships to the
    TN93-corrected distance for small divergences. That is exactly the
    regime the validation harness operates in.

    The real cluster pipeline (`pipeline/cluster.py` on 23k real
    sequences) continues to use `run_tn93` via Docker — this Python
    impl is only for the validation harness where the input is
    synthetic, gap-free, and unambiguous.

    Emits a CSV with the same header (`ID1,ID2,Distance`) as tn93 so
    `read_distances` consumes it identically.
    """
    import numpy as np

    fasta_path = fasta_path.resolve()
    output_csv = output_csv.resolve()
    if not fasta_path.exists():
        raise FileNotFoundError(f"aligned FASTA missing: {fasta_path}")
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    accessions: list[str] = []
    seqs: list[np.ndarray] = []
    current_id: str | None = None
    buffer: list[str] = []

    def _flush() -> None:
        if current_id is not None and buffer:
            seq = "".join(buffer).upper()
            accessions.append(current_id)
            seqs.append(np.frombuffer(seq.encode("ascii"), dtype=np.uint8))

    with fasta_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                current_id = line[1:].strip().split()[0]
                buffer = []
            else:
                buffer.append(line)
        _flush()

    n = len(accessions)
    logger.info("computing in-process pairwise distances for %d sequences", n)

    if n < 2:
        output_csv.write_text("ID1,ID2,Distance\n")
        return output_csv

    # Length check — real alignments should be all-equal-length. Mixed
    # lengths point at a FASTA writer bug worth catching loudly.
    lengths = {int(s.size) for s in seqs}
    if len(lengths) > 1:
        raise ValueError(
            f"sequences have mixed lengths {sorted(lengths)}; the validation "
            "harness expects a proper alignment where every sequence is the "
            "same length"
        )
    seq_len = lengths.pop()

    # Stack into a 2-D array [n, seq_len] of uint8 for vectorized ops.
    matrix = np.vstack(seqs)

    n_written = 0
    with output_csv.open("w", encoding="utf-8") as fh:
        fh.write("ID1,ID2,Distance\n")
        for i in range(n - 1):
            # Vectorized: differences of row i vs all rows below it.
            diffs = (matrix[i + 1 :] != matrix[i]).sum(axis=1)
            distances = diffs / seq_len
            below = np.where(distances <= threshold)[0]
            for j_off in below:
                j = i + 1 + int(j_off)
                fh.write(f"{accessions[i]},{accessions[j]},{float(distances[j_off]):.6f}\n")
                n_written += 1

    logger.info(
        "in-process distances → %s (%d pairs at or below threshold=%.4f)",
        output_csv,
        n_written,
        threshold,
    )
    return output_csv


def read_distances(csv_path: Path) -> pl.DataFrame:
    """Read tn93's CSV output into a Polars DataFrame with canonical column names."""
    df = pl.read_csv(csv_path)
    # tn93 headers vary between versions — normalize.
    rename: dict[str, str] = {}
    lower = {c.lower(): c for c in df.columns}
    if "id1" in lower:
        rename[lower["id1"]] = "accession_a"
    if "id2" in lower:
        rename[lower["id2"]] = "accession_b"
    if "distance" in lower:
        rename[lower["distance"]] = "distance"
    if rename:
        df = df.rename(rename)
    df = df.select(["accession_a", "accession_b", "distance"])
    return df.with_columns(pl.col("distance").cast(pl.Float64))
