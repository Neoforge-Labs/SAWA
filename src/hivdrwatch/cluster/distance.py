"""Pairwise nucleotide distance for cluster construction.

Two engines, same CSV shape (`ID1,ID2,Distance` header, then one row
per pair at or below the threshold — above-threshold pairs are
billions of useless zeros):

* `compute_pairwise_distances_python` — the production engine, used by
  `pipeline/cluster.py` and by the Fig 3 validation harness. Divides
  the mismatch count by the **pairwise unambiguous overlap** and
  refuses pairs that overlap in fewer than `DEFAULT_MIN_OVERLAP` sites.
* `run_tn93` — the pinned tn93 docker image, kept as the reference
  implementation for A/B comparison.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

DEFAULT_IMAGE = "hivdrwatch/tn93:1.0.17"
DEFAULT_THRESHOLD = 0.015  # HIV-TRACE convention (1.5%)
DEFAULT_AMBIGUITY_HANDLING = "resolve"  # tn93 flag

# Fewest shared unambiguous sites a pair must have before a distance is
# reported for it. Matches the `-l 100` default of tn93 itself.
DEFAULT_MIN_OVERLAP = 100

# The only characters that count toward a pairwise overlap. A gap is the
# absence of sequence and an IUPAC ambiguity code is a site the
# depositor could not call; neither is evidence that two sequences
# agree, so neither belongs in the denominator.
UNAMBIGUOUS_BASES = b"ACGT"


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

    # Rely on tn93's own defaults for overlap and output format:
    # `-t THRESHOLD -a AMBIGS -o OUTPUT INPUT`. This is exactly what
    # HIV-TRACE's tutorial invokes. Passing `-l 100 -f csv` explicitly
    # was measured to give the identical pair set on the filtered cohort
    # (only the row order moves, which is OpenMP scheduling), so the
    # flags buy nothing and the shorter line is the documented one.
    #
    # `-l` is the one place tn93's help and its code disagree: the help
    # says "at least OVERLAP nucleotides" but the comparison is
    # `totalNonGap <= min_overlap`, so 100 shared sites are refused and
    # 101 accepted. `DEFAULT_MIN_OVERLAP` uses the documented reading,
    # `overlap >= min_overlap`. The coverage gate leaves every surviving
    # record with >= 1,506 sites, so the one-site boundary cannot bind.
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


def _read_alignment(fasta_path: Path) -> tuple[list[str], np.ndarray]:
    """Parse an aligned FASTA into (accessions, [n, width] uint8 matrix)."""
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

    if not seqs:
        return accessions, np.empty((0, 0), dtype=np.uint8)

    # Length check — real alignments should be all-equal-length. Mixed
    # lengths point at a FASTA writer bug worth catching loudly.
    lengths = {int(s.size) for s in seqs}
    if len(lengths) > 1:
        raise ValueError(
            f"sequences have mixed lengths {sorted(lengths)}; the distance "
            "engine expects a proper alignment where every sequence is the "
            "same length"
        )
    return accessions, np.vstack(seqs)


def compute_pairwise_distances_python(
    fasta_path: Path,
    output_csv: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    min_overlap: int = DEFAULT_MIN_OVERLAP,
) -> Path:
    """In-process pairwise nucleotide distance over the pairwise overlap.

    This is the production distance engine: `pipeline/cluster.py`
    defaults `--distance-engine` to `python`, and `cluster_validate`
    (Fig 3) calls it directly. `run_tn93` is the reference path kept
    for A/B comparison. Earlier revisions of this docstring described
    the function as validation-harness-only; that stopped being true
    when the pipeline default was set, and the two now agree.

    ## The denominator is the overlap, never the alignment width

    A curated row is 3,012 columns wide whether or not any sequence was
    placed on them, so most rows are mostly dashes. Dividing the
    mismatch count by 3,012 treats a column where **neither** sequence
    has a base as evidence that the two agree there. Nobody sequenced
    it; it is not agreement, it is absence. Scored that way, two rows
    of pure gaps are distance 0.0 and link, and any two sequences are
    pulled toward 0 in proportion to how little of them overlaps — the
    metric rewards missing data.

    So: divide by the number of columns where both sequences carry an
    unambiguous base (A/C/G/T). Gap characters and IUPAC ambiguity
    codes (N, R, Y, …) leave the denominator — an uncalled site is not
    evidence of agreement either.

    ## Minimum overlap

    An overlap-aware distance from a handful of shared sites is a
    number with no information in it, so a pair whose overlap is below
    `min_overlap` gets no row at all rather than an unreliable one.
    The default matches the `-l 100` default of tn93 itself, whose
    semantics are the same: refuse the pair, do not clamp its distance.

    ## Relationship to TN93

    For sequences with no ambiguities, the Tamura-Nei 93 correction is
    a monotone function of the raw p-distance, and the clusterer only
    ever applies `distance <= threshold` to it, so cluster memberships
    are identical for the small divergences that matter here. What this
    function does **not** reproduce is tn93's `-a resolve` ambiguity
    handling: tn93 resolves an ambiguous site against its partner,
    where this drops the site. That is the conservative direction (it
    shrinks the denominator rather than inventing agreement) and it is
    recorded in `docs/evaluations/cluster-distance-denominator.md`.

    Emits a CSV with the same header (`ID1,ID2,Distance`) as tn93 so
    `read_distances` consumes it identically.
    """
    fasta_path = fasta_path.resolve()
    output_csv = output_csv.resolve()
    if not fasta_path.exists():
        raise FileNotFoundError(f"aligned FASTA missing: {fasta_path}")
    if min_overlap < 1:
        raise ValueError(f"min_overlap must be at least 1 site, got {min_overlap}")
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    accessions, matrix = _read_alignment(fasta_path)
    n = len(accessions)
    logger.info("computing in-process pairwise distances for %d sequences", n)

    if n < 2:
        output_csv.write_text("ID1,ID2,Distance\n")
        return output_csv

    unambiguous = np.isin(matrix, np.frombuffer(UNAMBIGUOUS_BASES, dtype=np.uint8))
    coverage = unambiguous.sum(axis=1)
    logger.info(
        "alignment width %d; unambiguous coverage per sequence: "
        "min=%d median=%d max=%d; %d sequences below min_overlap=%d",
        matrix.shape[1],
        int(coverage.min()),
        int(np.median(coverage)),
        int(coverage.max()),
        int((coverage < min_overlap).sum()),
        min_overlap,
    )

    n_written = 0
    n_refused = 0
    with output_csv.open("w", encoding="utf-8") as fh:
        fh.write("ID1,ID2,Distance\n")
        for i in range(n - 1):
            # An overlap can never exceed either sequence's own count of
            # unambiguous sites, so pairs that cannot possibly clear
            # min_overlap are dropped before any column-wise work.
            reachable = np.minimum(coverage[i + 1 :], coverage[i]) >= min_overlap
            n_refused += int((~reachable).sum())
            candidates = np.flatnonzero(reachable)
            if candidates.size == 0:
                continue

            shared = unambiguous[i + 1 :][candidates] & unambiguous[i]
            overlap = shared.sum(axis=1)
            keep = overlap >= min_overlap
            n_refused += int((~keep).sum())

            diffs = ((matrix[i + 1 :][candidates] != matrix[i]) & shared).sum(axis=1)
            distances = diffs / np.maximum(overlap, 1)
            for k in np.flatnonzero(keep & (distances <= threshold)):
                j = i + 1 + int(candidates[k])
                fh.write(f"{accessions[i]},{accessions[j]},{float(distances[k]):.6f}\n")
                n_written += 1

    logger.info(
        "in-process distances → %s (%d pairs at or below threshold=%.4f; "
        "%d pairs refused for overlap below %d sites)",
        output_csv,
        n_written,
        threshold,
        n_refused,
        min_overlap,
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
