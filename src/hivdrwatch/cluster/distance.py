"""Pairwise nucleotide distance for cluster construction.

Two engines, same CSV shape (`ID1,ID2,Distance` header, then one row
per pair at or below the threshold — above-threshold pairs are
billions of useless zeros):

* `compute_pairwise_distances_python` — the production engine, used by
  `pipeline/cluster.py` and by the Fig 3 validation harness. Scores
  every pair over the **pairwise unambiguous overlap** and refuses
  pairs that overlap in fewer than `DEFAULT_MIN_OVERLAP` sites. It can
  score that overlap two ways: the uncorrected mismatch fraction
  (`p-distance`, the default) or the Tamura-Nei 1993 correction
  (`tn93`), which is what HIV-TRACE thresholds.
* `run_tn93` — the pinned tn93 docker image, kept as the reference
  implementation for A/B comparison.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

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
#
# The order is tn93's own (`ValidChars` in `src/tn93_shared.cc` opens
# "ACGT"), and `_tn93_distances` indexes its base frequencies by it, so
# it is load-bearing rather than alphabetical by accident.
UNAMBIGUOUS_BASES = b"ACGT"

# How a pair's overlap is scored. `p-distance` is the uncorrected
# mismatch fraction and the only metric that has ever shipped; `tn93`
# applies the Tamura-Nei 1993 correction to the same overlap, which is
# the convention HIV-TRACE — and Paper 1's methods section — thresholds.
# Which of the two Paper 1 should claim is Kelyn's call, so the default
# is unchanged and the correction is reachable only by asking for it
# (`docs/evaluations/cluster-distance-denominator.md` §12, §13).
type DistanceMetric = Literal["p-distance", "tn93"]
DISTANCE_METRICS: tuple[DistanceMetric, ...] = ("p-distance", "tn93")
DEFAULT_DISTANCE_METRIC: DistanceMetric = "p-distance"

# Base code for a column that is a gap or an IUPAC ambiguity code.
_NOT_A_BASE = 255


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


def _base_codes(matrix: np.ndarray) -> np.ndarray:
    """Map an ASCII alignment onto base indices: A/C/G/T → 0/1/2/3, else 255."""
    codes = np.full(matrix.shape, _NOT_A_BASE, dtype=np.uint8)
    for index, base in enumerate(UNAMBIGUOUS_BASES):
        codes[matrix == base] = index
    return codes


def _substitution_cells(query: np.ndarray, block: np.ndarray) -> np.ndarray:
    """Per-pair 4×4 substitution counts over each pair's unambiguous overlap.

    `query` is one sequence's base codes and `block` is `[m, width]`
    codes for the m sequences it is compared against. Cell `[b1, b2, k]`
    counts the columns where `query` carries base `b1` and pair `k`
    carries base `b2`. A column where either side is a gap or an IUPAC
    code lands in no cell, which is tn93's `-a skip`.

    Slicing the block to the columns where `query` holds each base
    partitions the width, so all sixteen counts cost about four passes
    over the block rather than sixteen.
    """
    cells = np.zeros((4, 4, block.shape[0]), dtype=np.int64)
    for b1 in range(4):
        columns = block[:, query == b1]
        for b2 in range(4):
            cells[b1, b2] = (columns == b2).sum(axis=1)
    return cells


def _tn93_distances(cells: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Tamura-Nei 1993 distance per pair; returns `(overlap, distance)`.

    Follows `computeTN93` in tn93 v1.0.17 (`src/tn93_shared.cc`) term
    for term, including the two branches a vectorized rewrite is most
    likely to drop:

    * **Base frequencies are the pair's own, averaged over both of its
      sequences.** tn93 accumulates `nucFreq` from both the row and the
      column of every cell and divides by twice the overlap, so each
      pair is corrected against its own composition — not the cohort's,
      and not one sequence's.
    * **A pair missing one of the four bases falls back to Kimura
      2-parameter.** With a zero frequency the TN93 rate parameters
      divide by zero, so tn93 switches formula rather than producing a
      NaN. Short or low-complexity overlaps really do hit this.

    Where a logarithm's argument is not positive the pair is infinitely
    divergent under the model. tn93 returns its `TN93_MAX_DIST`
    sentinel (1000.0) and this returns `inf`; either is far above any
    surveillance threshold, so the pair gets no row — the same
    treatment a `min_overlap` refusal gets, rather than a number that
    would read as a distance.
    """
    overlap = cells.sum(axis=(0, 1))
    # Pairs with no overlap at all are refused on `min_overlap` by the
    # caller; keeping the divisor at 1 here just avoids dividing by it.
    sites = np.maximum(overlap, 1)

    base_totals = cells.sum(axis=1) + cells.sum(axis=0)
    frequency = base_totals / (2 * sites)

    purine = (cells[0, 2] + cells[2, 0]) / sites
    pyrimidine = (cells[1, 3] + cells[3, 1]) / sites
    identical = (cells[0, 0] + cells[1, 1] + cells[2, 2] + cells[3, 3]) / sites
    transversion = 1.0 - identical - purine - pyrimidine

    distance = np.full(overlap.shape, np.inf, dtype=np.float64)
    degenerate = (base_totals == 0).any(axis=0)

    full = ~degenerate
    if full.any():
        f_a, f_c, f_g, f_t = (frequency[b][full] for b in range(4))
        f_r = f_a + f_g
        f_y = f_c + f_t
        k1 = 2.0 * f_a * f_g / f_r
        k2 = 2.0 * f_c * f_t / f_y
        k3 = 2.0 * (f_r * f_y - f_a * f_g * f_y / f_r - f_c * f_t * f_r / f_y)
        tv = transversion[full]
        w_purine = 1.0 - purine[full] / k1 - 0.5 * tv / f_r
        w_pyrimidine = 1.0 - pyrimidine[full] / k2 - 0.5 * tv / f_y
        w_transversion = 1.0 - 0.5 * tv / (f_y * f_r)
        inside = (w_purine > 0.0) & (w_pyrimidine > 0.0) & (w_transversion > 0.0)
        corrected = np.full(tv.shape, np.inf)
        corrected[inside] = (
            -k1[inside] * np.log(w_purine[inside])
            - k2[inside] * np.log(w_pyrimidine[inside])
            - k3[inside] * np.log(w_transversion[inside])
        )
        distance[full] = corrected

    if degenerate.any():
        tv = transversion[degenerate]
        ti = purine[degenerate] + pyrimidine[degenerate]
        w_transition = 1.0 - 2.0 * ti - tv
        w_transversion = 1.0 - 2.0 * tv
        inside = (w_transition > 0.0) & (w_transversion > 0.0)
        corrected = np.full(tv.shape, np.inf)
        corrected[inside] = -0.5 * np.log(w_transition[inside]) - 0.25 * np.log(
            w_transversion[inside]
        )
        distance[degenerate] = corrected

    # A pair identical over its overlap sums to -0.0, every term being
    # -k x log(1.0), and would be written "-0". Adding +0.0 maps -0.0 to
    # +0.0 under IEEE 754 and leaves every other value untouched.
    distance += 0.0
    return overlap, distance


def compute_pairwise_distances_python(
    fasta_path: Path,
    output_csv: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    min_overlap: int = DEFAULT_MIN_OVERLAP,
    metric: DistanceMetric = DEFAULT_DISTANCE_METRIC,
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

    ## The two metrics

    `metric="p-distance"`, the default, divides the mismatch count by
    the overlap and thresholds that. `metric="tn93"` applies the
    Tamura-Nei 1993 correction to the same overlap: it counts purine
    (A↔G) and pyrimidine (C↔T) transitions and transversions
    separately, estimates the four base frequencies from the pair's own
    overlap, and thresholds the corrected distance. That is the
    HIV-TRACE convention, and it is what `run_tn93` computes.

    The correction is **not** a monotone transform of p across pairs,
    so one metric's pair set is not the other's reordered: two pairs
    with the same mismatch fraction get different corrected distances
    according to their substitution pattern and composition. On the
    filtered cohort the TN93-accepted and TN93-rejected pairs have
    overlapping p-ranges, so no p-threshold reproduces the TN93 pair
    set exactly (`docs/evaluations/cluster-distance-denominator.md`
    §12). At 0.015 the two cuts differ on 197 of 35,852 pairs and 3 of
    3,397 accessions, one of them joining or not joining a 58-member
    cluster.

    The default stays `p-distance` because the cut is a Paper 1 claim
    awaiting Kelyn (§12.4); `tn93` is verified against the pinned
    reference image pair for pair on the real cohort (§13) and is ready
    for the day that decision is made.

    What neither metric reproduces is tn93's `-a resolve` ambiguity
    handling: tn93 resolves an ambiguous site against its partner,
    where this drops the site. That is the conservative direction — it
    shrinks the denominator rather than inventing agreement, and it is
    tn93's `-a skip`. Pass `-a skip` for a like-for-like A/B.

    Emits a CSV with the same header (`ID1,ID2,Distance`) as tn93 so
    `read_distances` consumes it identically. The uncorrected path
    writes six decimal places, as it always has; the corrected path
    writes six significant figures, which is tn93's own `printf("%g")`,
    because the correction's spread across pairs at a fixed p-distance
    lives in the fifth decimal and the band that decides a link is
    narrower still (§12.2).
    """
    fasta_path = fasta_path.resolve()
    output_csv = output_csv.resolve()
    if not fasta_path.exists():
        raise FileNotFoundError(f"aligned FASTA missing: {fasta_path}")
    if min_overlap < 1:
        raise ValueError(f"min_overlap must be at least 1 site, got {min_overlap}")
    if metric not in DISTANCE_METRICS:
        raise ValueError(
            f"unknown distance metric {metric!r}; expected one of {DISTANCE_METRICS}"
        )
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    accessions, matrix = _read_alignment(fasta_path)
    n = len(accessions)
    logger.info("computing in-process %s distances for %d sequences", metric, n)

    if n < 2:
        output_csv.write_text("ID1,ID2,Distance\n")
        return output_csv

    codes = _base_codes(matrix)
    unambiguous = codes != _NOT_A_BASE
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
    precision = ".6g" if metric == "tn93" else ".6f"
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

            if metric == "tn93":
                overlap, distances = _tn93_distances(
                    _substitution_cells(codes[i], codes[i + 1 :][candidates])
                )
            else:
                shared = unambiguous[i + 1 :][candidates] & unambiguous[i]
                overlap = shared.sum(axis=1)
                diffs = ((matrix[i + 1 :][candidates] != matrix[i]) & shared).sum(axis=1)
                distances = diffs / np.maximum(overlap, 1)

            keep = overlap >= min_overlap
            n_refused += int((~keep).sum())

            for k in np.flatnonzero(keep & (distances <= threshold)):
                j = i + 1 + int(candidates[k])
                fh.write(
                    f"{accessions[i]},{accessions[j]},{float(distances[k]):{precision}}\n"
                )
                n_written += 1

    logger.info(
        "in-process %s distances → %s (%d pairs at or below threshold=%.4f; "
        "%d pairs refused for overlap below %d sites)",
        metric,
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
