"""TreeTime wrapper — molecular-clock time-scaling of an IQ-TREE topology.

Runs TreeTime in the pinned Docker image (same as IQ-TREE — both live in
`hivdrwatch/phylo:latest`). Inputs:

  * ML tree Newick (from `iqtree.run_iqtree`)
  * Alignment FASTA (same input MAFFT/mappy produced for IQ-TREE)
  * Tip dates as `name,date` CSV (we build this from the curated parquet)

Outputs are written under `<out_dir>/timed/` per TreeTime's convention:

  * `timetree.nexus`       — time-scaled tree
  * `divergence_tree.nexus` — same topology, divergence units
  * `dates.tsv`             — inferred dates for internal nodes
  * `timetree.log`          — run log
  * `ancestral_sequences.fasta` — optional; we don't request it

Determinism: TreeTime's clock estimation is deterministic given the
input dates + tree, so we don't need to pass a seed.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ._streaming import stream_subprocess

logger = logging.getLogger(__name__)

DEFAULT_IMAGE = "hivdrwatch/phylo:latest"
DEFAULT_MEMORY_LIMIT = "8g"


@dataclass(frozen=True)
class TreetimeResult:
    """Artifacts produced by a successful TreeTime run."""

    timetree_nexus: Path
    log_file: Path
    out_dir: Path
    wall_seconds: float
    clock_rate: float | None
    clock_rate_stddev: float | None
    tool_version: str


def build_tip_dates_csv(
    curated_parquet: pl.DataFrame,
    accessions: list[str],
    out_csv: Path,
) -> Path:
    """Emit TreeTime's expected `name,date` CSV for the given accessions.

    Uses `collection_year` from the curated parquet as the date (as a
    year-only decimal — TreeTime accepts either ISO date or decimal
    years). Rows with missing year are dropped; if a tree has too many
    date-less tips TreeTime falls back to strict clock which is fine.
    """
    required = {"accession", "collection_year"}
    missing = required - set(curated_parquet.columns)
    if missing:
        raise ValueError(f"curated parquet missing columns: {missing}")

    subset = curated_parquet.filter(
        pl.col("accession").is_in(accessions) & pl.col("collection_year").is_not_null()
    ).select(["accession", "collection_year"])

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8") as fh:
        fh.write("name,date\n")
        for row in subset.iter_rows(named=True):
            year = int(row["collection_year"])
            # TreeTime accepts decimal years — midyear (0.5) as the
            # coarse-grained anchor when we only know the year.
            fh.write(f"{row['accession']},{year + 0.5:.2f}\n")

    logger.info(
        "wrote %d tip dates → %s (from %d requested accessions)",
        subset.height,
        out_csv,
        len(accessions),
    )
    return out_csv


def run_treetime(
    treefile: Path,
    aln_fasta: Path,
    dates_csv: Path,
    out_dir: Path,
    *,
    image: str = DEFAULT_IMAGE,
    memory_limit: str = DEFAULT_MEMORY_LIMIT,
    stream_logger: logging.Logger | None = None,
) -> TreetimeResult:
    """Time-scale a Newick tree with TreeTime, given tip dates.

    All three input paths must share a parent directory so we can mount
    it once into `/work`. `out_dir` should be a subdirectory of that
    parent — TreeTime writes there.
    """
    treefile = treefile.resolve()
    aln_fasta = aln_fasta.resolve()
    dates_csv = dates_csv.resolve()
    out_dir = out_dir.resolve()

    for p in (treefile, aln_fasta, dates_csv):
        if not p.exists():
            raise FileNotFoundError(p)

    # All three inputs must be under a shared parent so /work maps cleanly.
    parents = {treefile.parent, aln_fasta.parent, dates_csv.parent}
    if len(parents) != 1:
        raise ValueError(f"treetime inputs must share a parent dir; got {parents}")
    work_dir = parents.pop()
    if not str(out_dir).startswith(str(work_dir)):
        raise ValueError(f"out_dir {out_dir} must be under work_dir {work_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_rel = out_dir.relative_to(work_dir)

    args = [
        "treetime",
        "--tree",
        f"/work/{treefile.name}",
        "--aln",
        f"/work/{aln_fasta.name}",
        "--dates",
        f"/work/{dates_csv.name}",
        "--outdir",
        f"/work/{out_rel}",
        # No reconstructed ancestral sequences — we don't need them and
        # they add non-trivial runtime.
        "--keep-root",
    ]
    docker_cmd = [
        "docker",
        "run",
        "--rm",
        "--memory",
        memory_limit,
        "-v",
        f"{work_dir}:/work",
        image,
        *args,
    ]
    logger.info("running treetime: %s", " ".join(docker_cmd))

    # Prefix the subtype dir name so per-subtype runs are distinguishable
    # in the log stream (e.g. `treetime/A1`).
    subtype_hint = out_dir.parent.name if out_dir.name == "timed" else out_dir.name

    t0 = time.monotonic()
    returncode, stdout_text, stderr_text = stream_subprocess(
        docker_cmd,
        prefix=f"treetime/{subtype_hint}",
        logger=stream_logger or logger,
    )
    elapsed = time.monotonic() - t0

    if returncode != 0:
        raise RuntimeError(
            f"treetime failed (exit {returncode}) in {elapsed:.1f}s\n"
            f"--- stdout tail ---\n{stdout_text[-3000:]}\n"
            f"--- stderr tail ---\n{stderr_text[-3000:]}"
        )

    timetree = out_dir / "timetree.nexus"
    log_file = out_dir / "timetree.log"
    if not timetree.exists():
        # Older TreeTime versions use different filename conventions
        # (e.g. `timetree.nex`); look for anything Nexus-y as a fallback.
        candidates = list(out_dir.glob("timetree.nex*"))
        if candidates:
            timetree = candidates[0]
        else:
            raise RuntimeError(
                f"treetime exit 0 but no timetree.nexus in {out_dir}; "
                f"stderr tail:\n{stderr_text[-500:]}"
            )

    clock_rate, clock_rate_stddev = _extract_clock_rate(stderr_text)
    version = _extract_version(stderr_text)

    logger.info(
        "treetime complete in %.1fs (clock=%s ± %s subs/site/year) → %s",
        elapsed,
        f"{clock_rate:.4g}" if clock_rate else "?",
        f"{clock_rate_stddev:.2g}" if clock_rate_stddev else "?",
        timetree.name,
    )

    return TreetimeResult(
        timetree_nexus=timetree,
        log_file=log_file if log_file.exists() else out_dir,
        out_dir=out_dir,
        wall_seconds=elapsed,
        clock_rate=clock_rate,
        clock_rate_stddev=clock_rate_stddev,
        tool_version=version,
    )


# ---- pure helpers ----


_CLOCK_RE = re.compile(
    r"clock\s+rate\s*:\s*(-?[\d.eE+-]+)\s*\+/-\s*(-?[\d.eE+-]+)",
    re.IGNORECASE,
)
_VERSION_RE = re.compile(r"TreeTime\s+([\d.]+)", re.IGNORECASE)


def _extract_clock_rate(stderr: str) -> tuple[float | None, float | None]:
    """Parse `clock rate: X +/- Y` from TreeTime stderr."""
    match = _CLOCK_RE.search(stderr)
    if not match:
        return (None, None)
    try:
        return (float(match.group(1)), float(match.group(2)))
    except ValueError:
        return (None, None)


def _extract_version(stderr: str) -> str:
    match = _VERSION_RE.search(stderr)
    return f"treetime {match.group(1)}" if match else "treetime"
