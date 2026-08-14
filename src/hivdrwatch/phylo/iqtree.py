"""IQ-TREE 2 wrapper — ML tree inference via the pinned Docker image.

Runs IQ-TREE 2 in the `hivdrwatch/phylo:latest` container with
GTR+F+I+G4 substitution model and a fixed random seed. Output artifacts
go to the caller-specified directory:

    <out_dir>/<prefix>.treefile   ← Newick (the ML tree)
    <out_dir>/<prefix>.log        ← run log
    <out_dir>/<prefix>.iqtree     ← consensus tree + stats

We deliberately do NOT run ultrafast bootstrap in the default path —
it doubles wall time and Paper 1 Figure 4 doesn't need branch support
values. Callers who want them pass `bootstrap=1000`.

Determinism: `--seed` is passed explicitly. Two runs with the same
input FASTA and same seed produce bit-identical `.treefile` output.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path

from ._streaming import stream_subprocess
from .schema import SubstitutionModel

logger = logging.getLogger(__name__)

DEFAULT_IMAGE = "hivdrwatch/phylo:latest"
DEFAULT_MODEL = SubstitutionModel.GTR_F_I_G4
DEFAULT_SEED = 42
DEFAULT_MEMORY_LIMIT = "8g"

# `-fast` skips the slow iterative-perturbation loop of IQ-TREE 2's
# default ML search and uses one round of NNI + one round of SPR
# instead. Order-of-magnitude speedup with alignment topologies that
# match the full search 95%+ of the time on HIV pol. This is what
# HIV-TRACE and Nextstrain use for interim builds. For final paper
# figures, run with `fast=False`.
DEFAULT_FAST = True


@dataclass(frozen=True)
class IqtreeResult:
    """What a successful IQ-TREE run produces."""

    treefile: Path
    log_file: Path
    prefix_dir: Path
    wall_seconds: float
    log_likelihood: float | None
    tool_version: str


def run_iqtree(
    input_fasta: Path,
    out_dir: Path,
    prefix: str,
    *,
    model: SubstitutionModel = DEFAULT_MODEL,
    seed: int = DEFAULT_SEED,
    thread_count: int | None = None,
    bootstrap: int = 0,
    fast: bool = DEFAULT_FAST,
    image: str = DEFAULT_IMAGE,
    memory_limit: str = DEFAULT_MEMORY_LIMIT,
    stream_logger: logging.Logger | None = None,
) -> IqtreeResult:
    """Run IQ-TREE 2 in Docker; return the artifact record.

    `thread_count=None` → IQ-TREE picks (`-T AUTO`). On big alignments
    this beats guessing; we let the tool do its thing.

    `fast=True` (default) enables IQ-TREE 2's `-fast` mode — order-of-
    magnitude faster ML search with minimal topology quality loss. Turn
    off for final paper-figure runs (`fast=False`) if you want the full
    iterative perturbation loop.
    """
    input_fasta = input_fasta.resolve()
    out_dir = out_dir.resolve()
    if not input_fasta.exists():
        raise FileNotFoundError(input_fasta)
    out_dir.mkdir(parents=True, exist_ok=True)

    threads = "AUTO" if thread_count is None else str(thread_count)

    args = [
        "iqtree2",
        "-s",
        f"/in/{input_fasta.name}",
        "-m",
        model.value,
        "--seed",
        str(seed),
        "-T",
        threads,
        "--prefix",
        f"/out/{prefix}",
        # `-redo` makes re-runs idempotent. Without it, IQ-TREE detects
        # a checkpoint from a prior successful run and refuses with
        # "previous run successfully finished". Prefect task caching
        # handles skip-when-unchanged at the outer layer.
        "-redo",
        # Deliberately NOT --quiet: we stream per-line to the logger
        # via stream_subprocess() so operators see progress in real
        # time instead of a black-box multi-hour wait.
    ]
    if fast:
        args.append("-fast")
    if bootstrap > 0:
        args.extend(["-B", str(bootstrap)])

    docker_cmd = [
        "docker",
        "run",
        "--rm",
        "--memory",
        memory_limit,
        "-v",
        f"{input_fasta.parent}:/in:ro",
        "-v",
        f"{out_dir}:/out",
        image,
        *args,
    ]
    logger.info("running iqtree2: %s", " ".join(docker_cmd))

    t0 = time.monotonic()
    returncode, stdout_text, stderr_text = stream_subprocess(
        docker_cmd,
        prefix=f"iqtree/{prefix}",
        logger=stream_logger or logger,
    )
    elapsed = time.monotonic() - t0

    if returncode != 0:
        # Stdout + stderr already streamed to the logger line-by-line;
        # include the tail in the exception too so it's visible in
        # Prefect task-failure summaries.
        raise RuntimeError(
            f"iqtree2 failed (exit {returncode}) in {elapsed:.1f}s\n"
            f"--- stdout tail ---\n{stdout_text[-5000:]}\n"
            f"--- stderr tail ---\n{stderr_text[-5000:]}"
        )

    treefile = out_dir / f"{prefix}.treefile"
    log_file = out_dir / f"{prefix}.log"
    if not treefile.exists():
        raise RuntimeError(f"iqtree2 exit 0 but no {treefile.name} produced; check {log_file}")

    logger.info(
        "iqtree2 complete in %.1fs → %s (%d bytes)",
        elapsed,
        treefile.name,
        treefile.stat().st_size,
    )

    ll = _extract_log_likelihood(log_file)
    version = _extract_version(log_file)

    return IqtreeResult(
        treefile=treefile,
        log_file=log_file,
        prefix_dir=out_dir,
        wall_seconds=elapsed,
        log_likelihood=ll,
        tool_version=version,
    )


# ---- pure helpers, unit-testable without Docker ----


_LL_RE = re.compile(r"BEST SCORE FOUND :\s*(-?[\d.]+)")
_VERSION_RE = re.compile(r"IQ-TREE (?:multicore )?version\s+([\d.\w-]+)")


def _extract_log_likelihood(log_file: Path) -> float | None:
    """Parse the final log-likelihood out of an IQ-TREE run log."""
    if not log_file.exists():
        return None
    try:
        text = log_file.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = _LL_RE.search(text)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _extract_version(log_file: Path) -> str:
    """Best-effort extraction of the IQ-TREE version string from the log."""
    if not log_file.exists():
        return "iqtree2"
    try:
        text = log_file.read_text(encoding="utf-8", errors="replace")[:2000]
    except OSError:
        return "iqtree2"
    match = _VERSION_RE.search(text)
    return f"iqtree2 {match.group(1)}" if match else "iqtree2"
