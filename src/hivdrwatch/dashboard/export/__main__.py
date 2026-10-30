"""Command-line entry point: `python -m hivdrwatch.dashboard.export`.

Reads the Sierra results parquet and its DVC pointer, builds the summary and
writes `static/summary.json`. This module and the CLI entry points are the only
places in the dashboard package allowed to `print`.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path

import polars as pl

from hivdrwatch.dashboard.export.drm import load_sdrm_panels
from hivdrwatch.dashboard.export.summary import build_summary, input_dvc_md5, write_summary

logger = logging.getLogger(__name__)

DEFAULT_INPUT = Path("data/interim/dr_results.parquet")
DEFAULT_REFERENCE = Path("data/reference/drm_positions.json")
DEFAULT_OUT = Path(__file__).resolve().parents[1] / "static" / "summary.json"

UNKNOWN_GIT_SHA = "unknown"


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="python -m hivdrwatch.dashboard.export",
        description="Write the resistance dashboard's summary.json from Sierra results.",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--pointer",
        type=Path,
        default=None,
        help="DVC pointer for --input (default: the input path plus .dvc)",
    )
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    # Required: the parquet holds only accessions that returned a result, so the ingested
    # total cannot be derived from it, and a default would go stale when it is regenerated.
    parser.add_argument("--ingested-count", type=int, required=True)
    return parser.parse_args(argv)


def pointer_path(args: argparse.Namespace) -> Path:
    pointer: Path | None = args.pointer
    if pointer is not None:
        return pointer
    input_path: Path = args.input
    return input_path.with_name(input_path.name + ".dvc")


def resolve_git_sha() -> str:
    """HEAD's SHA, or `unknown` when git is absent or this is not a repository."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
        )
    except OSError:
        logger.warning("git is unavailable; recording git_sha as %s", UNKNOWN_GIT_SHA)
        return UNKNOWN_GIT_SHA
    if completed.returncode != 0:
        logger.warning("git rev-parse failed; recording git_sha as %s", UNKNOWN_GIT_SHA)
        return UNKNOWN_GIT_SHA
    return completed.stdout.strip()


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)

    panels, sdrm_version = load_sdrm_panels(args.reference)
    summary = build_summary(
        pl.read_parquet(args.input),
        panels=panels,
        sdrm_version=sdrm_version,
        input_md5=input_dvc_md5(pointer_path(args)),
        git_sha=resolve_git_sha(),
        n_accessions_ingested=args.ingested_count,
    )
    write_summary(summary, args.out)

    meta = summary["meta"]
    print(
        f"wrote {args.out}: {meta['n_calls']} calls, "
        f"{meta['n_accessions_scored']} of {meta['n_accessions_ingested']} accessions scored"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
