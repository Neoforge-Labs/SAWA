"""Prefect flow: merged sequences → aligned + trimmed + QC'd + deduped parquet.

Reads `data/interim/east_africa_pol.parquet` (from the ingest flow) and
writes `data/interim/curated/aligned.parquet` — the silver layer that
feeds both the DR feature builder and the phylogenetic clustering stage.

Steps (each a Prefect task with retries + input-hashed caching):
    1. ensure the HXB2 pol reference exists (extract from LANL FASTA if not)
    2. run MAFFT --add via docker
    3. trim to PR / RT / IN
    4. run QC (ambiguity, APOBEC signature, premature stops, coverage)
    5. assign patient_group_id + is_earliest_per_patient
    6. write parquet

Run locally:
    prefect server start &
    python -m hivdrwatch.pipeline.curate \\
        --input  data/interim/east_africa_pol.parquet \\
        --lanl-dir data/raw/lanl \\
        --work-dir data/interim/curate \\
        --out data/interim/curated/aligned.parquet
"""

from __future__ import annotations

import argparse
import logging
from datetime import date
from pathlib import Path

import polars as pl
from prefect import flow, get_run_logger, task
from prefect.tasks import task_input_hash

from ..align import (
    CURATE_SCHEMA_VERSION,
    ensure_hxb2_reference,
    align_sequences,
)
from ..align.dedup import assign_patient_groups
from ..align.qc import (
    ambiguity_fraction,
    apobec_signature_ratio,
    has_premature_stop,
    passes_qc,
    qc_flags,
)
from ..align.trim import coverage_fraction, slice_all_regions
from ..resistance.cache import sequence_sha256


# ------------------------------------------------------------------
# Tasks
# ------------------------------------------------------------------


@task(cache_key_fn=task_input_hash)
def load_ingest_parquet(path: Path) -> pl.DataFrame:
    log = get_run_logger()
    df = pl.read_parquet(path)
    required = {"accession", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    log.info("loaded %d sequences from %s", df.height, path)
    return df


@task(cache_key_fn=task_input_hash)
def ensure_reference(reference_path: Path, lanl_dir: Path) -> Path:
    ref = ensure_hxb2_reference(reference_path, lanl_dir)
    return ref.source_path


@task(retries=1)
def align(
    df: pl.DataFrame,
    reference_fasta: Path,
    work_dir: Path,
    chunk_size: int = 1000,
    chunk_concurrency: int = 2,
    aligner: str = "mappy",
) -> dict[str, str]:
    log = get_run_logger()
    records = list(zip(df["accession"].to_list(), df["sequence"].to_list(), strict=True))
    log.info(
        "aligner=%s for %d input records (chunk_size=%d, chunk_concurrency=%d)",
        aligner, len(records), chunk_size, chunk_concurrency,
    )
    aligned = align_sequences(
        records,
        reference_fasta,
        work_dir,
        aligner=aligner,
        chunk_size=chunk_size,
        chunk_concurrency=chunk_concurrency,
    )
    log.info("aligned %d sequences (reference: %s)", len(aligned), reference_fasta.name)
    return aligned


@task
def curate_rows(
    df: pl.DataFrame, aligned: dict[str, str], reference_fasta: Path
) -> list[dict]:
    """Assemble the per-sequence curated row: trim + QC + provenance."""
    log = get_run_logger()

    ref_seq = _read_reference_ungapped_padded(reference_fasta, len(next(iter(aligned.values()), "")))
    today = date.today()
    rows: list[dict] = []

    for row in df.iter_rows(named=True):
        accession = row["accession"]
        aligned_seq = aligned.get(accession)
        if aligned_seq is None:
            continue  # MAFFT dropped it (empty / all-gap input)

        regions = slice_all_regions(aligned_seq)
        ambig = ambiguity_fraction(aligned_seq)
        apobec = apobec_signature_ratio(aligned_seq, ref_seq) if ref_seq else None
        coverage = {r: coverage_fraction(seq) for r, seq in regions.items()}
        stops = {r: has_premature_stop(seq) for r, seq in regions.items()}
        flags = qc_flags(
            ambiguity=ambig,
            apobec_ratio=apobec,
            region_coverage=coverage,
            region_stops=stops,
        )
        rows.append(
            {
                "accession": accession,
                "sequence_sha256": sequence_sha256(row["sequence"]),
                "aligned_pol": aligned_seq,
                "aligned_length": len(aligned_seq),
                "pr_nt": regions["PR"],
                "rt_nt": regions["RT"],
                "in_nt": regions["IN"],
                "ambiguity_fraction": ambig,
                "apobec_signature_ratio": apobec,
                "has_pr_stop": stops["PR"],
                "has_rt_stop": stops["RT"],
                "has_in_stop": stops["IN"],
                "qc_flags": [f.value for f in flags],
                "passes_qc": passes_qc(flags),
                "country": row.get("country"),
                "collection_year": row.get("collection_year"),
                "subtype": row.get("subtype"),
                "sequence_name": row.get("sequence_name"),
                "source_file": row.get("source_file"),
                "curated_at": today.isoformat(),
                "schema_version": CURATE_SCHEMA_VERSION,
            }
        )
    log.info("curated %d rows (dropped %d)", len(rows), df.height - len(rows))
    return rows


@task
def dedup(rows: list[dict]) -> list[dict]:
    return assign_patient_groups(rows)


@task
def write_curated_parquet(rows: list[dict], out_path: Path) -> Path:
    log = get_run_logger()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # `sequence_name` is provenance-only; drop before write to match schema.
    for row in rows:
        row.pop("sequence_name", None)
    df = pl.from_dicts(rows) if rows else pl.DataFrame()
    df.write_parquet(out_path, compression="zstd", compression_level=6)
    log.info("wrote %d curated rows → %s", df.height, out_path)

    if df.height:
        by_qc = df.group_by("passes_qc").len().sort("len", descending=True)
        log.info("QC breakdown:\n%s", by_qc)
        if "country" in df.columns:
            by_country = df.group_by("country").len().sort("len", descending=True)
            log.info("by country:\n%s", by_country)

    return out_path


# ------------------------------------------------------------------
# Flow
# ------------------------------------------------------------------


@flow(name="curate")
def curate_flow(
    input_parquet: Path,
    lanl_dir: Path,
    work_dir: Path,
    out_path: Path,
    reference_path: Path | None = None,
    *,
    chunk_size: int = 1000,
    chunk_concurrency: int = 2,
    aligner: str = "mappy",
) -> Path:
    reference_path = reference_path or (
        Path("data/reference") / "hxb2_pol.fasta"
    )
    df = load_ingest_parquet(input_parquet)
    ref_fasta = ensure_reference(reference_path, lanl_dir)
    aligned = align(
        df, ref_fasta, work_dir,
        chunk_size=chunk_size,
        chunk_concurrency=chunk_concurrency,
        aligner=aligner,
    )
    rows = curate_rows(df, aligned, ref_fasta)
    rows = dedup(rows)
    return write_curated_parquet(rows, out_path)


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------


def _read_reference_ungapped_padded(reference_fasta: Path, aligned_length: int) -> str:
    """Read reference and pad to aligned length if needed.

    After MAFFT --add --keeplength, the reference retains its own length.
    APOBEC scoring uses gap-free reference vs aligned query, so if lengths
    differ we bail (returning empty triggers apobec=None downstream).
    """
    from ..align.reference import _read_single_fasta

    seq = _read_single_fasta(reference_fasta)
    if aligned_length and len(seq) != aligned_length:
        return ""
    return seq


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.pipeline.curate")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/interim/east_africa_pol.parquet"),
        help="merged sequence parquet from hivdrwatch.ingest",
    )
    parser.add_argument(
        "--lanl-dir",
        type=Path,
        default=Path("data/raw/lanl"),
        help="where HXB2 reference is extracted from if not yet cached",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("data/interim/curate"),
        help="working directory for MAFFT input/output FASTA (mounted into docker)",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("data/reference/hxb2_pol.fasta"),
        help="HXB2 pol reference (auto-extracted from LANL if missing)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("data/interim/curated/aligned.parquet"),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1000,
        help="MAFFT --add chunk size — smaller = less memory, more chunks. "
             "1000 fits comfortably in 6 GB Docker; drop to 500 on tight boxes.",
    )
    parser.add_argument(
        "--chunk-concurrency",
        type=int,
        default=2,
        help="How many MAFFT chunks run in parallel. Each uses ~memory_limit "
             "memory. Default 2 assumes 16 GB Docker Desktop; bump to 4 with 32 GB.",
    )
    parser.add_argument(
        "--aligner",
        choices=["mappy", "biopython", "mafft"],
        default="mappy",
        help="'mappy' (default): in-process minimap2, ~700 seq/s. "
             "'biopython': pure Python fallback, ~5 seq/s. "
             "'mafft': Docker MAFFT --add for joint MSA (slow, rare).",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    curate_flow(
        input_parquet=args.input,
        lanl_dir=args.lanl_dir,
        work_dir=args.work_dir,
        out_path=args.out,
        reference_path=args.reference,
        chunk_size=args.chunk_size,
        chunk_concurrency=args.chunk_concurrency,
        aligner=args.aligner,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
