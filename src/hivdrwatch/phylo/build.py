"""Per-subtype phylogeny orchestrator.

Consumes the silver-layer curated parquet (`data/interim/curated/aligned.parquet`)
and produces, per subtype:

    data/interim/trees/<subtype>/
        <subtype>.fasta              — alignment subset (input to IQ-TREE)
        <subtype>.treefile           — ML tree (Newick)
        <subtype>.log                — IQ-TREE run log
        <subtype>.iqtree             — consensus tree + stats
        dates.csv                    — tip-date CSV for TreeTime
        timed/timetree.nexus         — time-scaled tree
        timed/timetree.log           — TreeTime run log

And a manifest parquet at `data/gold/trees/manifest.parquet` with one
`TreeMetadata` row per (subtype, kind) pair, so downstream stages can
look up "which tree, what vintage, built how" without re-parsing.

Subtype selection:
  * `MIN_TAXA_FOR_TREE = 20` — anything smaller is skipped (a tree of
    5 sequences is a caterpillar, not phylogenetics)
  * Only sequences that passed QC and have a non-null aligned_pol are
    included in the tree.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Literal

import polars as pl

from .iqtree import IqtreeResult, run_iqtree
from .schema import (
    PHYLO_SCHEMA_VERSION,
    SubstitutionModel,
    TreeKind,
    TreeMetadata,
)
from .treetime import TreetimeResult, build_tip_dates_csv, run_treetime

logger = logging.getLogger(__name__)

PartitionDecision = Literal[
    "built",
    "excluded_label",
    "below_min_taxa",
    "below_min_coverage",
]


PARTITION_REPORT_SCHEMA_VERSION = "0.1.0"


@dataclass(frozen=True)
class PartitionRow:
    """One row of the tree partition report — audit trail for a subtype label.

    Every candidate subtype label seen in curated data produces exactly
    one row, regardless of whether a tree was ultimately built. Used to
    reconstruct *why* the manifest contains the subtypes it does without
    reading IQ-TREE logs.
    """

    subtype: str
    n_available: int
    n_after_qc: int
    n_after_coverage: int
    decision: PartitionDecision
    decision_reason: str
    min_taxa_used: int
    min_coverage_used: float
    generated_at: str  # ISO-8601 UTC timestamp
    schema_version: str = PARTITION_REPORT_SCHEMA_VERSION

MIN_TAXA_FOR_TREE = 20

# Sequences with less than this fraction of non-gap positions in
# `aligned_pol` are dropped before tree inference. Most real GenBank
# HIV pol entries only cover a subset of pol (RT-only for genotyping,
# integrase-only for INSTI surveillance) — sending them to IQ-TREE
# alongside full-length sequences causes "problematic sequences" errors
# because they share no informative columns.
#
# 0.5 (50% non-gap) is the scientifically-defensible default: it drops
# single-gene submissions (RT-only ~44%, IN-only ~29%) and keeps
# sequences that span ≥ half of pol. Lower this if you know your cohort
# and want partial-coverage in; use per-gene trees for the cleaner
# alternative that recovers the single-gene entries (v0.2).
MIN_ALIGNED_COVERAGE = 0.5

# Subtype labels that shouldn't get their own tree.
#
#   * "unknown" — the subtype-caller couldn't classify. Those sequences
#     might be any subtype, so building a "subtype=unknown" tree mixes
#     A1 / C / D / etc. together and produces a phylogenetically
#     meaningless topology.
#   * "recombinant" — a heterogeneous label that lumps distinct CRFs
#     (CRF_01_AE, CRF_02_AG, CRF_35_AD, …) into one bucket. Those CRFs
#     do not share a common ancestor, so a single "recombinant" tree is
#     phylogenetically incoherent in the same way "unknown" is. If you
#     want a recombinant-focused tree, build per-CRF trees instead —
#     that's the scientifically defensible v0.3 direction.
#
# Both exclusions are recorded per-run in
# `data/gold/trees/partition_report.parquet` with `decision =
# "excluded_label"` so the reasoning stays visible to an auditor
# without requiring them to read the source.
EXCLUDED_SUBTYPES: frozenset[str] = frozenset({"unknown", "recombinant"})


def partition_by_subtype(
    curated: pl.DataFrame,
    *,
    min_taxa: int = MIN_TAXA_FOR_TREE,
    require_qc_pass: bool = True,
    excluded: frozenset[str] = EXCLUDED_SUBTYPES,
) -> dict[str, pl.DataFrame]:
    """Group curated rows by subtype; return only groups meeting `min_taxa`.

    Rows with null/empty subtype are dropped from consideration. Subtypes
    in `excluded` (default: `{"unknown"}`) are also dropped because
    mixing unclassified sequences into a "subtype" tree is
    phylogenetically incoherent. QC-failed rows are dropped by default;
    caller can override for exploration.
    """
    if "subtype" not in curated.columns:
        raise ValueError("curated frame missing 'subtype' column")

    df = curated.filter(
        pl.col("subtype").is_not_null()
        & (pl.col("subtype").str.strip_chars().str.len_chars() > 0)
    )
    if require_qc_pass and "passes_qc" in df.columns:
        df = df.filter(pl.col("passes_qc"))

    out: dict[str, pl.DataFrame] = {}
    for subtype in df["subtype"].unique().to_list():
        if subtype in excluded:
            logger.info("skipping excluded subtype: %s", subtype)
            continue
        subset = df.filter(pl.col("subtype") == subtype)
        if subset.height < min_taxa:
            logger.info(
                "skipping subtype %s: %d < %d min_taxa",
                subtype,
                subset.height,
                min_taxa,
            )
            continue
        out[subtype] = subset

    logger.info(
        "partition: %d subtypes eligible (>= %d taxa; excluded=%s): %s",
        len(out),
        min_taxa,
        sorted(excluded),
        sorted(out),
    )
    return out


def _non_gap_fraction(seq: str | None) -> float:
    """Fraction of positions that aren't gap characters (`-`)."""
    if not seq:
        return 0.0
    n = len(seq)
    if n == 0:
        return 0.0
    return sum(1 for c in seq if c != "-") / n


def filter_by_coverage(
    subset: pl.DataFrame,
    *,
    seq_col: str = "aligned_pol",
    min_coverage: float = MIN_ALIGNED_COVERAGE,
) -> pl.DataFrame:
    """Drop rows whose `seq_col` has less than `min_coverage` non-gap content.

    Logs the coverage distribution before + after so operators can see
    exactly how many sequences got dropped and by roughly how much.
    """
    if seq_col not in subset.columns:
        return subset

    coverages = [_non_gap_fraction(s) for s in subset[seq_col].to_list()]
    before = subset.height
    if before == 0:
        return subset

    # Compute a small histogram for the log — every 20 percentiles.
    import statistics

    q = statistics.quantiles(coverages, n=5) if len(coverages) >= 5 else []
    q_repr = ", ".join(f"{v:.2f}" for v in q) if q else "n/a"

    keep_mask = [c >= min_coverage for c in coverages]
    n_kept = sum(keep_mask)
    logger.info(
        "coverage filter (min=%.0f%%): %d/%d kept "
        "(before-quintiles: %s; dropped %d below threshold)",
        min_coverage * 100,
        n_kept,
        before,
        q_repr,
        before - n_kept,
    )

    return subset.filter(pl.Series("_keep", keep_mask))


def write_subtype_alignment(
    subset: pl.DataFrame,
    out_fasta: Path,
    *,
    seq_col: str = "aligned_pol",
    min_coverage: float | None = MIN_ALIGNED_COVERAGE,
) -> Path:
    """Write the aligned sequences for one subtype to FASTA.

    Uses `aligned_pol` by default (the full-length HXB2-aligned form
    from the curate stage). Falls back to unaligned `sequence` if
    aligned_pol isn't present.

    `min_coverage` (0..1) drops sequences whose non-gap fraction in
    `seq_col` is below the threshold — pass `None` to disable. See
    module constant `MIN_ALIGNED_COVERAGE` for the default rationale.
    """
    if seq_col not in subset.columns:
        seq_col = "sequence" if "sequence" in subset.columns else None
        if seq_col is None:
            raise ValueError("neither 'aligned_pol' nor 'sequence' in subset")

    if min_coverage is not None and seq_col == "aligned_pol":
        subset = filter_by_coverage(subset, seq_col=seq_col, min_coverage=min_coverage)

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as fh:
        for row in subset.iter_rows(named=True):
            seq = row.get(seq_col)
            if not seq:
                continue
            fh.write(f">{row['accession']}\n{seq}\n")
    logger.info(
        "wrote %d-taxa alignment for subtype (%s) → %s",
        subset.height,
        seq_col,
        out_fasta.name,
    )
    return out_fasta


def build_subtype_tree(
    subtype: str,
    subset: pl.DataFrame,
    trees_root: Path,
    *,
    seed: int = 42,
    bootstrap: int = 0,
    fast: bool = True,
    skip_treetime: bool = False,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
    min_taxa: int = MIN_TAXA_FOR_TREE,
    stream_logger: logging.Logger | None = None,
) -> tuple[list[TreeMetadata], IqtreeResult, TreetimeResult | None]:
    """End-to-end: alignment FASTA → IQ-TREE → TreeTime → manifest rows."""
    out_dir = trees_root / subtype
    out_dir.mkdir(parents=True, exist_ok=True)

    # Coverage filter first, so we can re-check min_taxa on the actual
    # tree-input set (not the pre-filter subtype-membership count).
    filtered = filter_by_coverage(subset, min_coverage=min_coverage)
    if filtered.height < min_taxa:
        raise ValueError(
            f"subtype {subtype}: only {filtered.height} sequences meet "
            f"min_coverage={min_coverage:.0%}; below min_taxa={min_taxa}. "
            f"Skipping. Either lower --min-coverage, or switch to per-gene "
            f"trees where more sequences will pass."
        )

    fasta_path = out_dir / f"{subtype}.fasta"
    write_subtype_alignment(filtered, fasta_path, min_coverage=None)  # already filtered

    ml_result = run_iqtree(
        fasta_path,
        out_dir,
        prefix=subtype,
        seed=seed,
        bootstrap=bootstrap,
        fast=fast,
        stream_logger=stream_logger,
    )

    ml_meta = TreeMetadata(
        subtype=subtype,
        kind=TreeKind.ML,
        n_taxa=filtered.height,
        newick_path=str(ml_result.treefile),
        log_path=str(ml_result.log_file),
        substitution_model=SubstitutionModel.GTR_F_I_G4,
        log_likelihood=ml_result.log_likelihood,
        bootstrap_replicates=bootstrap if bootstrap > 0 else None,
        random_seed=seed,
        tool_version=ml_result.tool_version,
        wall_seconds=ml_result.wall_seconds,
        built_at=date.today(),
    )

    if skip_treetime:
        return ([ml_meta], ml_result, None)

    dates_csv = out_dir / "dates.csv"
    build_tip_dates_csv(filtered, filtered["accession"].to_list(), dates_csv)

    tt_result = run_treetime(
        treefile=ml_result.treefile,
        aln_fasta=fasta_path,
        dates_csv=dates_csv,
        out_dir=out_dir / "timed",
        stream_logger=stream_logger,
    )

    # Compute observed time span from tip dates for the manifest.
    years = filtered["collection_year"].drop_nulls().to_list()
    time_span = (max(years) - min(years)) if years else 0.0

    tt_meta = TreeMetadata(
        subtype=subtype,
        kind=TreeKind.TIME_SCALED,
        n_taxa=filtered.height,
        newick_path=str(tt_result.timetree_nexus),
        log_path=str(tt_result.log_file),
        dates_csv_path=str(dates_csv),
        clock_rate=tt_result.clock_rate,
        clock_rate_stddev=tt_result.clock_rate_stddev,
        time_span_years=float(time_span),
        random_seed=seed,
        tool_version=tt_result.tool_version,
        wall_seconds=tt_result.wall_seconds,
        built_at=date.today(),
    )

    return ([ml_meta, tt_meta], ml_result, tt_result)


def build_partition_report(
    curated: pl.DataFrame,
    *,
    min_taxa: int = MIN_TAXA_FOR_TREE,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
    excluded: frozenset[str] = EXCLUDED_SUBTYPES,
) -> list[PartitionRow]:
    """One row per candidate subtype label — full accounting for the audit trail.

    Deterministic, does not invoke IQ-TREE. Safe to regenerate for old
    runs from a frozen curated parquet. Rows sorted by subtype so
    downstream diffs are stable.

    A row is emitted for every non-null, non-empty subtype label present
    in `curated`, so an auditor can see the fate of *every* label —
    including those excluded by policy and those below thresholds.
    """
    if "subtype" not in curated.columns:
        raise ValueError("curated frame missing 'subtype' column")

    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

    non_null = curated.filter(
        pl.col("subtype").is_not_null()
        & (pl.col("subtype").str.strip_chars().str.len_chars() > 0)
    )
    labels = sorted(non_null["subtype"].unique().to_list())

    rows: list[PartitionRow] = []
    for label in labels:
        subset = non_null.filter(pl.col("subtype") == label)
        n_available = subset.height

        if "passes_qc" in subset.columns:
            qc_pass = subset.filter(pl.col("passes_qc"))
        else:
            qc_pass = subset
        n_after_qc = qc_pass.height

        if "aligned_pol" in qc_pass.columns:
            cov_pass = filter_by_coverage(
                qc_pass, seq_col="aligned_pol", min_coverage=min_coverage
            )
        else:
            cov_pass = qc_pass
        n_after_coverage = cov_pass.height

        decision: PartitionDecision
        reason: str
        if label in excluded:
            decision = "excluded_label"
            reason = (
                f"label in EXCLUDED_SUBTYPES ({sorted(excluded)}) — "
                "phylogenetically incoherent to build a single tree "
                "over this label"
            )
        elif n_after_qc < min_taxa:
            decision = "below_min_taxa"
            reason = (
                f"n_after_qc={n_after_qc} < min_taxa={min_taxa}; "
                "increase cohort or lower --min-taxa (risks caterpillar tree)"
            )
        elif n_after_coverage < min_taxa:
            decision = "below_min_coverage"
            reason = (
                f"n_after_coverage={n_after_coverage} < min_taxa={min_taxa} "
                f"at min_coverage={min_coverage:.0%}; predominantly single-gene "
                "submissions — switch to per-gene trees to recover them"
            )
        else:
            decision = "built"
            reason = (
                f"passes all gates (n_after_coverage={n_after_coverage} "
                f"≥ min_taxa={min_taxa} at min_coverage={min_coverage:.0%})"
            )

        rows.append(
            PartitionRow(
                subtype=label,
                n_available=n_available,
                n_after_qc=n_after_qc,
                n_after_coverage=n_after_coverage,
                decision=decision,
                decision_reason=reason,
                min_taxa_used=min_taxa,
                min_coverage_used=min_coverage,
                generated_at=generated_at,
            )
        )
    return rows


def write_partition_report(rows: list[PartitionRow], out_path: Path) -> Path:
    """Serialize partition rows to a parquet audit trail."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        pl.DataFrame().write_parquet(out_path)
        return out_path
    df = pl.from_dicts([r.__dict__ for r in rows])
    df.write_parquet(out_path, compression="zstd", compression_level=6)
    logger.info(
        "wrote tree partition report: %d subtype rows → %s "
        "(decisions: %s)",
        len(rows),
        out_path,
        {r.subtype: r.decision for r in rows},
    )
    return out_path


def build_all_trees(
    curated: pl.DataFrame,
    trees_root: Path,
    *,
    subtypes: Iterable[str] | None = None,
    seed: int = 42,
    bootstrap: int = 0,
    fast: bool = True,
    skip_treetime: bool = False,
    min_taxa: int = MIN_TAXA_FOR_TREE,
    min_coverage: float = MIN_ALIGNED_COVERAGE,
    stream_logger: logging.Logger | None = None,
    partition_report_path: Path | None = None,
) -> list[TreeMetadata]:
    """Build ML + time-scaled trees for every eligible subtype.

    Returns the flat list of TreeMetadata rows across all subtypes;
    caller writes them to a manifest parquet.

    If `partition_report_path` is supplied, also emits a per-subtype
    audit-trail parquet at that path — one row per candidate subtype
    label (including excluded / below-threshold ones), so reviewers can
    reconstruct why the manifest contains what it does without reading
    IQ-TREE logs. The report reflects the same `min_taxa`, `min_coverage`,
    and `EXCLUDED_SUBTYPES` used for the run.
    """
    if partition_report_path is not None:
        report = build_partition_report(
            curated, min_taxa=min_taxa, min_coverage=min_coverage
        )
        write_partition_report(report, partition_report_path)

    partitions = partition_by_subtype(curated, min_taxa=min_taxa)
    if subtypes is not None:
        keep = set(subtypes)
        partitions = {k: v for k, v in partitions.items() if k in keep}
        logger.info("subtype filter: keeping %s", sorted(partitions))

    all_meta: list[TreeMetadata] = []
    for subtype, subset in partitions.items():
        try:
            meta_rows, _ml, _tt = build_subtype_tree(
                subtype=subtype,
                subset=subset,
                trees_root=trees_root,
                seed=seed,
                bootstrap=bootstrap,
                fast=fast,
                skip_treetime=skip_treetime,
                min_coverage=min_coverage,
                min_taxa=min_taxa,
                stream_logger=stream_logger,
            )
            all_meta.extend(meta_rows)
        except ValueError as e:
            # Controlled skip (insufficient coverage / below min_taxa after
            # filter). Not a bug — a data condition. Log at INFO, don't
            # dump a full traceback that suggests something went wrong.
            logger.info("subtype %s skipped: %s", subtype, e)
        except (RuntimeError, FileNotFoundError) as e:
            # Real failure (Docker error, IQ-TREE crash, missing input).
            # These deserve the full traceback.
            logger.exception("subtype %s failed: %s", subtype, e)
    return all_meta


def write_manifest(rows: list[TreeMetadata], out_path: Path) -> Path:
    """Serialize TreeMetadata rows to a parquet manifest."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        pl.DataFrame().write_parquet(out_path)
        return out_path
    df = pl.from_dicts([r.model_dump(mode="json") for r in rows])
    df.write_parquet(out_path, compression="zstd", compression_level=6)
    logger.info(
        "wrote tree manifest: %d rows (%d subtypes) → %s",
        len(rows),
        len({r.subtype for r in rows}),
        out_path,
    )
    return out_path


__all__ = [
    "EXCLUDED_SUBTYPES",
    "MIN_ALIGNED_COVERAGE",
    "MIN_TAXA_FOR_TREE",
    "PHYLO_SCHEMA_VERSION",
    "PartitionDecision",
    "PartitionRow",
    "build_all_trees",
    "build_partition_report",
    "build_subtype_tree",
    "partition_by_subtype",
    "write_manifest",
    "write_partition_report",
    "write_subtype_alignment",
]
