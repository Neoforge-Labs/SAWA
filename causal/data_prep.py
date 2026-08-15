"""Adapter: Paper 1 Parquet outputs → causal-analysis DataFrames.

`causal/` ships the DAGs, DoWhy identify+estimate+refute loop, feature
pruning, and HTML report. The one missing piece was a data adapter
converting HIV-DR-Watch's canonical outputs into the shape
`causal.estimate.estimate_dr_effect()` expects. This module is that
adapter.

## The shape estimate_dr_effect wants

    df: pl.DataFrame with columns
      - <mutation>          binary 0/1 per row (e.g. "RT_184_V")
      - label_binary        binary 0/1 outcome (resistant vs not)
      - subtype             DR-DAG covariate
      - country             DR-DAG covariate
      - collection_year     DR-DAG covariate

## What we produce

For each drug class (or drug, if `per_drug=True`), one row per
accession with:
  * mutation binary columns (from Sierra's per-drug `mutations` list,
    aggregated across all drugs in scope for that accession),
  * covariates joined from `curated/aligned.parquet`,
  * `label_binary = 1` iff the accession has ≥ intermediate resistance
    to *any* drug in scope — matches the Fig 2 / Fig 4 convention.

## Unit-of-analysis choice

We use one row per accession per scope (drug class or drug). The
alternative — patient-visit or clinic-quarter — needs clinical data
we don't have in the LANL/GenBank path. When KEMRI's DSA lands
(Paper 4), a Phase-2 module `causal/data_prep_kemri.py` will replace
this on the KEMRI branch; the DR-DAG stays identical.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# Sierra ships mutation strings like "M184V" (wild + position + variant).
# The DR-DAG uses "RT_184_V" (gene + position + variant). We ignore
# multi-variant strings like "M184MV" (mixed callouts) and insertions
# like "T69S_SS" — they aren't in the DAG's canonical mutation list.
_MUT_RE = re.compile(r"^[A-Z](\d+)([A-Z*])$")

# Sierra levels that count as "resistant" for the binary outcome.
# Matches paper1.fig02_aggregate.RESISTANT_LEVELS and
# paper1.fig04_dtg.RESISTANT_LEVELS — one definition across the paper.
RESISTANT_LEVELS: frozenset[str] = frozenset({"high", "intermediate"})


@dataclass(frozen=True)
class CausalDatasetCounts:
    """Compact "adapter built OK" summary for logs + Prefect UI."""

    n_accessions: int
    n_scopes: int  # drug classes if per_drug=False, else drugs
    n_dr_rows: int
    n_curated_rows: int
    n_mutation_features: int
    scope_kind: str  # "drug_class" or "drug"


def build_dr_causal_frames(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    *,
    mutations_of_interest: tuple[str, ...] | None = None,
    per_drug: bool = False,
) -> tuple[dict[str, pl.DataFrame], CausalDatasetCounts]:
    """Convert Paper 1 DR + curated outputs into causal-analysis frames.

    Returns `(frames, counts)` where `frames` is a dict keyed by scope
    value (drug class or drug name). Each frame is drop-in input for
    `causal.estimate.estimate_dr_effect(df, mutation="RT_184_V",
    data_scope=f"DR/drug_class={scope}")`.

    `mutations_of_interest` defaults to `DRDag().mutation_features` so
    the produced frames slot directly into the shipped DR-DAG.
    """
    # Lazy DAG import — this module is safe to import in a KEMRI branch
    # where the DAG may be swapped, without forcing DRDag construction.
    from .dag import DRDag

    if mutations_of_interest is None:
        mutations_of_interest = DRDag().mutation_features

    _require_cols(
        "dr_results",
        dr_results,
        {"accession", "gene", "mutations", "drug_class", "drug", "level"},
    )
    _require_cols(
        "curated",
        curated,
        {"accession", "subtype", "country", "collection_year"},
    )

    scope_col: str = "drug" if per_drug else "drug_class"

    # Per-accession binary matrix for the DAG's mutation features.
    mut_matrix = _extract_mutation_features(dr_results, mutations_of_interest)

    # Per (accession, scope) binary label: "any drug in scope is
    # ≥ intermediate resistant".
    label = (
        dr_results.with_columns(
            pl.col("level").is_in(list(RESISTANT_LEVELS)).alias("_r"),
        )
        .group_by(["accession", scope_col])
        .agg(pl.col("_r").any().alias("label_binary"))
        .with_columns(pl.col("label_binary").cast(pl.Int8))
    )

    # De-dup curated to one row per accession (some pipelines emit
    # duplicates when the same accession appears in multiple LANL folder
    # revisions).
    covs = curated.select(["accession", "subtype", "country", "collection_year"]).unique(
        subset=["accession"], keep="first"
    )

    frames: dict[str, pl.DataFrame] = {}
    for scope_value in sorted(
        str(v) for v in label[scope_col].unique().to_list() if v is not None
    ):
        sub_label = label.filter(pl.col(scope_col) == scope_value).select(
            ["accession", "label_binary"]
        )
        joined = (
            sub_label.join(mut_matrix, on="accession", how="inner")
            .join(covs, on="accession", how="inner")
            .filter(
                pl.col("subtype").is_not_null()
                & pl.col("country").is_not_null()
                & pl.col("collection_year").is_not_null()
            )
        )
        frames[scope_value] = joined

    counts = CausalDatasetCounts(
        n_accessions=int(mut_matrix.height),
        n_scopes=len(frames),
        n_dr_rows=int(dr_results.height),
        n_curated_rows=int(curated.height),
        n_mutation_features=len(mutations_of_interest),
        scope_kind=scope_col,
    )
    logger.info(
        "built %d causal frames (scope=%s): %d accessions × %d mutations "
        "from %d DR rows + %d curated rows",
        counts.n_scopes,
        counts.scope_kind,
        counts.n_accessions,
        counts.n_mutation_features,
        counts.n_dr_rows,
        counts.n_curated_rows,
    )
    return frames, counts


def load_dr_causal_dataset(
    dr_parquet: Path,
    curated_parquet: Path,
    *,
    per_drug: bool = False,
) -> tuple[dict[str, pl.DataFrame], CausalDatasetCounts]:
    """Convenience wrapper — read the Parquets and hand off to the builder."""
    dr = pl.read_parquet(dr_parquet)
    curated = pl.read_parquet(curated_parquet)
    return build_dr_causal_frames(dr, curated, per_drug=per_drug)


# ------------------------------------------------------------------
# internals
# ------------------------------------------------------------------


def _parse_mutation(gene: str | None, mut: str | None) -> str | None:
    """Sierra 'M184V' + gene 'RT' → DAG feature 'RT_184_V'.

    Returns None on any of:
      * missing gene or mutation string
      * insertions / deletions / mixed-variant strings
      * gene not in {PR, RT, IN}
    """
    if not gene or not mut:
        return None
    gene = gene.upper().strip()
    if gene not in {"PR", "RT", "IN"}:
        return None
    m = _MUT_RE.match(mut.strip())
    if not m:
        return None
    position, variant = m.group(1), m.group(2)
    return f"{gene}_{position}_{variant}"


def _extract_mutation_features(
    dr_results: pl.DataFrame,
    mutations_of_interest: tuple[str, ...],
) -> pl.DataFrame:
    """Per-accession binary matrix for the DAG's mutation features.

    A mutation counts as "present" for an accession if any DR call for
    that accession lists it. Sierra emits the same mutation redundantly
    across every drug it affects; we dedup by (accession, mutation).

    Every mutation in `mutations_of_interest` is guaranteed to appear
    as a column, even if no accession carries it — that keeps the
    downstream estimator interface stable across cohorts of varying
    resistance prevalence.
    """
    exploded = (
        dr_results.select(["accession", "gene", "mutations"])
        .filter(pl.col("mutations").list.len() > 0)
        .explode("mutations")
    )

    if exploded.height == 0:
        # No mutations at all — return an all-zero matrix with every
        # accession we saw in dr_results.
        accessions = dr_results.select("accession").unique(keep="first")
        for mut in mutations_of_interest:
            accessions = accessions.with_columns(pl.lit(0, dtype=pl.Int8).alias(mut))
        return accessions.select(["accession", *mutations_of_interest])

    # Parse Sierra strings into DAG format. Python loop is fine here —
    # exploded height is O(n_accessions × n_mutations_per_call), typically
    # ~10⁵ for a real East-Africa cohort.
    dag_muts: list[str | None] = [
        _parse_mutation(g, m)
        for g, m in zip(
            exploded["gene"].to_list(),
            exploded["mutations"].to_list(),
            strict=True,
        )
    ]
    exploded = exploded.with_columns(pl.Series("dag_mutation", dag_muts))
    interesting = exploded.filter(pl.col("dag_mutation").is_in(list(mutations_of_interest)))

    accessions = dr_results.select("accession").unique(keep="first")

    if interesting.height == 0:
        for mut in mutations_of_interest:
            accessions = accessions.with_columns(pl.lit(0, dtype=pl.Int8).alias(mut))
        return accessions.select(["accession", *mutations_of_interest])

    # (accession, mutation) presence → wide matrix.
    matrix = (
        interesting.select(["accession", "dag_mutation"])
        .unique()
        .with_columns(pl.lit(1, dtype=pl.Int8).alias("_present"))
        .pivot(values="_present", index="accession", on="dag_mutation")
    )

    # Backfill accessions that had mutations but none in our interest set,
    # and ensure every mutation column exists (fill missing as 0).
    matrix = accessions.join(matrix, on="accession", how="left").fill_null(0)
    for mut in mutations_of_interest:
        if mut not in matrix.columns:
            matrix = matrix.with_columns(pl.lit(0, dtype=pl.Int8).alias(mut))
        else:
            matrix = matrix.with_columns(pl.col(mut).cast(pl.Int8))

    return matrix.select(["accession", *mutations_of_interest])


def _require_cols(name: str, df: pl.DataFrame, cols: set[str]) -> None:
    missing = cols - set(df.columns)
    if missing:
        raise ValueError(f"{name} missing required columns: {sorted(missing)}")


__all__ = [
    "RESISTANT_LEVELS",
    "CausalDatasetCounts",
    "build_dr_causal_frames",
    "load_dr_causal_dataset",
]
