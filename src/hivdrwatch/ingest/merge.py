"""Join GenBank sequences with LANL curation into the canonical dataset.

Output: data/interim/east_africa_pol.parquet
Contract: LANL subtype wins when present, GenBank fills gaps.

When GenBank is unavailable (empty frame, e.g. `--skip-genbank`) we build the
canonical dataset from LANL alone — the Table export ships raw sequences and
the FASTA has aligned ones, so LANL is now self-sufficient as a source.

The Parquet includes the schema_version and ingestion date so downstream
analyses can pin against a specific vintage.
"""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path

import polars as pl

from .schema import SCHEMA_VERSION

logger = logging.getLogger(__name__)


def merge_sources(
    genbank_df: pl.DataFrame,
    lanl_df: pl.DataFrame,
    out_path: Path,
) -> pl.DataFrame:
    """LEFT join GenBank ⟵ LANL on accession, coalescing metadata.

    If genbank_df has no sequence-bearing rows, promote lanl_df to the
    canonical output directly (LANL-only pipeline).
    """
    if "accession" not in lanl_df.columns:
        raise ValueError("lanl_df missing 'accession' column")

    genbank_empty = genbank_df.height == 0 or "sequence" not in genbank_df.columns
    if genbank_empty:
        logger.info("merge_sources: no GenBank sequences → using LANL as base frame")
        joined = _promote_lanl(lanl_df)
    else:
        joined = _join_genbank_with_lanl(genbank_df, lanl_df)

    joined = joined.with_columns(
        merged_at=pl.lit(date.today()),
        merged_schema_version=pl.lit(SCHEMA_VERSION),
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    joined.write_parquet(out_path, compression="zstd", compression_level=6)
    logger.info("wrote %d records → %s", joined.height, out_path)

    _log_summary(joined)
    return joined


# ---- internals ----


def _join_genbank_with_lanl(
    genbank_df: pl.DataFrame, lanl_df: pl.DataFrame
) -> pl.DataFrame:
    if "accession" not in genbank_df.columns:
        raise ValueError("genbank_df missing 'accession' column")

    lanl_slim_cols = [
        c
        for c in [
            "accession",
            "subtype_lanl",
            "country_code",
            "collection_year",
            "is_drug_naive",
            "risk_group",
            "pubmed_id",
            "lanl_se_id",
        ]
        if c in lanl_df.columns
    ]
    lanl_slim = lanl_df.select(lanl_slim_cols)
    if "pubmed_id" in lanl_slim.columns:
        lanl_slim = lanl_slim.rename({"pubmed_id": "pubmed_id_lanl"})

    joined = genbank_df.join(lanl_slim, on="accession", how="left")

    # Coalesce: LANL curation preferred, GenBank fallback.
    coalesce_exprs = []

    if "subtype_lanl" in joined.columns:
        coalesce_exprs.append(
            pl.coalesce([pl.col("subtype_lanl"), pl.col("subtype")]).alias("subtype")
        )
        coalesce_exprs.append(
            pl.when(pl.col("subtype_lanl").is_not_null())
            .then(pl.lit("lanl_curated"))
            .otherwise(pl.col("subtype_source"))
            .alias("subtype_source")
        )

    if "country_code" in joined.columns and "country" in joined.columns:
        coalesce_exprs.append(
            pl.coalesce([pl.col("country_code"), pl.col("country")]).alias("country")
        )

    if coalesce_exprs:
        joined = joined.with_columns(coalesce_exprs)

    joined = joined.drop(
        [
            c
            for c in ["subtype_lanl", "country_code", "pubmed_id_lanl"]
            if c in joined.columns
        ]
    )
    return joined


def _promote_lanl(lanl_df: pl.DataFrame) -> pl.DataFrame:
    """Reshape a LANL-only DataFrame into the canonical column layout.

    Expects lanl_df to already carry a `sequence` column (either from the
    TSV Sequence field or hydrated from FASTA upstream). Rows without a
    sequence are dropped.
    """
    if "sequence" not in lanl_df.columns:
        raise ValueError(
            "LANL-only merge requires a 'sequence' column. Either include the "
            "Sequence field in the LANL Table export or hydrate from FASTA "
            "via LANLIngester.to_sequence_records()."
        )

    df = lanl_df.filter(pl.col("sequence").is_not_null())

    exprs = []
    # Only add these when not already present — the SequenceRecord path
    # already writes canonical columns; the raw-DataFrame path does not.
    if "source" not in df.columns:
        exprs.append(pl.lit("lanl").alias("source"))
    if "gene_region" not in df.columns:
        exprs.append(pl.lit("pol").alias("gene_region"))
    if "schema_version" not in df.columns:
        exprs.append(pl.lit(SCHEMA_VERSION).alias("schema_version"))

    if "subtype_lanl" in df.columns and "subtype" not in df.columns:
        exprs.append(pl.col("subtype_lanl").alias("subtype"))
        exprs.append(
            pl.when(pl.col("subtype_lanl").is_not_null())
            .then(pl.lit("lanl_curated"))
            .otherwise(pl.lit(None))
            .alias("subtype_source")
        )

    if "country_code" in df.columns and "country" not in df.columns:
        exprs.append(pl.col("country_code").alias("country"))

    if exprs:
        df = df.with_columns(exprs)

    drop = [c for c in ["subtype_lanl", "country_code"] if c in df.columns]
    if drop:
        df = df.drop(drop)

    return df


def _log_summary(df: pl.DataFrame) -> None:
    """Print a quick sanity check that shows up in CI logs & the paper appendix."""
    logger.info("=== merge summary ===")
    logger.info("total records: %d", df.height)
    if "country" in df.columns:
        by_country = df.group_by("country").len().sort("len", descending=True)
        logger.info("by country:\n%s", by_country)
    if "subtype" in df.columns:
        by_subtype = df.group_by("subtype").len().sort("len", descending=True)
        logger.info("by subtype:\n%s", by_subtype)
    if "collection_year" in df.columns:
        yr = df.select(
            pl.col("collection_year").min().alias("min"),
            pl.col("collection_year").max().alias("max"),
            pl.col("collection_year").is_not_null().sum().alias("with_year"),
        )
        logger.info("year coverage:\n%s", yr)
