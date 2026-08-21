"""Country × year outbreak-feature aggregation from cluster output.

Emits the frame consumed by `causal/dag.py::OutbreakDag` and written to
`data/gold/outbreak_features.parquet`. Columns (one row per
country × year):

  * country
  * year
  * n_sequences_deposited      — accessions in this (country, year)
  * n_new_cluster_members      — accessions joining any cluster this year
  * n_active_clusters          — clusters with any deposit this year
  * top_cluster_size           — largest cluster with any deposit this year
  * mean_growth_rate           — mean per-cluster growth rate this year
                                 (from `cluster_growth_rate`)
  * subtype_a1_fraction        — fraction of this year's deposits typed A1
  * subtype_c_fraction         — fraction typed C
  * subtype_d_fraction         — fraction typed D
  * mean_subtype_entropy       — Shannon entropy of subtype composition
  * deposit_rate               — deposits-per-year normalised by country's
                                 max-year deposit count (proxy for
                                 sequencing effort)

This module deliberately does no cluster detection or metadata joining
of its own — it consumes what `pipeline/cluster.py` has already
produced. That keeps the outbreak-feature schema stable across cluster
threshold changes and separates the "what does an outbreak look like
in the data" question from the "what is a cluster" question.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable

import polars as pl

logger = logging.getLogger(__name__)

# Subtype fractions we report by name. Everything else rolls into
# `mean_subtype_entropy` and does not get its own column, per the
# OutbreakDag schema in causal/dag.py.
REPORTED_SUBTYPES: tuple[str, ...] = ("A1", "C", "D")


def build_outbreak_features(
    *,
    sequence_metadata: pl.DataFrame,
    cluster_assignments: pl.DataFrame,
    cluster_growth: pl.DataFrame,
    cluster_subtypes: pl.DataFrame,
) -> pl.DataFrame:
    """Aggregate cluster output into a country × year outbreak-feature frame.

    All inputs must share `accession` (except `cluster_growth` and
    `cluster_subtypes`, which are keyed by `cluster_id`). Rows with
    null `country` or `collection_year` in `sequence_metadata` are
    dropped; downstream columns handle the resulting empty strata
    by producing 0 counts and null rates rather than raising.
    """
    _require_cols(
        "sequence_metadata", sequence_metadata, {"accession", "country", "collection_year"}
    )
    _require_cols("cluster_assignments", cluster_assignments, {"accession", "cluster_id"})

    seq = sequence_metadata.filter(
        pl.col("country").is_not_null() & pl.col("collection_year").is_not_null()
    ).with_columns(pl.col("collection_year").cast(pl.Int32).alias("year"))

    seq_cluster = seq.join(cluster_assignments, on="accession", how="left")

    # ---- per-(country, year) counts ------------------------------------

    deposits = seq.group_by(["country", "year"]).agg(pl.len().alias("n_sequences_deposited"))

    new_members = (
        seq_cluster.filter(pl.col("cluster_id").is_not_null())
        .group_by(["country", "year"])
        .agg(pl.len().alias("n_new_cluster_members"))
    )

    active_clusters = (
        seq_cluster.filter(pl.col("cluster_id").is_not_null())
        .group_by(["country", "year"])
        .agg(pl.col("cluster_id").n_unique().alias("n_active_clusters"))
    )

    # ---- top-cluster size per (country, year) --------------------------

    cluster_country_year_size = (
        seq_cluster.filter(pl.col("cluster_id").is_not_null())
        .group_by(["country", "year", "cluster_id"])
        .agg(pl.len().alias("cluster_size_this_year"))
    )
    top_cluster = cluster_country_year_size.group_by(["country", "year"]).agg(
        pl.col("cluster_size_this_year").max().alias("top_cluster_size")
    )

    # ---- mean per-cluster growth this year -----------------------------

    growth_col = _first_present_col(
        cluster_growth, ("growth_rate", "mean_growth_rate", "rate")
    )
    if growth_col and "cluster_id" in cluster_growth.columns:
        year_col = _first_present_col(cluster_growth, ("year", "collection_year"))
        # If cluster_growth carries a per-cluster × per-year row, join through
        # the cluster_country_year_size frame; if per-cluster only, replicate.
        growth_join_keys = ["cluster_id", year_col] if year_col else ["cluster_id"]
        _growth = cluster_growth.select([*growth_join_keys, growth_col]).rename(
            {growth_col: "_g"}
        )
        if year_col and year_col != "year":
            _growth = _growth.rename({year_col: "year"})
        joined = cluster_country_year_size.join(
            _growth, on=[c if c != year_col else "year" for c in growth_join_keys], how="left"
        )
        mean_growth = joined.group_by(["country", "year"]).agg(
            pl.col("_g").mean().alias("mean_growth_rate")
        )
    else:
        mean_growth = deposits.select(["country", "year"]).with_columns(
            pl.lit(None, dtype=pl.Float64).alias("mean_growth_rate")
        )

    # ---- subtype fractions + Shannon entropy ---------------------------

    subtype_counts = (
        seq.filter(pl.col("subtype").is_not_null())
        .group_by(["country", "year", "subtype"])
        .agg(pl.len().alias("n"))
    )
    subtype_totals = subtype_counts.group_by(["country", "year"]).agg(
        pl.col("n").sum().alias("n_typed")
    )
    subtype_wide = subtype_counts.join(subtype_totals, on=["country", "year"]).with_columns(
        (pl.col("n") / pl.col("n_typed")).alias("frac"),
    )

    fractions_df = subtype_totals.select(["country", "year"]).join(
        _pivot_subtype_fractions(subtype_wide),
        on=["country", "year"],
        how="left",
    )

    # Shannon entropy via native Polars aggregation: -Σ frac · ln(frac).
    # `when-then` guards frac == 0 so log(0) doesn't produce NaN.
    entropy_df = subtype_wide.group_by(["country", "year"]).agg(
        (
            -pl.when(pl.col("frac") > 0)
            .then(pl.col("frac") * pl.col("frac").log())
            .otherwise(0.0)
            .sum()
        ).alias("mean_subtype_entropy")
    )

    # ---- deposit rate --------------------------------------------------

    per_country_max = deposits.group_by("country").agg(
        pl.col("n_sequences_deposited").max().alias("_country_max")
    )
    deposit_rate = (
        deposits.join(per_country_max, on="country")
        .with_columns(
            (pl.col("n_sequences_deposited") / pl.col("_country_max"))
            .fill_nan(0.0)
            .alias("deposit_rate")
        )
        .select(["country", "year", "deposit_rate"])
    )

    # ---- outer-join everything on (country, year) ----------------------

    out = deposits
    for extra in (
        new_members,
        active_clusters,
        top_cluster,
        mean_growth,
        fractions_df,
        entropy_df,
        deposit_rate,
    ):
        out = out.join(extra, on=["country", "year"], how="left")

    # Fill NA count columns with 0 (a stratum with no clusters has 0
    # active clusters, not "unknown"); leave rate columns null.
    zero_fill = [
        "n_new_cluster_members",
        "n_active_clusters",
        "top_cluster_size",
    ]
    out = out.with_columns([pl.col(c).fill_null(0).cast(pl.Int64) for c in zero_fill])

    # Ensure every column named in the OutbreakDag schema exists, even
    # if the underlying data didn't have any (e.g., a subtype absent
    # from a country gets a 0.0 fraction, not a missing column).
    _expected = {
        "country",
        "year",
        "n_sequences_deposited",
        "n_new_cluster_members",
        "n_active_clusters",
        "top_cluster_size",
        "mean_growth_rate",
        *(f"subtype_{s.lower()}_fraction" for s in REPORTED_SUBTYPES),
        "mean_subtype_entropy",
        "deposit_rate",
    }
    for c in _expected - set(out.columns):
        out = out.with_columns(pl.lit(0.0, dtype=pl.Float64).alias(c))

    out = out.select(sorted(_expected)).sort(["country", "year"])
    logger.info(
        "built outbreak features: %d (country × year) rows across %d countries",
        out.height,
        out["country"].n_unique(),
    )
    return out


# ------------------------------------------------------------------
# internals
# ------------------------------------------------------------------


def _require_cols(name: str, df: pl.DataFrame, cols: set[str]) -> None:
    missing = cols - set(df.columns)
    if missing:
        raise ValueError(f"{name} missing required columns: {sorted(missing)}")


def _first_present_col(df: pl.DataFrame, candidates: Iterable[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _pivot_subtype_fractions(subtype_wide: pl.DataFrame) -> pl.DataFrame:
    """Wide-cast fractional-subtype rows to one column per REPORTED_SUBTYPES.

    Subtypes outside REPORTED_SUBTYPES contribute to entropy but not to
    named columns.
    """
    parts = [subtype_wide.select(["country", "year"]).unique()]
    for s in REPORTED_SUBTYPES:
        _slice = subtype_wide.filter(pl.col("subtype") == s).select(
            ["country", "year", pl.col("frac").alias(f"subtype_{s.lower()}_fraction")]
        )
        parts.append(_slice)
    out = parts[0]
    for p in parts[1:]:
        out = out.join(p, on=["country", "year"], how="left")
    for s in REPORTED_SUBTYPES:
        col = f"subtype_{s.lower()}_fraction"
        out = out.with_columns(pl.col(col).fill_null(0.0))
    return out


__all__ = ["REPORTED_SUBTYPES", "build_outbreak_features"]
