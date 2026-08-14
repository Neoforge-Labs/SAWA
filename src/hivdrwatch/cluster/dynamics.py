"""Per-cluster dynamics: growth over time, geography, subtype mix.

We deliberately avoid dated-tree phylogenetics here (that's the
BEAST/TreeTime path — expensive and a full engineering effort of its
own). Instead we use *deposit-year growth* as a cheap proxy for cluster
activity: how many new members entered the cluster in a given year.

That's a biased proxy — deposit year ≠ infection year, and sequencing
rates dominate the signal. We correct for that downstream in the outbreak
features by explicitly including per-country-year deposit rate as a
covariate. That's the biggest methodological trap in this problem;
callers should keep it front-of-mind.
"""

from __future__ import annotations

import logging
from collections import Counter

import polars as pl

logger = logging.getLogger(__name__)


def cluster_time_series(
    cluster_assignments: pl.DataFrame,
    sequence_metadata: pl.DataFrame,
) -> pl.DataFrame:
    """Return per-(cluster, country, year) new-member counts.

    Joins `cluster_assignments` (cols: accession, cluster_id, cluster_size)
    with `sequence_metadata` (cols: accession, country, collection_year).
    Rows without a collection year are dropped; those without a country
    are grouped under 'unknown'.
    """
    required = {"accession", "cluster_id", "cluster_size"}
    if not required.issubset(cluster_assignments.columns):
        raise ValueError(f"cluster_assignments missing {required}")
    if not {"accession", "country", "collection_year"}.issubset(
        sequence_metadata.columns
    ):
        raise ValueError(
            "sequence_metadata must have accession, country, collection_year"
        )

    joined = cluster_assignments.join(
        sequence_metadata.select("accession", "country", "collection_year"),
        on="accession",
        how="inner",
    ).filter(pl.col("collection_year").is_not_null())

    joined = joined.with_columns(
        pl.col("country").fill_null("unknown").alias("country")
    )

    return (
        joined.group_by("cluster_id", "country", "collection_year")
        .agg(
            n_new_members=pl.col("accession").n_unique(),
            cluster_size=pl.col("cluster_size").first(),
        )
        .sort("cluster_id", "country", "collection_year")
    )


def cluster_growth_rate(
    time_series: pl.DataFrame, *, window_years: int = 3
) -> pl.DataFrame:
    """Compute the rolling growth rate per cluster per country.

    growth_rate = (n_new in last `window_years`) / max(1, n_new in prior window)

    A cluster with 10 new members in 2021-2023 and 2 in 2018-2020 has a
    growth rate of 5.0. We cap the denominator at 1 to keep new
    fast-growing clusters (no prior activity) from dividing by zero.
    """
    if window_years < 1:
        raise ValueError("window_years must be >= 1")

    # For every (cluster, country), pivot to year → n_new, compute sums over windows.
    rows: list[dict] = []
    for (cluster_id, country), group in time_series.group_by(
        ["cluster_id", "country"], maintain_order=True
    ):
        year_to_n = {
            int(y): int(n)
            for y, n in zip(group["collection_year"], group["n_new_members"], strict=True)
        }
        min_year, max_year = min(year_to_n), max(year_to_n)
        for anchor in range(min_year, max_year + 1):
            recent_lo = anchor - window_years + 1
            prior_hi = recent_lo - 1
            prior_lo = prior_hi - window_years + 1
            recent = sum(
                year_to_n.get(y, 0) for y in range(recent_lo, anchor + 1)
            )
            prior = sum(year_to_n.get(y, 0) for y in range(prior_lo, prior_hi + 1))
            rows.append(
                {
                    "cluster_id": cluster_id,
                    "country": country,
                    "year": anchor,
                    "recent_n": recent,
                    "prior_n": prior,
                    "growth_rate": recent / max(1, prior),
                }
            )

    return pl.from_dicts(rows) if rows else pl.DataFrame()


def cluster_subtype_mix(
    cluster_assignments: pl.DataFrame, sequence_metadata: pl.DataFrame
) -> pl.DataFrame:
    """Return per-cluster majority subtype + subtype entropy (mixedness)."""
    if "subtype" not in sequence_metadata.columns:
        raise ValueError("sequence_metadata missing 'subtype'")

    joined = cluster_assignments.join(
        sequence_metadata.select("accession", "subtype"),
        on="accession",
        how="inner",
    )

    rows: list[dict] = []
    for (cluster_id,), group in joined.group_by(["cluster_id"], maintain_order=True):
        subtypes = group["subtype"].to_list()
        counter = Counter([s for s in subtypes if s])
        total = sum(counter.values())
        if total == 0:
            rows.append(
                {
                    "cluster_id": cluster_id,
                    "majority_subtype": None,
                    "majority_fraction": None,
                    "subtype_entropy": None,
                }
            )
            continue
        majority_subtype, majority_n = counter.most_common(1)[0]
        probs = [c / total for c in counter.values()]
        # Shannon entropy in bits — 0 = pure, higher = mixed.
        import math

        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
        rows.append(
            {
                "cluster_id": cluster_id,
                "majority_subtype": majority_subtype,
                "majority_fraction": majority_n / total,
                "subtype_entropy": entropy,
            }
        )
    return pl.from_dicts(rows) if rows else pl.DataFrame()


def cluster_geographic_spread(
    cluster_assignments: pl.DataFrame, sequence_metadata: pl.DataFrame
) -> pl.DataFrame:
    """Per-cluster: how many countries represented, majority country, fraction."""
    if "country" not in sequence_metadata.columns:
        raise ValueError("sequence_metadata missing 'country'")

    joined = cluster_assignments.join(
        sequence_metadata.select("accession", "country"),
        on="accession",
        how="inner",
    ).filter(pl.col("country").is_not_null())

    rows: list[dict] = []
    for (cluster_id,), group in joined.group_by(["cluster_id"], maintain_order=True):
        countries = group["country"].to_list()
        counter = Counter(countries)
        total = sum(counter.values())
        if total == 0:
            continue
        majority_country, majority_n = counter.most_common(1)[0]
        rows.append(
            {
                "cluster_id": cluster_id,
                "n_countries": len(counter),
                "majority_country": majority_country,
                "majority_country_fraction": majority_n / total,
            }
        )
    return pl.from_dicts(rows) if rows else pl.DataFrame()
