"""Unit tests for clustering: assembly, dynamics, and graph construction.

The `hivdrwatch.features.outbreak_features` module and its Naomi
label-join helpers were deleted in the Phase 0 cleanup (see MODULES.md
Phase 0 alignment note) — they were out-of-scope for Paper 1. Tests
that covered those symbols were removed with them.
"""

from __future__ import annotations

import polars as pl
import pytest

from hivdrwatch.cluster.assembly import (
    assign_clusters,
    build_cluster_graph,
    cluster_size_distribution,
)
from hivdrwatch.cluster.dynamics import (
    cluster_geographic_spread,
    cluster_growth_rate,
    cluster_subtype_mix,
    cluster_time_series,
)

# ---------- assembly ----------


class TestClusterAssembly:
    def test_connected_components(self) -> None:
        distances = pl.DataFrame(
            {
                "accession_a": ["A", "B", "D"],
                "accession_b": ["B", "C", "E"],
                "distance": [0.01, 0.01, 0.01],
            }
        )
        all_accessions = ["A", "B", "C", "D", "E", "F"]
        result = assign_clusters(distances, all_accessions, threshold=0.015)

        # Groups: {A,B,C}, {D,E}, {F singleton} → 3 clusters
        assert result.n_clusters == 3
        assert result.n_singletons == 1

        by_acc = {r["accession"]: r for r in result.df.iter_rows(named=True)}
        # A/B/C share a cluster_id
        assert by_acc["A"]["cluster_id"] == by_acc["B"]["cluster_id"]
        assert by_acc["A"]["cluster_id"] == by_acc["C"]["cluster_id"]
        assert by_acc["A"]["cluster_size"] == 3
        # D/E share
        assert by_acc["D"]["cluster_id"] == by_acc["E"]["cluster_id"]
        # F alone
        assert by_acc["F"]["cluster_size"] == 1

    def test_deterministic_id_ordering(self) -> None:
        distances = pl.DataFrame(
            {"accession_a": ["A"], "accession_b": ["B"], "distance": [0.01]}
        )
        r1 = assign_clusters(distances, ["A", "B", "C"], threshold=0.015)
        r2 = assign_clusters(distances, ["A", "B", "C"], threshold=0.015)
        assert r1.df.to_dicts() == r2.df.to_dicts()

    def test_size_distribution(self) -> None:
        assignments = pl.DataFrame(
            {
                "accession": ["A", "B", "C", "D"],
                "cluster_id": ["C-1", "C-1", "C-2", "C-3"],
                "cluster_size": [2, 2, 1, 1],
                "tn93_threshold": [0.015] * 4,
            }
        )
        dist = cluster_size_distribution(assignments)
        d = {
            int(row["cluster_size"]): int(row["n_clusters"])
            for row in dist.iter_rows(named=True)
        }
        assert d == {1: 2, 2: 1}


# ---------- dynamics ----------


@pytest.fixture
def toy_data() -> tuple[pl.DataFrame, pl.DataFrame]:
    metadata = pl.DataFrame(
        {
            "accession": ["A1", "A2", "A3", "B1", "B2", "C1"],
            "country": ["KE", "KE", "KE", "UG", "UG", "TZ"],
            "collection_year": [2018, 2019, 2020, 2019, 2020, 2020],
            "subtype": ["A1", "A1", "A1", "C", "C", "D"],
        }
    )
    cluster_df = pl.DataFrame(
        {
            "accession": ["A1", "A2", "A3", "B1", "B2", "C1"],
            "cluster_id": ["C-1", "C-1", "C-1", "C-2", "C-2", "C-3"],
            "cluster_size": [3, 3, 3, 2, 2, 1],
            "tn93_threshold": [0.015] * 6,
        }
    )
    return metadata, cluster_df


class TestClusterTimeSeries:
    def test_counts_new_members_per_year(self, toy_data) -> None:
        metadata, cluster_df = toy_data
        ts = cluster_time_series(cluster_df, metadata)
        # C-1 KE 2018 → 1 new member; 2019 → 1; 2020 → 1
        c1_ke = ts.filter((pl.col("cluster_id") == "C-1") & (pl.col("country") == "KE")).sort(
            "collection_year"
        )
        assert c1_ke["n_new_members"].to_list() == [1, 1, 1]


class TestClusterGrowthRate:
    def test_growth_over_window(self, toy_data) -> None:
        metadata, cluster_df = toy_data
        ts = cluster_time_series(cluster_df, metadata)
        growth = cluster_growth_rate(ts, window_years=2)
        assert growth.height > 0
        # Growth rate is always finite and non-negative.
        assert (growth["growth_rate"] >= 0).all()


class TestClusterSubtypeMix:
    def test_majority_and_entropy(self, toy_data) -> None:
        metadata, cluster_df = toy_data
        mix = cluster_subtype_mix(cluster_df, metadata)
        by_id = {r["cluster_id"]: r for r in mix.iter_rows(named=True)}
        assert by_id["C-1"]["majority_subtype"] == "A1"
        assert by_id["C-1"]["majority_fraction"] == 1.0
        assert by_id["C-1"]["subtype_entropy"] == 0.0  # pure


class TestClusterGeographicSpread:
    def test_single_country(self, toy_data) -> None:
        metadata, cluster_df = toy_data
        spread = cluster_geographic_spread(cluster_df, metadata)
        by_id = {r["cluster_id"]: r for r in spread.iter_rows(named=True)}
        assert by_id["C-1"]["n_countries"] == 1
        assert by_id["C-1"]["majority_country"] == "KE"


# ---------- graph builder ----------


class TestBuildClusterGraph:
    def test_isolated_nodes_preserved(self) -> None:
        distances = pl.DataFrame(
            {"accession_a": [], "accession_b": [], "distance": []},
            schema={"accession_a": pl.Utf8, "accession_b": pl.Utf8, "distance": pl.Float64},
        )
        g = build_cluster_graph(distances, ["A", "B", "C"])
        assert set(g.nodes) == {"A", "B", "C"}
        assert g.number_of_edges() == 0
