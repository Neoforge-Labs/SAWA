"""Assemble sequences into clusters via connected components on the
below-threshold pair graph.

This is the HIV-TRACE algorithm in its purest form: an edge is added
between two sequences iff their TN93 distance is ≤ threshold. Clusters
are connected components of that graph. Singletons (nodes with no
below-threshold neighbors) get their own single-member cluster so
downstream joins don't have to deal with missing rows.

Cluster IDs are assigned deterministically by descending cluster size,
then lexicographically by the smallest accession — same input, same IDs
across runs.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass

import networkx as nx
import polars as pl

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ClusterAssignments:
    """A cluster assignment result — one row per accession."""

    df: pl.DataFrame  # cols: accession, cluster_id, cluster_size
    threshold: float
    n_clusters: int
    n_singletons: int


def build_cluster_graph(distances: pl.DataFrame, all_accessions: Iterable[str]) -> nx.Graph:
    """Return the below-threshold adjacency graph with every accession as a node.

    All accessions are added as nodes even if they have no edges — this
    keeps singletons in the cluster set.
    """
    g = nx.Graph()
    g.add_nodes_from(all_accessions)
    for row in distances.iter_rows(named=True):
        g.add_edge(row["accession_a"], row["accession_b"], distance=row["distance"])
    logger.info("cluster graph: %d nodes, %d edges", g.number_of_nodes(), g.number_of_edges())
    return g


def assign_clusters(
    distances: pl.DataFrame,
    all_accessions: Iterable[str],
    *,
    threshold: float,
) -> ClusterAssignments:
    """Compute connected-component clusters and return assignments."""
    g = build_cluster_graph(distances, all_accessions)
    components = list(nx.connected_components(g))

    # Deterministic ordering: descending size, then ascending min accession.
    components.sort(key=lambda c: (-len(c), min(c)))

    rows: list[dict] = []
    singletons = 0
    for i, component in enumerate(components):
        cluster_id = f"C-{i:05d}"
        size = len(component)
        if size == 1:
            singletons += 1
        for accession in component:
            rows.append(
                {
                    "accession": accession,
                    "cluster_id": cluster_id,
                    "cluster_size": size,
                    "tn93_threshold": threshold,
                }
            )

    df = pl.from_dicts(rows) if rows else pl.DataFrame()
    logger.info(
        "assigned %d accessions → %d clusters (%d singletons) at threshold=%.4f",
        df.height,
        len(components),
        singletons,
        threshold,
    )
    return ClusterAssignments(
        df=df,
        threshold=threshold,
        n_clusters=len(components),
        n_singletons=singletons,
    )


def cluster_size_distribution(assignments: pl.DataFrame) -> pl.DataFrame:
    """Return a size-histogram: cluster_size → count of clusters that size."""
    per_cluster = assignments.select("cluster_id", "cluster_size").unique()
    return (
        per_cluster.group_by("cluster_size")
        .len()
        .rename({"len": "n_clusters"})
        .sort("cluster_size")
    )
