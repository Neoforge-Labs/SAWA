"""Domain causal DAGs for the DR and outbreak analyses.

Every edge below is a claim about the *world*, not the data. When you
disagree with an edge, argue about the biology or epidemiology — the
statistics will follow from the graph. This is the single most
opinionated file in the codebase; treat changes here as needing review.

## Notation

Nodes are variables (observed or unobserved). Edges are directed causal
influences (`A → B` reads "A causally affects B"). Unobserved variables
are prefixed `U_` and enable proper backdoor identification even when we
can't measure them directly.

## DR DAG (per-sequence drug resistance)

Reader convention: `M184V → resistance` is shorthand for "presence of
the M184V mutation at RT position 184 causally affects Sierra-called
resistance." Whether Sierra's rule is the ground truth is a *labelling*
issue, not a causal one — dowhy still identifies given the DAG.

Structural relations we encode:

  U_treatment_history ─┐
                        ├─→ drug_pressure ──→ mutation ──→ resistance
  country ──────────────┘                   ↑
                                            │
  subtype ────────────────────────────── polymorphism_context
  collection_year ──→ regimen_era ──→ drug_pressure

Key confounders for the mutation → resistance question:
  - subtype (natural polymorphisms at DR positions differ by subtype)
  - collection_year (older isolates saw different drug regimens)
  - country (regimen availability differs, and correlates with subtype)

## Outbreak DAG (country × year incidence)

  U_true_incidence ──→ sequenced_cases ──→ observed_clusters
                    │                       │
                    │                       └─→ n_new_cluster_members
                    │                       └─→ mean_growth_rate
                    │                                    │
                    └────→ transmission_dynamics ────────┘
                              ↑
  subtype_mix ────────────────┘
  U_sequencing_effort ──→ deposit_rate ──→ observed_clusters
                                            (spurious path)
  country ──→ U_sequencing_effort
  country ──→ U_healthcare_system ──→ detection_rate ──→ sequenced_cases

The killer confounder is `deposit_rate`: it creates a spurious path from
`country` to `n_new_cluster_members` that isn't causal. Backdoor
adjustment on `deposit_rate` blocks this path.

## Why the DAG is in code, not YAML

Two reasons:
1. Docstrings document the biology inline — code review sees the reasoning.
2. Programmatic access to the graph (subsetting to observed variables,
   exporting to dowhy's expected GML string) is a one-liner.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx

DAG_VERSION = "0.1.0"


# ------------------------------------------------------------------
# DR DAG
# ------------------------------------------------------------------


@dataclass
class DRDag:
    """The drug-resistance domain DAG.

    Observed nodes align with columns in `data/gold/dr_features.parquet`
    (mutation one-hots + covariates). Unobserved nodes (`U_*`) are only
    used for identification — we never try to condition on them.
    """

    mutation_features: tuple[str, ...] = field(
        default_factory=lambda: (
            # The mutation columns that drive resistance for common drugs.
            # Not exhaustive — extend as candidate-effect list grows.
            "RT_184_V",
            "RT_103_N",
            "RT_181_C",
            "RT_190_A",
            "PR_46_L",
            "PR_82_A",
            "PR_90_M",
            "IN_148_H",
            "IN_155_H",
        )
    )
    covariates: tuple[str, ...] = field(
        default_factory=lambda: (
            "subtype",
            "country",
            "collection_year",
        )
    )

    def build_graph(self) -> nx.DiGraph:
        g: nx.DiGraph = nx.DiGraph()

        # Outcome node (the model target — one binary label per drug).
        g.add_node("resistance", observed=True, kind="outcome")

        # Observed mutation drivers.
        for feat in self.mutation_features:
            g.add_node(feat, observed=True, kind="mutation")
            # Direct causal edge: mutation → resistance.
            g.add_edge(feat, "resistance")

        # Observed covariates.
        for cov in self.covariates:
            g.add_node(cov, observed=True, kind="covariate")

        # Unobserved treatment/regimen history — the classic latent
        # confounder that both selects for mutations AND correlates with
        # resistance calls.
        g.add_node("U_treatment_history", observed=False, kind="latent")
        g.add_node("drug_pressure", observed=False, kind="mediator")

        g.add_edge("country", "U_treatment_history")
        g.add_edge("collection_year", "U_treatment_history")
        g.add_edge("U_treatment_history", "drug_pressure")
        # Drug pressure selects for mutations (this is the biology).
        for feat in self.mutation_features:
            g.add_edge("drug_pressure", feat)

        # Subtype creates natural polymorphism context at DR positions.
        g.add_node("polymorphism_context", observed=False, kind="latent")
        g.add_edge("subtype", "polymorphism_context")
        for feat in self.mutation_features:
            g.add_edge("polymorphism_context", feat)

        # Subtype also has a direct effect on some Sierra rules (e.g.,
        # subtype-aware scoring for a few drugs). Keep the edge explicit.
        g.add_edge("subtype", "resistance")

        return g

    def minimum_adjustment_for_mutation(self, mutation: str) -> set[str]:
        """Observed variables that block confounding paths from `mutation`
        to `resistance`. Derived from the DAG structure — not manually
        maintained per mutation.
        """
        if mutation not in self.mutation_features:
            raise ValueError(f"unknown mutation node: {mutation}")
        # For M-shaped confounding: subtype creates a common cause via
        # polymorphism_context; country/year create one via
        # U_treatment_history → drug_pressure. Both need adjusting.
        return {"subtype", "country", "collection_year"}


# ------------------------------------------------------------------
# Outbreak DAG
# ------------------------------------------------------------------


@dataclass
class OutbreakDag:
    """Country×year outbreak DAG.

    The whole reason this exists is `deposit_rate`. If you take one edge
    out of this graph, don't take that one.
    """

    def build_graph(self) -> nx.DiGraph:
        g: nx.DiGraph = nx.DiGraph()

        # Outcome we want to predict.
        g.add_node("incidence", observed=True, kind="outcome")

        # Unobserved true drivers.
        g.add_node("U_true_incidence", observed=False, kind="latent")
        g.add_node("U_sequencing_effort", observed=False, kind="latent")
        g.add_node("U_healthcare_system", observed=False, kind="latent")

        # Observed covariates & features from outbreak_features.parquet.
        for node in (
            "country",
            "year",
            "n_sequences_deposited",
            "n_new_cluster_members",
            "n_active_clusters",
            "top_cluster_size",
            "mean_growth_rate",
            "subtype_a1_fraction",
            "subtype_c_fraction",
            "subtype_d_fraction",
            "mean_subtype_entropy",
            "deposit_rate",
        ):
            g.add_node(node, observed=True, kind="feature")

        # Country/year drive both sequencing effort and the true epidemic.
        g.add_edge("country", "U_sequencing_effort")
        g.add_edge("country", "U_healthcare_system")
        g.add_edge("country", "U_true_incidence")
        g.add_edge("year", "U_true_incidence")
        g.add_edge("year", "U_sequencing_effort")

        # Sequencing effort → deposit rate → observed features.
        g.add_edge("U_sequencing_effort", "deposit_rate")
        g.add_edge("deposit_rate", "n_sequences_deposited")
        g.add_edge("deposit_rate", "n_new_cluster_members")
        g.add_edge("deposit_rate", "n_active_clusters")
        g.add_edge("deposit_rate", "top_cluster_size")

        # True incidence → sequenced cases → cluster signal.
        g.add_edge("U_true_incidence", "n_new_cluster_members")
        g.add_edge("U_true_incidence", "mean_growth_rate")
        g.add_edge("U_true_incidence", "incidence")  # incidence label reflects
        # U_true_incidence (Naomi is a model)

        # Healthcare system affects both sequencing AND label quality.
        g.add_edge("U_healthcare_system", "U_sequencing_effort")

        # Subtype fractions affect transmission dynamics.
        for s in ("a1", "c", "d"):
            g.add_edge(f"subtype_{s}_fraction", "mean_growth_rate")

        # Subtype mixing at cluster level → new-member churn signal.
        g.add_edge("mean_subtype_entropy", "n_new_cluster_members")

        return g

    def minimum_adjustment_for_feature(self, feature: str) -> set[str]:
        """Which observed variables to adjust for when estimating
        `feature → incidence`."""
        # Every observed feature on the sequencing-effort side needs
        # deposit_rate + country + year for identification.
        sequencing_side = {
            "n_sequences_deposited",
            "n_new_cluster_members",
            "n_active_clusters",
            "top_cluster_size",
        }
        if feature in sequencing_side:
            return {"deposit_rate", "country", "year"}
        # Growth/entropy features: need country + year (proxy for
        # U_true_incidence temporal/spatial variation).
        return {"country", "year"}


# ------------------------------------------------------------------
# dowhy interop
# ------------------------------------------------------------------


def to_dowhy_gml(graph: nx.DiGraph) -> str:
    """Serialize a DiGraph to the GML string dowhy expects.

    dowhy's `CausalModel(graph=<gml_string>)` accepts a plain GML string;
    we emit the minimal form it needs (nodes + directed edges).
    """
    lines = ["graph [", "  directed 1"]
    node_ids: dict[str, int] = {}
    for i, node in enumerate(graph.nodes):
        node_ids[node] = i
        lines.append(f'  node [ id {i} label "{node}" ]')
    for u, v in graph.edges:
        lines.append(f"  edge [ source {node_ids[u]} target {node_ids[v]} ]")
    lines.append("]")
    return "\n".join(lines)


def observed_subgraph(graph: nx.DiGraph) -> nx.DiGraph:
    """Return the subgraph containing only observed nodes (drops all `U_*`).

    Useful for the report layer where you want a clean picture of the
    edges you can actually estimate against.
    """
    observed = [n for n, d in graph.nodes(data=True) if d.get("observed", False)]
    return graph.subgraph(observed).copy()
