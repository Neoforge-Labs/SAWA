"""Unit tests for the causal analysis stack.

Full dowhy round-trip is *not* exercised in unit tests — dowhy is heavy
and its output is stochastic (bootstrap refuters). Instead we test:
  * DAG construction (nodes, edges, adjustment sets)
  * GML serialization is valid enough for dowhy
  * Estimate wrapper handles missing/insufficient data cleanly
  * Report generation shape
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import polars as pl
import pytest

from causal import (
    DRDag,
    OutbreakDag,
    effects_to_frame,
    estimate_dr_effect,
    estimate_outbreak_effect,
    observed_subgraph,
    to_dowhy_gml,
    write_html_report,
)
from causal.schema import (
    CausalEffect,
    EstimationMethod,
    RefutationResult,
    RefutationVerdict,
)


# ---------- DR DAG ----------


class TestDRDag:
    def test_mutation_features_are_nodes(self) -> None:
        g = DRDag().build_graph()
        assert "RT_184_V" in g.nodes
        assert "resistance" in g.nodes
        assert g.has_edge("RT_184_V", "resistance")

    def test_covariates_present(self) -> None:
        g = DRDag().build_graph()
        for cov in ("subtype", "country", "collection_year"):
            assert cov in g.nodes

    def test_confounder_edges_exist(self) -> None:
        g = DRDag().build_graph()
        # subtype → polymorphism_context → mutation (backdoor)
        assert g.has_edge("subtype", "polymorphism_context")
        assert g.has_edge("polymorphism_context", "RT_184_V")

    def test_minimum_adjustment_set(self) -> None:
        dag = DRDag()
        adj = dag.minimum_adjustment_for_mutation("RT_184_V")
        assert adj == {"subtype", "country", "collection_year"}

    def test_unknown_mutation_raises(self) -> None:
        with pytest.raises(ValueError, match="unknown mutation"):
            DRDag().minimum_adjustment_for_mutation("RT_999_Z")


# ---------- Outbreak DAG ----------


class TestOutbreakDag:
    def test_deposit_rate_confounds_cluster_signal(self) -> None:
        g = OutbreakDag().build_graph()
        assert g.has_edge("deposit_rate", "n_new_cluster_members")
        assert g.has_edge("U_true_incidence", "n_new_cluster_members")

    def test_adjustment_for_sequencing_side_includes_deposit_rate(self) -> None:
        dag = OutbreakDag()
        adj = dag.minimum_adjustment_for_feature("n_new_cluster_members")
        assert "deposit_rate" in adj
        assert "country" in adj
        assert "year" in adj

    def test_adjustment_for_dynamics_features(self) -> None:
        dag = OutbreakDag()
        adj = dag.minimum_adjustment_for_feature("mean_growth_rate")
        # Growth rate is on the U_true_incidence path — country + year only
        assert adj == {"country", "year"}


# ---------- GML serialization ----------


class TestGML:
    def test_serialization_shape(self) -> None:
        g = DRDag().build_graph()
        gml = to_dowhy_gml(g)
        assert gml.startswith("graph [")
        assert "directed 1" in gml
        assert 'label "resistance"' in gml
        # One node line per node, one edge line per edge
        assert gml.count("node [") == len(g.nodes)
        assert gml.count("edge [") == len(g.edges)

    def test_observed_subgraph_drops_latents(self) -> None:
        g = DRDag().build_graph()
        obs = observed_subgraph(g)
        assert all(not n.startswith("U_") for n in obs.nodes)
        assert "polymorphism_context" not in obs.nodes
        assert "resistance" in obs.nodes


# ---------- estimate wrappers ----------


class TestEstimateDR:
    def test_missing_column_returns_unidentified(self) -> None:
        df = pl.DataFrame(
            {"label_binary": [0, 1] * 20, "subtype": ["A1"] * 40}
        )
        eff = estimate_dr_effect(df, "RT_184_V")
        assert eff.identified is False
        assert eff.n_observations == 0
        assert "skipped" in eff.data_scope

    def test_insufficient_rows_returns_unidentified(self) -> None:
        df = pl.DataFrame(
            {
                "RT_184_V": [1, 0, 1, 0],
                "label_binary": [1, 0, 1, 0],
                "subtype": ["A1"] * 4,
                "country": ["KE"] * 4,
                "collection_year": [2020] * 4,
            }
        )
        eff = estimate_dr_effect(df, "RT_184_V")
        assert eff.identified is False
        assert eff.n_observations == 0


class TestEstimateOutbreak:
    def test_missing_feature_returns_unidentified(self) -> None:
        df = pl.DataFrame(
            {"label_incidence": [1.0] * 40, "country": ["KE"] * 40}
        )
        eff = estimate_outbreak_effect(df, "mean_growth_rate")
        assert eff.identified is False


# ---------- report ----------


class TestReport:
    def _fake_effect(self, treatment: str = "RT_184_V") -> CausalEffect:
        return CausalEffect(
            treatment=treatment,
            outcome="label_binary",
            adjustment_set=["country", "subtype"],
            data_scope="DR/drug=3TC",
            n_observations=1000,
            method=EstimationMethod.LOGISTIC_REGRESSION,
            estimate=0.87,
            identified=True,
            identification_criterion="backdoor",
            refutations=[
                RefutationResult(
                    method="placebo_treatment_refuter",
                    original_estimate=0.87,
                    new_estimate=0.02,
                    verdict=RefutationVerdict.PASSES,
                ),
                RefutationResult(
                    method="random_common_cause",
                    original_estimate=0.87,
                    new_estimate=0.85,
                    verdict=RefutationVerdict.PASSES,
                ),
            ],
            dag_version="0.1.0",
            estimated_at=date(2026, 8, 12),
        )

    def test_frame_shape(self) -> None:
        df = effects_to_frame([self._fake_effect(), self._fake_effect("PR_46_L")])
        assert df.height == 2
        assert "treatment" in df.columns
        assert "ref_pass" in df.columns
        assert df["ref_pass"].to_list() == [2, 2]
        assert df["ref_fail"].to_list() == [0, 0]

    def test_html_report_writes(self, tmp_path: Path) -> None:
        out = tmp_path / "report.html"
        write_html_report([self._fake_effect()], out)
        content = out.read_text()
        assert "<html>" in content
        assert "RT_184_V" in content
        assert "DR/drug=3TC" in content
