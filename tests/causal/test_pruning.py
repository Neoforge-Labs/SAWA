"""Unit tests for causal-derived feature pruning."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from causal.pruning import (
    PruneConfig,
    apply_prune_config,
    derive_prune_config,
    load_prune_config,
)


def _effects_frame(rows: list[dict]) -> pl.DataFrame:
    """Match the shape emitted by causal/report.py:effects_to_frame."""
    base_defaults = {
        "outcome": "label_binary",
        "adjustment_set": "country,subtype",
        "method": "backdoor.generalized_linear_model",
        "identified": True,
        "criterion": "backdoor",
        "n": 500,
        "ci_low": None,
        "ci_high": None,
        "p_value": None,
        "dag_version": "0.1.0",
        "estimated_at": "2026-08-12",
        "ref_pass": 3,
        "ref_fail": 0,
        "ref_inconclusive": 0,
    }
    filled = [{**base_defaults, **r} for r in rows]
    return pl.from_dicts(filled)


class TestDeriveConfig:
    def test_drops_ci_spans_zero_small_estimate(self) -> None:
        df = _effects_frame(
            [
                {
                    "treatment": "RT_184_V",
                    "data_scope": "DR/drug=3TC",
                    "estimate": 0.9,
                    "ci_low": 0.5,
                    "ci_high": 1.3,
                },
                {
                    "treatment": "PR_82_A",
                    "data_scope": "DR/drug=3TC",
                    "estimate": 0.05,
                    "ci_low": -0.4,
                    "ci_high": 0.5,  # width 0.9, |est| 0.05 << 0.45
                },
            ]
        )
        cfg = derive_prune_config(df)
        assert "PR_82_A" in cfg.drop
        assert "RT_184_V" in cfg.keep
        assert "spans 0" in cfg.prune_reasons["PR_82_A"]

    def test_drops_repeated_refutation_failures(self) -> None:
        df = _effects_frame(
            [
                {
                    "treatment": "n_active_clusters",
                    "data_scope": "outbreak",
                    "estimate": 0.3,
                    "ref_pass": 1,
                    "ref_fail": 2,
                    "ref_inconclusive": 0,
                },
            ]
        )
        cfg = derive_prune_config(df)
        assert "n_active_clusters" in cfg.drop
        assert "fails 2/3" in cfg.prune_reasons["n_active_clusters"]

    def test_keeps_features_with_solid_refutations(self) -> None:
        df = _effects_frame(
            [
                {
                    "treatment": "RT_184_V",
                    "data_scope": "DR/drug=3TC",
                    "estimate": 0.9,
                },
            ]
        )
        cfg = derive_prune_config(df)
        assert "RT_184_V" in cfg.keep
        assert cfg.drop == ()

    def test_keep_wins_across_scopes(self) -> None:
        """If dropped in one scope but kept in another → keep."""
        df = _effects_frame(
            [
                {
                    "treatment": "RT_184_V",
                    "data_scope": "DR/drug=3TC",
                    "estimate": 0.05,
                    "ci_low": -0.4,
                    "ci_high": 0.5,
                },
                {
                    "treatment": "RT_184_V",
                    "data_scope": "DR/drug=EFV",
                    "estimate": 0.9,
                },
            ]
        )
        cfg = derive_prune_config(df)
        assert "RT_184_V" in cfg.keep
        assert "RT_184_V" not in cfg.drop

    def test_scope_prefix_filter(self) -> None:
        df = _effects_frame(
            [
                {
                    "treatment": "PR_82_A",
                    "data_scope": "DR/drug=3TC",
                    "estimate": 0.02,
                    "ci_low": -0.5,
                    "ci_high": 0.6,
                },
                {
                    "treatment": "PR_82_A",
                    "data_scope": "outbreak",  # irrelevant scope
                    "estimate": 0.9,
                },
            ]
        )
        cfg = derive_prune_config(df, scope_prefix="DR")
        assert "PR_82_A" in cfg.drop

    def test_unidentified_rows_ignored(self) -> None:
        df = _effects_frame(
            [
                {
                    "treatment": "PR_82_A",
                    "data_scope": "DR/drug=3TC",
                    "estimate": 0.0,
                    "identified": False,
                    "ref_fail": 3,
                },
            ]
        )
        cfg = derive_prune_config(df)
        # Unidentified rows shouldn't drive drop OR keep.
        assert "PR_82_A" not in cfg.drop
        assert "PR_82_A" not in cfg.keep


class TestApplyConfig:
    def test_prefix_match_drops_related_columns(self) -> None:
        cfg = PruneConfig(
            dag_version="0.1.0",
            generated_at=__import__("datetime").date(2026, 8, 12),
            keep=("RT_184",),
            drop=("PR_46",),
            prune_reasons={"PR_46": "test"},
        )
        cols = ["PR_46_L", "PR_46_I", "PR_46_covered", "RT_184_V", "IN_148_H"]
        kept, dropped = apply_prune_config(cols, cfg)
        assert set(dropped) == {"PR_46_L", "PR_46_I", "PR_46_covered"}
        assert set(kept) == {"RT_184_V", "IN_148_H"}

    def test_no_config_keeps_everything(self) -> None:
        cols = ["PR_46_L", "RT_184_V"]
        kept, dropped = apply_prune_config(cols, None)
        assert kept == cols
        assert dropped == []


class TestLoadRoundtrip:
    def test_load_missing_returns_none(self, tmp_path: Path) -> None:
        assert load_prune_config(tmp_path / "does_not_exist.json") is None

    def test_write_load_roundtrip(self, tmp_path: Path) -> None:
        cfg = PruneConfig(
            dag_version="0.1.0",
            generated_at=__import__("datetime").date(2026, 8, 12),
            keep=("A", "B"),
            drop=("C",),
            prune_reasons={"C": "spurious"},
        )
        p = cfg.write(tmp_path / "prune.json")
        loaded = load_prune_config(p)
        assert loaded == cfg

        payload = json.loads(p.read_text())
        assert payload["drop"] == ["C"]
        assert payload["prune_reasons"] == {"C": "spurious"}
