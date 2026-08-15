"""Unit tests for foundation_models.evaluate."""

from __future__ import annotations

import numpy as np
import polars as pl

from foundation_models.evaluate import (
    MIN_POSITIVES_TO_SCORE,
    compute_metrics,
    per_subtype_metrics,
)


class TestComputeMetrics:
    def test_perfect_classifier(self) -> None:
        y_true = [0, 0, 1, 1]
        y_prob = [0.1, 0.2, 0.9, 0.95]
        m = compute_metrics(y_true, y_prob)
        assert m.auroc == 1.0
        assert m.pr_auc == 1.0
        assert m.n == 4
        assert m.n_positive == 2

    def test_random_classifier(self) -> None:
        rng = np.random.default_rng(42)
        y_true = rng.integers(0, 2, size=1000)
        y_prob = rng.random(1000)
        m = compute_metrics(y_true, y_prob)
        assert m.auroc is not None
        assert 0.4 < m.auroc < 0.6  # noise floor

    def test_single_class_returns_none_for_ranking_metrics(self) -> None:
        y_true = [0] * 10
        y_prob = [0.3] * 10
        m = compute_metrics(y_true, y_prob)
        assert m.auroc is None
        assert m.pr_auc is None
        # But Brier + ECE are still computed
        assert m.brier is not None
        assert m.ece is not None

    def test_empty_returns_zeroed_row(self) -> None:
        m = compute_metrics([], [])
        assert m.n == 0
        assert m.n_positive == 0
        assert m.auroc is None

    def test_calibration_captured(self) -> None:
        # Well-calibrated: probs match label rate
        y_true = [0, 0, 0, 1, 1] * 20
        y_prob = [0.4] * 100  # constant → ECE ≈ |0.4 - 0.4| = 0
        m = compute_metrics(y_true, y_prob)
        assert m.brier is not None
        # brier for constant p=0.4 on 40% positives: (0.4 - actual)²
        assert m.brier < 0.3

    def test_shape_mismatch_raises(self) -> None:
        import pytest as _pytest

        with _pytest.raises(ValueError, match="shape"):
            compute_metrics([0, 1], [0.5])


class TestPerSubtypeMetrics:
    def test_row_per_drug_subtype_cell(self) -> None:
        df = pl.DataFrame(
            {
                "drug": ["3TC"] * 20 + ["EFV"] * 20,
                "subtype": (["A1"] * 10 + ["C"] * 10) * 2,
                "label_binary": [0, 1] * 20,
                "prob_test": [0.2, 0.8] * 20,  # near-perfect classifier
            }
        )
        metrics = per_subtype_metrics(df, prob_col="prob_test")
        assert set(
            zip(metrics["drug"].to_list(), metrics["subtype"].to_list(), strict=True)
        ) == {
            ("3TC", "A1"),
            ("3TC", "C"),
            ("EFV", "A1"),
            ("EFV", "C"),
        }
        assert all(auc > 0.95 for auc in metrics["auroc"].to_list())

    def test_cells_below_min_positives_dropped(self) -> None:
        df = pl.DataFrame(
            {
                "drug": ["3TC"] * 10 + ["EFV"] * 10,
                "subtype": ["A1"] * 20,
                "label_binary": [0, 1] * 5 + [0] * 10,  # 3TC has 5 pos, EFV has 0
                "prob_test": [0.5] * 20,
            }
        )
        metrics = per_subtype_metrics(df, prob_col="prob_test", min_positives=5)
        drugs = metrics["drug"].to_list()
        assert "3TC" in drugs
        assert "EFV" not in drugs  # dropped for insufficient positives

    def test_empty_input_returns_empty_frame_with_schema(self) -> None:
        df = pl.DataFrame(
            schema={
                "drug": pl.Utf8,
                "subtype": pl.Utf8,
                "label_binary": pl.Int8,
                "prob_test": pl.Float64,
            }
        )
        metrics = per_subtype_metrics(df, prob_col="prob_test")
        assert metrics.height == 0
        # Schema is still correct so downstream aggregations don't break
        assert "auroc" in metrics.columns


class TestConstants:
    def test_min_positives_default(self) -> None:
        assert MIN_POSITIVES_TO_SCORE == 5
