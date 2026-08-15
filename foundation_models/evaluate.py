"""Evaluation metrics for the Paper 3 benchmark.

Metrics per drug × subtype cell:
  * AUROC       — area under ROC curve
  * PR-AUC      — area under precision-recall curve (better for
                  imbalanced labels, common for rare mutations)
  * Brier score — probabilistic calibration
  * ECE         — expected calibration error (10-bin default)
  * n_positive  — positive-class count in the cell (for the
                  "insufficient signal" filter)

`per_subtype_metrics(df, prob_col, ...)` groups a labeled DataFrame
by drug × subtype and returns a Polars metrics frame ready for
Fig 2 of Paper 3.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

# Cells below this positive-class count aren't scored — the metric
# would be dominated by noise.
MIN_POSITIVES_TO_SCORE: int = 5


@dataclass(frozen=True)
class ClassificationMetrics:
    """One row of metrics for a (drug × subtype × model) cell."""

    n: int
    n_positive: int
    auroc: float | None
    pr_auc: float | None
    brier: float | None
    ece: float | None

    def as_dict(self) -> dict:
        return asdict(self)


def compute_metrics(
    y_true: np.ndarray | list[int],
    y_prob: np.ndarray | list[float],
    *,
    n_ece_bins: int = 10,
) -> ClassificationMetrics:
    """Point-estimate metrics for one cell.

    Returns None for the ranking metrics (AUROC, PR-AUC) when the cell
    has fewer than 2 classes present, since they're undefined. Brier
    + ECE are still computed as long as any predictions exist.
    """
    y_true_arr = np.asarray(y_true, dtype=np.int8)
    y_prob_arr = np.asarray(y_prob, dtype=np.float64)
    if y_true_arr.shape != y_prob_arr.shape:
        raise ValueError(f"y_true shape {y_true_arr.shape} != y_prob shape {y_prob_arr.shape}")

    n = int(y_true_arr.size)
    n_positive = int(y_true_arr.sum())

    if n == 0:
        return ClassificationMetrics(0, 0, None, None, None, None)

    try:
        from sklearn.metrics import (
            average_precision_score,
            brier_score_loss,
            roc_auc_score,
        )
    except ImportError as e:  # pragma: no cover
        raise ImportError("scikit-learn required for compute_metrics") from e

    n_classes_present = len(np.unique(y_true_arr))
    auroc = pr_auc = None
    if n_classes_present >= 2:
        auroc = float(roc_auc_score(y_true_arr, y_prob_arr))
        pr_auc = float(average_precision_score(y_true_arr, y_prob_arr))

    brier = float(brier_score_loss(y_true_arr, y_prob_arr))
    ece = _expected_calibration_error(y_true_arr, y_prob_arr, n_bins=n_ece_bins)

    return ClassificationMetrics(
        n=n,
        n_positive=n_positive,
        auroc=auroc,
        pr_auc=pr_auc,
        brier=brier,
        ece=ece,
    )


def per_subtype_metrics(
    df: pl.DataFrame,
    *,
    prob_col: str,
    label_col: str = "label_binary",
    subtype_col: str = "subtype",
    drug_col: str = "drug",
    min_positives: int = MIN_POSITIVES_TO_SCORE,
) -> pl.DataFrame:
    """Group a labeled DataFrame by (drug, subtype) → metrics rows.

    Cells with fewer than `min_positives` positive labels are omitted.
    Result columns:
        drug, subtype, n, n_positive, auroc, pr_auc, brier, ece
    """
    rows: list[dict] = []
    for (drug, subtype), group in df.group_by([drug_col, subtype_col]):
        y_true = group[label_col].to_numpy()
        y_prob = group[prob_col].to_numpy()
        m = compute_metrics(y_true, y_prob)
        if m.n_positive < min_positives:
            continue
        rows.append({"drug": drug, "subtype": subtype, **m.as_dict()})

    if not rows:
        return pl.DataFrame(
            schema={
                "drug": pl.Utf8,
                "subtype": pl.Utf8,
                "n": pl.Int64,
                "n_positive": pl.Int64,
                "auroc": pl.Float64,
                "pr_auc": pl.Float64,
                "brier": pl.Float64,
                "ece": pl.Float64,
            }
        )
    return pl.from_dicts(rows).sort(["drug", "subtype"])


# ------------------------------------------------------------------
# internals
# ------------------------------------------------------------------


def _expected_calibration_error(
    y_true: np.ndarray, y_prob: np.ndarray, *, n_bins: int
) -> float:
    """Expected calibration error with uniform-width bins."""
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    total = 0.0
    n = y_true.size
    for i in range(n_bins):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        # Last bin includes right edge so 1.0 doesn't fall out.
        mask = (y_prob >= lo) & (y_prob < hi if i < n_bins - 1 else y_prob <= hi)
        if not mask.any():
            continue
        bin_confidence = float(y_prob[mask].mean())
        bin_accuracy = float(y_true[mask].mean())
        total += (mask.sum() / n) * abs(bin_confidence - bin_accuracy)
    return float(total)


__all__ = [
    "MIN_POSITIVES_TO_SCORE",
    "ClassificationMetrics",
    "compute_metrics",
    "per_subtype_metrics",
]
