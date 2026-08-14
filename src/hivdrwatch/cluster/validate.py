"""Pair-based cluster-quality metrics.

The standard metric for evaluating transmission-cluster detection
against a ground-truth membership (from a simulator or a validated
epi curve) is the **same-cluster-pair** confusion matrix:

  For every unordered pair (i, j), i != j:
    same_pred = predicted[i] == predicted[j]
    same_true = truth[i] == truth[j]

    TP: same_pred AND same_true
    FP: same_pred AND NOT same_true
    FN: NOT same_pred AND same_true
    TN: NOT same_pred AND NOT same_true

Precision = TP / (TP + FP)
Recall    = TP / (TP + FN)
F1        = 2 * P * R / (P + R)

This is what HIV-TRACE, PhyClip, and every other transmission-cluster
paper reports. Also emitted: adjusted Rand index for a summary
robustness statistic (unaffected by different cluster label
conventions).

## Why not just Rand / NMI

They're useful sanity numbers but Fig 3 wants the pair-based confusion
because it maps directly to public-health decisions ("if we act on the
predicted clusters, what fraction of our alerts are real and what
fraction of real outbreaks did we miss?"). Precision/recall are the
right axes for a public-health surveillance tool.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from itertools import combinations

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ClusterValidationResult:
    """One evaluation of a predicted cluster labeling against ground truth."""

    n_sequences: int
    n_pairs: int

    tp: int
    fp: int
    fn: int
    tn: int

    precision: float
    recall: float
    f1: float
    adjusted_rand_index: float | None

    # Slice for diagnostics
    n_predicted_clusters: int
    n_truth_clusters: int

    def as_dict(self) -> dict:
        return {
            "n_sequences": self.n_sequences,
            "n_pairs": self.n_pairs,
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
            "tn": self.tn,
            "precision": self.precision,
            "recall": self.recall,
            "f1": self.f1,
            "adjusted_rand_index": self.adjusted_rand_index,
            "n_predicted_clusters": self.n_predicted_clusters,
            "n_truth_clusters": self.n_truth_clusters,
        }


def score_clusters(
    predicted: dict[str, str],
    truth: dict[str, str],
) -> ClusterValidationResult:
    """Score a predicted `{accession: cluster_id}` against a truth mapping.

    Accessions present in `truth` but missing from `predicted` are
    treated as singletons in the prediction (each in its own cluster).
    Accessions present in `predicted` but missing from `truth` are
    ignored — the ground truth doesn't know about them, so scoring them
    would be meaningless.
    """
    scored_accessions = sorted(set(predicted) & set(truth))
    if not scored_accessions:
        raise ValueError("no overlap between predicted and truth accessions")

    if len(scored_accessions) < len(predicted):
        n_dropped = len(predicted) - len(scored_accessions)
        logger.info(
            "dropped %d predicted accessions not in ground truth", n_dropped
        )

    tp = fp = fn = tn = 0
    for a, b in combinations(scored_accessions, 2):
        same_pred = predicted[a] == predicted[b]
        same_true = truth[a] == truth[b]
        if same_pred and same_true:
            tp += 1
        elif same_pred and not same_true:
            fp += 1
        elif not same_pred and same_true:
            fn += 1
        else:
            tn += 1

    n_pairs = tp + fp + fn + tn
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    ari = _adjusted_rand_index(predicted, truth, scored_accessions)

    n_pred = len({predicted[a] for a in scored_accessions})
    n_true = len({truth[a] for a in scored_accessions})

    result = ClusterValidationResult(
        n_sequences=len(scored_accessions),
        n_pairs=n_pairs,
        tp=tp, fp=fp, fn=fn, tn=tn,
        precision=precision, recall=recall, f1=f1,
        adjusted_rand_index=ari,
        n_predicted_clusters=n_pred,
        n_truth_clusters=n_true,
    )
    logger.info(
        "cluster validation: n=%d, P=%.3f, R=%.3f, F1=%.3f, ARI=%.3f "
        "(pred %d clusters, truth %d clusters)",
        result.n_sequences, precision, recall, f1, ari or float("nan"),
        n_pred, n_true,
    )
    return result


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------


def _adjusted_rand_index(
    predicted: dict[str, str],
    truth: dict[str, str],
    accessions: list[str],
) -> float | None:
    """Compute Adjusted Rand Index directly from mappings.

    Same formula sklearn uses; reimplemented here to avoid pulling in
    sklearn just for this. Returns None if degenerate (all same or all
    different in either partition, which would divide by zero).
    """
    try:
        # Use a light dependency-free implementation via scipy comb —
        # but scipy might not be installed. Roll our own binomial:
        def comb2(n: int) -> int:
            return n * (n - 1) // 2 if n >= 2 else 0

        # Build the contingency table: pred_labels × truth_labels
        pred_labels = sorted({predicted[a] for a in accessions})
        truth_labels = sorted({truth[a] for a in accessions})
        pi = {p: i for i, p in enumerate(pred_labels)}
        ti = {t: j for j, t in enumerate(truth_labels)}
        table = [[0] * len(truth_labels) for _ in pred_labels]
        for a in accessions:
            table[pi[predicted[a]]][ti[truth[a]]] += 1

        n = len(accessions)
        row_sums = [sum(row) for row in table]
        col_sums = [sum(table[i][j] for i in range(len(table))) for j in range(len(truth_labels))]

        sum_comb_cells = sum(comb2(cell) for row in table for cell in row)
        sum_comb_rows = sum(comb2(rs) for rs in row_sums)
        sum_comb_cols = sum(comb2(cs) for cs in col_sums)
        total = comb2(n)
        if total == 0:
            return None

        expected = sum_comb_rows * sum_comb_cols / total
        max_index = (sum_comb_rows + sum_comb_cols) / 2
        denom = max_index - expected
        if denom == 0:
            return None
        return (sum_comb_cells - expected) / denom
    except Exception as e:  # noqa: BLE001
        logger.warning("ARI computation failed: %s", e)
        return None


__all__ = [
    "ClusterValidationResult",
    "score_clusters",
]
