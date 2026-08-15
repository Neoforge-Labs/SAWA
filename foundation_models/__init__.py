"""Foundation-model resistance benchmark — the Paper 3 substrate.

Paper 3's thesis: fine-tuned foundation models (ESM-2, ESM-C, Evo)
underperform non-B HIV-1 subtypes when trained on subtype-B-dominant
data. Subtype-stratified fine-tuning closes the gap. The benchmark
compares FMs against classical baselines (k-mer logistic regression,
one-hot random forest) across per-drug × per-subtype cells.

## Zero-GPU foundations (shipped)

  * `data.py`      — subtype-stratified train/val/test splits from
                     SAWA's Paper 1 Parquet outputs (dr_results +
                     curated). Reproducible via seeded RNG.
  * `kmers.py`     — protein/nucleotide k-mer featurizer. Sparse or
                     dense output; supports arbitrary k.
  * `baselines.py` — scikit-learn LR + RF baselines. Run in seconds
                     on a laptop. These are what the FMs must beat.
  * `evaluate.py`  — AUROC, PR-AUC, Brier score, expected calibration
                     error. Per-subtype breakdown for Paper 3's Fig 2.

## GPU-required (skeleton only)

  * `models/esm2.py`   — ESM-2 fine-tuning head (placeholder).
  * `models/evo.py`    — Evo nucleotide FM wrapper (placeholder).
  * `train.py`         — subtype-stratified curriculum trainer
                         (placeholder).
  * `failure_analysis.py` — "where do FMs fail on non-B" analysis
                            (placeholder; Fig 2 of Paper 3).

Populate the GPU modules when compute lands. The zero-GPU components
produce a full baseline benchmark table that the fine-tuned models
must beat — de-risks the paper's thesis before you spend a single
A100-hour.
"""

from .data import (
    FMDataset,
    SplitCounts,
    build_fm_dataset,
    stratified_train_val_test_split,
)
from .evaluate import (
    ClassificationMetrics,
    compute_metrics,
    per_subtype_metrics,
)
from .kmers import kmer_counts, kmer_matrix

__all__ = [
    "ClassificationMetrics",
    "FMDataset",
    "SplitCounts",
    "build_fm_dataset",
    "compute_metrics",
    "kmer_counts",
    "kmer_matrix",
    "per_subtype_metrics",
    "stratified_train_val_test_split",
]
