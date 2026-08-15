# `foundation_models/` — Paper 3 substrate

Foundation-model resistance benchmark for **Paper 3**. Compares fine-tuned foundation models (ESM-2, ESM-C, Evo) against non-FM baselines on per-drug × per-subtype resistance prediction.

## Shipped (zero-GPU)

- **`data.py`** — subtype-stratified train/val/test splits from SAWA Paper 1 outputs. Two split strategies:
  - `stratified` — 70/15/15 stratified on (subtype × drug × label). Fair baseline.
  - `subtype_holdout` — train + val on subtype B; test on all non-B. **This is the paper's thesis-testing split.**
- **`kmers.py`** — nucleotide/protein k-mer featurizer. Sparse CSR matrix, deterministic vocabulary, `normalize=True` default.
- **`baselines.py`** — the two non-FM baselines FMs must beat:
  - `KmerLogisticRegression` (L2, liblinear solver)
  - `OneHotRandomForest` (positional one-hot + RF)
- **`evaluate.py`** — AUROC, PR-AUC, Brier, ECE per (drug × subtype) cell. `min_positives=5` filter.

## GPU-required (skeleton only — populate when compute lands)

- `models/esm2.py` — ESM-2 fine-tuning head
- `models/evo.py` — Evo nucleotide-FM wrapper
- `train.py` — subtype-stratified curriculum trainer
- `failure_analysis.py` — "where do FMs fail on non-B" analysis (Paper 3 Fig 2)

## Reproduce the baseline benchmark

```python
from pathlib import Path
import polars as pl
from foundation_models import build_fm_dataset
from foundation_models.baselines import KmerLogisticRegression, OneHotRandomForest
from foundation_models.evaluate import per_subtype_metrics

# 1. Load Paper 1 outputs + split
dr = pl.read_parquet("data/interim/dr_results.parquet")
curated = pl.read_parquet("data/interim/curated/aligned.parquet")
ds = build_fm_dataset(dr, curated, strategy="subtype_holdout", seed=42)

# 2. Fit + predict per drug
per_drug_test = ds.test.partition_by("drug", as_dict=True)
metric_rows = []
for drug, test_df in per_drug_test.items():
    train_df = ds.train.filter(pl.col("drug") == drug)
    if train_df.is_empty():
        continue
    lr = KmerLogisticRegression(k=3, seed=42).fit(
        train_df["sequence"].to_list(),
        train_df["label_binary"].to_list(),
    )
    probs = lr.predict_proba(test_df["sequence"].to_list())
    scored = test_df.with_columns(pl.Series("prob_lr", probs))
    metric_rows.append(per_subtype_metrics(scored, prob_col="prob_lr"))

metrics = pl.concat(metric_rows)
metrics.write_parquet("data/gold/paper3_baseline_metrics.parquet")
```

This is what the FMs need to beat, per subtype. If a fine-tuned ESM-2 doesn't beat `KmerLogisticRegression` on subtype A1, that's a headline result.

## Tests

```bash
pytest tests/foundation_models -v
```
