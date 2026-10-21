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

## GPU workflow (Paper 3 headline)

Four commands, one Paper 3 result:

1. **`make fm-prepare`** — turns curated + Sierra Parquets into per-drug-class train/val/test Parquets under `data/gold/fm/<class>/<strategy>/`. Uses **the primary gene per drug class** — PR for PI, RT for NRTI/NNRTI, IN for INSTI — because that's the region Sierra actually scores. This is scientifically correct (no cross-gene leakage) AND ~8× faster than concatenated pol on attention.

2. **`make fm-esm-smoke`** — 1 epoch on 500 rows to verify GPU + tokenizer + training loop. Expected ~20-40s on an M-series Mac. If AUROC comes back as a real number (not NaN), your plumbing works.

3. **`make fm-esm-train`** — full 3-epoch fine-tune on INSTI. Real numbers on a first run (INSTI stratified, 22,950 train / 5,002 test): **test AUROC 0.71, AUPRC 0.086** in ~32 min on Apple Silicon MPS at BS=16. Loss trajectory 0.69 → 0.2, val monotonically improving epoch-over-epoch → try `--max-epochs 5-10` for a real paper run.

4. **`make fm-compare CLASS=INSTI`** — the essential head-to-head. Runs kmer-LR + one-hot-RF on the same test split as ESM-2 and emits `data/gold/fm/comparison_insti_stratified.{parquet,json}`. Without this, the ESM-2 number is unmoored. If ESM-2 beats kmer-LR by ≥0.03 AUROC, that's your Paper 3 headline. If it doesn't, that's the more interesting finding — write it up.

**Install path.** `uv sync --extra fm` installs torch/transformers/accelerate on your GPU rig. Deferred imports throughout keep the rest of SAWA importable without `[fm]` deps.

**Model choice.** Default backbone `facebook/esm2_t6_8M_UR50D` (~8M params, <2GB VRAM at BS=16). Scale-up (t12-35M / t30-150M) is a one-line `--model-name` change with the same tokenizer. Only go bigger after t6 clears the kmer-LR baseline — if it doesn't, more parameters won't rescue it.

**MPS gotchas.** fp16 autocast on MPS silently falls back to CPU for several ops in PyTorch 2.x. `esm.py` detects MPS and forces fp32 (with a warning), which is ~500× faster in practice on Apple Silicon. Override with `--no-fp16` on CUDA to silence the warning explicitly.

**Sequence-length trap.** Prior versions concatenated PR+RT+IN → ~830 AA sequences → measured 150s/step on MPS (broken). The per-drug-class primary gene (INSTI → IN → 289 AA) drops attention cost by (289/830)² ≈ 8× → 0.3-1s/step at steady state. Keep this default unless you have a specific reason to concatenate.

**Split strategies.**
- `stratified` — 70/15/15 stratified on (subtype × drug × label). Fair baseline; FMs look strong here regardless of subtype-generalization.
- `subtype_holdout` — train on `--holdout-train-subtypes` (default `A1 C D` — dominant East African), test on all others (B, unknown, recombinant, minor subtypes). **This is the paper's thesis-testing split.** Run with `make fm-esm-holdout`. If test AUROC drops meaningfully vs stratified, the FM learned subtype composition, not biology. If it holds up, that's the *stronger* result.

**Checkpoints** land at `data/gold/fm/checkpoints/<drug_class>/<strategy>/best.pt` alongside `run_metadata.json` (test AUROC/AUPRC + full config for reproducibility). Both are DVC-trackable — add them to `scripts/setup_dvc.sh` after your first successful run.

**Ambiguity filter.** `prepare_training_data.py --max-aa-x-fraction` drops AA sequences whose fraction of `X` (unknown residues) exceeds the threshold. Default 1.0 = disabled. On the current East African cohort, 0.50 drops ~57% of rows — verify positive rate before enabling.

## Reproduce the baseline benchmark (CPU-only, no GPU)

## Reproduce the baseline benchmark

```python
from pathlib import Path
import polars as pl
from foundation_models import build_fm_dataset
from foundation_models.baselines import KmerLogisticRegression, OneHotRandomForest
from foundation_models.evaluate import per_subtype_metrics

# 1. Load Paper 1 outputs + split
dr = pl.read_parquet("data/interim/dr_results.parquet")
# curated schema at v0.2.0+ exposes `aligned_pol`; build_fm_dataset
# auto-detects and gap-strips. Pass sequence_column= to override.
curated = pl.read_parquet("data/interim/curated/aligned.parquet").filter(pl.col("passes_qc"))
ds = build_fm_dataset(
    dr, curated, strategy="subtype_holdout", seed=42, holdout_train_subtypes=("A1", "C", "D")
)

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
