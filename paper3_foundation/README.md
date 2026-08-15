# `paper3_foundation/` — foundation-model resistance benchmark manuscript

Substrate for **Paper 3** — subtype-stratified evaluation of ESM-2 / Evo vs classical baselines on HIV-1 drug-resistance prediction.

## Contents

- `manuscript.qmd` — Quarto skeleton with 4 figure placeholders.
- `refs.bib` — starter bibliography (ESM-2, Evo, Stanford HIVdb, SAWA v0.1.0).
- `figures/` — build artifacts (PDF + PNG); DVC-tracked.

## Code lives in `foundation_models/`

See [`foundation_models/README.md`](../foundation_models/README.md) for the shipped zero-GPU components (data, k-mers, baselines, evaluate) and the GPU-required skeletons (ESM-2, Evo, train, failure_analysis).

## Reproduce the baseline benchmark

```bash
# 1. Ensure Paper 1 outputs exist
make ingest sierra curate

# 2. Fit baselines + evaluate (adapter → LR + RF → per-subtype metrics)
python -m foundation_models.baselines  # once wired to a CLI

# 3. Render Fig 1
python notebooks/fig01_paper3_baselines.py  # TBA
```

## Preregistration

`paper3_foundation/preregistration.md` — TBA. Locks: split strategy (subtype_holdout), baseline hyperparameters, primary metric (AUROC gap between train-subtype and non-train subtypes), minimum positives per cell (5).
