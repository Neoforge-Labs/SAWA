# `paper3_foundation/` — HIV-1 FM benchmark manuscript, prereg, and figures

Substrate for **Paper 3** — a protocol / benchmark manuscript on
subtype-stratified evaluation of foundation models (ESM-2, Evo) vs
non-FM baselines for HIV-1 drug-resistance prediction. Frozen and
Zenodo-DOI-able at tag `paper3-v0.1.0`; the full FM comparison ships
as Paper 3 v2 after the SAWA v0.1.0 real-data run and fine-tune.

## Contents

- `manuscript.qmd` — full v1 Quarto manuscript (protocol + non-FM
  baselines + synthetic-ground-truth results + preregistered FM
  comparison).
- `preregistration.md` — locked plan for the v2 FM comparison.
  Filed on OSF at tag SHA.
- `refs.bib` — bibliography.
- `figures/` — build artifacts (PDF + PNG); DVC-tracked.

## Regenerate

```bash
# 1. Run the synthetic-ground-truth baseline benchmark that
#    generates the numbers cited in §Results.
make paper3-bench

# 2. Render the Quarto manuscript.
make paper3-render
# → paper3_foundation/manuscript.pdf and .html
```

## Cite

Once tagged and Zenodo-issued:

```
Njeri K. (2026). A subtype-stratified benchmark for HIV-1
drug-resistance prediction: protocol, non-FM baselines, and
preregistered foundation-model comparison [v0.1.0]. Zenodo.
https://doi.org/10.5281/zenodo.XXXXXXX
```

## Code lives in `foundation_models/`

See [`foundation_models/README.md`](../foundation_models/README.md)
for the shipped zero-GPU components (data, k-mers, baselines,
evaluate) and the GPU-required v2 skeletons (ESM-2, Evo, train,
failure_analysis).

## v2 FM comparison

After tagging v0.1.0:

```bash
# Real-data baseline run
make fm-baselines STRATEGY=subtype_holdout

# GPU-required FM fine-tune (v2)
python -m foundation_models.train --model esm2 --split subtype_holdout ...
python -m foundation_models.train --model evo  --split subtype_holdout ...
```
