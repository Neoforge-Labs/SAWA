# `paper2_causal/` — CausalPhylo manuscript, prereg, and figures

Substrate for **Paper 2** — a protocol / methodology manuscript on
causal identification of HIV-1 drug-resistance drivers using SAWA's
Paper 1 outputs as the data source. Frozen and Zenodo-DOI-able at
tag `paper2-v0.1.0`; the full real-data application ships as
CausalPhylo v2 once the SAWA v0.1.0 real-data run completes.

## Contents

- `manuscript.qmd` — full v1 Quarto manuscript (protocol +
  synthetic-ground-truth validation + preregistered real-data
  analysis).
- `preregistration.md` — locked analysis plan for the v2 real-data
  run. Filed on OSF at tag SHA.
- `refs.bib` — bibliography.
- `figures/` — build artifacts (PDF + PNG); DVC-tracked.

## Regenerate

```bash
# 1. Run the synthetic-ground-truth benchmark that generates the
#    numbers cited in §Results.
make paper2-bench

# 2. Render the Quarto manuscript.
make paper2-render
# → paper2_causal/manuscript.pdf and .html
```

## Cite

Once tagged and Zenodo-issued:

```
Njeri K. (2026). CausalPhylo: a causal-identification framework for
HIV-1 drug-resistance phylodynamics in sub-Saharan Africa (protocol,
software, synthetic validation) [v0.1.0]. Zenodo.
https://doi.org/10.5281/zenodo.XXXXXXX
```

## Code lives elsewhere

Everything the paper cites is code in `causal/` at the repo root:

- `causal/dag.py` — `DRDag`, `OutbreakDag`
- `causal/estimate.py` — DoWhy identify → estimate → 3-refuter loop
- `causal/pruning.py` — feature pruning from the effect frame
- `causal/report.py` — HTML report + Polars effect frame
- `causal/data_prep.py` — adapter turning Paper 1 Parquets into
  the frames DoWhy consumes
- `causal/flows/causal.py` — Prefect flow chaining everything
  end-to-end
- `benchmarks/paper2_synth_recovery.py` — synthetic-ground-truth
  benchmark

## v2 real-data application

After tagging v0.1.0:

```bash
# Ensure Paper 1 outputs exist (or `dvc pull`).
make ingest sierra curate

# Full identify → estimate → refute against real SAWA v0.1.0 data.
make causal PER_SUBTYPE=1
```
