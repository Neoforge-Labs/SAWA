# `paper2_causal/` — CausalPhylo manuscript + figures

Substrate for **Paper 2** — causal identification of drug-resistance drivers using SAWA's Paper 1 outputs as the data source.

## Contents

- `manuscript.qmd` — Quarto manuscript scaffold with 4 figure placeholders (DAG, effect forest, subtype interaction, refutation landscape). Renders to PDF + HTML.
- `refs.bib` — starter bibliography (DoWhy, Pearl, Sierra, WHO DTG, SAWA v0.1.0 self-citation).
- `figures/` — build artifacts (PDF + PNG); DVC-tracked.

## Code lives elsewhere

Everything the paper cites is code in `causal/` at the repo root:

- `causal/dag.py` — `DRDag`, `OutbreakDag`
- `causal/estimate.py` — DoWhy identify → estimate → 3-refuter loop
- `causal/pruning.py` — feature pruning from the effect frame
- `causal/report.py` — HTML report + Polars effect frame
- `causal/data_prep.py` — **adapter** turning Paper 1 Parquets into the frames DoWhy consumes
- `causal/flows/causal.py` — Prefect flow chaining everything end-to-end

## Regenerate the paper

```bash
# 1. Run Paper 1's ingest → sierra → curate first (or `dvc pull`).
make ingest sierra curate

# 2. Run the causal flow — reads dr_results + curated, writes
#    causal_effects.parquet + causal_report.html to data/gold/.
uv run python -m causal.flows.causal \
    --dr-parquet     data/interim/dr_results.parquet \
    --curated-parquet data/interim/curated/aligned.parquet \
    --out-dir        data/gold \
    --bootstrap-samples 200 \
    --per-subtype

# 3. Regenerate figures — TBA notebooks under paper2_causal/notebooks/.

# 4. Render the manuscript.
quarto render paper2_causal/manuscript.qmd
```

## Preregistration

`paper2_causal/preregistration.md` — TBA. Analysis plan (DAG, primary hypotheses, refutation thresholds, sensitivity analyses) is locked before running the estimator on the frozen data vintage, per CLAUDE.md's confirmatory-analysis rule.
