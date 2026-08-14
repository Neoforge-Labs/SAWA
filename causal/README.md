# `causal` — DAG-based identification for HIV phylodynamics (Paper 2)

Substrate for Paper 2 (CausalPhylo). Builds domain DAGs for HIV drug-resistance and outbreak-dynamics questions, runs DoWhy identify → estimate → refute, prunes features by causal criteria, and emits a self-contained HTML report.

Sits at the repo root — not under `src/hivdrwatch/` — because it reads from the Paper 1 outputs but doesn't ship as part of the surveillance library.

## Modules

| File | Purpose |
|---|---|
| `dag.py` | `DRDag` (mutation → resistance) and `OutbreakDag` (deposit-rate confounder baked in). Minimum-adjustment queries. |
| `estimate.py` | DoWhy identify + estimate + 3-refuter loop (placebo / random-common-cause / subset). Bootstrap CIs. Per-subtype stratified estimator. Batch drivers for DR and outbreak. |
| `pruning.py` | Rule-based feature pruning from the effect frame (drop when CI spans 0 with small effect, or fails ≥ 2/3 refutations); keep-wins-across-scopes. |
| `report.py` | Polars effect-frame serializer + self-contained HTML report grouped by data scope. |
| `schema.py` | `CausalEffect`, `RefutationResult`, `EstimationMethod`, `RefutationVerdict` — the shared types every module speaks. |
| `flows/causal.py` | Prefect flow: DR gold + outbreak gold → effects Parquet + HTML + prune configs. MLflow-tracked. |

## Run it

```bash
python -m causal.flows.causal \
    --dr-features data/gold/dr_features.parquet \
    --outbreak-features data/gold/outbreak_features.parquet \
    --out-dir data/gold/causal
```

Output: `effects.parquet`, `report.html`, `prune_dr.json`, `prune_outbreak.json`.

## Tests

```bash
pytest tests/causal -v
```

`test_causal.py` covers DAG construction, GML serialization, and the estimator wrappers' error paths. Full DoWhy round-trip isn't exercised in unit tests — DoWhy is heavy and its output is stochastic. Integration coverage comes from Paper 2's regeneration pass.
