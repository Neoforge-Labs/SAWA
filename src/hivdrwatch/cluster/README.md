# `cluster` — TN93 distances + connected-components clustering + validation

Detects transmission clusters using the HIV-TRACE convention: TN93 pairwise distances below a 0.015 threshold, then connected components on the resulting graph. Also ships the ground-truth validation harness that anchors Paper 1 Figure 3.

## Method

- Distance: `tn93` (Docker image `hivdrwatch/tn93:1.0.13`) — the CDC HIV-TRACE distance tool, TN93 with ambiguity resolution.
- Threshold: 0.015 (1.5 %) — HIV-TRACE convention.
- Graph: undirected, one node per accession, one edge per below-threshold pair.
- Membership: connected components. One accession belongs to exactly one cluster. Singletons get their own component ID.
- Determinism: cluster IDs are hashes of the sorted member accessions, so a re-run on the same input yields identical IDs.

## Dynamics

`dynamics.py` adds per-cluster derived features consumed by causal analyses:

- `cluster_time_series` — new-member count per (cluster, year)
- `cluster_subtype_mix` — entropy of subtype composition
- `cluster_geographic_spread` — country diversity
- `cluster_growth_rate` — rolling window on time-series

## Validation

Paper 1 Fig 3 requires proving our clusterer recovers ground-truth clusters and matches HIV-TRACE on the same input.

- `simulate.simulate_clusters(...)` — JC69-based synthetic generator producing FASTA + `truth_membership()` mapping. Fully deterministic given seed.
- `validate.score_clusters(predicted, truth)` — pair-based P/R/F1 + Adjusted Rand Index. This is the standard metric HIV-TRACE, PhyClip and every other transmission-cluster paper reports.
- `hivtrace.py` — thin wrapper around the HIV-TRACE Docker image that returns membership in the same `{accession: cluster_id}` shape as our clusterer, so both can be scored against the same ground truth.

## Run it

```bash
# Real data
python -m hivdrwatch.pipeline.cluster \
    --aligned-fasta data/interim/curate/aligned.fasta \
    --curated data/interim/curated/aligned.parquet \
    --work-dir data/interim/cluster \
    --gold-dir data/gold

# Validation harness (Paper 1 Fig 3)
make cluster-validate
```

## Tests

```bash
pytest tests/cluster -v
```
