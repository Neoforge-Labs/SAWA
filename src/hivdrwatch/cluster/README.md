# `cluster` — TN93 distances + connected-components clustering + validation

Detects transmission clusters using the HIV-TRACE convention: TN93 pairwise distances below a 0.015 threshold, then connected components on the resulting graph. Also ships the ground-truth validation harness that anchors Paper 1 Figure 3.

## Method

- Distance: `distance.compute_pairwise_distances_python`, the in-process engine both `pipeline/cluster.py` and the Fig 3 harness use. `distance.run_tn93` wraps `tn93` (Docker image `hivdrwatch/tn93:1.0.17`), the CDC HIV-TRACE distance tool, and is kept as the reference implementation for A/B comparison.
- Denominator: the **pairwise overlap** — the columns where both sequences carry an unambiguous base. Not the alignment width. A curated row is 3,012 columns wide whether or not anything was placed on them, so dividing by the width counts a column where neither sequence has a base as agreement, and two all-gap rows come out as distance 0.0. See `docs/evaluations/cluster-distance-denominator.md`.
- Minimum overlap: 100 shared unambiguous sites (`DEFAULT_MIN_OVERLAP`, matching tn93's own `-l` default). A pair below it gets no row rather than a distance estimated from too few sites.
- Input: records that pass QC and carry at least `MIN_ALIGNED_COVERAGE` (50 %) non-gap content, the same two gates `phylo/build.py` applies to tree input.
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
    --gold-dir data/gold \
    --min-overlap 100 \
    --min-coverage 0.5

# Validation harness (Paper 1 Fig 3)
make cluster-validate
```

## Tests

```bash
pytest tests/cluster -v
```
