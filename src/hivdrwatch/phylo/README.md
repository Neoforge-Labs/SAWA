# `phylo` — per-subtype phylogenetics

Builds per-subtype maximum-likelihood trees with IQ-TREE 2 and rescales them on a molecular clock with TreeTime. Subtype-stratification is deliberate: a joint tree over A1 + C + D is dominated by subtype topology and hides transmission signal.

## Method

- Model: `GTR+F+I+G4` (default; overridable via `--model`).
- Options: `-fast -redo --seed 42` for reproducible fast runs; `--no-fast` for the ML-quality path used before submission.
- Coverage filter: only sequences aligning ≥ `MIN_ALIGNED_COVERAGE` (default 0.5) of the pol reference enter the tree.
- Skip: subtypes with < `MIN_TAXA_FOR_TREE` (default 20) taxa are skipped with an INFO log — not a failure.
- Excluded: `EXCLUDED_SUBTYPES = {"unknown", "recombinant"}`. Both are heterogeneous labels — an "unknown" tree mixes A1/C/D indiscriminately; a "recombinant" tree lumps CRF_01_AE, CRF_02_AG, CRF_35_AD and friends into one bucket. A single tree over either label is phylogenetically incoherent. Extend as needed; per-CRF trees are the scientifically clean way to recover the recombinant signal.

## Audit trail

Every `make phylo` run emits `data/gold/trees/partition_report.parquet` alongside the manifest — one row per candidate subtype label seen in curated data. Columns: `subtype, n_available, n_after_qc, n_after_coverage, decision, decision_reason, min_taxa_used, min_coverage_used, generated_at`.

Decision codes:

- `built` — passed all gates, tree exists in manifest.
- `excluded_label` — subtype label is in `EXCLUDED_SUBTYPES`.
- `below_min_taxa` — fewer than `min_taxa` QC-passing sequences.
- `below_min_coverage` — QC passed, but fewer than `min_taxa` sequences met `min_coverage` (predominantly single-gene submissions; switch to per-gene trees to recover).

The EDA notebook §E cross-checks the partition report against the manifest's `n_taxa` and flags a "stale manifest" if on-disk trees don't reflect the current curated data — catches the case where curated has been re-ingested but `make phylo` hasn't been re-run.

## Live logging

Both `iqtree.py` and `treetime.py` shell out through `_streaming.py`, which pipes stdout+stderr line-by-line into a supplied logger. When called from a Prefect flow the logger is `get_run_logger()`, so IQ-TREE progress appears live in the Prefect UI.

## Run it

```bash
python -m hivdrwatch.pipeline.phylo \
    --curated data/interim/curated/aligned.parquet \
    --trees-root data/interim/trees \
    --manifest data/gold/trees/manifest.parquet \
    --min-coverage 0.5 \
    --min-taxa 20
```

Smoke a single subtype:

```bash
SUBTYPES=C make phylo-smoke
```

Output:

- one Newick per subtype under `data/interim/trees/<subtype>/`,
- `data/gold/trees/manifest.parquet` — `TreeMetadata` for every built tree,
- `data/gold/trees/partition_report.parquet` — per-subtype accounting for every candidate label (built, excluded, or below threshold), keyed to the thresholds used.

## Tests

```bash
pytest tests/phylo -v
```
