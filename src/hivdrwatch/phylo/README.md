# `phylo` — per-subtype phylogenetics

Builds per-subtype maximum-likelihood trees with IQ-TREE 2 and rescales them on a molecular clock with TreeTime. Subtype-stratification is deliberate: a joint tree over A1 + C + D is dominated by subtype topology and hides transmission signal.

## Method

- Model: `GTR+F+I+G4` (default; overridable via `--model`).
- Options: `-fast -redo --seed 42` for reproducible fast runs; `--no-fast` for the ML-quality path used before submission.
- Coverage filter: only sequences aligning ≥ `MIN_ALIGNED_COVERAGE` (default 0.5) of the pol reference enter the tree.
- Skip: subtypes with < `MIN_TAXA_FOR_TREE` (default 20) taxa are skipped with an INFO log — not a failure.
- Excluded: `EXCLUDED_SUBTYPES = {"unknown"}`; extend as needed.

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

Output: one Newick per subtype under `data/interim/trees/<subtype>/`, plus a `manifest.parquet` listing what got built and what got skipped.

## Tests

```bash
pytest tests/phylo -v
```
