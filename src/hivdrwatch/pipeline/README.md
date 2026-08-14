# `pipeline` — Prefect 3 orchestration

Every stage of HIV-DR-Watch runs as a Prefect flow so a fresh checkout can `make up && make ingest && make curate && make sierra && make cluster && make phylo` and see the whole thing progress live in the Prefect UI at http://localhost:4200.

## Flows

| Module | Purpose |
|---|---|
| `ingest.py` | Merged LANL + GenBank ingest → bronze Parquet. |
| `dr_ingest.py` | Stanford HIVdb Sierra calls with sha256 cache. |
| `curate.py` | Align (`mappy`), trim, QC, dedup → silver Parquet + FASTA. |
| `cluster.py` | tn93 distances + connected components + dynamics → gold cluster assignments. |
| `phylo.py` | Per-subtype IQ-TREE 2 + TreeTime → trees under `data/interim/trees/`. |
| `cluster_validate.py` | Synthetic ground truth → our clusterer → P/R/F1/ARI. Extended in Batch B to also run HIV-TRACE side-by-side. |

Each flow is also a plain Python entry point — `python -m hivdrwatch.pipeline.<name>` works without a Prefect server, in which case it runs sequentially and logs to stderr.

## Storage tail (Batch G)

`ingest`, `dr_ingest`, and `cluster` each end with an optional Postgres load via `pipeline/_storage_tail.py`. Enable per-invocation:

```bash
# CLI flag
python -m hivdrwatch.pipeline.ingest --to-postgres --skip-genbank ...
# Env var (applies to every flow)
HIVDRWATCH_LOAD_TO_POSTGRES=1 make ingest sierra cluster
# Makefile shortcut
TO_POSTGRES=1 make ingest
```

The tail is **best-effort and never crashes the flow.** If the DB is unreachable, the tail logs a WARNING and continues; the Parquet on disk remains the source of truth. Idempotent — reruns upsert-or-skip via `ON CONFLICT DO NOTHING` (Postgres) or `INSERT OR IGNORE` (SQLite in tests).

## Worker

The Prefect worker imports the flow module at startup, so **restart the worker after any code change** to the flows or the modules they import:

```bash
make restart-worker
```

## Observability

- **Prefect UI** (`localhost:4200`) — flow runs, live logs, retries. Set `PREFECT_LOGGING_EXTRA_LOGGERS=hivdrwatch,causal` (already in `.env.example`) to forward stdlib logs.
- **MLflow** (`localhost:5000`) — run parameters, metrics, artifacts. Every task that produces a result JSON also logs to MLflow through `hivdrwatch.tracking.tracked_run`.
- **Live subprocess streaming** — IQ-TREE / TreeTime / MAFFT stdout+stderr stream line-by-line into the Prefect run logger via `phylo/_streaming.py`.

## Tests

```bash
pytest tests/pipeline -v
```
