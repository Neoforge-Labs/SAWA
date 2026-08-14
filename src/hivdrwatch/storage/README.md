# `storage` — Postgres + DuckDB persistence layer

The system of record for the live surveillance service. Every downstream row (resistance call, cluster assignment) traces back through `pipeline_run_id` to the canonical ingest run that produced its sequence.

## Backends

- **Postgres 16** — write-side system of record. Provisioned by the compose stack (`postgres:16-alpine` service). Optional PostGIS for country centroids (fall back to plain lat/lon otherwise).
- **DuckDB** — read-side analytics. Reads Parquet directly for ad-hoc queries and `postgres_scan`s the hosted schema when needed.

## Tables

| Table | Grain | Uniqueness |
|---|---|---|
| `pipeline_runs` | one canonical run | primary key `id` |
| `sequences` | one HIV-1 pol sequence | `(schema_version, accession)` |
| `dr_calls` | one Sierra call per (sequence, drug) | `(accession, drug_class, drug, algorithm_version)` |
| `cluster_memberships` | one sequence's cluster in one run | `(accession, pipeline_run_id)` |

Every table carries `pipeline_run_id` FK — you can always trace a value back to the exact ingest / sierra / cluster run that produced it.

## Migrations

Alembic-managed under `migrations/`. The initial migration creates all four tables + unique constraints. Subsequent schema changes bump `STORAGE_SCHEMA_VERSION`.

```bash
# Apply migrations (bootstrap or upgrade)
alembic -c src/hivdrwatch/storage/alembic.ini upgrade head
# Generate a new revision
alembic -c src/hivdrwatch/storage/alembic.ini revision -m "add cluster_dynamics_snapshot"
```

## Loader

```python
from hivdrwatch.storage import get_session
from hivdrwatch.storage.loader import load_sequences

with get_session() as sess:
    counts = load_sequences(Path("data/interim/east_africa_pol.parquet"), session=sess)
    print(f"inserted {counts.inserted}, skipped {counts.skipped} (dup)")
```

Idempotent — reruns upsert-or-skip via `ON CONFLICT DO NOTHING`.

## Analytics

```python
from hivdrwatch.storage.duckdb import open_analytics
from pathlib import Path

with open_analytics() as db:
    db.open_parquet_view(Path("data/interim/dr_results.parquet"), "dr")
    db.open_parquet_view(Path("data/interim/curated/aligned.parquet"), "cur")
    df = db.sql("""
        SELECT country, subtype, COUNT(*) AS n
        FROM dr JOIN cur USING (accession)
        WHERE level IN ('high', 'intermediate')
        GROUP BY 1, 2
        ORDER BY n DESC
    """)
```

## Configuration

Env vars (defaults match the compose stack):

- `HIVDRWATCH_DATABASE_URL` — full override, e.g. `sqlite:///./test.db` for tests
- `POSTGRES_USER`, `POSTGRES_PASSWORD`, `POSTGRES_HOST`, `POSTGRES_PORT` — component overrides
- `HIVDRWATCH_DATABASE_NAME` — defaults to `hivdrwatch`

## Tests

```bash
pytest tests/storage -v
```

Tests run against SQLite via `HIVDRWATCH_DATABASE_URL=sqlite:///:memory:` — no Postgres required.
