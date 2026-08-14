# Deployment

Two audiences: (1) a **reviewer** who just wants to reproduce a paper, and (2) an **operator** running HIV-DR-Watch as a live surveillance service.

## Prereqs

- Docker + Docker Compose v2
- `uv` (Python dep manager) — https://docs.astral.sh/uv/
- `dvc[s3]` — installed via `uv sync` or `pip install "dvc[s3]"`

That's it. Everything else (Postgres, MinIO, MLflow, Prefect, tn93, IQ-TREE, MAFFT, HIV-TRACE) runs in containers.

---

## Reviewer path — reproduce a paper

```bash
git clone https://github.com/Neoforge-Labs/SAWA
cd hiv-dr-watch
git checkout paper1-v1.0.0             # pinned Zenodo-linked tag
make setup                              # .env + uv sync
make up                                 # postgres + minio + prefect + mlflow
make images                             # tn93 + MAFFT + IQ-TREE + HIV-TRACE
dvc pull                                # download frozen data vintage from Zenodo/MinIO
make paper1                             # regenerate every Paper 1 figure
```

Output figures land under `paper1_tools/figures/*.{pdf,png}` and should be **bit-identical** to the versions in the published preprint.

---

## Operator path — running HIV-DR-Watch as a live service

### 1. First-time infrastructure

```bash
cp .env.example .env                    # rotate every password before exposing beyond localhost
make up
make images
```

### 2. First-time DVC remote

The bundled MinIO in `docker-compose.yml` is dev-only. For production, point DVC at a durable S3 bucket:

```bash
dvc remote modify --local minio url s3://<your-prod-bucket>/
dvc remote modify --local minio endpointurl <your-s3-endpoint>
dvc remote modify --local minio access_key_id <your-key>
dvc remote modify --local minio --password secret_access_key <your-secret>
```

The `--local` flag stores credentials in `.dvc/config.local` which is gitignored — never commit real credentials.

### 3. Tracking raw data with DVC

Any new download of raw LANL or GenBank data:

```bash
# Drop the files under data/raw/... as normal (see ingest/README.md).
dvc add data/raw/lanl/ke/results.tsv
git add data/raw/lanl/ke/results.tsv.dvc data/raw/lanl/ke/.gitignore
git commit -m "data(ingest): refresh Kenya LANL export (2026-08 vintage)"
dvc push                                # ship the file to the remote
```

Downstream: run the pipeline. Its Parquet outputs go under `data/interim/`, `data/processed/`, `data/gold/` and should also be `dvc add`ed after each canonical run.

### 4. Scheduled surveillance (Prefect)

```bash
make worker                             # host-side Prefect worker
prefect deployment run 'hivdrwatch/ingest-incremental'
prefect deployment run 'hivdrwatch/sierra-refresh'
```

Deployments are declared under `pipeline/deployments.py` (Phase 2 — not yet shipped). Until then, run the stages manually via `hivdrwatch <subcommand>`.

### 5. Public API

```bash
docker compose up -d api                # FastAPI on http://localhost:8080
curl http://localhost:8080/health
curl http://localhost:8080/docs         # auto-generated OpenAPI docs
```

Read-only endpoints. Write endpoints (Phase 2) require the DSA sign-off with KEMRI.

---

## Backups

- **Postgres** — `pg_dumpall` nightly; retention per your policy.
- **MinIO** — versioning enabled by default in the compose config; bucket-level replication to a second region for disaster recovery.
- **DVC remote** — treat like any other S3 bucket; versioning + lifecycle rules apply.

## Monitoring

- **Prefect UI** — flow-level runs and failures.
- **MLflow** — pipeline-stage metrics per run (row counts, wall time, throughput).
- **Any external APM** (Sentry, Datadog) — set the standard `SENTRY_DSN` / `DD_*` env vars and the FastAPI + Prefect processes will pick them up.

## Rolling back

Every published preprint has a `paperN-vX.Y.Z` git tag AND a Zenodo DOI. Rolling back a bad release:

```bash
git checkout paper1-v1.0.0
dvc checkout                            # snap data/ back to the tagged vintage
docker compose up -d
```

If a schema-breaking change slipped through, `SCHEMA_VERSION` in `src/hivdrwatch/ingest/schema.py` will differ; re-ingest is the safe path.

## Security

- Never commit `.env`, `.dvc/config.local`, or KEMRI-related credentials.
- Env vars only for API keys (NCBI, Sierra, Zenodo).
- KEMRI data is untouchable until the DSA is signed (CLAUDE.md); code paths that could ever run against it must have a `require_dsa()` guard at entry.
