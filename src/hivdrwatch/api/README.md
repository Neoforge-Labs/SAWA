# `api` — public read-only HTTP service

FastAPI application exposing the storage layer as JSON. Read-only in v1; the write path is gated behind CLAUDE.md's `require_dsa()` check pending KEMRI's Data Sharing Agreement.

## Endpoints

| Method | Path | Purpose |
|---|---|---|
| GET | `/health` | Liveness probe. |
| GET | `/version` | Package + schema-version dump. |
| GET | `/sequences` | Filtered sequence listing (country, subtype, year, source). |
| GET | `/dr-calls` | Sierra resistance calls (accession, drug_class, level). |
| GET | `/clusters` | Transmission-cluster memberships (cluster_id, accession, min size). |
| GET | `/docs` | Interactive OpenAPI docs (Swagger UI). |
| GET | `/redoc` | ReDoc-rendered API docs. |
| GET | `/openapi.json` | Raw OpenAPI 3.1 spec. |

Every response includes an `X-Request-Id` header (echoed from the caller when provided, else UUID4) and an `X-Response-Time-ms` timing header.

## Filters + pagination

All list endpoints accept:

- `?limit=N` (default 100, max 1000)
- `?offset=N` (default 0)

Plus endpoint-specific filters — see `/docs` for the current schema.

## Configuration

Env vars prefixed `HIVDRWATCH_API_`:

- `HIVDRWATCH_API_HOST` (default `0.0.0.0`)
- `HIVDRWATCH_API_PORT` (default `8080`)
- `HIVDRWATCH_API_CORS_ORIGINS` (default `["*"]`; JSON list)
- `HIVDRWATCH_API_RATE_LIMIT_DEFAULT` (default `100/minute` per IP)
- `HIVDRWATCH_API_ENABLE_DOCS` (default `true`)

Database connection is inherited from `hivdrwatch.storage.session` — same env vars as everything else.

## Run it

```bash
# Local dev
uvicorn hivdrwatch.api.main:app --reload --host 0.0.0.0 --port 8080

# Compose service (see docker-compose.yml)
docker compose up -d api
curl http://localhost:8080/health
curl http://localhost:8080/docs
```

## Tests

```bash
pytest tests/api -v
```

Tests use SQLite via `dependency_overrides[get_db]` — no live Postgres required.
