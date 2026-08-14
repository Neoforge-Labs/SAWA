# `resistance` — Stanford HIVdb Sierra client

Wraps the Stanford HIVdb Sierra GraphQL endpoint. Returns per-drug resistance levels (susceptible / potential-low / low / intermediate / high) and per-mutation scores against the pinned algorithm version.

## Design decisions

- **Sierra, not Sierra-local.** The hosted endpoint is the algorithm source of truth. We pin the returned `algorithmVersion` in every cache key so a Stanford upgrade never silently changes past results.
- **SHA256-keyed cache.** Cache key = `sha256(sequence + algorithm_version)`. Same sequence + same algorithm = zero network. Cache lives under `data/raw/hivdb/` (DVC-tracked, not git).
- **Batching to 100 sequences/request.** Sierra's public quota is generous; we stay well under it. 429s trigger tenacity retry with exponential backoff.
- **Offline mode.** With `HIVDRWATCH_OFFLINE=1`, the client only reads the cache; a miss raises rather than calling out. Used in CI.

## Run it

```bash
python -m hivdrwatch.pipeline.dr_ingest \
    --input data/interim/east_africa_pol.parquet \
    --cache-dir data/raw/hivdb \
    --out-dir data/interim
```

Output: `data/interim/dr_results.parquet` with `(accession, drug_class, drug, level, mutations, algorithm_version, called_at)`.

## Tests

```bash
pytest tests/resistance -v
```

`pytest-httpx` mocks the GraphQL endpoint; no real network in tests. The cache round-trip test uses a synthetic sequence under `tests/fixtures/`.
