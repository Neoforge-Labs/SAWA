# `ingest` — LANL and GenBank ingestion

Two sources, one canonical `SequenceRecord` output. Every downstream stage reads only that record; sources beyond LANL and GenBank plug in by producing the same shape.

## Sources

**GenBank** (programmatic, incremental). Pulled via Biopython Entrez. Set `NCBI_API_KEY` in `.env` to raise the rate limit from 3 → 10 req/s. Resume-safe: partial batches are cached and validated against the terminator (`//`); on repeated failure the batch bisects until the poison record is isolated. Per-country watermarks record the last successful MDAT so incremental runs don't refetch history.

**LANL** (manual, curated). LANL rate-blocks scrapers, so downloads are manual and version-controlled via DVC. Procedure:

1. https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html
2. Virus **HIV-1** · Subtype **A1, C, D** (ctrl-click) + tick "include recombinants"
3. Genomic region **Pol** · Min length **800 bp**
4. Sampling country: one country per search
5. Download → format **Table** (TSV) · check all metadata columns
6. Save as `data/raw/lanl/<cc>/results.tsv` (also accepts `results.txt`)
7. Optional: also download the alignment FASTA and save alongside as `hiv-db.fasta` for the fallback sequence source
8. Repeat for KE, UG, TZ, RW, BI, ET

## Contract

- `subtype` — LANL curation wins when present; GenBank note field is the fallback.
- `country` — LANL `Country` (mapped to ISO 3166-1 alpha-2) wins; folder name is the fallback; requested country is the last-resort default.
- Duplicate accessions across LANL files — later mtime wins.
- Sequences shorter than `--min-length` (default 800 bp) are dropped; phylogenetic signal is too weak below that.

## Run it

```bash
# LANL only (offline)
python -m hivdrwatch.ingest \
    --lanl-dir data/raw/lanl \
    --out data/interim/east_africa_pol.parquet \
    --skip-genbank

# LANL + GenBank (needs NCBI_EMAIL in .env)
python -m hivdrwatch.pipeline.ingest \
    --email $NCBI_EMAIL \
    --countries KE UG TZ RW BI ET \
    --lanl-dir data/raw/lanl \
    --genbank-dir data/raw/genbank \
    --out data/interim/east_africa_pol.parquet
```

Output: Parquet with `schema_version` and `merged_at` columns for downstream pinning.

## Schema versioning

The output schema is versioned by `SCHEMA_VERSION` in `schema.py`. Bump per the table in `CLAUDE.md` ("Data and schemas"). Every Parquet output stamps the version; downstream code should read the stamp before assuming column presence.

## Tests

```bash
pytest tests/ingest -v
```

Golden-file integration tests (frozen GenBank batches under `tests/fixtures/`) go in `tests/ingest/test_integration.py`; add after the first real pull whose behavior you want to freeze.
