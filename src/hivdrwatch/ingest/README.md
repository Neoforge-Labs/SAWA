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

- `subtype` — LANL curation wins when present; the GenBank `/note=` field is the fallback. Both go through `schema.subtype_from_designation`, so the two sources cannot disagree about the same string. Two or more named parents, or a mosaic marker, resolve to `recombinant`; a named CRF resolves to itself. `A` and `F` are parent clades for records that do not say A1/A2 or F1/F2 — ingest does not promote `A` to `A1`.
- `country` — **the country the record itself states, never the search query's.** GenBank: the `source` feature's `/geo_loc_name=` (NCBI's replacement for `/country=`), mapped via `schema.COUNTRY_NAME_TO_CODE`. LANL: the `Country` column, folder name as fallback. A record that states no country gets `country = None`; there is no fallback, because the country search over-matches and returns records from outside the requested country. `country_raw` keeps the stated value verbatim including any subdivision, and `country_source` says where `country` came from.
- `gene_region` — derived from the sequence's extent: `"pol"` up to `MAX_POL_LENGTH` (4,000 nt), `"genome"` above it. A pol subregion is still pol; a gag→nef deposit is not.
- One row per accession. `merge_sources` collapses the copies a record picks up from being returned by several country searches.
- Duplicate accessions across LANL files — later mtime wins.
- Sequences shorter than `--min-length` (default 800 bp) are dropped; phylogenetic signal is too weak below that.

`docs/evaluations/ingest-provenance.md` records what the pre-0.2.0 contract did to the shipped corpus, and `scripts/ingest_provenance_probe.py` re-measures it.

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

The output schema is versioned by `SCHEMA_VERSION` in `schema.py` (currently **0.2.0**). Bump per the table in `CLAUDE.md` ("Data and schemas"); the migration notes live in the `schema.py` module docstring and `CHANGELOG.md`. Every Parquet output stamps the version; downstream code should read the stamp before assuming column presence.

0.1.0 vintages **cannot be migrated in place** — the correct `subtype` and `country` are only recoverable from the source records, so re-run ingest against `data/raw/genbank` (the flat-file cache means this needs no network).

## Tests

```bash
pytest tests/ingest -v
```

Golden-file integration tests (frozen GenBank batches under `tests/fixtures/`) go in `tests/ingest/test_integration.py`; add after the first real pull whose behavior you want to freeze.
