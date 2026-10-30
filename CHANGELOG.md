# Changelog

Schema-version bumps and anything that changes the meaning of a shipped
column. Format loosely follows [Keep a Changelog]; versions here are the
per-stage `*_SCHEMA_VERSION` constants, not a package version.

[Keep a Changelog]: https://keepachangelog.com/en/1.1.0/

## dashboard summary `SUMMARY_SCHEMA_VERSION` 1.2.0 — 2026-10-30

**The page's schema gate refuses any `summary.json` not at 1.2.0; a 1.1.0
file exists only in this branch's own history, never in a release.**

### Changed

- Summary schema 1.1.0 → 1.2.0 (minor: new required field). `meta` gains
  `n_distinct_sequences`, the number of distinct `sequence_sha256` values
  among the scored accessions. Every count the page publishes is of
  accessions, and several accessions carry one sequence, so the page needs
  this number to state the unit and how far it differs from distinct
  sequences. Any 1.1.0 file lacks it and must be re-exported.
- `DedupReport` gains `n_distinct_sequences`.

## dashboard summary `SUMMARY_SCHEMA_VERSION` 1.1.0 — 2026-10-30

**New module. The page's schema gate refuses any `summary.json` not at
1.1.0 instead of half-drawing it; a 1.0.0 file exists only in this branch's
own history, never in a release.**

### Added

- `src/hivdrwatch/dashboard/`: `export/` deduplicates Sierra results,
  matches PR and RT surveillance mutations and aggregates four views into
  a deterministic `summary.json`; a dependency-free static page draws it.
  Registered as `hivdrwatch dashboard` and `make dashboard`.
- `.github/workflows/pages.yml` publishes the static page to GitHub Pages
  once Pages is enabled on the repository (it is not yet), and the
  `summary-input-drift` job in `test.yml` fails when the DVC pointer has
  moved since the summary was exported. The guard compares
  `meta.input_dvc_md5` against the committed pointer only; it does not
  verify the aggregation or the parquet.

### Changed

- Summary schema 1.0.0 → 1.1.0 (minor: new required fields). `meta` gains
  `n_sdrm_observations` and `n_sdrm_observations_multistate`, so the page
  states the measured size of the single-state attribution rule rather
  than a typed-in figure, and `input_file_md5`, the digest of the file
  actually read, kept apart from `input_dvc_md5`, which is only what the
  DVC pointer claims. The page requires all three, so this is a
  compatibility break in the JSON contract and not a cosmetic version
  move. Any 1.0.0 file lacks them and must be re-exported.

### Not done

- Geographic and subtype-stratified views: the shipped vintage's country
  and subtype are unreliable (`docs/evaluations/ingest-provenance.md`),
  so they return after a provenance rebuild.
- Integrase surveillance mutations: the vendored SDRM list has no
  dolutegravir-era integrase block.
- Flow 1's API-backed Auspice + deck.gl dashboard (`MODULES.md`).

## ingest `SCHEMA_VERSION` 0.2.0 — 2026-10-21

**`subtype` and `country` change meaning. 0.1.0 and 0.2.0 records must not
be pooled, and a 0.1.0 vintage cannot be migrated in place — re-run ingest
against `data/raw/genbank`.**

### Fixed

- `genbank._normalize_subtype` matched the declared subtype by substring,
  so `"SUBTYPE: C"` matched inside `"SUBTYPE: CRF01_AE"` and every CRF was
  stored as subtype C. Matching is now token-based via the shared
  `schema.subtype_from_designation`. 10,069 of the 42,106 shipped curated
  rows are relabelled; 1,952 rows shipped as C are CRFs.
- `genbank._normalize_country` mapped only the eight East African country
  names and `_parse` fell back to the search query's country for
  everything else, so a record stating Thailand was stored as Kenyan,
  Ugandan or Tanzanian. Compounded by `_extract_source_features` reading
  only `/country=`, which NCBI replaced with `/geo_loc_name=` — not one
  record in the 36,510-accession cache carries `/country=`, so the
  fallback fired for every record in the corpus. 8,702 shipped rows carry
  a country the record contradicts; 4,600 state a country outside East
  Africa. There is no longer any fallback to the query.
- `merge_sources` wrote one row per country search per record, so an
  accession returned by three searches became three rows. 5,596 of the
  42,106 shipped rows are duplicates of 3,156 accessions, and every
  country-stratified count double- or triple-counted them.
- `align.dedup.assign_patient_groups` keyed on a `sequence_name` column
  that GenBank rows do not have, collapsing the whole corpus into 6
  patient groups with `is_earliest_per_patient` true for 6 rows. Unnamed
  sequences now get a per-accession group, as the module docstring already
  promised.
- `genbank._parse` labelled every record `gene_region="pol"`, so 5,905
  multi-gene deposits (3,559 of them over 8 kb) entered the pipeline as
  "pol". Now derived from the sequence's extent.

### Added

- `Subtype.A` and `Subtype.F` — parent clades, for the ~2,000 records
  whose note says only `subtype: A` or `subtype: F` and which were stored
  as `unknown`. Promoting A → A1 would fabricate a clade split the
  depositor did not make.
- `SequenceRecord.country_raw` — the country verbatim as the source states
  it, subdivision included (`"Kenya: Nairobi"`), which was previously
  discarded.
- `SequenceRecord.country_source` — `"genbank_source_feature"`,
  `"lanl_curated"`, or `None`. Never the search query.
- `schema.COUNTRY_NAME_TO_CODE` / `country_name_to_code` — one country map
  for both ingesters; `lanl._country_name_to_code` now delegates to it.
- `scripts/ingest_provenance_probe.py` — re-derives subtype and country
  for the whole corpus from the local GenBank cache and reports the delta.

### Changed

- `subtype_source` is `"genbank"` only when a subtype was actually
  resolved from the note; previously any note at all claimed provenance,
  including notes that resolved to `unknown`.
- `lanl._coerce_subtype` delegates to the shared classifier, so a LANL
  `A1/D` and a GenBank note saying `subtype: A1/D` cannot disagree. Side
  effect: LANL mosaic shorthands (`AD`, `AC`, `ACD`) now resolve to
  `recombinant` instead of `unknown`, and `A`/`F` to themselves.

### Not done

- The `data/interim/` sample regeneration the bump protocol asks for.
  Rebuilding the vintage is a pipeline run, and the subtype relabelling
  needs Kelyn's review first. Stale artifacts are listed in
  `docs/evaluations/ingest-provenance.md` §7.
