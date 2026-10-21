# Changelog

Schema-version bumps and anything that changes the meaning of a shipped
column. Format loosely follows [Keep a Changelog]; versions here are the
per-stage `*_SCHEMA_VERSION` constants, not a package version.

[Keep a Changelog]: https://keepachangelog.com/en/1.1.0/

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
