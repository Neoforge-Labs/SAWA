# `dashboard` — static resistance summary

A single static page that shows what Stanford HIVdb Sierra made of one data vintage, with the denominator stated beside every figure. An export step reduces `data/interim/dr_results.parquet` to one committed file, `static/summary.json`; the page (`static/index.html`, `app.js`, `styles.css`, no build step, no dependencies) draws only from that file and makes no live calls. CI publishes `static/` to GitHub Pages.

This is not the API-backed Auspice + deck.gl dashboard that `MODULES.md` Flow 1 specifies. That one is not started. This module is the resistance summary only.

## What it computes

`export/` turns the parquet into four views, written under `views` in `summary.json`:

- `prevalence_by_drug_class` — sequences scored for a class, and how many score intermediate or high-level resistance to any drug in it.
- `levels_by_drug` — every drug across Sierra's five levels, zero-filled so a missing level is a zero and not an absent bar.
- `gene_coverage` — distinct accessions with a Sierra score for each gene (PR, RT, IN). This sets the denominators above.
- `sdrm_frequency` — surveillance drug-resistance mutations (SDRMs) per panel (PR, RT NRTI list, RT NNRTI list), counted per sequence.

Before any aggregation, `dedup.py` drops exact duplicate rows (the parquet holds up to three copies of each call) and refuses to proceed if that is not lossless: it raises if rows sharing an (accession, drug) key differ in any column, or if one accession carries more than one `sequence_sha256`. `meta` additionally requires a single `algorithm_version`, `algorithm_publish_date` and `called_at` across the frame.

`summary.json` is byte-deterministic: aggregations end in explicit sorts, keys are serialised sorted, floats are rounded to four places. Two runs over the same input produce the same file, which is what lets it be committed and reviewed in a diff.

## Run it

```bash
make dashboard
# or
python -m hivdrwatch.dashboard.export --ingested-count 36510
hivdrwatch dashboard --ingested-count 36510
```

`--ingested-count` is required and has no default. The parquet holds only accessions that returned a Sierra result, so the number that entered ingest cannot be derived from it, and a default would silently go stale when the vintage is regenerated. 36510 is the current vintage's total; `make dashboard INGESTED_COUNT=<n>` overrides it.

Other flags: `--input` (default `data/interim/dr_results.parquet`), `--pointer` (default the input path plus `.dvc`), `--reference` (default `data/reference/drm_positions.json`), `--out` (default `static/summary.json`).

Preview the page locally with any static file server, for example `python -m http.server -d src/hivdrwatch/dashboard/static`; opening `index.html` from disk usually fails because browsers block `fetch` over `file://`.

## Inputs and outputs

| | Path | Notes |
|---|---|---|
| In | `data/interim/dr_results.parquet` | Sierra results, DVC-tracked. Not in git. |
| In | `data/interim/dr_results.parquet.dvc` | The pointer; its `md5` is recorded in `meta.input_dvc_md5`. |
| In | `data/reference/drm_positions.json` | Surveillance positions and resistant states (`SDRM-2019`). |
| Out | `static/summary.json` | `meta` (provenance and counts) and `views` (the four lists above). Schema 1.1.0. |

`meta.schema_version` is `SUMMARY_SCHEMA_VERSION` in `export/summary.py`. The page declares the version it was built for in `<body data-summary-schema>` and refuses to draw charts from a file that declares another, showing a banner instead. See `CHANGELOG.md` for the 1.1.0 compatibility break.

## What is deliberately not here

**No geographic or subtype axis.** In this data vintage the country and subtype fields were written from the wrong source. Every country value reflects the search query rather than the record, and the same evaluation counts 3,108 accessions that state a country outside East Africa; the subtype matcher mislabelled CRFs as subtype C. Any country- or subtype-stratified figure would present the search design as a finding. The page says so and links `docs/evaluations/ingest-provenance.md`. These views return after a provenance rebuild (ingest `SCHEMA_VERSION` 0.2.0 fixes the cause; the shipped vintage has not been regenerated).

**Surveillance mutations cover PR and RT only, not integrase.** The vendored list predates the dolutegravir era. Bennett et al. 2009, the cited source, has no integrase component, and the `IN` block in `data/reference/drm_positions.json` is the raltegravir/elvitegravir-era list: it covers positions 66, 92, 121, 140, 143, 147, 148 and 155, and is missing R263K, G118R and E138K. Showing it would attribute to Bennett a list it does not contain. Integrase appears in the prevalence and level views, which rest on Sierra's own scoring and not on this list, with a thinner denominator (9,365 of 35,777 scored accessions have integrase coverage in the current vintage).

**Mixtures are attributed to one state.** A mixture counts as carrying an SDRM when any of its states is on the list. Where it matches more than one state, `match_sdrm` attributes it to the first in the panel's sorted order. `meta.n_sdrm_observations_multistate` of `meta.n_sdrm_observations` measures how often that happens (427 of 27,391 in the current vintage), and the page states the figure from those fields.

## What the drift guard guarantees

The `summary-input-drift` job in `.github/workflows/test.yml` runs `check_input_not_stale`, which compares `meta.input_dvc_md5` in the committed `summary.json` against the `md5` in the committed `.dvc` pointer. That is one pointer value against another.

- It **detects** `dvc add` moving the pointer after an export was made: the summary then records a stale hash.
- It does **not** verify the aggregation. CI never re-derives a count.
- It does **not** detect a parquet that disagrees with its own pointer. The only DVC remote is a MinIO on `localhost` (`.dvc/config`), which no hosted runner can reach, so CI never holds the data.
- `meta.input_file_md5` is the digest of the file actually hashed at export time. It cannot be checked in CI at all, and the guard does not consult it.

On the current vintage `input_dvc_md5` and `input_file_md5` differ: the file the summary was computed from does not hash to the md5 its pointer records. The page footer shows both digests and discloses the mismatch.

The deploy is not gated on that job either: `pages.yml` triggers independently of `test.yml`, so a stale `summary.json` that reaches `main` is published unless branch protection requires the test workflow.

## Tests

```bash
pytest tests/dashboard -v
pytest tests/cli -v
```

`tests/dashboard/` covers deduplication invariants, SDRM parsing and matching, each aggregate, summary determinism and provenance, the drift check, and the static page (bindings resolve against the real `summary.json`, the schema gate, the digest note). The CLI registration is picked up by the contract tests in `tests/cli/test_app.py`.
