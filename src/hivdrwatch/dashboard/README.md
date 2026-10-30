# `dashboard` — static resistance summary

A single static page that shows what Stanford HIVdb Sierra made of one data vintage, with the denominator stated beside every figure. An export step reduces `data/interim/dr_results.parquet` to one committed file, `static/summary.json`; the page (`static/index.html`, `app.js`, `styles.css`, no build step, no dependencies) draws only from that file and makes no live calls. The `pages.yml` workflow publishes `static/` to GitHub Pages once Pages is enabled on the repository, which it is not yet.

This is not the API-backed Auspice + deck.gl dashboard that `MODULES.md` Flow 1 specifies. That one is not started. This module is the resistance summary only.

## What it computes

`export/` turns the parquet into four views, written under `views` in `summary.json`:

- `prevalence_by_drug_class` — sequences scored for a class, and how many score intermediate or high-level resistance to any drug in it.
- `levels_by_drug` — every drug across Sierra's five levels, zero-filled so a missing level is a zero and not an absent bar.
- `gene_coverage` — distinct accessions with a Sierra score for each gene (PR, RT, IN). This sets the denominators above.
- `sdrm_frequency` — surveillance drug-resistance mutations (SDRMs) per panel (PR, RT NRTI list, RT NNRTI list), counted per sequence.

Before any aggregation, `dedup.py` drops exact duplicate rows (the parquet holds up to three copies of each call) and refuses to proceed if that is not lossless: it raises if rows sharing an (accession, drug) key differ in any column, or if one accession carries more than one `sequence_sha256`. It also reports `n_distinct_sequences`, because the converse does not hold: several accessions can carry one sequence. `meta` additionally requires a single `algorithm_version`, `algorithm_publish_date`, `called_at` and DR `schema_version` across the frame (the publish date alone may be null), and an ingested count no smaller than the scored count.

`summary.json` is byte-deterministic for a given input at a given commit: aggregations end in explicit sorts, keys are serialised sorted, floats are rounded to four places. Two runs over the same input **at the same commit** produce the same file, which is what lets it be committed and reviewed in a diff. Across commits the file differs by `meta.git_sha`, which is the commit that produced it; a `-dirty` suffix means the working tree had uncommitted changes (or could not be read) when it was built, and a summary carrying one is not reproducible from the commit it names. A test fails if the committed file does.

## Run it

The input is DVC-tracked and is not in git, so pull it first (`dvc pull data/interim/dr_results.parquet`); the commands below fail without it.

```bash
make dashboard INGESTED_COUNT=36510
# or
python -m hivdrwatch.dashboard.export --ingested-count 36510
hivdrwatch dashboard --ingested-count 36510
```

`--ingested-count` is required and has no default, here or in `make dashboard` (as `INGESTED_COUNT`). The parquet holds only accessions that returned a Sierra result, so the number that entered ingest cannot be derived from it, and a default would silently go stale when the vintage is regenerated. 36510 is the current vintage's total. Commit or discard other changes first, or `meta.git_sha` is recorded `-dirty`.

Other flags: `--input` (default `data/interim/dr_results.parquet`), `--pointer` (default the input path plus `.dvc`), `--reference` (default `data/reference/drm_positions.json`), `--out` (default `static/summary.json`).

Preview the page locally with any static file server, for example `python -m http.server -d src/hivdrwatch/dashboard/static`; opening `index.html` from disk usually fails because browsers block `fetch` over `file://`.

## Inputs and outputs

| | Path | Notes |
|---|---|---|
| In | `data/interim/dr_results.parquet` | Sierra results, DVC-tracked. Not in git. |
| In | `data/interim/dr_results.parquet.dvc` | The pointer; its `md5` is recorded in `meta.input_dvc_md5`. |
| In | `data/reference/drm_positions.json` | Surveillance positions and resistant states (`SDRM-2019`). |
| Out | `static/summary.json` | `meta` (provenance and counts) and `views` (the four lists above). Schema 1.2.0. |

`meta.schema_version` is `SUMMARY_SCHEMA_VERSION` in `export/summary.py`. The page declares the version it was built for in `<body data-summary-schema>` and refuses to draw charts from a file that declares another, showing a banner instead. See `CHANGELOG.md` for the 1.1.0 and 1.2.0 compatibility breaks. `meta.dr_schema_version` is read from the data frame, not from the code constant.

## What is deliberately not here

**No geographic or subtype axis.** In this data vintage the country and subtype fields are both wrong, for different reasons. Country was written from the search query rather than the record, so every country value reflects the query, and the evaluation counts 3,108 accessions that state a country outside East Africa. Subtype was misparsed: the matcher compared substrings, so every CRF was stored as subtype C. Any country- or subtype-stratified figure would present the search design as a finding. The page says so and links `docs/evaluations/ingest-provenance.md`. These views return after a provenance rebuild (ingest `SCHEMA_VERSION` 0.2.0 fixes the cause; the shipped vintage has not been regenerated).

**Surveillance mutations cover PR and RT only, not integrase.** The vendored list predates the dolutegravir era. It is the Stanford/WHO SDRM list (`SDRM-2019`, cited in the data file to the WHO resistance-list page), whose PR and RT content follows Bennett et al. 2009. That paper has no integrase component, and the `IN` block in `data/reference/drm_positions.json` is the raltegravir/elvitegravir-era list: it covers positions 66, 92, 121, 140, 143, 147, 148 and 155, and is missing R263K, G118R and E138K. Showing it would attribute to Bennett a list it does not contain. Integrase appears in the prevalence and level views, which rest on Sierra's own scoring and not on this list, with a thinner denominator (9,365 of 35,777 scored accessions have integrase coverage in the current vintage).

**Counts are of accessions, not patients.** The unit throughout is the accession. `meta.n_distinct_sequences` counts the distinct sequences among them (34,962 against 35,777 accessions in the current vintage), so 815 accessions repeat a sequence already present under another, clustered inside single submission series. They are not removed: `patient_group_id` is degenerate in this vintage (6 groups for 36,510 accessions), so de-duplicating by patient is not available, and de-duplicating by sequence would be a different denominator. The page states the unit and the repeat count.

**No treatment history.** The Sierra results carry no treatment-history field, so the prevalence figures are not pretreatment drug-resistance estimates: the corpus mixes surveillance, treatment-failure and trial sequences in unknown proportion. The page says so beside the figures.

**Mixtures are attributed to one state, the alphabetically first.** A mixture counts as carrying an SDRM when any of its states is on the list. Where it matches more than one state, `match_sdrm` attributes it to the alphabetically first matching state, so `M46IL` is counted as `M46I` and `V82AF` as `V82A`. The effect is directional, not noise: at a position with several surveillance states the earlier label is inflated at the later one's expense. `meta.n_sdrm_observations_multistate` of `meta.n_sdrm_observations` measures how often that happens (427 of 27,391 in the current vintage), and the page states the rule and the figure from those fields.

## What the drift guard guarantees

The `summary-input-drift` job in `.github/workflows/test.yml` runs `check_input_not_stale`, which compares `meta.input_dvc_md5` in the committed `summary.json` against the `md5` in the committed `.dvc` pointer. That is one pointer value against another. The same comparison also runs as a test in `tests/dashboard/test_summary.py`, so `pytest tests/ -q` catches drift before a push and not only after it.

- It **detects** `dvc add` moving the pointer after an export was made: the summary then records a stale hash.
- It does **not** verify the aggregation. CI never re-derives a count.
- It does **not** detect a parquet that disagrees with its own pointer. The only DVC remote is a MinIO on `localhost` (`.dvc/config`), which no hosted runner can reach, so CI never holds the data.
- `meta.input_file_md5` is the digest of the file actually hashed at export time. It cannot be checked in CI at all, and the guard does not consult it.

## What the digests do and do not establish

On the current vintage the DVC pointer, the DVC cache and the file the summary was computed from all agree: `input_dvc_md5` and `input_file_md5` are both `a8279a4d…`. That was not always so. The pointer had recorded `e51cbb39…`, which the cache still held, while the working copy was a re-run on 2026-08-17 that was never `dvc add`ed. The two frames are identical in every resistance call (`level`, `score`, `gene`, `algorithm_version`, `sequence_sha256`, over the same 728,087 accession-drug pairs); only `called_at` differs, 2026-08-16 in the old vintage and 2026-08-17 in the re-run. That was a provenance failure, not a data failure. The summary was re-exported once the pointer was updated.

Agreement is not a reproducibility guarantee, and the page's footer notice, which shows only when the two digests differ, stays in place to catch the next divergence.

- `summary.json` is byte-deterministic for a given input at a given commit (see above).
- The upstream ingest is **not** yet byte-reproducible. Re-running the pipeline that produced this parquet can emit byte-different output even when the calls are identical (row ordering is one known cause, addressed separately), so a matched digest says that this file is the one that was summarised, not that it can be regenerated byte for byte.

The deploy is not gated on that job either. `pages.yml` triggers independently of `test.yml`, and branch protection requiring the test workflow is not currently enabled on `main`, so a stale `summary.json` that reaches `main` is published even while the drift job is red. The guard is advisory today; enabling that protection is what would make it binding.

## Tests

```bash
pytest tests/dashboard -v
pytest tests/cli -v
```

`tests/dashboard/` covers deduplication invariants, SDRM parsing and matching, each aggregate, summary determinism and provenance, the drift check (including a check of the committed summary against the committed DVC pointer, so drift fails locally and not only in CI), `make dashboard`, and the static page (bindings resolve against the real `summary.json`, the schema gate, the digest note, null handling run under node when it is installed).

`tests/fixtures/dashboard/sdrm_panels_golden.json` pins the PR, RT NRTI and RT NNRTI positions and states exactly as `data/reference/drm_positions.json` holds them. The page cites those panels to Bennett et al. 2009 (updated 2019), while the data file's own `notes` call them candidate model features, and nothing has yet checked them against a published list. The golden makes any edit to a position or state a deliberate, reviewable change; it is intentionally not regenerated by `--update-golden`. The CLI registration is picked up by the contract tests in `tests/cli/test_app.py`.
