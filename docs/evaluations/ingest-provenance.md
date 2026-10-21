# Three ingest provenance bugs, and what they did to the corpus

**Scope.** `src/hivdrwatch/ingest/` only. Nothing here is about resistance
scoring or alignment. The three bugs were found incidentally while
investigating alignment coverage (`align-coverage.md` §7, branch
`work/cool-burnell-e07285`); the alignment fix is independent of this one.

**One-line summary.** The two metadata fields that Paper 1 figure 2 and
Paper 3's benchmark split on — `subtype` and `country` — were both being
written from the wrong source. 10,069 of the 42,106 shipped curated rows
carry a subtype the records do not support, 8,702 carry a country the
records contradict, and 5,596 rows are duplicate copies of records
already in the table.

> **Kelyn's call before anything is regenerated.** Subtype labels feed two
> papers. The relabelling below moves 1,952 rows out of subtype C, creates
> a `recombinant` stratum of 5,482 rows where there were 186, and
> introduces two new enum members (`A`, `F`). Review the mapping before
> the vintage is rebuilt.

---

## 1. Method

- Source metadata re-derived from the local GenBank flat-file cache
  (`data/raw/genbank`, 1,450 batches, 85,917 cache entries, **36,510
  unique accessions**) — all of it, no sampling.
- Every number is measured against the shipped
  `data/interim/curated/aligned.parquet` (42,106 rows), joined to the
  cache by accession. **42,106 / 42,106 rows join.**
- The pre-fix implementations are reproduced verbatim in
  `scripts/ingest_provenance_probe.py` (`legacy_normalize_subtype`,
  `legacy_normalize_country`), so "the shipped value is the buggy one" is
  a measurement rather than an assumption. For **all 42,106 rows** the
  shipped `subtype` equals what the buggy parser produces from that
  record's own note — the curated vintage was labelled entirely by the
  GenBank path, with no LANL curation overriding it.

Reproduce:

```bash
python scripts/ingest_provenance_probe.py scan-cache \
    --genbank-dir data/raw/genbank --out /tmp/gb_meta.jsonl
python scripts/ingest_provenance_probe.py report \
    --cache /tmp/gb_meta.jsonl \
    --curated data/interim/curated/aligned.parquet \
    --fm-dir data/gold/fm --out /tmp/provenance.json
```

---

## 2. Bug 1 — every CRF was stored as subtype C

`_normalize_subtype` tested `f"SUBTYPE: {token}" in upper` for each enum
member in declaration order. `Subtype` declares `C` before the CRF
members, and `"SUBTYPE: C"` is a substring of `"SUBTYPE: CRF01_AE"`, so
the loop returned C for every CRF. The comment said "Order matters: check
CRFs before bare letters" — the CRF check sat *after* the loop, where it
could never be reached for a note that named a CRF.

Substring matching failed in three further ways: a bare letter with no
enum member (`subtype: A`) fell through to `unknown`; a mosaic that
*began* with a known subtype was read as that subtype (`subtype: A1/D
recombinant` → `A1`, `subtype: B/CRF01_AE` → `B`); and any spelling other
than `subtype: ` / `subtype ` was invisible (`subtype:A1`, `genotype: D`,
`clade A1`).

### What changed

Matching is now token-based (`schema.subtype_from_designation`, shared
with the LANL ingester so the two cannot disagree about the same string).
The declaration is isolated first, then scanned left to right, stopping at
the first token that is not part of a designation — which is what keeps
trailing prose out of the answer. Two or more named parents, or any mosaic
marker, resolve to `recombinant`; a named CRF resolves to itself.

The free-text forms tested are the ones that actually occur: all 2,511
distinct `/note=` values containing "subtype" were surveyed across the
cache, and the ~60 representative forms are parametrized in
`tests/ingest/test_provenance.py`. 24 of 25 of the headline cases were
wrong before the fix.

### Corpus effect

**10,069 rows (7,913 accessions) are relabelled.**

| shipped | corrected | rows |
|---|---|---|
| A1 | recombinant | 3,407 |
| unknown | A | 1,981 |
| C | **CRF01_AE** | 1,890 |
| unknown | recombinant | 1,424 |
| unknown | D | 512 |
| C | recombinant | 243 |
| unknown | A1 | 149 |
| B | recombinant | 138 |
| unknown | C | 101 |
| D | recombinant | 51 |
| unknown | F | 48 |
| A2 | recombinant | 38 |
| C | CRF10_CD | 33 |
| C | CRF02_AG | 23 |
| C | CRF16_A2D / CRF35_AD / CRF21_A2D | 6 |
| recombinant | CRF10_CD | 6 |
| C | unknown | 8 |
| other (≤ 4 rows each) | | 11 |

**1,952 rows shipped as subtype C are CRFs.** Subtype C is 7,109 shipped
rows, so **27.5% of the C stratum is not subtype C.**

Distribution over the 42,106 shipped rows:

| subtype | shipped | corrected |
|---|---|---|
| unknown | 19,540 | 15,326 |
| A1 | 10,749 | 7,491 |
| C | 7,109 | 5,007 |
| D | 3,459 | 3,918 |
| B | 932 | 796 |
| recombinant | 186 | **5,482** |
| A | — | 1,982 |
| CRF01_AE | — | 1,890 |
| F | — | 48 |
| G / A2 / J / H | 131 | 97 |
| CRF10_CD / CRF02_AG / CRF16_A2D / CRF35_AD / CRF21_A2D / F1 | — | 69 |

### Judgement call: bare `A` and bare `F`

1,981 rows say only `subtype: A` and 48 say only `subtype: F`. The enum
had A1/A2 and F1/F2 but no parent clade, so all of them were stored as
`unknown`.

**Decision: add `A` and `F` to the `Subtype` enum rather than promote A →
A1 or keep them as `unknown`.**

- Promoting A → A1 fabricates a distinction the depositor did not make.
  A1 vs A2 is a real clade split, and Paper 3's `subtype_holdout` split
  trains on A1 — a silently promoted A would put ~2,000 sequences of
  unasserted sub-clade into the training stratum.
- Keeping them `unknown` throws away information that is certainly
  present: these sequences are definitely not C, D or B.
- A parent-clade member says exactly what the record says. Rolling A1/A2/A
  up to "A" for a figure is a downstream decision, and ingest no longer
  has to make it.

`subtype: A (A1)` resolves to A1, not to `recombinant`: a parent clade
alongside one of its own sub-clades is the depositor being more specific,
not a list of mosaic parents.

### Residual: 9 malformed notes

Nine accessions have a note the fix cannot parse and the buggy parser got
wrong anyway: eight say `subtype: CFR21,A2` (a transposition of CRF21) and
one says `subtype: D10_CDB`. They shipped as C and D; they are now
`unknown`. That is a strict improvement — both former labels were wrong —
and typo tolerance is not worth encoding in the parser.

---

## 3. Bug 2 — the country was the search query's, not the record's

`_normalize_country` mapped only the eight `EAST_AFRICA_COUNTRIES` names
and returned `None` for anything else, and `_parse` then did
`country = _normalize_country(...) or requested_country`. A record whose
source feature said Thailand was stored as Kenyan, Ugandan or Tanzanian
depending on which country search returned it.

**It was worse than that.** NCBI renamed the source qualifier `/country=`
to `/geo_loc_name=`, and `/geo_loc_name=` is the **only** form present in
the cache:

| source qualifier | cache entries |
|---|---|
| `/country=` | **0** |
| `/geo_loc_name=` | 85,611 |
| neither | 306 |

`_extract_source_features` only ever looked for `country`, so
`_normalize_country` was called with `None` for every record in the
corpus and the fallback fired **every time**. Every country value in the
shipped table is query provenance. It happens to be right for most rows
only because the search usually returns records from the country asked
for.

### What changed

- `_stated_country` reads `/geo_loc_name=` first, then legacy `/country=`.
- `_normalize_country` maps **every** country, not just East African ones,
  via a shared `COUNTRY_NAME_TO_CODE` in `ingest/schema.py` (the LANL
  ingester's smaller private map now delegates to it, so the two sources
  cannot disagree). It covers every name present in the corpus — all 30
  map, none is unmapped — plus the sub-Saharan neighbours a future pull
  will reach.
- **There is no fallback to the requested country.**
- Two new optional fields carry the provenance: `country_raw` (verbatim,
  subdivision kept — `"Kenya: Nairobi"`, which the old code discarded
  entirely) and `country_source` (`"genbank_source_feature"`,
  `"lanl_curated"`, or `None`).
- A record whose stated country is not the requested one is logged at INFO
  and an unmappable name at WARNING, so the search's over-matching is
  visible instead of silent.

### Judgement call: drop the record, or keep it with its true country?

**Decision: keep it, with its own country and an explicit
`country_source`. Drop nothing, and drop the query fallback entirely —
including for records with no country field at all.**

- A record stating Thailand is not bad data; it is correct data about
  Thailand. Ingest's job is faithful provenance. Cohort definition belongs
  to the analysis layer, which already filters on country — so a correctly
  labelled Thai record falls out of every East-Africa-stratified analysis
  by itself, at no cost, while a dropped one destroys the audit trail that
  shows the search over-matching.
- Dropping at ingest would also be irreversible from the Parquet: the next
  person asking "why does the cohort not match the NCBI hit count?" would
  have nothing to look at.
- The fallback is **not** retained for the 277 records with no country
  field. The 3,108 non-East-African records are proof that the query's
  country field is not a statement about the record, so inheriting it for
  the silent 277 would fabricate geography with no better justification
  than for the Thai ones. They get `country = None` and are simply absent
  from country-stratified analysis.

### Corpus effect

**8,702 rows (5,516 accessions) carry a country the record contradicts.**
**4,600 of them state a country outside East Africa.**

| shipped | stated | rows |
|---|---|---|
| UG | KE | 1,105 |
| KE | UG | 847 |
| TZ | KE | 838 |
| TZ | **TH** | 767 |
| UG | **TH** | 702 |
| KE | **TH** | 681 |
| TZ | UG | 522 |
| UG | TZ | 348 |
| KE | TZ | 334 |
| UG | **MW** | 316 |
| TZ | **BR** | 300 |
| UG | **ZA** | 188 |
| UG | **ZW** | 183 |
| UG | **IN** | 181 |
| … 16 further destinations | | 1,190 |

Note the first three rows: the wrong country is often another East
African country, which is why this never looked obviously broken.

Per accession, the corrected geography of the corpus is:

| | accessions |
|---|---|
| in `EAST_AFRICA_COUNTRIES` | 33,125 |
| stated, elsewhere (TH 788, ZA 441, BR 368, MW 343, IN 183, ZW 183, ZM 154, IL 136, + 22 more) | 3,108 |
| no country stated | 277 |

Burundi is the sharpest case: **247 shipped rows claim BI; exactly 3
accessions in the whole cache state Burundi.**

### Why the search over-matches

`_build_query` filters on `"<name>"[Country]`. Whatever that field
resolves to at E-utilities, it returns records whose source feature says
Thailand, Brazil and India. The RV217 records are the clearest example —
a multi-site study (Kenya, Uganda, Tanzania **and** Thailand) whose Thai
deposits match an East African country search. Narrowing the query is a
separate change and is **not** made here; the fix makes the over-matching
visible and harmless instead of silent and corrupting.

---

## 4. Bug 3 — one row per country search, and six patient groups

`fetch_country` runs one search per country and `pipeline/ingest.py`
appends every country's records to one list; `merge_sources` wrote them
all. A record returned by three searches became three rows, and because
bug 2 gave each copy a different country, nothing downstream could see
they were the same record.

| | |
|---|---|
| curated rows | 42,106 |
| unique accessions | 36,510 |
| unique `sequence_sha256` | 35,616 |
| **surplus rows** | **5,596** |
| accessions appearing 2× | 716 |
| accessions appearing 3× | 2,440 |
| accessions where the copies share a `patient_group_id` | **0** |

Every country-stratified count in the shipped corpus double- or
triple-counts those 3,156 accessions.

In the cache the multiplicity is higher still — 29,282 accessions were
returned by two or more country searches (up to 8) — so the 5,596 surplus
rows are what survived the length filter and schema validation, not the
ceiling.

### A worse finding in the same place

`patient_group_id` has **6 distinct values across 42,106 rows** — one per
country — and `is_earliest_per_patient` is `True` for exactly **6 rows**.

`assign_patient_groups` keys on `(country, cohort, patient_token)` and
takes the token from a `sequence_name` column. **GenBank rows have no
`sequence_name` at all**, so `extract_patient_token(None)` returned
`("unknown", "no_name")` for every row and the key collapsed to
`(country, "unknown", "no_name")`. The module docstring promises the
opposite — "Sequences that can't be tied to a patient get a per-accession
group_id … Better false-independent than false-shared" — so the code
contradicted its own stated contract.

Nothing in production reads `is_earliest_per_patient` yet (only
`tests/align/test_curate.py` and the schema), so no published figure is
affected. Any analysis that had started using it would have seen six
sequences.

### What changed, and the judgement call on `sequence_sha256`

1. `merge_sources` collapses to one row per accession
   (`_collapse_duplicate_accessions`, `keep="first"`, order-stable) and
   logs how many it dropped. An accession is the record's identity; three
   copies are not three observations.
2. `assign_patient_groups` uses the accession as the patient token when a
   row has no `sequence_name`, honouring the docstring. Named
   longitudinal sequences (`KER2003_1998`, `KER2003_2001`) still share a
   group.

**Decision: do not collapse patient groups on `sequence_sha256`.**
Identical pol sequences from different accessions are not evidence of the
same patient — that is exactly the signal `cluster/` is built to detect,
and a short PR-only fragment can be identical between unrelated patients.
Merging on sequence identity would delete the epidemiology that Paper 1
figure 3 reports. Deduplicating on **accession**, upstream at merge, is
the correct and sufficient fix: with one row per accession, the `country`
component of the patient key can no longer split a record from itself.

`country` is deliberately **left in** `PatientKey`. It is the module's
documented conservatism ("better false-independent than false-shared"),
and with the accession fix in place it no longer causes the collapse.

---

## 5. `gene_region` was hardcoded to "pol"

`_parse` set `gene_region="pol"` unconditionally, so a gag→nef deposit
entered the pipeline as 8.7 kb labelled "pol".

| | accessions | shipped rows |
|---|---|---|
| > 4,000 nt | 5,905 | 9,087 |
| > 8,000 nt | 3,559 | — |

**Fixed, minimally:** `_infer_gene_region(length)` returns `"genome"`
above `MAX_POL_LENGTH = 4000` and `"pol"` otherwise. HXB2 pol is 3,012 nt;
the 3,012–4,000 band (78 accessions) is pol plus flanking sequence, and a
pol *subregion* is still pol. Safe to change: nothing in the codebase
filters on `gene_region` — it is carried to Parquet, loaded to Postgres as
`String(8)`, and reported by the API, nothing more.

**Not done, deliberately:**

- **Deriving the label from the record's gene/CDS features.** More correct
  in principle, but GenBank gene annotation is inconsistent across
  depositors — which is why `_build_query` has to OR across gene-name
  variants in the first place — and `/product=` is free text. Length holds
  for every record; feature-derived labels would hold for some. Worth
  doing as its own change, with a vocabulary (`pr`, `rt`, `in`, `pr-rt`,
  …) chosen against the 8-character Postgres column.
- **Trimming whole-genome deposits to pol at ingest.** This would be a
  regression: `pipeline/dr_ingest.py` sends Sierra the **raw** sequence
  from this record, and Sierra handles whole genomes correctly. Trimming
  at ingest would narrow what the resistance caller sees in order to fix a
  labelling problem. The aligner already restricts the sequence to pol for
  `phylo/` and `cluster/`, which is the right place for it.

---

## 6. Schema version

`SCHEMA_VERSION` **0.1.0 → 0.2.0** (minor), per the CLAUDE.md table:
`subtype` and `country` change semantics, and `country_raw` /
`country_source` are new optional fields with defaults. The migration note
is in the `ingest/schema.py` module docstring.

A 0.1.0 vintage **cannot be migrated in place** — the correct subtype and
country are only recoverable from the source records — so the fix is a
re-run of ingest against `data/raw/genbank`, not a column transform.

The bump protocol also asks for a regenerated `data/interim/` sample in
DVC. That is a pipeline run and is **not** done here; see §7.

---

## 7. Downstream artifacts that are now stale

Nothing below has been regenerated. Every one of these was built from the
0.1.0 vintage.

| artifact | why it is stale |
|---|---|
| `data/interim/east_africa_pol.parquet` | the ingest output itself. Also note the copy on disk is LANL-only (32,115 rows, all `source=lanl`) and is **not** the vintage `curated/aligned.parquet` was built from. |
| `data/interim/curated/aligned.parquet` | 42,106 rows: 5,596 duplicates, 10,069 wrong subtypes, 8,702 wrong countries, 6 patient groups. Also stale for the independent alignment fix on `work/cool-burnell-e07285`. |
| `data/interim/trees/{A1,C,D}` + `data/gold/trees/partition_report.parquet` | per-subtype trees are partitioned on the field that changed. The C tree drew from 1,077 post-QC C sequences, of which a share are CRFs (including Thai CRF01_AE); `recombinant` had 25 available sequences and should have had thousands. |
| **Paper 1 figure 2** (resistance heatmap, country × year × subtype) | both axes moved. Of 394 populated country × year × subtype cells, **43 disappear and 324 are new** (675 cells after the fix); only 351 are unchanged. Every cell count also carries the duplicate rows. |
| **Paper 1 figure 3** (cluster parity) + `data/gold/cluster_assignments.parquet`, `data/gold/outbreak_features.parquet` | duplicate rows entered the distance matrix as independent sequences, and `cluster/features.py` aggregates on `country` and `subtype`. (Figure 3 is also gated on the distance-denominator question in `align-coverage.md` §8 — that one is still Kelyn's call.) |
| **Paper 3 `subtype_holdout` splits** (`data/gold/fm/*/subtype_holdout/*.parquet`) | the split *is* the subtype label. The INSTI train stratum is A1 7,064 / C 4,626 / D 1,980 — the C arm is contaminated with CRFs — and the test stratum is 15,345 `unknown` rows, of which 1,981 are now `A` and 1,424 `recombinant`. The stratified splits are affected too, less severely. |
| `data/gold/fm/comparison_*.{parquet,json}` + `data/gold/fm/checkpoints/` | fitted on the above. |
| `data/interim/sierra_subtypes.parquet`, `data/interim/dr_results.parquet` | **not stale for resistance calls** — `dr_ingest` sends Sierra the raw sequence, and the duplicate rows hash to the same `sequence_sha256`, so the cache collapses them. Stale only in that they carry per-row country/subtype columns. |
| **Paper 1 figure 4** (DTG case study) | stratified by country and year. Affected by the country fix and the duplicate rows, not by the alignment fix. |
| DVC pointers for all of the above | still addressing the 0.1.0 outputs. |

**Unaffected:** Paper 2 (`causal/`) has no data yet — `causal/data_prep.py`
is still `·` to build. `benchmarks/results/*.json` measure throughput, not
labels.

---

## 8. What is left

1. **Kelyn reviews the subtype relabelling** (§2) before the vintage is
   rebuilt. The `A`/`F` enum members and the 5,482-row `recombinant`
   stratum are the two decisions worth a second opinion.
2. **Re-run ingest, then curate**, and refresh DVC. This also needs the
   alignment fix from `work/cool-burnell-e07285` merged first, or the
   vintage will be rebuilt twice.
3. **Narrow the GenBank country search** (§3). The fix makes over-matching
   visible; it does not stop 3,108 off-cohort records being fetched and
   parsed on every run.
4. **CI lint is already red on `main`.** `ruff check .` has three
   pre-existing auto-fixable findings (`notebooks/eda_v0.2.0.py` F541 × 2,
   `src/hivdrwatch/phylo/build.py` UP017), and `ruff format --check .`
   would reformat 10 pre-existing files (`foundation_models/` × 5,
   `notebooks/` × 2, `src/hivdrwatch/{phylo/build,pipeline/cluster}.py`,
   `tests/phylo/test_phylo.py`). Untouched here per the "don't reformat
   files you didn't change" rule, but `.github/workflows/test.yml` runs
   both and is failing on them.
5. **mypy is not actually run in CI**, contrary to CLAUDE.md and
   MODULES.md: `.github/workflows/test.yml` runs `ruff check` and
   `ruff format --check` only. `mypy src/hivdrwatch/{ingest,align}` has a
   pre-existing backlog of 22 errors (23 before this change), mostly
   missing Biopython stubs.
6. **Feature-derived `gene_region`** (§5), as its own change.
