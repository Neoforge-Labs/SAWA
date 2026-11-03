# Design — resistance dashboard (clean axes)

**Date:** 2026-10-30 · **Owner:** Kelyn · **Status:** awaiting approval · §14 resolved 2026-10-30
**Target paper:** Paper 1 (Flow 1, sequencing item 4) · **Module:** `dashboard/`

---

## 1. Intent and the why

SAWA's resistance calls are the most finished artifact in the repo — 35,777
sequences scored against Stanford HIVdb — and they are currently reachable
only by running `make fig2` with a DVC pull in front of it. Kelyn wants a
**URL he can put on his portfolio**, so a reader who will never clone the
repo can see that the pipeline produces real resistance surveillance
output.

The audience is therefore a **portfolio visitor**: prospective PhD
supervisors (Weisse, Atkins, Leigh Brown), Dr. Chohan, and peers. They
skim. They will not install anything. Several of them are domain experts
who will spot a methodological error instantly, and the recruitment goal
in CONTEXT.md means a wrong number in front of this audience costs more
than a missing feature.

This is **not** the Flow 1 reviewer dashboard and **not** the KEMRI
deployment UI. It is the first slice of `dashboard/`, built so those can
grow from it rather than replace it.

## 2. Goals

- A static page at a stable public URL, linkable from a portfolio.
- Shows resistance findings that are **defensible against today's data**.
- Shows its own provenance: algorithm version, cohort size, data vintage,
  per-class denominators.
- Loads fast on a poor connection and on a phone (CONTEXT.md: the dashboard
  audience is African clinical and public-health researchers; even though
  this slice targets a portfolio reader, the constraint is cheap to honour
  now and expensive to retrofit).
- Grows into the stratified version (§5, Approach B) without a rewrite.

## 3. Non-goals

- No country, subtype, or year stratification. No map, no heatmap. (§4)
- No phylogenetic tree. CONTEXT.md commits SAWA to embedding Auspice, and
  `phylo/auspice.py` already exports Auspice v2 JSON; drawing trees here
  would duplicate it.
- No live API calls. The page does not depend on the FastAPI service being
  up, so a dead container can never break the portfolio link.
- No per-sequence record browser. Interesting to a reviewer, noise to a
  portfolio reader.
- No write path, no auth, no KEMRI or PHI code path.

## 4. Anti-goals — what would make this a failure even if it ships

- **Rendering the poisoned strata.** `docs/evaluations/ingest-provenance.md`
  documents that the shipped vintage's `subtype` and `country` were written
  from the wrong source: 10,069 rows carry an unsupported subtype, 8,702 a
  contradicted country, and 27.5% of the subtype-C stratum is not C. The
  doc explicitly reserves the relabelling decision for Kelyn. A map or a
  country × subtype heatmap on this vintage would publish known-wrong
  science to the exact audience being courted. **This is the single thing
  the design exists to prevent.**
- **Reporting a mutation chart of natural polymorphisms.** See §7.2.
- **Reporting inflated counts.** See §7.1.
- A page that looks like a product demo rather than research
  infrastructure. CONTEXT.md's voice note applies: calm, competent,
  infrastructure-flavoured, not an ML demo.

## 5. Approaches considered

**A — Clean axes now, static page. CHOSEN (Kelyn, 2026-10-30).**
Render only axes that come from Sierra's own output (drug, drug class,
level, score, gene, mutations, algorithm version), which key on
`accession + sequence_sha256` and never touch the poisoned ingest
metadata. Ships without any pipeline decision. Strict subset of B's page.

**B — Settle the subtype relabelling, re-ingest, then ship the full
landscape. DEFERRED, next.**
The visually stronger version (map + country × year × drug-class heatmap).
Compute is cheap: `data/raw/` holds the LANL dir, 1,453 GenBank batches
(85,917 cache entries) and the HIVdb Sierra cache, so a rebuild needs no
network and `mappy` aligns ~1,900 seq/s. The cost is Kelyn's review of a
mapping that moves 1,952 rows out of subtype C, grows `recombinant` from
186 to 5,482 rows, and adds two enum members (`A`, `F`) — labels that feed
Paper 1 Figure 2 and Paper 3's benchmark split. Sequenced after A so the
link exists while that decision gets unhurried attention.

**C — Ship stratified now, labelled "provisional". DROPPED.**
A caveat does not undo a wrong map in front of a prospective supervisor,
and the label would itself advertise that the pipeline's provenance was
not trusted. Rejected on audience grounds, not technical ones.

## 6. Constraints

- **Repo conventions (CLAUDE.md):** Python 3.12+, `from __future__ import
  annotations`, mypy strict, no `print()` in library code, no `# TODO`,
  specific exceptions, module docstring required, per-module README
  required, MODULES.md row added in the same PR, Conventional Commits,
  feature branch (`feat/dashboard/...`), never commit to `main`, coverage
  ≥ 70% on touched modules.
- **No raw data in git.** Aggregate JSON is not raw data and follows the
  existing precedent of committed `benchmarks/results/*.json`; nothing
  lands under `data/`.
- **Determinism.** Same input parquet → bit-identical JSON, so the page is
  reviewable in diffs.
- **No new heavy dependencies.** Polars is already in the tree and does the
  aggregation; the page should need no build toolchain (§9).
- The source parquet is DVC-tracked and absent in CI, so CI cannot
  regenerate the JSON. The committed JSON is the deploy input.

## 7. Two data findings that shape the design

### 7.1 `dr_results.parquet` contains exact duplicate rows

848,590 rows reduce to 728,087 fully distinct rows — **120,503 redundant
(14.2%)**, up to 3 identical copies of the same call. Verified: every one
of the 35,777 accessions has exactly one `sequence_sha256`, and
`distinct(accession, drug)` equals `distinct(accession, sha, drug)` equals
`n_unique()` over all columns, so the duplicates carry no information and
`.unique()` is lossless.

Consequence: any count taken from the raw parquet is inflated ~14%,
unevenly across levels. **The export deduplicates before aggregating, and
asserts the dedup is lossless rather than assuming it.** The headline call
count is 728,087, not 848,590.

This duplication is in the resistance output and is distinct from the 5,596
duplicate rows the provenance doc found in the curated table. Its root
cause is not diagnosed here and is out of scope for this work — it should
be raised as a GitHub issue against `resistance/` or `pipeline/dr_ingest`.

### 7.2 Sierra returns every variant, not just resistance mutations

The `mutations` column holds 16,357 distinct strings, and the most frequent
are natural polymorphisms, not DRMs: PR R41K (27,694 sequences), RT V35T
(27,439), PR M36I (26,641), RT K122E (25,275). A "top mutations" chart off
this column shows background variation and would read as a domain error to
a phylodynamics reviewer.

**The mutation view filters to the Stanford surveillance DRM list already
vendored at `data/reference/drm_positions.json`** (pinned `SDRM-2019`,
Bennett et al. 2009 upd. 2019), and names that list and its pin on the
page. Filtering keys on (gene, position) with the SDRM-defined states, per
that file's own `notes`.

**Counting convention (my call, flagged for override).** Sierra's strings
are messier than a plain substitution list: of the 16,357 distinct values,
~8,800 are plain substitutions, ~4,900 are two-state mixtures (`E298EQ`),
~1,700 are three-to-six-state (`T215TFILMS`), 455 are deletions
(`T39del`), 387 are stop codons (`I220*`) and 76+ are insertions in an
underscore form (`T69S_SS`). The rule:

- A mixture counts as carrying an SDRM if **any** of its states is an SDRM
  state at that position. This is the surveillance convention and it
  roughly triples hit volume versus strict single-state matching, so it
  moves every number in §8.4 — override if you read it differently.
- Stop codons and frameshift-style strings are excluded.
- Deletions and insertions count only where the SDRM file names a `del` or
  `ins` state at that position (RT 69 `ins` is the one such case), matched
  against the underscore form explicitly.
- The page states the convention in words next to the panel.

### 7.3 Denominators differ per drug class

Gene coverage is uneven: PR 31,116 · RT 30,881 · **IN 9,365**. Only 26% of
the cohort has integrase coverage, because most deposited *pol* sequences
cover PR+RT only. So there is no single cohort N, and every prevalence
figure carries its own denominator. This is also the honest framing of why
the DTG/INSTI question in Paper 1 Figure 4 is open rather than answered.

## 8. What the page shows

Deduped, from the 2026-08-17 `HIVDB_10.2` vintage (`called_at` is a
single value across the parquet; `algorithm_publish_date` is 2026-04-26).

**What is asserted, and what is not.** The export hard-fails only on
*structural* invariants, which hold for any vintage: that `.unique()` is
lossless (§7.1), that each accession has exactly one `sequence_sha256`,
that every `drug_class` is in `style.DRUG_CLASS_ORDER`, and that every
`level` is a `ResistanceLevel` member. The measured quantities below are
**recorded into `meta`, never asserted** — asserting 728,087 or 24.0%
would break the moment the parquet is regenerated, which is the next piece
of work (Approach B). The golden-file test pins exact numbers for the
synthetic fixture only, never for the real vintage.

1. **Header / provenance band.** What SAWA is, in two sentences. Cohort:
   **35,777 of 36,510 ingested accessions** returned a Sierra result
   (733 did not), giving 728,087 resistance calls across 27 drugs in 4
   classes. Algorithm `HIVDB_10.2`, called 2026-08-17, algorithm published
   2026-04-26. Per §14.1 the cohort is described by its retrieval process
   — "public HIV-1 *pol* sequences retrieved by an East-Africa-targeted
   LANL/GenBank query" — and the band carries **no region claim** about
   the sequences themselves. A one-line note
   that geographic and subtype breakdowns are deliberately absent pending a
   provenance rebuild, linking the evaluation doc — turning the main
   limitation into evidence of rigour.
2. **Resistance prevalence by drug class.** Sequence counted resistant to a
   class if any drug in it scores `high` or `intermediate` (same rule as
   `paper1/fig02_aggregate.RESISTANT_LEVELS`, reused not reimplemented).
   The denominator is "Sierra returned a score for this gene", which is
   **not** the same as "the SDRM codons were sequenced" — 10,605 of the
   30,881 RT-scored sequences are under 1,100nt, and the parquet carries
   no `firstAA`/`lastAA` because the GraphQL query does not request
   `alignedGeneSequences`, so true codon coverage cannot be computed
   without re-calling Sierra. The page defines the denominator in words
   rather than letting a reader assume the stronger one. Measured:
   NNRTI 24.0% (7,400/30,881) · NRTI 17.5% (5,413/30,881) · INSTI 3.2%
   (300/9,365) · PI 2.6% (799/31,116). Each bar shows its denominator.
3. **Per-drug level breakdown.** All 27 drugs, stacked across Sierra's five
   levels, grouped by class in `DRUG_CLASS_ORDER`. Leads with NVP 22.3% and
   EFV 17.8% high-level — the historic NNRTI first-line signal — and shows
   DTG and the other INSTIs on a visibly thinner denominator.
4. **Surveillance DRM frequency.** Top SDRM-list mutations by the number of
   sequences carrying them, split by gene (PR/RT/IN), per §7.2.
5. **Gene coverage.** PR/RT/IN coverage of the cohort, which explains the
   denominators in (2) and (3).
6. **Footer.** Repo link, CITATION.cff / Zenodo DOI, algorithm pin, data
   vintage date, and the generating commit SHA.

## 9. How it is built

**Python export** — `src/hivdrwatch/dashboard/`:

- `aggregate.py` — pure Polars functions, one per view in §8, each taking a
  DataFrame and returning a DataFrame. Mirrors the shape of
  `paper1/fig02_aggregate.py` (pure, unit-testable, no plotting import) and
  reuses its `RESISTANT_LEVELS` and `style.DRUG_CLASS_ORDER` rather than
  restating the thresholds.
- `drm.py` — loads and applies the SDRM filter from
  `data/reference/drm_positions.json`.
- `export.py` — dedups, runs the aggregations, and writes one
  `summary.json` with a `meta` block
  (algorithm version, vintage date, cohort counts, SDRM pin, git SHA,
  `DR_SCHEMA_VERSION`) plus one key per view. Sorted keys, fixed float
  precision, trailing newline, so the output is deterministic and diffable.
- `__main__.py` — CLI entry, wired into `cli/` alongside the other stages
  and given a `make dashboard` target.

**Layout.** `src/hivdrwatch/dashboard/` splits into `export/` (the Python
above) and `static/` (the page), so the API-backed Auspice/deck.gl
dashboard MODULES.md Flow 1 specifies can later land beside them as
`service/` without this slice having claimed its name. This slice does not
touch Flow 1's acceptance criteria.

**Static page** — `src/hivdrwatch/dashboard/static/`: one `index.html`,
one stylesheet, one small ES-module script, and the committed
`summary.json`. Charts drawn as inline SVG from the JSON. No framework, no
bundler, no npm — which keeps the page dependency-free, fast on a weak
connection, and reviewable as plain text. Colours come from
`paper1/style.COLORBLIND_PALETTE` (Okabe-Ito), so the page matches the
manuscript figures and stays colourblind-safe; the five resistance levels
get an ordered ramp from that palette.

**Deploy** — a GitHub Actions workflow publishing the static directory to
GitHub Pages on push to `main`. CI does not run the export (no DVC data,
§6); it publishes the committed JSON.

**Drift guard — a weaker guarantee than first drafted.** The only DVC
remote is `s3://hivdrwatch-dvc/` at `http://localhost:9000` with dev
credentials, so no GitHub-hosted runner can fetch the parquet and
recompute. The guard is therefore: `summary.json`'s `meta` records the
md5 of the input parquet as the `.dvc` pointer states it, and CI fails if
that recorded hash does not match the committed pointer. This detects a
*stale input* with no data present, but it does **not** verify the
aggregation. A self-hosted runner or a reachable remote would allow the
stronger check; until then the doc should not claim it.

## 10. Testing

- Unit tests per aggregation in `tests/dashboard/`, on a small synthetic
  fixture under `tests/fixtures/` named `_synthetic_*` per CLAUDE.md.
- A test asserting the dedup is lossless on a fixture built to contain
  duplicates, so §7.1's invariant is enforced rather than documented.
- A golden-file test on `summary.json` for the synthetic fixture, with
  `--update-golden`, matching the repo's existing parser-test convention.
- A test that the SDRM filter rejects the known polymorphisms in §7.2
  (R41K, M36I, V35T, K122E) and admits known DRMs (e.g. RT M184V, K103N),
  so §7.2's trap cannot regress.
- Coverage ≥ 70% on `dashboard/`.

## 11. Documentation obligations

- `src/hivdrwatch/dashboard/README.md` — what it does, how to run it,
  inputs/outputs, tests (CLAUDE.md requires this per module).
- A `dashboard/` row in MODULES.md §1.1, status `~`, blocking Paper 1,
  worded to cover **only** the static resistance summary and to state
  that Flow 1's API-backed Auspice + deck.gl dashboard is not started and
  that the stratified views await the provenance rebuild. The row adds to
  Flow 1's contract; it does not redefine its acceptance criteria.
- A CHANGELOG entry.
- No CONTEXT.md change: scope, venue and authorship are unaffected.

## 12. What is left to the builder

- Exact SVG chart geometry, spacing, type scale, and responsive
  breakpoints, within §9's palette and no-framework constraint.
- The copy in the header's two-sentence description and the limitation
  note, in CONTEXT.md's voice.
- How many DRMs the §8.4 view shows before truncating, and whether it
  truncates by rank or by a prevalence floor.
- JSON key naming and whether views nest under a `views` object.
- Whether the staleness check in §9 runs as a separate workflow or a job
  within `test.yml`.
- Whether `make dashboard` writes into the static directory in place or via
  a staging path.

## 13. Open decisions for Kelyn

Recorded as decided unless he says otherwise:

- **Hosting:** GitHub Pages off this repo, not a personal domain — free,
  versioned with the code, and the repo is already JOSS-shaped and
  outward-facing. Prerequisite **verified**: the repo is public and Pages
  is available; it needs enabling once (§14.3).
- **Module name:** `dashboard/` under `src/hivdrwatch/`, matching
  MODULES.md's wording in sequencing item 4 and Flow 1.
- The duplicate-row root cause (§7.1) is filed as an issue, not fixed here.

---

## 14. Open questions — Kelyn's call, blocking the build

Surfaced by the builder check (2026-10-30). Each one changes what gets
built; none can be settled from the repo.

### 14.1 Can the page say "East Africa"?

§8.1 made "East Africa" a property of the cohort. It cannot be, on this
vintage. `ingest-provenance.md` §3 establishes that `/country=` appears in
**zero** cache entries (NCBI renamed it `/geo_loc_name=`), so
`_normalize_country` received `None` for every record and the query
fallback fired every time: **every country value in the shipped table is
query provenance, not a statement about the record.** Corrected against
the records, **3,108 accessions state a country outside East Africa** —
Thailand ~2,150 rows (the RV217 Thai arm), plus MW, BR, ZA, ZW, IN, ZM, IL
and 16 further destinations — and 277 state none.

So the one geographic claim the page makes would be sourced from exactly
the field §4 exists to keep off the page, and ~8.5% of the cohort
contradicts it. It also colours how the prevalence numbers read: 24.0%
NNRTI is currently a figure over a corpus including ~3,100 non-East-African
sequences.

Three paths:

- **(a) Keep the region word, disclose the over-match** in the limitation
  note.
- **(b) Drop the region word**, describe the cohort by process: "public
  HIV-1 *pol* sequences retrieved by an East-Africa-targeted LANL/GenBank
  query". Accurate on this vintage without needing corrected geography.
- **(c) Exclude the 3,108** — which requires the corrected country, i.e.
  Approach B, i.e. this slice stops shipping first.

**My recommendation: (b).** It is true today, it costs nothing, it keeps
the sequencing decision intact, and describing a cohort by its retrieval
process reads as careful rather than evasive to the audience in §1.

> **RESOLVED (Kelyn, 2026-10-30): (b).** No region claim anywhere on the
> page. The cohort is described by its retrieval process. The limitation
> note links `ingest-provenance.md` so a reader can see why, and the
> geographic view arrives with Approach B.

### 14.2 Does an integrase DRM panel ship in this slice?

`data/reference/drm_positions.json` self-describes as *"Positions used as
candidate features for the DR model"* — a feature-engineering artifact, not
a publication reference. Two problems with citing it on the page as
"Stanford SDRM list, Bennett et al. 2009 upd. 2019":

- **Bennett 2009 has no integrase component.** The file's IN block is 8
  positions (66, 92, 121, 140, 143, 147, 148, 155) — the
  raltegravir/elvitegravir-era major list — and it is **missing the
  DTG-era mutations** (R263K, G118R, E138K/A). Citing Bennett for the IN
  panel misattributes it to a list that does not contain it.
- The IN panel it produces is thin: T66A 31, N155H 27, G140A 25, G148K 19,
  S147G 15, down to Y143H 1, out of 9,365 integrase-covered sequences.

PR and RT look like a faithful SDRM list (M184V 4,383, K103N 3,025, G190A
1,503, Y181C 1,367 — which reads exactly right).

**My recommendation: ship PR and RT only**, cite Bennett 2009 for those
two, and state plainly that integrase DRMs are not shown because the
vendored list predates the DTG era. That turns the gap into the same
honest-limitation move as §8.1 and avoids putting a misattributed list in
front of Leigh Brown. Replacing the IN block with a current list
(e.g. the Stanford INSTI major list) is its own piece of reference-data
work with its own citation, and it touches Paper 1 Figure 4's story.

> **RESOLVED (Kelyn, 2026-10-30): PR and RT only.** Bennett et al. 2009
> cited for those two genes. The page states that integrase DRMs are not
> shown because the vendored list predates the DTG era. Replacing the IN
> block is separate work, not part of this slice.

### 14.3 Is `Neoforge-Labs/SAWA` going public?

The whole deliverable rests on this and it is unverified. Origin is
`git@github.com:Neoforge-Labs/SAWA.git`; GitHub Pages needs the repo
public or Pages enabled on the org plan. If neither, the hosting decision
in §13 reopens.

Going public also publishes, in the same push:

- `paper1_tools/preregistration_fig04.md`, which MODULES.md says is
  **awaiting OSF submission**. Publishing a preregistration before filing
  it is a research-integrity call, not a hosting detail.
- The unsubmitted Paper 1 manuscript skeleton.
- `.dvc/config` with `minioadmin/minioadmin` — dev-only credentials for a
  localhost endpoint, so low risk, but it should not go out unreviewed.

> **RESOLVED (verified 2026-10-30, no decision needed).**
> `Neoforge-Labs/SAWA` is **already public** (`visibility: PUBLIC`), and
> the org is on an enterprise plan, so GitHub Pages is available at no
> cost. Pages is not yet configured (`GET /repos/.../pages` → 404), so the
> only remaining step is enabling it — a repository settings change for
> Kelyn to make or approve, not a design question.
>
> The exposure concerns above do not survive checking:
>
> - `paper1_tools/preregistration_fig04.md` is already public **by
>   design** — its own header reads "pending OSF submission — *locked here
>   in git until then*", i.e. git is deliberately being used as the
>   timestamp ahead of OSF filing. My raising it was over-cautious.
> - No `.env`, key, or credential file is tracked on `main`.
> - `.dvc/config` carries `minioadmin/minioadmin` against
>   `http://localhost:9000` — docker-compose dev defaults for an endpoint
>   that is unreachable from anywhere but the developer's own machine. Not
>   a secret; worth tidying someday, not a blocker and not this work.
>
> **No scrub pass is needed.**
