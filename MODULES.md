# HIV-DR-Watch — Modules & Flows: Definition of Done

*What must exist in this repo — and what must actually run end-to-end —
for each of the four preprints (Papers 1–4) plus the JOSS citation to be
publishable. Read alongside CONTEXT.md.*

Status legend: **✓** shipped · **~** scaffolded · **·** to build

---

## Phase 0 alignment (2026-08-13)

Repo layout re-aligned to this document's contract via
`scripts/phase0_cleanup.sh`. Deletions + renames + moves recorded as
one commit:

- **Deleted (out-of-scope for Paper 1):** `src/hivdrwatch/models/`,
  `src/hivdrwatch/features/`, `src/hivdrwatch/naomi/`, and the four
  pipeline flows that fed them (`train_dr.py`, `train_outbreak.py`,
  `dr_features.py`, `naomi_ingest.py`). Also `notebooks/dashboard.py`
  and `data/gold/` + `data/models/`.
- **Moved:** `src/hivdrwatch/causal/` → `causal/` at repo root (per
  monorepo layout in CONTEXT.md); associated flow moved to
  `causal/flows/causal.py`.
- **Renamed** to match this document's naming: `hivdb/` → `resistance/`,
  `curate/` → `align/`, `flows/` → `pipeline/`.
- **Scaffolded** empty `workflows/` for future Nextflow paper analyses.
- **Updated** every import + `pyproject.toml` packages + Makefile targets.

---

## Part 1 — Required modules

### 1.1 Core platform (`src/hivdrwatch/`) — anchors Paper 1

| Module | Status | What it does | Blocks paper |
|---|---|---|---|
| `ingest/` | ✓ | LANL TSV + FASTA parsing; GenBank Entrez fetcher with resume-safe batching + bisection + watermarks; canonical `SequenceRecord` schema; Parquet outputs with version stamps. Tested end-to-end against real Kenya + full East Africa cohort. | 1, 2, 3, 4 |
| `resistance/` | ~ | Stanford HIVdb Sierra GraphQL client with sha256 cache + algorithm-version pinning; batching to 100/req; offline replay mode. **Missing:** concordance benchmark vs held-out annotated set. | 1, 2, 3 |
| `align/` | ~ | `mappy` (minimap2 bindings) in-process pairwise-to-HXB2 with `--keeplength` semantics; ~1,900 seq/s single-thread at the `asm20` default. Keeps every forward-strand alignment block, not just the largest. Docker MAFFT `--add` retained as opt-in for MSA cases; Biopython PairwiseAligner as portable fallback. `NO_ALIGNMENT` blocks rows with no PR/RT/IN site covered. **Missing:** incremental `--add` cache for the live path; per-subtype alignment cache. **Note:** the shipped `curated/aligned.parquet` vintage predates the `asm20` fix — half its rows are all-gap; regenerate before using it for figures. See `docs/evaluations/align-coverage.md`. | 1, 2, 4 |
| `phylo/` | · | IQ-TREE 2 wrapper (GTR+F+I+G4); TreeTime time-scaling; per-subtype reference trees; incremental grafting for new sequences | 1, 2, 4 |
| `cluster/` | ✓ | Distance clustering via tn93 (Docker) + networkx connected components at 0.015 threshold; deterministic cluster IDs; time-series + subtype-entropy + geographic-spread dynamics. HIV-TRACE parity via `hivtrace.py` (bundled `hivnetworkcsv` in `docker/hivtrace_full/`); `cluster_validate` runs both clusterers and scores pair-agreement. JC69 synthetic ground-truth simulator in `simulate.py`. Real-FAVITES upgrade is Phase 2 polish. | 1 |
| `storage/` | ✓ | SQLAlchemy 2.x tables (`pipeline_runs`, `sequences`, `dr_calls`, `cluster_memberships`) with per-table uniqueness constraints; Alembic-managed migrations under `migrations/`; idempotent Parquet→Postgres loader with `ON CONFLICT DO NOTHING`; DuckDB analytics wrapper (`open_analytics()`) reading Parquet directly or `postgres_scan` when hosted. SQLite fallback for tests. | 1, 2, 3, 4 |
| `api/` | ✓ | FastAPI service with read-only `/sequences`, `/dr-calls`, `/clusters` + `/health`, `/version`. Auto-OpenAPI at `/docs`. `X-Request-Id` middleware, `slowapi` per-IP rate limiting, CORS. `create_app()` factory; dep-injection via `get_db`. Docker Compose service exposed on `API_PORT` (default 8080). Writes gated behind DSA. | 1, 4 |
| `pipeline/` | ~ | Prefect 3 flows for `ingest`, `dr_ingest` (Sierra), `curate` (align + trim + QC + dedup), `cluster` — all runnable + MLflow-tracked. **Missing:** daily/weekly scheduled deployments; 30-day unattended run proven. | 1, 4 |
| `cli/` | ✓ | Unified `hivdrwatch` Typer entry point wrapping every stage: `ingest`, `sierra`, `curate`, `cluster`, `validate-clusters`, `phylo`, `bench-*`, `version`, `commands`. Thin dispatch — extra args forward verbatim to each module's `main(argv)`. Tests cover discovery + forwarding + contract. | 1, 4 |

**Definition of "core platform done":** given a fresh clone,
`docker compose up` starts every service; `hivdrwatch pipeline run --once`
produces a Parquet dataset, a resistance report, a subtype-stratified
tree per subtype, and a cluster set — all with schema versions stamped.

### 1.2 Causal work (`causal/`) — anchors Paper 2

| Module | Status | What it does |
|---|---|---|
| `causal/dag.py` | ✓ | The domain DAG: exposure(s), outcome, covariates, unobserved confounders. Two DAGs: DRDag (mutation → resistance) + OutbreakDag (deposit-rate confounder baked in). |
| `causal/estimate.py` | ✓ | DoWhy identify + estimate + 3-refuter loop (placebo / random-common-cause / subset); bootstrap percentile CIs; per-subtype stratified estimator; batch drivers for DR and outbreak. **Merges `identify.py` + `estimate.py`** from the original spec — split if reviewers want the separation. |
| `causal/pruning.py` | ✓ | Rule-based feature pruning from causal effect frame (drop when CI spans 0 with small effect, or fails ≥2/3 refutations); keep-wins-across-scopes. Named `pruning.py` not `feature_pruning.py`. |
| `causal/report.py` | ✓ | Effect-frame Polars serializer + self-contained HTML report grouped by data scope. |
| `causal/flows/causal.py` | ✓ | Prefect flow: DR gold + outbreak gold → effects parquet + HTML + prune configs; MLflow-tracked. |
| `causal/sensitivity.py` | · | Rosenbaum bounds, e-values, tipping-point analysis (beyond the 3 refutations already in estimate.py). |
| `causal/counterfactual.py` | · | Counterfactual queries — the "what if DTG had rolled out in 2017?" story. |
| `causal/data_prep.py` | · | Adapters converting HIV-DR-Watch Parquet outputs into the causal dataset (unit of analysis: patient-visit or clinic-quarter). |

**Definition of "causal work done":** given the Paper 1 Parquet outputs,
`python -m causal.paper2 run-all` produces the DAG SVG, identification
result printout, point estimates + CIs for the primary and secondary
exposures, a sensitivity analysis table, and one counterfactual figure.

### 1.3 Foundation-model benchmark (`foundation_models/`) — anchors Paper 3

| Module | Status | What it does |
|---|---|---|
| `fm/data.py` | ✓ | Subtype-stratified train/val/test splits from SAWA Paper 1 outputs. Two strategies: `stratified` (70/15/15 on subtype×drug×label) and `subtype_holdout` (train on B, test on non-B — Paper 3's thesis). Deterministic given seed. |
| `fm/kmers.py` | ✓ | Nucleotide/protein k-mer featurizer → scipy sparse CSR matrix; explicit vocab support for fit/transform. |
| `fm/tokenize.py` | · | Protein tokenization for ESM-family; nucleotide tokenization for Evo (GPU path) |
| `fm/models/esm2.py` | · | ESM-2 / ESM-C fine-tuning heads for resistance prediction (GPU) |
| `fm/models/evo.py` | · | Evo fine-tuning wrapper (GPU) |
| `fm/baselines.py` | ✓ | Non-FM baselines: `KmerLogisticRegression` (L2, liblinear) + `OneHotRandomForest`. Runs in seconds. These are what the FMs must beat. |
| `fm/train.py` | · | Fine-tuning loop with subtype-stratified curriculum (GPU) |
| `fm/evaluate.py` | ✓ | AUROC, PR-AUC, Brier, ECE (10-bin default). `per_subtype_metrics` returns a Polars metrics frame; cells with < 5 positives dropped. |
| `fm/failure_analysis.py` | · | The "where do foundation models fail on non-B subtypes" analysis — Figure 2 of Paper 3 (GPU) |

**Definition of "foundation-model work done":** `python -m foundation_models.paper3 benchmark`
reproduces the full benchmark table (models × subtypes × drug classes)
and emits the four Paper 3 figures.

### 1.4 KEMRI deployment (`deployments/kemri/`) — anchors Paper 4

| Module | Status | What it does |
|---|---|---|
| `deployments/kemri/docker-compose.kemri.yml` | · | Compose overrides for low-bandwidth clinic hardware |
| `deployments/kemri/offline_sync.py` | · | Offline-first sync module; conflict-free replicated data type (CRDT) or last-write-wins queue |
| `deployments/kemri/data_adapters/` | · | KEMRI-specific data ingestion adapters (post-DSA only) |
| `deployments/kemri/audit_log.py` | · | Deployment audit trail — what ran, when, on what data version. Becomes Paper 4's methods section |
| `deployments/kemri/README.md` | · | The deployment runbook Chohan's team actually follows |

**Definition of "KEMRI deployment done":** the stack runs on a mid-range
Linux box at KEMRI Nairobi with intermittent connectivity; one full
surveillance cycle completes offline; results sync when connectivity
returns; audit log captures every step.

### 1.5 Cross-cutting infrastructure (repo-wide) — required for all papers

| Concern | Requirement | Status |
|---|---|---|
| Environment | Pixi or uv lockfile pinning every dep including MAFFT / IQ-TREE via conda-forge/bioconda | ~ (uv + `uv.lock` present; MAFFT/IQ-TREE not yet in the env — mappy in-process for now; MAFFT is Docker-only if invoked) |
| Containerization | Docker Compose stack that runs from a fresh clone with no host deps beyond Docker | ~ (docker-compose.yml stands up Postgres + MinIO + Prefect server + MLflow. Analysis containers — MAFFT, tn93 — separate. No unified `docker compose up` producing a paper artifact yet.) |
| Data versioning | DVC tracking `data/raw/`, `data/interim/`, `data/processed/`; hashes cited in every manuscript | ~ (DVC configured with MinIO remote in `.dvc/config`; per-directory `.gitignore` allows `*.dvc` pointers; `DEPLOYMENT.md` documents the flow. Awaiting first `dvc add` on a real vintage.) |
| Reproducible pipelines | Nextflow pipelines under `workflows/` following nf-core conventions for every paper's analysis | ~ (empty `workflows/` scaffolded with a README explaining the Prefect-vs-Nextflow split. No `.nf` files yet.) |
| CI | GitHub Actions running: unit tests, integration smoke test, full benchmark suite, figure regeneration | ✓ (`.github/workflows/{test,benchmarks,figures}.yml` — ruff+mypy+pytest on PR; weekly benchmark refresh with auto-commit; nightly figure regeneration.) |
| Software citation | `CITATION.cff` on day 1; Zenodo–GitHub integration for auto-DOI on tagged release | ✓ (`CITATION.cff` + `.zenodo.json`; DOI issues automatically on `paperN-vX.Y.Z` GitHub Release once the Zenodo–GitHub integration is enabled in the repo settings.) |
| Reproducibility artifact | Every release tagged with `paperN-vX.Y.Z` for the frozen state that matches the preprint | · |
| Manuscripts | Quarto (`manuscript.qmd`) per paper under `paperN_*/`, pulling figures from notebooks | ~ (`paper1_tools/manuscript.qmd` + `refs.bib` skeleton with all 5 figure placeholders wired. Paper 2-4 skeletons pending.) |
| Preregistration | OSF preregistration for the Paper 1 case-study analysis; filed before running the analysis | ~ (`paper1_tools/preregistration_fig04.md` locked in git — cohort, primary hypothesis, statistical plan, sensitivity analyses, BH-FDR. Awaiting OSF submission before Fig 4 runs on real data.) |
| Benchmarks | `benchmarks/` folder with results JSON checked in and updated by CI; feeds paper tables directly | ✓ (`benchmarks/{scale_ingest,scale_align,scale_cluster}.py` + `BenchmarkResult` schema + `measure_block` context manager; results under `benchmarks/results/*.json`; weekly CI refresh; Fig 5 consumes.) |
| Documentation | README (JOSS-shaped), CONTRIBUTING.md, DEPLOYMENT.md, per-module READMEs | ✓ (JOSS-shaped root `README.md`, `DEPLOYMENT.md`, per-module READMEs for `ingest/`, `resistance/`, `align/`, `phylo/`, `cluster/`, `pipeline/`, `paper1/`, `causal/`, `cli/`. `CONTRIBUTING.md` still to write.) |
| Experiment tracking | MLflow server (via docker-compose) with per-run parameters, metrics, artifacts | ✓ (`src/hivdrwatch/tracking.py` no-op wrapper; every pipeline flow logs) |

---

## Part 2 — End-to-end flows that must work

A module existing in isolation isn't enough. These are the integration
flows that must actually run — each one is what makes a specific claim
in a specific paper defensible.

### Flow 1 — Cold-start ingestion → published dashboard *(Paper 1)*

```
LANL TSV/FASTA (manual)  ─┐
                          ├─→ ingest ─→ Parquet ─→ resistance ─→ Postgres
GenBank Entrez ───────────┘                    ─→ align ─→ phylo ─→ cluster
                                                                     │
                                                                     ↓
                                                        api ─→ dashboard
                                                                (Auspice + deck.gl)
```

**Acceptance:** on a fresh clone, `docker compose up && make cold-start`
produces a working dashboard at `localhost:8080` with real Kenya data
inside 60 minutes on a laptop. **This is Paper 1 Figure 1.**

### Flow 2 — Incremental daily update *(Paper 1, Paper 4)*

```
Prefect scheduled job (04:00 UTC daily)
  → GenBank fetch (only accessions with MDAT > last_ingest_date)
  → dedupe against LANL + prior GenBank pulls
  → resistance-call (Sierra with local cache)
  → incremental alignment (MAFFT --add against reference)
  → tree grafting (TreeTime incremental)
  → cluster re-detection (only affected components)
  → dashboard refresh
```

**Acceptance:** the flow runs unattended for 30 consecutive days on a
staging server without human intervention. Failure recovery: any single
step can fail and be retried without corrupting downstream state.

### Flow 3 — Full paper regeneration *(Papers 1, 2, 3, 4)*

```
make paper1  # regenerates every Paper 1 figure from raw data
make paper2  # same for Paper 2 (depends on Paper 1 Parquet outputs)
make paper3  # same for Paper 3
make paper4  # same for Paper 4
make all     # regenerates everything, in dependency order
```

**Acceptance:** on a machine with the pinned env, `make all` completes
in under 8 hours, produces bit-identical figures across runs given the
same seeds, and re-renders each `manuscript.qmd` to PDF.

### Flow 4 — Cluster detection with ground-truth validation *(Paper 1 Figure 3)*

```
FAVITES / VirusPopGen simulator
  → synthetic transmission network with known infector→infectee edges
  → export as FASTA + ground-truth edge list
  → run HIV-TRACE on FASTA
  → run our custom clusterer on FASTA
  → compare both against ground truth: precision/recall/F1
  → emit results.json + figure
```

**Acceptance:** the harness runs in CI on every PR touching `cluster/`;
regressions block merge.

### Flow 5 — DTG-INSTI resistance case study *(Paper 1 Figure 4 — the biological headline)*

```
Preregister analysis on OSF (before touching data)
  → select DTG-rollout cohort (post-2019 sequences from countries where
    DTG became first-line)
  → compute INSTI major-mutation prevalence by year × country × subtype
  → statistical test: DiD or ITS against pre-rollout baseline
  → sensitivity analysis: alternate cohort definitions
  → figure + table with 95% CIs
```

**Acceptance:** the preregistration is public, the analysis code is a
single notebook (`notebooks/fig04_dtg_case_study.py`), and re-running
it against the frozen Parquet vintage produces the identical figure.

### Flow 6 — Causal identification pipeline *(Paper 2)*

```
Paper 1 Parquet outputs
  → causal/data_prep: aggregate to unit of analysis
  → causal/dag: load DAG, validate against data columns
  → causal/identify: DoWhy backdoor / IV identification
  → causal/estimate: point estimate + CI for each exposure
  → causal/feature_pruning: reduce covariates, re-estimate
  → causal/sensitivity: refutation tests + Rosenbaum bounds
  → causal/counterfactual: intervention query
  → paper2_causal/figures/: emit all 4 figures
```

**Acceptance:** every step is deterministic given seed; every estimate
comes with a refutation-test result; the DAG is committed as a `.dot`
file that renders identically for reviewers.

### Flow 7 — Foundation-model benchmark *(Paper 3)*

```
Paper 1 Parquet
  → fm/data: subtype-stratified splits (train B-dominant, test A1/C/D)
  → fm/tokenize: prep for ESM-2, ESM-C, Evo
  → fm/train: fine-tune each model (checkpoints to shared storage)
  → fm/baselines: fit logistic regression + RF baselines
  → fm/evaluate: AUROC/PR-AUC/calibration per model × subtype × drug class
  → fm/failure_analysis: identify which mutations FMs miss on non-B
  → paper3_foundation/figures/: emit all 4 figures + benchmark table
```

**Acceptance:** benchmark reproducible on a single A100 (or equivalent)
in under 24 hours; results JSON is checked into the repo and cited in
the paper.

### Flow 8 — KEMRI offline deployment cycle *(Paper 4)*

```
Deploy stack on KEMRI hardware (docker-compose.kemri.yml)
  → ingest local sequence batch (KEMRI FASTA, offline)
  → run full pipeline offline
  → audit_log records every step with data hashes
  → when connectivity returns: offline_sync pushes results to central instance
  → central instance merges without duplicates
```

**Acceptance:** the deployment runs a real weekly cycle at KEMRI
Nairobi for at least 8 weeks before Paper 4 submission. Audit logs
become the paper's methods section.

---

## Part 3 — Publishability checklist (per paper)

Every paper must clear this checklist before preprint submission.

- [ ] Reviewer can `git clone && docker compose up` and get a working stack
- [ ] `make paperN` regenerates every figure from raw data
- [ ] Every headline claim has a benchmark script in `benchmarks/`
- [ ] All datasets are DVC-tracked; data hashes cited in the manuscript
- [ ] `CITATION.cff` updated with authors, ORCIDs, current version
- [ ] Zenodo DOI issued for the tagged release matching the preprint
- [ ] Preregistration (OSF) filed for any confirmatory analysis
- [ ] FAIR data statement in the manuscript
- [ ] Code coverage ≥ 70% on modules cited in the paper
- [ ] Docs updated: per-module README, DEPLOYMENT.md if applicable
- [ ] Manuscript rendered to PDF from Quarto with zero warnings
- [ ] Independent reviewer (Chohan, Weisse, or a peer) has run the flow end-to-end without hand-holding

---

## Part 4 — Sequencing (what to build in what order)

Build order is dictated by paper dependency, not by module aesthetic
preference:

1. **Finish `ingest/`** (real integration test on Kelyn's LANL Kenya download)
2. **Build `resistance/`** — unblocks Paper 1's headline figure and Paper 3's benchmark target
3. **Build `align/` + `phylo/` + `cluster/`** — completes Paper 1's technical spine
4. **Stand up `storage/` + `api/` + minimal `dashboard/`** — the artifact reviewers will click
5. **Wire the CI + reproducibility scaffold** — before writing manuscript, not after
6. **Write Paper 1 manuscript in `paper1_tools/`** — first preprint out
7. **Wire `causal/data_prep`** so the existing `causal/` code runs against Paper 1 outputs
8. **Write Paper 2 manuscript** — second preprint
9. **Build `foundation_models/`** — third preprint
10. **Deploy at KEMRI (with Chohan)** — fourth preprint, co-authored

Anything not on this list is Not Yet.
