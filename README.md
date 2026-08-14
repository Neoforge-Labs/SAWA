# SAWA — Sub-Saharan Africa Watch for Antiretroviral resistance

*(also: HIV-DR-Watch — the Python package name)*

**An open real-time antiretroviral-resistance surveillance platform, starting with HIV-1 in East Africa.**

SAWA ingests public HIV-1 *pol* sequences from LANL and GenBank, calls resistance mutations against Stanford HIVdb, builds subtype-stratified phylogenies, detects transmission clusters, and stratifies resistance prevalence across country, year, and subtype — with a design bias toward the non-B subtypes (A1, C, D) that dominate East African epidemics but that most existing tools underserve.

The name is a backronym and a word: **sawa** means "OK / fine / it's good" in Swahili. The scope name — *Antiretroviral resistance* — is deliberately broader than HIV so the same platform can extend to HBV and TB drug resistance in later releases.

Everything reproduces from `git clone && docker compose up && make paper1`.

---

## Statement of need

Existing HIV drug-resistance surveillance tooling is one of three things: (a) country-specific and closed, (b) built for the sequencing-rich, subtype-B-dominant global north, or (c) not maintained. Sub-Saharan Africa carries most of the global HIV-1 burden and almost none of it is subtype B. Pre-treatment NNRTI resistance in East Africa crossed the 10 % WHO threshold years before the 2019 dolutegravir rollout; INSTI resistance under DTG monotherapy is now emerging in the same populations but is largely uncharted in the public literature. Public health teams and reference laboratories need an open pipeline that runs on the sequences they have, on the hardware they have, in the subtypes they see.

HIV-DR-Watch fills that gap by:

- Treating LANL Table exports and GenBank Entrez pulls as first-class ingest sources with a shared canonical schema.
- Calling resistance via the Stanford HIVdb Sierra GraphQL API, algorithm-version-pinned and locally cached so runs are reproducible offline.
- Aligning to HXB2 K03455 with `mappy` (~700 seq/s in-process) so a full East-Africa cohort aligns in under a minute on a laptop.
- Building per-subtype phylogenies with IQ-TREE 2 + TreeTime and detecting transmission clusters with tn93 + connected components, using the same conventions HIV-TRACE reviewers expect.
- Emitting figures from Marimo notebooks under `notebooks/` so every claim in the paper maps to a code path, and every benchmark in `benchmarks/` regenerates in CI.

The immediate deliverable is a tools paper (Paper 1, Virus Evolution). Three sibling preprints (causal identification, foundation-model benchmark, KEMRI deployment) live in the same repo under `causal/`, `foundation_models/`, and `deployments/kemri/` and cite the shared library.

---

## Installation

Prereqs: **Docker** (for tn93, IQ-TREE, MAFFT images) and **uv** (Python dep manager).

```bash
git clone https://github.com/neoforge-labs/hiv-dr-watch
cd hiv-dr-watch
make setup                    # copy .env.example → .env; uv sync
make images                   # build tn93 + MAFFT + IQ-TREE images
make up                       # docker compose up (postgres + minio + prefect + mlflow)
```

That's it. `make help` lists everything else.

---

## Quick start

```bash
# 1. Ingest — LANL Tables + GenBank incremental pull
make ingest                                # LANL only (offline)
NCBI_EMAIL=you@example.com make ingest-genbank    # LANL + GenBank

# 2. Call resistance — Stanford HIVdb Sierra with local cache
make sierra

# 3. Align to HXB2 pol
make curate

# 4. Cluster (tn93 + connected components)
make cluster

# 5. Per-subtype phylogenies (IQ-TREE 2, fast mode)
make phylo

# 6. Regenerate paper figures
make fig2       # resistance landscape heatmap
make fig5       # scalability plot
```

The Prefect UI is at http://localhost:4200, MLflow at http://localhost:5000, MinIO at http://localhost:9001.

---

## Functionality

The core platform (`src/hivdrwatch/`) is organized by pipeline stage. Each subpackage has its own README explaining inputs, outputs, and tests.

| Module | Purpose |
|---|---|
| `ingest/` | LANL TSV + FASTA parser, GenBank Entrez fetcher with resume-safe batching and per-country watermarks. Canonical `SequenceRecord` schema. |
| `resistance/` | Stanford HIVdb Sierra GraphQL client. SHA256-keyed local cache. Algorithm version pinned in every response. |
| `align/` | Pairwise-to-HXB2 via `mappy` (minimap2 bindings); Docker MAFFT `--add --keeplength` and Biopython fallbacks. Trim, QC, dedup. |
| `phylo/` | IQ-TREE 2 (GTR+F+I+G4) wrapper + TreeTime molecular-clock scaling. Per-subtype orchestrator. |
| `cluster/` | tn93 pairwise distances + networkx connected components at the standard 0.015 threshold. Deterministic cluster IDs. JC69 ground-truth simulator + pair-based P/R/F1/ARI validator. |
| `pipeline/` | Prefect 3 flows chaining every stage. Live subprocess streaming to the Prefect UI. |
| `paper1/` | Pure-Python aggregators consumed by figure notebooks. Kept out of the notebooks so figures are unit-testable. |
| `causal/` (repo root) | DAGs, DoWhy identify + estimate + refute, feature pruning, HTML report — foundations for Paper 2. |

All headline claims are backed by scripts in `benchmarks/` whose JSON results are checked in and refreshed weekly by CI. Anything missing that link doesn't ship.

---

## Target audience

- **Phylodynamics and public-health-genomics researchers** who need an open, reproducible pipeline over the specific subtypes they work with.
- **African reference laboratories** (starting with KEMRI Nairobi) that need something they can deploy on modest hardware, run offline, and hand off to non-bioinformatics staff.
- **Method developers** who want a substrate they can plug new resistance predictors, subtypers, or causal-identification methods into without rebuilding the surveillance stack.

The dashboard (Phase 2) is written for the third audience and its bandwidth constraints. The pipeline is written for the first two.

---

## Reproducibility guarantees

- **Every figure comes from a Marimo notebook under `notebooks/`.** The PDF/PNG is a build artifact.
- **Every stochastic step sets a seed derived from its config hash.** Two runs with the same config produce bit-identical outputs.
- **Every dependency is pinned** via `uv.lock`; bio tools (MAFFT, IQ-TREE, tn93) via bioconda inside Docker.
- **Data isn't in git.** DVC tracks `data/raw/`, `data/interim/`, `data/processed/`.
- **Every release** is tagged `paperN-vX.Y.Z` and gets a Zenodo DOI. Preprints cite the DOI, not the branch.

---

## Contributing

Read `CLAUDE.md` first — the working conventions apply to human contributors too (Conventional Commits, small atomic commits, no `# TODO`s in code, KEMRI data gated behind DSA sign-off). PRs run `ruff`, `mypy --strict`, and `pytest` in CI.

---

## Papers

| # | Title | Venue | Timeline |
|---|---|---|---|
| 1 | HIV-DR-Watch — tools paper | Virus Evolution | Q1 2027 |
| 2 | CausalPhylo — DAG + DoWhy for HIV phylodynamics | PLOS CompBio | Q3 2027 |
| 3 | Foundation-model resistance benchmark | Bioinformatics / NeurIPS D&B | Q4 2027 |
| 4 | KEMRI deployment case study | Bioinformatics App Note | Q1 2028 |
| 5 | JOSS software citation | Journal of Open Source Software | Continuous |

See `CONTEXT.md` for the narrative arc and `MODULES.md` for the module-level definition of done.

---

## Citation

If you use HIV-DR-Watch in research, please cite via the `CITATION.cff` file in this repository (GitHub renders a "Cite this repository" button in the sidebar). On tagged releases a Zenodo DOI is issued and should be preferred.

---

## License

MIT. See `LICENSE`.
