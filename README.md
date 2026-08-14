# SAWA — Sub-Saharan Africa Watch for Antiretroviral resistance

[![test](https://github.com/Neoforge-Labs/SAWA/actions/workflows/test.yml/badge.svg)](https://github.com/Neoforge-Labs/SAWA/actions/workflows/test.yml)
[![benchmarks](https://github.com/Neoforge-Labs/SAWA/actions/workflows/benchmarks.yml/badge.svg)](https://github.com/Neoforge-Labs/SAWA/actions/workflows/benchmarks.yml)
[![figures](https://github.com/Neoforge-Labs/SAWA/actions/workflows/figures.yml/badge.svg)](https://github.com/Neoforge-Labs/SAWA/actions/workflows/figures.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/downloads/)

**An open, reproducible platform for real-time antiretroviral-resistance surveillance in sub-Saharan Africa** — starting with HIV-1 in East Africa. SAWA ingests public HIV-1 *pol* sequences (LANL, GenBank), calls resistance mutations against Stanford HIVdb Sierra, aligns to HXB2 pol with `mappy`, builds per-subtype phylogenies with IQ-TREE 2, detects transmission clusters at the HIV-TRACE convention, and renders publication-quality figures — all from a single `docker compose up` on a laptop.

**Sawa** means *OK / fine / it's good* in Swahili. The Python package is `hivdrwatch`; the public/brand name is SAWA.

> **Reproduce every paper figure from raw data:** `git clone && make setup && make images && make up && make paper1`

---

## Statement of need

Existing HIV drug-resistance surveillance tooling is one of three things: **country-specific and closed**, **built for the sequencing-rich subtype-B-dominant global north**, or **unmaintained**. Sub-Saharan Africa carries most of the global HIV-1 burden, and almost none of it is subtype B. Pre-treatment NNRTI resistance in East Africa crossed the 10 % WHO threshold years before the 2019 dolutegravir rollout; INSTI resistance under DTG-based first-line is now emerging in the same populations but is largely uncharted in the public literature.

Public-health teams, reference laboratories, and phylodynamics researchers need an open pipeline that runs on the sequences they have, on the hardware they have, in the subtypes they see. **SAWA fills that gap** by:

- Treating **LANL Table exports and GenBank Entrez pulls as first-class ingest sources** with a shared canonical schema (versioned, DVC-tracked).
- Calling resistance via the **Stanford HIVdb Sierra GraphQL API**, algorithm-version-pinned and locally cached so runs are fully reproducible — including offline.
- Aligning to **HXB2 K03455 pol with `mappy`** (~700 seq/s in-process) so a full East-Africa cohort of ~25 k sequences aligns in under a minute on a laptop. Docker MAFFT `--add --keeplength` and Biopython fallbacks are available for cases requiring MSA-quality alignment.
- Building **per-subtype phylogenies with IQ-TREE 2 (GTR+F+I+G4) and TreeTime** molecular-clock scaling, plus **transmission-cluster detection with tn93 + connected-components at the standard 0.015 threshold** — the same conventions HIV-TRACE reviewers expect. Cluster parity vs HIV-TRACE is validated on simulated ground truth.
- Emitting **every figure from a Marimo notebook** under `notebooks/`, so every claim in every paper maps to a code path.
- Backing **every scalability claim with a benchmark script** under `benchmarks/` whose JSON results are checked into git and refreshed weekly by CI.

The immediate deliverable is a tools paper (*Paper 1, Virus Evolution*, in preparation). Three sibling preprints — causal identification (`causal/`), foundation-model resistance benchmark (`foundation_models/`), KEMRI deployment case study (`deployments/kemri/`) — live in the same repository and cite the shared library.

---

## Installation

**Prerequisites**

- [Docker](https://docs.docker.com/get-docker/) + Docker Compose v2 (for tn93, IQ-TREE, HIV-TRACE, MAFFT containers, and the compose stack)
- [uv](https://docs.astral.sh/uv/getting-started/installation/) (Python dependency manager)
- Python 3.12+ (installed by `uv sync`)

**Bootstrap**

```bash
git clone https://github.com/Neoforge-Labs/SAWA
cd SAWA
make setup           # copies .env.example → .env, runs `uv sync`
make images          # builds tn93, MAFFT, IQ-TREE, HIV-TRACE images (~5 min first time)
make up              # docker compose up -d — starts postgres, minio, prefect, mlflow
```

At the end you should see:

| Service | URL | Purpose |
|---|---|---|
| Prefect UI | http://localhost:4200 | Flow runs, live task logs, retries |
| MLflow UI | http://localhost:5000 | Per-run parameters, metrics, artifacts |
| MinIO console | http://localhost:9001 | S3-compatible artifact + DVC remote |
| SAWA API | http://localhost:8080 | Read-only public endpoints + `/docs` |

Verify:

```bash
uv run pytest -q                # 337 passing
uv run ruff check .             # clean
curl http://localhost:8080/health
curl http://localhost:8080/version
```

**Uninstall / reset:** `make down` stops services; `docker compose down -v` also drops volumes.

---

## Quick start — reproduce a paper

```bash
# 1. Ingest (LANL only, offline — no NCBI credentials needed)
make ingest

#    With GenBank incremental pull (needs NCBI_EMAIL in .env)
NCBI_EMAIL=you@example.com make ingest-genbank

# 2. Stanford HIVdb Sierra resistance calls (cached; second run is instant)
make sierra

# 3. Align to HXB2 pol, trim, QC, dedup
make curate

# 4. Transmission clusters (tn93 Docker + connected components)
make cluster

# 5. Per-subtype maximum-likelihood + time-scaled trees
make phylo

# 6. Cluster-detection validation vs HIV-TRACE on simulated ground truth
make cluster-validate WITH_HIVTRACE=1

# 7. Populate performance benchmarks
make benchmarks

# 8. Regenerate all five Paper 1 figures from raw data
make paper1
open paper1_tools/figures/*.pdf
```

To also load the pipeline outputs into Postgres (optional; the API then returns real data):

```bash
TO_POSTGRES=1 make ingest sierra cluster
curl "http://localhost:8080/sequences?country=KE&subtype=A1&limit=5"
```

---

## Overview

The core Python package (`src/hivdrwatch/`) is organized by pipeline stage. Every module has its own `README.md` explaining inputs, outputs, and tests.

| Module | Purpose |
|---|---|
| [`ingest/`](src/hivdrwatch/ingest/README.md) | LANL Table + FASTA parser, GenBank Entrez fetcher with resume-safe batching + per-country watermarks. Canonical `SequenceRecord` schema. |
| [`resistance/`](src/hivdrwatch/resistance/README.md) | Stanford HIVdb Sierra GraphQL client. SHA256-keyed local cache. Algorithm version pinned in every response. |
| [`align/`](src/hivdrwatch/align/README.md) | Pairwise-to-HXB2 via `mappy` (default); Docker MAFFT `--add --keeplength` and Biopython fallbacks. Trim, QC, patient-hash dedup. |
| [`phylo/`](src/hivdrwatch/phylo/README.md) | IQ-TREE 2 (GTR+F+I+G4) + TreeTime. Per-subtype orchestrator. |
| [`cluster/`](src/hivdrwatch/cluster/README.md) | tn93 + connected-components at 0.015 threshold. JC69 ground-truth simulator + pair-based P/R/F1/ARI validator + HIV-TRACE parity via `hivnetworkcsv`. |
| [`pipeline/`](src/hivdrwatch/pipeline/README.md) | Prefect 3 flows chaining every stage. Live subprocess streaming to the Prefect UI. Optional Postgres storage tail. |
| [`storage/`](src/hivdrwatch/storage/README.md) | SQLAlchemy 2.x schema + Alembic migrations + dialect-aware idempotent loader + DuckDB analytics wrapper. |
| [`api/`](src/hivdrwatch/api/README.md) | Read-only FastAPI service — `/sequences`, `/dr-calls`, `/clusters`, `/health`, `/version`. Auto-OpenAPI at `/docs`. Slowapi rate-limiting, CORS, request-id middleware. |
| [`cli/`](src/hivdrwatch/cli/README.md) | Unified `hivdrwatch` Typer entry point wrapping every stage. |
| [`paper1/`](src/hivdrwatch/paper1/README.md) | Pure-Python aggregators consumed by Fig 1-5 notebooks. Kept out of the notebooks so figures are unit-testable. |
| [`causal/`](causal/README.md) | Substrate for Paper 2: `DRDag`, `OutbreakDag`, DoWhy identify → estimate → refute, feature pruning, HTML report. |

**Cross-cutting:**

- `benchmarks/` — scale harness (`scale_ingest`, `scale_align`, `scale_cluster`). JSON checked into git, refreshed weekly by CI, feeds Fig 5.
- `.github/workflows/` — `test.yml` (ruff + mypy + pytest on every PR), `benchmarks.yml` (weekly refresh with auto-commit), `figures.yml` (nightly figure regeneration).
- `docker/` — pinned batch container images for MAFFT, tn93, IQ-TREE + TreeTime, HIV-TRACE.
- `data/reference/` — canonical fixtures (HXB2 pol, Stanford DRM position table).
- `data/{raw,interim,processed,gold}/` — DVC-tracked; `.dvc` pointers in git, files on MinIO/S3 remote.

See [`DEPLOYMENT.md`](DEPLOYMENT.md) for the operator runbook (reviewer path + production path, secrets handling, backups, rolling back).

---

## Documentation

- **Reviewer/reproducer path:** this README + [`DEPLOYMENT.md`](DEPLOYMENT.md).
- **Per-module contracts:** each subpackage has a `README.md` (linked in the table above).
- **Project scope + roadmap:** [`CONTEXT.md`](CONTEXT.md) (why we're building this), [`MODULES.md`](MODULES.md) (module-by-module definition of done).
- **Working conventions:** [`CLAUDE.md`](CLAUDE.md) — commits, code style, secrets, DSA gating.
- **Preregistration:** [`paper1_tools/preregistration_fig04.md`](paper1_tools/preregistration_fig04.md) — DTG-era INSTI case study, locked before running the analysis.
- **API reference:** live at `http://localhost:8080/docs` once `make up` is running; OpenAPI 3.1 spec at `/openapi.json`.

---

## Testing

**337 tests, one skip, ruff-clean, format-clean.**

```bash
uv run pytest -q                # full suite with coverage
uv run pytest tests/ingest -q   # single module
uv run ruff check .             # linter
uv run ruff format --check .    # formatter
```

Test types (per `CLAUDE.md`):

- **Unit tests** — pure helpers, run on every commit.
- **Integration tests** — cached fixtures for LANL/GenBank/Sierra under `tests/fixtures/`; frozen responses so runs are deterministic and offline-friendly.
- **Golden-file tests** — parser outputs checked into the repo; regenerate deliberately with `--update-golden`.
- **Benchmarks** — under `benchmarks/`, results JSON committed, regressions block merge.

Coverage on any module touched by a PR is required to stay above 70 %.

---

## Community guidelines

- **Report a bug or request a feature:** open an issue at https://github.com/Neoforge-Labs/SAWA/issues. Please include the SAWA version (`hivdrwatch version`), Python + OS, and a minimal reproduction.
- **Contribute code:** fork, branch (`feat/<module>/<short-desc>` or `paper/<N>/<section>`), open a draft PR early. Rebase on `main` before requesting review. One reviewer minimum for anything under `src/hivdrwatch/` core or `causal/` core; docs-only PRs can self-merge.
- **Code style:** Conventional Commits, type hints everywhere (`mypy --strict`), `ruff` clean. See `CLAUDE.md` for the full working conventions.
- **Security:** report vulnerabilities to `security@neoforge.io` — please **do not** open a public issue for security matters.
- **KEMRI data / PHI:** completely untouchable until the Data Sharing Agreement is signed. Code paths that could ever run against KEMRI data have a `require_dsa()` guard at entry.

---

## Papers

SAWA is the substrate for four planned first-author preprints plus a JOSS software citation. Papers 2–4 live in sibling packages (`causal/`, `foundation_models/`, `deployments/kemri/`) and import from the shared library.

| # | Working title | Venue | Status |
|---|---|---|---|
| 1 | **SAWA** — tools paper | Virus Evolution | code-complete, drafting |
| 2 | **CausalPhylo** — DAG + DoWhy for HIV phylodynamics | PLOS Comp Biol | Q3 2027 |
| 3 | **Foundation-model resistance benchmark** on non-B subtypes | Bioinformatics / NeurIPS D&B | Q4 2027 |
| 4 | **KEMRI deployment case study** (co-authored) | Bioinformatics App Note | Q1 2028 |
| 5 | **SAWA** — software citation | Journal of Open Source Software | continuous |

---

## Citation

If you use SAWA in research, please cite it. On tagged releases a Zenodo DOI is issued automatically — please cite the DOI. GitHub renders a "Cite this repository" button (top-right of the repo page) that reads `CITATION.cff`.

```bibtex
@software{sawa_2026,
  author       = {Njeri, Kelyn},
  title        = {{SAWA}: an open real-time antiretroviral-resistance
                   surveillance platform for sub-Saharan Africa},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.1.0-sawa},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://github.com/Neoforge-Labs/SAWA}
}
```

## Acknowledgements

SAWA is built with and depends on: [Stanford HIVdb Sierra](https://hivdb.stanford.edu/), [LANL HIV Sequence Database](https://www.hiv.lanl.gov/), [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/), [IQ-TREE 2](http://www.iqtree.org/), [TreeTime](https://treetime.readthedocs.io/), [HIV-TRACE](https://github.com/veg/hivtrace), [tn93](https://github.com/veg/tn93), [DoWhy](https://www.pywhy.org/dowhy/), [Prefect 3](https://www.prefect.io/), [FastAPI](https://fastapi.tiangolo.com/), [SQLAlchemy](https://www.sqlalchemy.org/), [Polars](https://pola.rs/), [Marimo](https://marimo.io/), and many more open-source projects. Domain framing shaped by early conversations with the KEMRI-Wellcome Trust Research Programme.

## Authors

- **Kelyn Njeri** — NeoForge Labs, Nairobi, Kenya. Lead maintainer, corresponding author. `kelyn@neoforge.io` · [ORCID pending]

Co-authors will be added on manuscript submission per the authorship contract in `MODULES.md`.

---

## License

MIT — see [`LICENSE`](LICENSE). Third-party licenses of bundled/dependency components are the responsibility of the user; SAWA does not redistribute their code.
