# SAWA (HIV-DR-Watch) — Project Context

*SAWA = Sub-Saharan Africa Watch for Antiretroviral resistance. The
Python package + repo name is `hivdrwatch`; the public/brand name is
SAWA. Use SAWA in manuscripts, talks, and outward-facing docs; use
`hivdrwatch` in code, imports, and CLI invocations.*

*A briefing document for Claude Code and any other collaborators joining
the build. Read this before writing code, opening a PR, or making any
architectural decision.*

---

## What we're building

**HIV-DR-Watch** is an open, real-time genomic surveillance platform for
HIV-1 drug resistance in sub-Saharan Africa, starting with East Africa
(Kenya first, then Uganda, Tanzania, Rwanda, Burundi, Ethiopia). It
ingests public HIV *pol* sequences (LANL, GenBank) plus — later —
partner clinical data (KEMRI), calls resistance mutations against the
Stanford HIVdb gold standard, builds subtype-stratified phylogenies,
detects transmission clusters, and exposes everything through a public
dashboard that African researchers and public health teams can actually
use.

The v1 unit of value is a **tools paper + running service**: someone can
`git clone` the repo, `docker compose up`, and reproduce every figure in
the paper on their laptop; separately they can visit the hosted dashboard
and see the current state of HIV drug resistance in East Africa.

This repo is also the **substrate for four related first-author preprints
plus a JOSS software citation** (see "The papers" below). The core
platform ships first; subsequent papers live in sibling subpackages that
import from the shared library and cite the tools paper.

## The papers

Five publishable artifacts come out of this repo, sequenced so each
depends on the one before. Paper 1 must ship before any other is
publicly staged, but implementation can proceed in parallel where the
substrate is ready.

### Paper 1 — HIV-DR-Watch (tools paper) — **Q1 2027**

**Venue:** Virus Evolution (Oxford). Fallback: Bioinformatics.
**Working title:** *"HIV-DR-Watch: an open real-time drug-resistance
surveillance platform for sub-Saharan Africa."*
**Thesis:** existing HIV drug-resistance surveillance is either
country-specific, closed-source, or optimized for subtype B — this tool
is open, sub-Saharan-Africa-first, and subtype-A1/C/D-native.

Five anchor figures:
1. Pipeline overview (ingest → subtype → resistance-call → align → tree → cluster → dashboard)
2. Resistance landscape heatmap (country × year × drug class, subtype-stratified)
3. Cluster validation vs HIV-TRACE on simulated ground truth
4. **Case study: DTG resistance signals** — first detection of emerging INSTI resistance in East African A1/D following the 2019 WHO dolutegravir rollout. Biological headline.
5. Scalability (wall-clock vs sequence count, log-log)

Everything in `src/hivdrwatch/` should ultimately feed one of these five.

### Paper 2 — CausalPhylo (methods paper) — **Q3 2027**

**Venue:** PLOS Computational Biology or Epidemics.
**Thesis:** causal identification (DAG + DoWhy + feature pruning) of
which interventions — ART regimen changes, DTG rollout, adherence
programs — actually *caused* observed resistance-frequency shifts, vs
what's explained by sampling bias and subtype confounding.

Lives in `causal/`. Code is already implemented — the work here is
applying it rigorously to HIV phylodynamics data from Paper 1,
constructing a defensible DAG with domain input from Chohan and Leigh
Brown, and writing the manuscript. Primary target audience: Leigh
Brown.

### Paper 3 — Foundation-model resistance benchmark — **Q4 2027**

**Venue:** Bioinformatics, or NeurIPS/ICML datasets & benchmarks track.
**Thesis:** foundation models (ESM-2, ESM-C, Evo) fine-tuned on HIV
*pol* underperform on non-B subtypes; subtype-stratified fine-tuning
closes the gap. This is the AI4BI-native paper — reviewers there care
exactly about this failure mode.

Lives in `foundation_models/`. Uses the HIV-DR-Watch pipeline as the
data substrate (that's why Paper 1 ships first).

### Paper 4 — KEMRI deployment (co-authored) — **Q1 2028**

**Venue:** Bioinformatics Applications Note, or JAMIA Open.
**Thesis:** deploying HIV-DR-Watch inside a real African reference lab
— what changed in the workflow, what broke, what we'd redesign.

Co-authored with Chohan; target senior author is one of the Edinburgh
supervisors. Technical novelty is modest; the collaboration is the
point. This is the paper that most cheaply converts you from applicant
to recruit. Lives in `deployments/kemri/`.

### Paper 5 — JOSS software citation — **continuous**

Submit to the Journal of Open Source Software once Paper 1 is on
bioRxiv. Two-page peer-reviewed software citation. Cheap, high-value,
adds a second DOI without competing for the same reviewer attention as
Paper 1.

## Repo shape (monorepo, one Zenodo DOI lineage)

```
hiv-dr-watch/
├── src/hivdrwatch/          # core platform — Paper 1
├── causal/                  # DAG + DoWhy work — Paper 2
├── foundation_models/       # ESM-2/Evo benchmark — Paper 3
├── deployments/kemri/       # KEMRI-specific configs — Paper 4
├── paper1_tools/            # manuscript, figures, refs
├── paper2_causal/
├── paper3_foundation/
├── paper4_kemri/
└── benchmarks/              # CI-tracked evaluation harness (shared)
```

Each `paperN_*/manuscript.qmd` pulls figures from notebooks and imports
from the shared library. One repo, one Zenodo DOI per release, four
preprints that all cite the same software.

## What "done" looks like

A reviewer at Virus Evolution can:

- Clone the repo and run `docker compose up` on a laptop and get a working stack
- Run `make paper` and regenerate every figure from raw data
- Read the manuscript in `paper/` alongside the code that produced it
- Cite the software via the Zenodo DOI (auto-issued on each GitHub release)
- Trust the numbers because every claim has a benchmark script in `benchmarks/` that runs in CI

The repo layout, tooling choices (Pixi/uv, DVC, Nextflow, GitHub Actions,
Zenodo integration, JOSS-shaped README), and the "no figure lives outside
`notebooks/`" rule all exist to make this bar reachable.

## Data reality (updated 2026-08-11 after real LANL download)

- **LANL is a *full* sequence source**, not curation-only. The Table TSV
  download contains sequences directly, plus HXB2 positional coords, plus
  curated subtype. The FASTA download is pre-aligned.
- **GenBank is a supplementary top-up** — used for accessions deposited
  after the last LANL refresh (LANL updates monthly), or as cross-validation.
- **KEMRI (via Dr. Chohan)** is planned for Phase 2, only after a data-sharing agreement is signed. Do not build anything that assumes KEMRI data is available today.
- Kenya alone gives us ~9,200 pol sequences 1994–2023 — enough for a real phylogenetic analysis and the v1 preprint.

## Boundaries — what NOT to build

- **Don't reinvent subtypers.** Use COMET + SCUEAL (peer-reviewed, containerized). Rolling our own is a two-year rabbit hole.
- **Don't reinvent tree visualization.** Embed Auspice (Nextstrain). Using it signals community awareness.
- **Keep `src/hivdrwatch/` classical** — Stanford HIVdb rules + phylogenetics. Foundation-model work belongs in `foundation_models/` (Paper 3), causal identification in `causal/` (Paper 2). The core platform must remain defensible without invoking either.
- **Don't add features that don't map to a paper figure or a benchmark.** Scope creep is the enemy of shipping preprints. Every new module names its target paper.
- **Don't touch KEMRI data or PHI in code paths that could ever run without the DSA.** Assume everything in the repo will be scrutinized.
- **Don't spread beyond the five papers listed.** Ebola/Marburg, mirror-life, embodied cognition, and general causal-RL work belong in *different repos* — they'd dilute the AI4BI narrative and confuse the Zenodo lineage.

## Why this matters (the PhD application angle)

The real goal is to be *recruited* into the UKRI AI4BI CDT at Edinburgh
(September 2029 entry) rather than merely applying. Target supervisors:
Dr. Andrea Weisse, Prof. Katie Atkins, Prof. Andrew Leigh Brown
(Edinburgh) + Dr. Bhavna Chohan (KEMRI Nairobi, wet-lab co-supervisor).

This repo produces **four of the five planned preprints** (Papers 1–4
above, plus the JOSS citation). Together they read as *"this candidate
is already doing the thesis, we just need to fund him."* The fifth
external preprint — an Ebola/Marburg spillover model reusing the phylo
stack — lives in a separate repo (`viral-spillover-watch`) so the
narrative here stays focused on HIV.

**The single highest-leverage move** is a co-authored paper with Chohan
*and* one Edinburgh supervisor before January 2028. HIV-DR-Watch is the
substrate that makes that collaboration natural to propose.

## Voice and audience

The paper is written for a **phylodynamics / public-health-genomics
reader who is skeptical of AI hype**. It's a tools paper first, methods
paper second. It should read as calm, competent, and infrastructure-flavored
— not as an ML demo. The dashboard is written for **African clinical and
public-health researchers**, so the UI is bandwidth-conscious, works
offline where possible, and never assumes desktop-class hardware or
familiarity with Nextstrain's conventions.

The lead author is a Kenyan engineer building tooling for the region he
lives in. That framing is not decoration — it changes what "usable"
means and what the case study is *about*.

---

**Owner:** Kelyn (CTO, NeoForge Labs) · **Preprint roadmap:** Paper 1 Q1 2027 · Paper 2 Q3 2027 · Paper 3 Q4 2027 · Paper 4 Q1 2028 · JOSS continuous · **Schema version:** see `src/hivdrwatch/ingest/schema.py::SCHEMA_VERSION`
