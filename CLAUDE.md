# CLAUDE.md — Working Conventions

*Instructions for Claude Code (or any AI collaborator) working in this
repo. Read at the start of every session. Applies alongside CONTEXT.md
(why we're building) and MODULES.md (what specifically must exist).*

---

## Session start

Every session, in order:

1. Read `CONTEXT.md` — check whether the papers list, venues, or
   supervisors have changed
2. Read `MODULES.md` — check which modules are shipped (✓) vs scaffolded
   (~) vs to-build (·) so you know what already exists
3. Run `git log --oneline -20` to see recent activity
4. Run `git status` to see uncommitted work
5. Only then start on the requested task

If CONTEXT.md and MODULES.md contradict each other, **MODULES.md wins on
scope, CONTEXT.md wins on framing.** Flag the drift to Kelyn.

## Commits and branches

**Commit style:** Conventional Commits.

```
feat(ingest): add LANL FASTA parser with header cross-validation
fix(resistance): handle Sierra timeout with exponential backoff
docs(context): update Paper 2 venue to PLOS CompBio
test(cluster): add FAVITES ground-truth harness
chore(deps): pin biopython to 1.84
refactor(phylo): extract subtype-reference cache into own module
```

Scopes match top-level module names (`ingest`, `resistance`, `causal`,
`fm`, `deployments`, `paper1`, etc.).

**Commit size:** small, atomic, individually revertable. If a commit
touches more than 3 modules, split it. If the message needs "and" more
than once, split it.

**Branches:**
- `main` is always releasable — never commit directly
- Feature branches: `feat/<module>/<short-desc>` (e.g. `feat/resistance/sierra-client`)
- Paper branches: `paper/<N>/<section>` (e.g. `paper/1/methods-section`)
- Never force-push shared branches; force-push your own feature branches freely

**PRs:** open early as draft. Rebase, don't merge, before marking ready.
One reviewer minimum (Kelyn) for anything touching `src/hivdrwatch/`
core or `causal/` core. Docs-only PRs can self-merge.

## Code style

**Python 3.12+ only.** Use `from __future__ import annotations` in every
module. Type hints are non-optional — mypy strict mode is enforced in CI.

**Naming:**
- Modules: `snake_case`, singular (`ingest`, not `ingesters`)
- Classes: `PascalCase`, no `I` prefixes for interfaces
- Functions: `snake_case`, verbs (`fetch_country`, not `country_fetcher`)
- Private helpers: prefix with `_`
- Constants: `SCREAMING_SNAKE_CASE`, module-level
- Type aliases: `PascalCase`, prefer `type X = ...` over old-style

**Error handling:**
- Raise specific exceptions, not `Exception` or `RuntimeError`
- Catch narrowly — `except OSError:` not `except:`
- One `except Exception as e:` per module maximum, and only in
  ingestion loops where one bad record must not stop the run — always
  log with `logger.exception(...)` in that case
- Never `pass` in an except block

**Logging:**
- `logger = logging.getLogger(__name__)` at module top
- Never `print()` in library code — only in CLI entry points
- INFO for lifecycle events (started, completed, N records processed)
- WARNING for degraded but recoverable state (schema fallback, rate
  limit approached)
- ERROR for anything that lost data or failed a claim
- Never log PHI or KEMRI patient data at any level

**Docstrings:** module-level docstring required; function docstrings
required when the signature doesn't explain the intent. Prefer clear
type hints over prose. No boilerplate `Args:/Returns:/Raises:` blocks
unless the function is public API.

**Comments:** explain *why*, never *what*. If the code isn't
self-explanatory, rewrite the code before adding a comment. Never write
`# TODO`, `# FIXME`, or `# HACK` — open a GitHub issue instead. Never
apologize in comments (`# this is ugly but...`).

## Data and schemas

**Raw data never enters git.** DVC-tracked only. `.gitignore` blocks
`data/raw/`, `data/interim/`, `data/processed/` — do not remove those
lines.

**KEMRI data and PHI:** untouchable until the Data Sharing Agreement is
signed and Kelyn confirms in writing. Code paths that could ever run
against KEMRI data must have a `require_dsa()` guard at entry. If in
doubt, don't touch it.

**Bumping `SCHEMA_VERSION`** (in `src/hivdrwatch/ingest/schema.py`):

| Change type | Version bump | Example |
|---|---|---|
| New optional field with default | Patch (0.2.0 → 0.2.1) | Add `pubmed_id: str \| None = None` |
| New required field, or changed field semantics | Minor (0.2.0 → 0.3.0) | Add `is_aligned: bool` field with no default |
| Removed field, renamed field, changed type in incompatible way | Major (0.2.0 → 1.0.0) | Rename `patient_id_hashed` → `patient_code` |

Every bump requires: (1) a migration note in the schema module
docstring, (2) a regenerated `data/interim/` sample checked into DVC,
(3) a CHANGELOG entry.

**Never invent data.** Placeholder sequences, mock patients, or
synthetic "example" records must live under `tests/fixtures/` and be
named `_synthetic_*` so they can never be mistaken for real data.

## Testing

**Test before PR.** `pytest tests/ -q` must pass. Coverage on any module
you touched must be ≥ 70%.

**Test types:**
- Unit tests for pure helpers — fast, run on every commit
- Integration tests using cached fixtures (real LANL/GenBank responses
  frozen under `tests/fixtures/`) — run in CI, not locally by default
- Golden-file tests for parsers — the fixture's expected output is
  checked into the repo; regenerate deliberately with a `--update-golden` flag
- Benchmarks under `benchmarks/` — run in CI, results JSON committed,
  regressions block merge

**No test-only branches in production code.** If a test needs a
back-door, that's a signal to refactor the code, not to add
`if os.getenv("TESTING"):`.

**Determinism:** every test that uses randomness must set an explicit
seed. Every benchmark that uses randomness must record the seed in its
results JSON.

## Reproducibility

**No figure lives outside `notebooks/`.** Every figure in every paper
is generated by a Marimo notebook (not Jupyter). The notebook is the
source of truth; the PDF/PNG is a build artifact.

**Deterministic seeds always.** Any stochastic pipeline step sets a
seed derived from the pipeline's config hash. Two runs with the same
config produce bit-identical outputs.

**Choosing an orchestrator:**
- **Nextflow** — paper analysis pipelines (`workflows/`). Reviewers
  expect it. Follow nf-core conventions.
- **Prefect 3** — the live surveillance service (`src/hivdrwatch/pipeline/`).
  Scheduled, monitored, has a UI.
- **Plain Python script** — one-off exploratory work in `notebooks/` or
  `scripts/`. Never for anything that runs in production.

Don't mix them. Don't wrap a Prefect flow in a Nextflow process, or
vice versa.

## Documentation

**When adding a module, update MODULES.md in the same PR.** Add the
row, set the status flag, name the paper it blocks. PRs that add code
without updating MODULES.md get bounced.

**When changing scope, venue, or authorship, update CONTEXT.md.** Rare
event. Kelyn approves.

**Per-module README required.** Every top-level module directory
(`src/hivdrwatch/ingest/`, `causal/`, `foundation_models/`, etc.) has a
README.md explaining: what the module does, how to run it, its inputs
and outputs, its tests.

**Never regenerate docs from scratch.** Edit in place. A regenerated
README loses the domain knowledge that accumulates in prose.

## Dependencies

**Pin everything.** Use `pixi` or `uv` — lockfiles are checked in.
Never `pip install` at runtime.

**No new dependencies without justification.** Before adding a
top-level dep:
1. Can this be done with stdlib? (Usually yes for small utilities.)
2. Is there already a dep in the tree that does this?
3. Is the new dep actively maintained (last release < 12 months)?
4. Does it have a license compatible with MIT? (No GPL, no AGPL.)

If yes to (1) or (2), don't add it. If you must add it, note the
justification in the commit message.

**Bio tools (MAFFT, IQ-TREE, HIV-TRACE) are installed via
conda-forge/bioconda in the Pixi env**, never system-installed. The
container reproduces the env; no host deps beyond Docker.

## Manuscripts

**Don't touch `paperN_*/manuscript.qmd` once the preprint is
submitted** — a `paperN-vX.Y.Z` git tag marks the frozen state. Fixes
go on a `paper/N/revisions` branch and merge back after acceptance.

**Tag releases per paper.** When Paper 1 goes to bioRxiv, tag
`paper1-v1.0.0`. The Zenodo DOI is issued off that tag and cited in
the manuscript.

**Every claim in every paper has a code path.** A statement like *"the
pipeline handles 100k sequences in under 4 hours"* must map to a
benchmark in `benchmarks/` whose result JSON matches the number in the
paper. No numbers in manuscripts that aren't in a results file.

## Secrets and configuration

**Env vars only.** `.env.example` is checked in with dummy values.
`.env` is gitignored. Never commit real credentials, API keys, or
KEMRI connection strings.

**Configuration lives in `config/` as YAML.** Not JSON, not TOML at
the config level. Pydantic-Settings loads and validates.

**NCBI, Sierra, and Zenodo API keys** go in env vars, never in code
or config files.

## When to stop and ask Kelyn

Do not silently proceed on any of these — surface and wait:

- **Any code path that could touch KEMRI data** before Kelyn confirms the DSA is signed
- **Schema-breaking changes** (major `SCHEMA_VERSION` bump)
- **Adding a new preprint** beyond the four in MODULES.md
- **Changing a paper's target venue** or authorship
- **Anything that costs money** beyond negligible cloud compute (large training runs, paid APIs, hosted services)
- **Any change that would touch code inside `paperN_*/` for a paper already submitted**
- **Bringing on a new collaborator or dependency** on a specific human's availability

## Don'ts (quick reference)

- Don't refactor code unrelated to the task
- Don't reformat files you didn't otherwise change
- Don't add `# TODO` / `# FIXME` comments — open issues
- Don't apologize in code, commits, or PRs
- Don't ship placeholder or mock data outside `tests/fixtures/`
- Don't rebuild something that already exists in `causal/` (Kelyn's already implemented the core)
- Don't reinvent COMET, SCUEAL, HIV-TRACE, IQ-TREE, TreeTime, or Auspice
- Don't add a subtype-B-only benchmark or default (the whole point is non-B)
- Don't commit raw data, secrets, or `.env` files
- Don't force-push shared branches
- Don't merge your own PRs on core modules
- Don't touch `main` directly

---

**Owner:** Kelyn · **Last updated:** 2026-08-11 · **When in doubt, ask.**
