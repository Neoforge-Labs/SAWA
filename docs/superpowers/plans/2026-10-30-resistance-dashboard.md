# Resistance Dashboard (Clean Axes) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A static, publicly hosted page rendering SAWA's Stanford HIVdb resistance calls on axes that are defensible against the current data vintage.

**Architecture:** A pure-Polars export package reads `data/interim/dr_results.parquet`, deduplicates it, aggregates four views, and writes a deterministic `summary.json`. A dependency-free HTML/CSS/ES-module page renders inline SVG from that committed JSON. GitHub Actions publishes the static directory to Pages and fails if the recorded input hash drifts from the committed DVC pointer.

**Tech Stack:** Python 3.12, Polars, stdlib `json`/`re`/`hashlib`, Typer (CLI dispatch only), pytest. No JS framework, no bundler, no npm.

**Spec:** `docs/superpowers/specs/2026-10-30-resistance-dashboard-design.md`

## Global Constraints

- Python 3.12+; `from __future__ import annotations` at the top of every module.
- mypy strict passes; `ruff check` and `ruff format --check` pass. Line length 95, lint rules `E,F,W,I,N,UP,B,SIM,RUF`.
- `logger = logging.getLogger(__name__)` per module; no `print()` outside `__main__`.
- Raise specific exceptions, never bare `Exception`/`RuntimeError`. Never `pass` in an except block.
- No `# TODO`, `# FIXME`, `# HACK` comments. No apologies in code or commits.
- No new third-party dependencies. Where stdlib suffices, use stdlib (the DVC pointer is parsed with `re`, not PyYAML, which is not in the tree).
- Module docstring required on every module; per-module `README.md` required.
- Synthetic test data lives under `tests/fixtures/` and is named `_synthetic_*`.
- Coverage ≥ 70% on `src/hivdrwatch/dashboard/`.
- Conventional Commits, scope `dashboard`. Work on the current feature branch; never commit to `main`.
- **No region claim anywhere in page copy or JSON** (spec §14.1). The cohort is described as "public HIV-1 *pol* sequences retrieved by an East-Africa-targeted LANL/GenBank query".
- **DRM panels cover PR, RT_NRTI and RT_NNRTI only; the `IN` block is excluded** (spec §14.2), and the page states that integrase DRMs are not shown because the vendored list predates the DTG era.
- Reuse `paper1.fig02_aggregate.RESISTANT_LEVELS` and `paper1.style.{COLORBLIND_PALETTE,DRUG_CLASS_ORDER}`; do not restate thresholds, palettes or orderings.

## Review Focus

Input classes the spec implies but which no task's happy-path tests exercise. Each has a test added to the owning task.

1. **A vintage missing a drug class entirely** (no INSTI rows at all) — prevalence must omit the class, not divide by zero. Owner: Task 3.
2. **A sequence scored for a gene but carrying an empty mutation list** (fully susceptible) — must still count in the DRM denominator, not vanish via `explode`. Owner: Task 3.
3. **Mutation strings that are pure ambiguity or malformed** (`X`, `""`, `T215`) — the matcher must reject them without raising. Owner: Task 2.
4. **Float formatting drifting across platforms** — two runs on different machines must produce byte-identical JSON. Owner: Task 4.
5. **A `summary.json` whose `schema_version` the page does not recognise** — the page must say so visibly, not render blank charts. Owner: Task 5.

---

### Task 1: Deduplication and structural invariants

**Files:**
- Create: `src/hivdrwatch/dashboard/__init__.py`, `src/hivdrwatch/dashboard/export/__init__.py`, `src/hivdrwatch/dashboard/export/dedup.py`
- Test: `tests/dashboard/test_dedup.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `class DedupReport` (frozen dataclass): `n_rows_in: int`, `n_rows_out: int`, `n_redundant: int`
  - `def dedup_dr_results(df: pl.DataFrame) -> tuple[pl.DataFrame, DedupReport]`
  - `class NonRedundantDuplicateError(ValueError)`
  - `class AmbiguousSequenceHashError(ValueError)`

- [ ] **Step 1: Write the failing tests**

```python
def test_dedup_removes_exact_duplicate_rows():
    df = _synthetic_dr_frame(rows=[_row("AB1", "3TC"), _row("AB1", "3TC"), _row("AB1", "AZT")])
    out, report = dedup_dr_results(df)
    assert out.height == 2
    assert report == DedupReport(n_rows_in=3, n_rows_out=2, n_redundant=1)


def test_dedup_raises_when_duplicates_carry_different_values():
    rows = [_row("AB1", "3TC", level="high"), _row("AB1", "3TC", level="susceptible")]
    with pytest.raises(NonRedundantDuplicateError):
        dedup_dr_results(_synthetic_dr_frame(rows=rows))


def test_dedup_raises_when_accession_has_two_sequence_hashes():
    rows = [_row("AB1", "3TC", sha="a" * 64), _row("AB1", "AZT", sha="b" * 64)]
    with pytest.raises(AmbiguousSequenceHashError):
        dedup_dr_results(_synthetic_dr_frame(rows=rows))


def test_dedup_is_idempotent():
    df = _synthetic_dr_frame(rows=[_row("AB1", "3TC"), _row("AB1", "3TC")])
    once, _ = dedup_dr_results(df)
    twice, report = dedup_dr_results(once)
    assert twice.equals(once)
    assert report.n_redundant == 0
```

Define `_row(accession, drug, *, level="susceptible", sha="a"*64, ...)` and `_synthetic_dr_frame(rows)` as local helpers building a frame with the `DRResult` columns.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/dashboard/test_dedup.py -v`
Expected: FAIL — `ModuleNotFoundError: hivdrwatch.dashboard`

- [ ] **Step 3: Implement `dedup_dr_results` in `src/hivdrwatch/dashboard/export/dedup.py`**

Drop fully identical rows with `df.unique()`. Then assert losslessness: if `df.select(["accession", "drug"]).n_unique()` differs from the deduped height, duplicates disagreed on some column — raise `NonRedundantDuplicateError` naming the count. Separately, if any accession maps to more than one `sequence_sha256`, raise `AmbiguousSequenceHashError`. Log the report at INFO.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/dashboard/test_dedup.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add src/hivdrwatch/dashboard tests/dashboard/test_dedup.py
git commit -m "feat(dashboard): dedup Sierra results and enforce structural invariants"
```

---

### Task 2: Surveillance DRM matcher

**Files:**
- Create: `src/hivdrwatch/dashboard/export/drm.py`
- Test: `tests/dashboard/test_drm.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `SDRM_PANELS: tuple[str, ...] = ("PR", "RT_NRTI", "RT_NNRTI")`
  - `class SdrmPanel` (frozen dataclass): `name: str`, `gene: str`, `positions: dict[int, frozenset[str]]`
  - `def load_sdrm_panels(path: Path) -> tuple[tuple[SdrmPanel, ...], str]` — returns the panels and the pinned list version
  - `def parse_mutation(raw: str) -> tuple[int, frozenset[str]] | None` — position and observed states, or `None` if unparseable
  - `def match_sdrm(raw: str, panel: SdrmPanel) -> str | None` — canonical `"<ref><pos><state>"` label if `raw` carries an SDRM in `panel`, else `None`

- [ ] **Step 1: Write the failing tests**

The reference file keys genes as `PR`, `RT_NRTI`, `RT_NNRTI`, `IN`; `IN` is excluded per the global constraints. `RT_NRTI` position 69 carries an `ins` state.

```python
def test_loads_only_pr_and_rt_panels():
    panels, version = load_sdrm_panels(REFERENCE_PATH)
    assert tuple(p.name for p in panels) == ("PR", "RT_NRTI", "RT_NNRTI")
    assert all(p.gene in {"PR", "RT"} for p in panels)
    assert version == "SDRM-2019"


@pytest.mark.parametrize("raw", ["R41K", "M36I", "V35T", "K122E", "D177E"])
def test_known_polymorphisms_are_not_sdrms(raw):
    panels, _ = load_sdrm_panels(REFERENCE_PATH)
    assert all(match_sdrm(raw, p) is None for p in panels)


@pytest.mark.parametrize(
    ("raw", "panel_name", "expected"),
    [("M184V", "RT_NRTI", "M184V"), ("K103N", "RT_NNRTI", "K103N"), ("M46I", "PR", "M46I")],
)
def test_known_drms_match_their_panel(raw, panel_name, expected):
    panel = _panel(panel_name)
    assert match_sdrm(raw, panel) == expected


def test_mixture_counts_as_carrying_any_matching_state():
    # Surveillance convention (spec §7.2): any state in the mixture counts.
    assert match_sdrm("M184MV", _panel("RT_NRTI")) == "M184V"
    assert match_sdrm("T215TFILMS", _panel("RT_NRTI")) is not None


def test_mixture_with_no_matching_state_does_not_match():
    # RT_NRTI 184 states are {V, I}; T is neither.
    assert match_sdrm("M184MT", _panel("RT_NRTI")) is None


def test_both_sdrm_states_at_a_position_match():
    # 184 lists V and I: M184I is an SDRM in its own right, not a near-miss.
    assert match_sdrm("M184I", _panel("RT_NRTI")) == "M184I"


def test_stop_codons_and_deletions_are_excluded():
    # Both sit at real SDRM positions (RT_NRTI 184 and 69), so these assert the
    # exclusion rule rather than a position miss.
    assert match_sdrm("M184*", _panel("RT_NRTI")) is None
    assert match_sdrm("T69del", _panel("RT_NRTI")) is None


def test_rt69_insertion_matches_the_ins_state():
    assert match_sdrm("T69S_SS", _panel("RT_NRTI")) == "T69ins"


@pytest.mark.parametrize("raw", ["X", "", "T215", "...", "69"])
def test_malformed_strings_are_rejected_without_raising(raw):
    assert parse_mutation(raw) is None
    assert match_sdrm(raw, _panel("RT_NRTI")) is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/dashboard/test_drm.py -v`
Expected: FAIL — `ImportError: cannot import name 'load_sdrm_panels'`

- [ ] **Step 3: Implement the matcher in `src/hivdrwatch/dashboard/export/drm.py`**

`load_sdrm_panels` reads the JSON, keeps only `SDRM_PANELS`, maps each panel to its parquet `gene` (`PR`→`PR`, both `RT_*`→`RT`), and returns `algorithm_version_pinned` alongside. Raise a specific error if a named panel is absent from the file.

`parse_mutation` matches `^([A-Z])(\d+)([A-Z]+)$` for substitutions, and the underscore insertion form `^([A-Z])(\d+)[A-Z]_[A-Z]+$` → states `{"ins"}`. Anything containing `*` or `del`, and anything not matching, returns `None`. The states set is every residue after the position, so `T215TFILMS` yields `{T,F,I,L,M,S}`.

`match_sdrm` returns the canonical label built from the reference residue, position, and the first matching state in the panel's sorted state order, so the label is deterministic when a mixture matches more than one. **The reference residue comes from the input string's leading letter** — the JSON stores only positions and their resistant states, no consensus residues, so do not look for them there.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/dashboard/test_drm.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/hivdrwatch/dashboard/export/drm.py tests/dashboard/test_drm.py
git commit -m "feat(dashboard): match surveillance DRMs with mixture-aware parsing"
```

---

### Task 3: View aggregations

**Files:**
- Create: `src/hivdrwatch/dashboard/export/aggregate.py`
- Test: `tests/dashboard/test_aggregate.py`

**Interfaces:**
- Consumes: `dedup_dr_results` (Task 1); `load_sdrm_panels`, `match_sdrm`, `SdrmPanel` (Task 2).
- Produces, each taking a deduped frame and returning a `pl.DataFrame`:
  - `def prevalence_by_drug_class(df) -> pl.DataFrame` → `[drug_class, n_sequences, n_resistant, frac_resistant]`
  - `def levels_by_drug(df) -> pl.DataFrame` → `[drug_class, drug, level, n_sequences, frac]`
  - `def gene_coverage(df) -> pl.DataFrame` → `[gene, n_sequences]`
  - `def sdrm_frequency(df, panels: tuple[SdrmPanel, ...]) -> pl.DataFrame` → `[panel, mutation, n_sequences, n_panel_denominator, frac]`

- [ ] **Step 1: Write the failing tests**

```python
def test_prevalence_counts_sequence_resistant_if_any_drug_in_class_is():
    # One accession: high to NVP, susceptible to EFV -> NNRTI-resistant once.
    out = prevalence_by_drug_class(_frame_two_nnrti_calls())
    row = out.filter(pl.col("drug_class") == "NNRTI").row(0, named=True)
    assert row["n_sequences"] == 1 and row["n_resistant"] == 1


def test_prevalence_uses_resistant_levels_from_fig02():
    from hivdrwatch.paper1.fig02_aggregate import RESISTANT_LEVELS

    assert RESISTANT_LEVELS == frozenset({"high", "intermediate"})
    out = prevalence_by_drug_class(_frame_with_level("low"))
    assert out.filter(pl.col("drug_class") == "NNRTI").row(0, named=True)["n_resistant"] == 0


def test_prevalence_omits_a_drug_class_absent_from_the_vintage():
    out = prevalence_by_drug_class(_frame_without_insti())
    assert "INSTI" not in out["drug_class"].to_list()
    assert out["frac_resistant"].is_not_nan().all()


def test_denominators_are_per_class_not_global():
    # PR scored for 2 accessions, INSTI for 1.
    out = prevalence_by_drug_class(_frame_uneven_gene_coverage())
    by = {r["drug_class"]: r["n_sequences"] for r in out.iter_rows(named=True)}
    assert by["PI"] == 2 and by["INSTI"] == 1


def test_levels_by_drug_covers_all_five_levels_per_drug():
    out = levels_by_drug(_frame_mixed_levels())
    assert set(out.filter(pl.col("drug") == "NVP")["level"]) == {
        "susceptible",
        "potential_low",
        "low",
        "intermediate",
        "high",
    }


def test_gene_coverage_counts_distinct_accessions_per_gene():
    assert dict(gene_coverage(_frame_uneven_gene_coverage()).iter_rows()) == {"PR": 2, "IN": 1}


def test_sdrm_frequency_counts_each_accession_once_per_mutation():
    # Same accession, M184V listed on both its RT rows.
    out = sdrm_frequency(_frame_repeated_mutation(), _panels())
    assert out.filter(pl.col("mutation") == "M184V").row(0, named=True)["n_sequences"] == 1


def test_sdrm_denominator_includes_sequences_with_no_mutations():
    # Two RT accessions: one carries M184V, one has an empty mutation list.
    out = sdrm_frequency(_frame_one_drm_one_clean(), _panels())
    row = out.filter(pl.col("mutation") == "M184V").row(0, named=True)
    assert row["n_sequences"] == 1
    assert row["n_panel_denominator"] == 2
    assert row["frac"] == pytest.approx(0.5)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/dashboard/test_aggregate.py -v`
Expected: FAIL — `ImportError: cannot import name 'prevalence_by_drug_class'`

- [ ] **Step 3: Implement the four aggregations in `src/hivdrwatch/dashboard/export/aggregate.py`**

Follow the shape of `paper1/fig02_aggregate.py`: validate required columns up front and raise `ValueError` naming the missing ones; no plotting imports. `prevalence_by_drug_class` reduces to per-`(accession, drug_class)` "resistant to any drug in class" exactly as fig02 does, importing `RESISTANT_LEVELS` rather than restating it. `sdrm_frequency` must take its denominator from the distinct accessions scored for the panel's gene **before** exploding the mutation lists, so sequences with no mutations stay in it. Order output rows deterministically (`drug_class` by `style.DRUG_CLASS_ORDER`, then the frame's natural keys).

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/dashboard/test_aggregate.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/hivdrwatch/dashboard/export/aggregate.py tests/dashboard/test_aggregate.py
git commit -m "feat(dashboard): aggregate prevalence, per-drug levels, coverage and DRMs"
```

---

### Task 4: Deterministic `summary.json` and export CLI

**Files:**
- Create: `src/hivdrwatch/dashboard/export/summary.py`, `src/hivdrwatch/dashboard/export/__main__.py`
- Test: `tests/dashboard/test_summary.py`
- Create fixture: `tests/fixtures/dashboard/_synthetic_dr_results.parquet` (built by the test module's own helper and committed)
- Create golden: `tests/fixtures/dashboard/_synthetic_summary.json`

**Interfaces:**
- Consumes: Tasks 1-3.
- Produces:
  - `SUMMARY_SCHEMA_VERSION: str = "1.0.0"`
  - `def input_dvc_md5(pointer: Path) -> str`
  - `def build_summary(df: pl.DataFrame, *, panels, sdrm_version: str, input_md5: str, git_sha: str, n_accessions_ingested: int) -> dict`
  - `def write_summary(summary: dict, out: Path) -> None`
  - `def main(argv: list[str] | None = None) -> int`

- [ ] **Step 1: Write the failing tests**

```python
def test_input_dvc_md5_reads_the_pointer():
    assert input_dvc_md5(POINTER) == "e51cbb39d5ebf3bdcaff05d4df52b515"


def test_meta_records_provenance_without_asserting_data_values():
    s = build_summary(
        _synthetic(),
        panels=_panels(),
        sdrm_version="SDRM-2019",
        input_md5="abc",
        git_sha="deadbee",
        n_accessions_ingested=4,
    )
    meta = s["meta"]
    assert meta["schema_version"] == SUMMARY_SCHEMA_VERSION
    assert meta["algorithm_version"] == "HIVDB_10.2"
    assert meta["algorithm_publish_date"] == "2026-04-26"
    assert meta["called_at"] == "2026-08-17"
    assert meta["sdrm_list_version"] == "SDRM-2019"
    assert meta["input_dvc_md5"] == "abc"
    assert meta["n_accessions_scored"] == 3
    assert meta["n_accessions_ingested"] == 4
    assert meta["n_calls"] == 7


def test_summary_carries_no_region_claim():
    blob = json.dumps(build_summary(_synthetic(), ...))
    for banned in ("East Africa", "east africa", "Kenya", "Uganda", "Tanzania"):
        assert banned not in blob


def test_summary_excludes_integrase_drms():
    s = build_summary(_synthetic(), ...)
    assert {v["panel"] for v in s["views"]["sdrm_frequency"]} <= {"PR", "RT_NRTI", "RT_NNRTI"}


def test_write_summary_is_byte_identical_across_runs(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    s = build_summary(_synthetic(), ...)
    write_summary(s, a)
    write_summary(build_summary(_synthetic(), ...), b)
    assert a.read_bytes() == b.read_bytes()


def test_fractions_are_rounded_to_four_places():
    s = build_summary(_synthetic(), ...)
    for v in s["views"]["prevalence_by_drug_class"]:
        assert v["frac_resistant"] == round(v["frac_resistant"], 4)


def test_golden_summary_matches(golden_update):
    s = build_summary(_load_fixture_parquet(), ...)
    _assert_golden(s, GOLDEN_PATH, update=golden_update)
```

Add a `--update-golden` pytest flag via a `conftest.py` fixture in `tests/dashboard/`, matching the repo's parser-test convention.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/dashboard/test_summary.py -v`
Expected: FAIL — `ImportError: cannot import name 'build_summary'`

- [ ] **Step 3: Implement `summary.py`**

`input_dvc_md5` extracts the first `md5:` value from the pointer text with `re`; raise a specific error if absent. `build_summary` runs dedup, then the four aggregations, and returns `{"meta": {...}, "views": {...}}` where every view is a list of row dicts. `meta` records — never asserts — the measured counts, plus the algorithm version, `algorithm_publish_date`, `called_at`, SDRM list version, input md5, git SHA, `DR_SCHEMA_VERSION` and `SUMMARY_SCHEMA_VERSION`. Raise a specific error if `algorithm_version` or `called_at` is not single-valued across the frame. `write_summary` uses `json.dumps(..., sort_keys=True, indent=2, ensure_ascii=False)` with a trailing newline, and rounds every float to 4 places before serialising.

- [ ] **Step 4: Implement `__main__.py`**

`main(argv)` takes `--input` (default `data/interim/dr_results.parquet`), `--pointer` (default the input path plus `.dvc`), `--reference` (default `data/reference/drm_positions.json`), `--out` (default the package's `static/summary.json`), and `--ingested-count`. Resolve the git SHA with `subprocess.run(["git", "rev-parse", "HEAD"])`, falling back to `"unknown"` on `OSError` or a non-zero return. Return `0` on success. `print()` is allowed here and nowhere else.

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/dashboard -v`
Expected: PASS (all tasks' tests)

- [ ] **Step 6: Generate the real summary and eyeball it**

Run: `uv run python -m hivdrwatch.dashboard.export --ingested-count 36510`
Expected: writes `src/hivdrwatch/dashboard/static/summary.json`; `meta.n_accessions_scored` is 35777, `meta.n_calls` is 728087.

- [ ] **Step 7: Commit**

```bash
git add src/hivdrwatch/dashboard/export tests/dashboard tests/fixtures/dashboard \
        src/hivdrwatch/dashboard/static/summary.json
git commit -m "feat(dashboard): emit a deterministic resistance summary JSON"
```

---

### Task 5: Static page

**Files:**
- Create: `src/hivdrwatch/dashboard/static/index.html`, `styles.css`, `app.js`
- Test: `tests/dashboard/test_static.py`

**Interfaces:**
- Consumes: `summary.json` as written by Task 4; `paper1.style.COLORBLIND_PALETTE` and `DRUG_CLASS_ORDER` (values copied into `styles.css` as custom properties and into `app.js` as a frozen array — the page cannot import Python).
- Produces: nothing for later tasks.

- [ ] **Step 1: Write the failing tests**

Python-side structural tests; no browser, no JS runtime.

```python
def test_page_declares_the_schema_version_it_expects():
    assert f'data-summary-schema="{SUMMARY_SCHEMA_VERSION}"' in INDEX.read_text()


def test_page_handles_an_unrecognised_schema_version_visibly():
    js = APP_JS.read_text()
    assert "data-summary-schema" in js
    assert "schemaMismatch" in js


def test_page_makes_no_region_claim():
    text = INDEX.read_text()
    for banned in ("East Africa", "Kenya", "Uganda", "Tanzania", "sub-Saharan"):
        assert banned not in text


def test_page_describes_the_cohort_by_retrieval_process():
    assert "East-Africa-targeted" in INDEX.read_text()


def test_page_states_why_integrase_drms_are_absent():
    text = INDEX.read_text()
    assert "integrase" in text.lower() and "dolutegravir" in text.lower()


def test_page_states_the_mixture_counting_convention():
    assert "mixture" in INDEX.read_text().lower()


def test_page_links_the_provenance_evaluation():
    assert "ingest-provenance.md" in INDEX.read_text()


def test_page_loads_no_external_scripts_or_styles():
    text = INDEX.read_text()
    assert '<script src="http' not in text and "cdn" not in text.lower()


def test_footer_carries_the_citation_and_provenance_pins():
    text = INDEX.read_text()
    assert "CITATION.cff" in text
    assert "github.com/Neoforge-Labs/SAWA" in text
    assert "HIVDB_10.2" in text


def test_palette_matches_paper1_style():
    from hivdrwatch.paper1.style import COLORBLIND_PALETTE

    css = STYLES.read_text()
    for colour in COLORBLIND_PALETTE[:4]:
        assert colour.lower() in css.lower()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/dashboard/test_static.py -v`
Expected: FAIL — `FileNotFoundError: index.html`

- [ ] **Step 3: Build the page**

`index.html` carries the §8 sections in order: provenance band, prevalence by drug class, per-drug level breakdown, SDRM frequency, gene coverage, footer. It sets `data-summary-schema` on `<body>`. Copy in §8's exact wording decisions: the cohort line "35,777 of 36,510 ingested accessions", the retrieval-process description, the denominator definition in words, the integrase absence note, the mixture convention, and the limitation note linking `docs/evaluations/ingest-provenance.md`.

`app.js` is one ES module: `fetch("./summary.json")`, compare `meta.schema_version` against the `data-summary-schema` attribute and render a visible banner instead of charts when they differ, then draw each view as inline SVG. Bars are `<rect>` with `<text>` labels; every prevalence bar renders its own denominator. No external requests.

`styles.css` defines the palette as custom properties on `:root`, a dark-mode block under `@media (prefers-color-scheme: dark)`, and a single-column layout at phone width with a 16px gutter and no horizontal scroll.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/dashboard/test_static.py -v`
Expected: PASS

- [ ] **Step 5: Look at it**

Run: `uv run python -m http.server -d src/hivdrwatch/dashboard/static 8099`
Open `http://localhost:8099` and confirm all five sections render with real numbers, then check it at 375px wide.

- [ ] **Step 6: Commit**

```bash
git add src/hivdrwatch/dashboard/static tests/dashboard/test_static.py
git commit -m "feat(dashboard): render the resistance summary as a static page"
```

---

### Task 6: Pages deploy and input-drift guard

**Files:**
- Create: `.github/workflows/pages.yml`
- Modify: `.github/workflows/test.yml` (add one job), `src/hivdrwatch/dashboard/export/summary.py` (add `check_input_not_stale`)
- Test: `tests/dashboard/test_summary.py` (add the drift cases)

**Interfaces:**
- Consumes: `meta.input_dvc_md5` in `summary.json` (Task 4); `data/interim/dr_results.parquet.dvc`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Write the drift check as a script with a test**

```python
def test_drift_check_passes_when_hashes_agree(tmp_path): ...
def test_drift_check_fails_when_the_pointer_moved(tmp_path): ...
```

Implement `check_input_not_stale(summary: Path, pointer: Path) -> int` in `src/hivdrwatch/dashboard/export/summary.py`, returning `0` when `meta.input_dvc_md5` equals the pointer's md5 and `1` with a logged ERROR otherwise.

- [ ] **Step 2: Run tests to verify they fail, then pass**

Run: `uv run pytest tests/dashboard/test_summary.py -k drift -v`

- [ ] **Step 3: Add the drift job to `test.yml`**

A job running `uv run python -c "...check_input_not_stale..."` against the committed files. It needs no DVC data. Add a comment stating plainly that this detects a stale input only and does not verify the aggregation, because the sole DVC remote is a localhost MinIO unreachable from CI.

- [ ] **Step 4: Write `pages.yml`**

Trigger on push to `main` touching `src/hivdrwatch/dashboard/static/**`, plus `workflow_dispatch`. Permissions `pages: write`, `id-token: write`. Steps: `actions/checkout@v4`, `actions/configure-pages@v5`, `actions/upload-pages-artifact@v3` with `path: src/hivdrwatch/dashboard/static`, `actions/deploy-pages@v4`. No Python, no uv — it publishes committed files.

- [ ] **Step 5: Verify the workflows parse**

Run: `uv run ruff check . && uv run pytest tests/dashboard -q`
Expected: PASS. Note in the commit that Pages must be enabled once in repository settings (spec §14.3) before the first deploy succeeds.

- [ ] **Step 6: Commit**

```bash
git add .github/workflows/pages.yml .github/workflows/test.yml \
        src/hivdrwatch/dashboard/export/summary.py tests/dashboard/test_summary.py
git commit -m "ci(dashboard): publish the static page and guard against input drift"
```

---

### Task 7: Wiring and documentation

**Files:**
- Create: `src/hivdrwatch/dashboard/README.md`
- Modify: `src/hivdrwatch/cli/app.py` (`_SUBCOMMANDS` table), `Makefile`, `MODULES.md` (§1.1 table), `CHANGELOG.md`
- Test: `tests/cli/` (existing contract test picks up the new subcommand)

**Interfaces:**
- Consumes: `hivdrwatch.dashboard.export.__main__:main` (Task 4).
- Produces: nothing.

- [ ] **Step 1: Add the CLI entry**

One row in `_SUBCOMMANDS`: `"dashboard": ("hivdrwatch.dashboard.export.__main__", "Aggregate Sierra resistance calls into the static dashboard's summary JSON.")`.

- [ ] **Step 2: Run the CLI contract tests**

Run: `uv run pytest tests/cli -v`
Expected: PASS, with `dashboard` discovered.

- [ ] **Step 3: Add the Makefile target**

```make
.PHONY: dashboard
dashboard:
	python -m hivdrwatch.dashboard.export --input $(INTERIM_DIR)/dr_results.parquet --ingested-count $${INGESTED_COUNT:-36510}
```

Add the line to the `help` target.

- [ ] **Step 4: Write `src/hivdrwatch/dashboard/README.md`**

What it does, how to run it, inputs and outputs, tests — per CLAUDE.md. State the three scope facts a reader will ask about: no geographic or subtype axis (and why, linking the provenance evaluation), PR/RT DRMs only (and why), and the drift guard's real strength.

- [ ] **Step 5: Add the MODULES.md row**

In §1.1, status `~`, blocking Paper 1. Word it to cover only the static resistance summary and to state that Flow 1's API-backed Auspice + deck.gl dashboard is not started and that stratified views await the provenance rebuild. Do not edit Flow 1's acceptance criteria.

- [ ] **Step 6: Add the CHANGELOG entry**

- [ ] **Step 7: Run the full suite**

Run: `uv run ruff check . && uv run ruff format --check . && uv run mypy src/hivdrwatch/dashboard && uv run pytest tests/ -q`
Expected: PASS. Then `uv run pytest tests/dashboard --cov=src/hivdrwatch/dashboard --cov-report=term-missing` and confirm ≥ 70%.

- [ ] **Step 8: Commit**

```bash
git add src/hivdrwatch/dashboard/README.md src/hivdrwatch/cli/app.py Makefile \
        MODULES.md CHANGELOG.md
git commit -m "docs(dashboard): document the module and register it in the CLI and MODULES"
```
