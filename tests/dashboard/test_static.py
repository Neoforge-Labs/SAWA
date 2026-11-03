"""Structural tests for the static resistance page.

There is no browser and no JS runtime here, so these tests read the three
shipped files and pin what a reader would see. Three habits keep them honest:

- Prose claims are asserted against the text of the element that makes the claim
  (found by id), never against the whole file, so a word that survives in a comment
  or in another paragraph cannot satisfy a test about a different one.
- Sentences that carry vintage numbers are checked as the page renders them, by
  substituting the real `summary.json` into the `data-bind` slots the way `app.js`
  does. A mistyped binding key or a wrong wording fails; a regenerated vintage does not.
- Every list the page re-states from Python (class order, level order, the palette)
  is compared against the Python source rather than a second hand-typed copy.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass, field
from html.parser import HTMLParser
from itertools import pairwise
from pathlib import Path
from typing import Any

import pytest

from hivdrwatch.dashboard.export.aggregate import LEVEL_ORDER
from hivdrwatch.dashboard.export.summary import SUMMARY_SCHEMA_VERSION
from hivdrwatch.paper1.style import COLORBLIND_PALETTE, DRUG_CLASS_ORDER

REPO = Path(__file__).resolve().parents[2]
STATIC = REPO / "src" / "hivdrwatch" / "dashboard" / "static"
INDEX = STATIC / "index.html"
STYLES = STATIC / "styles.css"
APP_JS = STATIC / "app.js"
SUMMARY: dict[str, Any] = json.loads((STATIC / "summary.json").read_text())

_VOID = {"meta", "link", "br", "hr", "img", "input"}
_SVG_NAMESPACE = "http://www.w3.org/2000/svg"
_RETRIEVAL_PHRASE = "East-Africa-targeted"


# --- a minimal HTML tree, enough to ask "what text does this element show?" ---------


@dataclass
class _Node:
    tag: str
    attrs: dict[str, str]
    children: list[_Node | str] = field(default_factory=list)


class _TreeBuilder(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.root = _Node("#root", {})
        self._stack = [self.root]

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        node = _Node(tag, {k: v or "" for k, v in attrs})
        self._stack[-1].children.append(node)
        if tag not in _VOID:
            self._stack.append(node)

    def handle_endtag(self, tag: str) -> None:
        for i in range(len(self._stack) - 1, 0, -1):
            if self._stack[i].tag == tag:
                del self._stack[i:]
                return

    def handle_data(self, data: str) -> None:
        self._stack[-1].children.append(data)


def _parse(html: str) -> _Node:
    builder = _TreeBuilder()
    builder.feed(html)
    return builder.root


def _walk(node: _Node) -> Iterator[_Node]:
    yield node
    for child in node.children:
        if isinstance(child, _Node):
            yield from _walk(child)


def _by_id(root: _Node, element_id: str) -> _Node:
    matches = [n for n in _walk(root) if n.attrs.get("id") == element_id]
    assert len(matches) == 1, f"expected exactly one #{element_id}, found {len(matches)}"
    return matches[0]


def _share(part: int, whole: int) -> str:
    """Python twin of `fmtShare` in app.js."""
    return f"{part / whole * 100:.1f}%" if whole else "0.0%"


def _landscape_cell(cls: str, country: str, year: str) -> dict[str, Any]:
    (row,) = [
        r
        for r in SUMMARY["views"]["resistance_landscape"]
        if (r["drug_class"], r["country"], r["collection_year"]) == (cls, country, int(year))
    ]
    return row


def _names(names: list[str]) -> str:
    """Python twin of `listNames` (Intl.ListFormat, English, long conjunction)."""
    if len(names) <= 1:
        return "".join(names)
    if len(names) == 2:
        return f"{names[0]} and {names[1]}"
    return f"{', '.join(names[:-1])}, and {names[-1]}"


def _drawn_cell_sizes() -> list[int]:
    rows = SUMMARY["views"]["resistance_landscape"]
    return [r["n_sequences"] for r in rows if r["frac_resistant"] is not None]


def _median(values: list[int]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    return ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2


def _bound_value(path: str, attrs: dict[str, str]) -> str:
    """Python twin of `resolveBinding` + `applyBindings` in app.js."""
    scope, *rest = path.split(".")
    meta = SUMMARY["meta"]
    views = SUMMARY["views"]
    value: Any
    if scope == "meta":
        value = meta[rest[0]]
    elif scope == "gene":
        value = next(r["n_sequences"] for r in views["gene_coverage"] if r["gene"] == rest[0])
    elif scope == "cell":
        cls, country, year, field = rest
        row = _landscape_cell(cls, country, year)
        value = row["n_sequences"] if field == "n" else f"{row['frac_resistant'] * 100:.1f}%"
    elif scope == "ratio":
        cls, country_a, year_a, country_b, year_b = rest
        a = _landscape_cell(cls, country_a, year_a)["frac_resistant"]
        b = _landscape_cell(cls, country_b, year_b)["frac_resistant"]
        value = round(a / b)
    elif scope == "derived":
        absent = [
            c["name"]
            for c in meta["landscape_countries"]
            if not any(
                r["country"] == c["code"] and r["frac_resistant"] is not None
                for r in views["resistance_landscape"]
            )
        ]
        value = {
            "n_unscored": meta["n_accessions_ingested"] - meta["n_accessions_scored"],
            "n_drugs": len({r["drug"] for r in views["levels_by_drug"]}),
            "n_classes": len({r["drug_class"] for r in views["levels_by_drug"]}),
            "n_repeated_sequences": meta["n_accessions_scored"] - meta["n_distinct_sequences"],
            "pct_sdrm_multistate": _share(
                meta["n_sdrm_observations_multistate"], meta["n_sdrm_observations"]
            ),
            "landscape_max_cell_n": max(_drawn_cell_sizes()),
            "landscape_median_cell_n": _median(_drawn_cell_sizes()),
            "landscape_region_names": _names([c["name"] for c in meta["landscape_countries"]]),
            "landscape_absent_countries": _names(absent),
            "pct_landscape_undated": _share(
                meta["landscape_n_excluded_no_year"],
                meta["landscape_n_accessions"] + meta["landscape_n_excluded_no_year"],
            ),
        }[rest[0]]
    else:
        raise AssertionError(f"unknown binding scope {scope!r}")
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.0f}" if value == int(value) else f"{value:,}"
    text = str(value)
    return text[: int(attrs["data-trunc"])] if "data-trunc" in attrs else text


def _rendered_text(node: _Node) -> str:
    """The element's text as a reader sees it once app.js has filled the bindings."""
    if node.tag in {"script", "style", "noscript"}:
        return ""
    if "data-bind" in node.attrs:
        return _bound_value(node.attrs["data-bind"], node.attrs)
    return "".join(_rendered_text(c) if isinstance(c, _Node) else c for c in node.children)


def _raw_text(node: _Node) -> str:
    """Text as written in the file, before any binding is applied."""
    return "".join(_raw_text(c) if isinstance(c, _Node) else c for c in node.children)


def _squash(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _note(element_id: str) -> str:
    return _squash(_rendered_text(_by_id(PAGE, element_id)))


def _js_function(name: str) -> str:
    js = APP_JS.read_text()
    start = js.index(f"function {name}(")
    depth = 0
    for i in range(js.index("{", start), len(js)):
        depth += {"{": 1, "}": -1}.get(js[i], 0)
        if depth == 0:
            return js[start : i + 1]
    raise AssertionError(f"unbalanced braces in {name}")


def _js_frozen_array(name: str) -> list[str]:
    match = re.search(rf"const {name} = Object\.freeze\(\[([^\]]*)\]\)", APP_JS.read_text())
    assert match, f"{name} is not declared as a frozen array"
    return re.findall(r'"([^"]+)"', match.group(1))


def _js_frozen_map(name: str) -> dict[str, str]:
    match = re.search(rf"const {name} = Object\.freeze\(\{{([^}}]*)\}}\)", APP_JS.read_text())
    assert match, f"{name} is not declared as a frozen object"
    return dict(re.findall(r'(\w+):\s*"([^"]+)"', match.group(1)))


PAGE = _parse(INDEX.read_text())


# --- schema gate ----------------------------------------------------------------------


def test_body_declares_the_schema_version_summary_json_carries() -> None:
    body = next(n for n in _walk(PAGE) if n.tag == "body")
    assert body.attrs.get("data-summary-schema") == SUMMARY_SCHEMA_VERSION
    assert SUMMARY["meta"]["schema_version"] == SUMMARY_SCHEMA_VERSION


# The expected schema must come from the page's own attribute, read either way.
_ATTRIBUTE_READ = (
    r'document\.body\.(?:getAttribute\(\s*"data-summary-schema"\s*\)|dataset\.summarySchema)'
)


def _expected_identifier(main: str) -> str:
    match = re.search(rf"(?:const|let)\s+(\w+)\s*=\s*{_ATTRIBUTE_READ}\s*;", main)
    assert match, "the expected schema is not read from the body's data-summary-schema"
    return match.group(1)


def test_schema_gate_runs_before_anything_is_drawn_and_shows_a_banner() -> None:
    main = _js_function("main")
    gate = main.index("schemaMismatch(")
    assert gate < main.index("applyBindings(") < main.index("drawAll(")
    mismatch_branch = main[gate : main.index("return;", gate)]
    assert "showBanner(" in mismatch_branch
    _expected_identifier(main)


def test_schema_gate_decision_logic_is_pinned() -> None:
    js = APP_JS.read_text()
    # The decision: a version that is not identical to the expected one is a mismatch.
    # `return false` here would draw charts against any schema. Parameter names, and extra
    # guards joined by `||`, are free to change; the inequality of the two inputs is not.
    decision = re.search(
        r"function schemaMismatch\((\w+),\s*(\w+)\)\s*\{\s*"
        r"return\s+(?:[^;]*?\|\|\s*)?(\w+)\s*!==\s*(\w+)\s*;\s*\}",
        js,
    )
    assert decision, "schemaMismatch no longer ends in an inequality of its two arguments"
    params, compared = decision.group(1, 2), decision.group(3, 4)
    assert sorted(params) == sorted(compared) and params[0] != params[1]

    main = _js_function("main")
    # Polarity at the call site: a leading `!` would show the banner on a match. The
    # version compared must be the summary's `schema_version`, not another version field.
    expected = _expected_identifier(main)
    call = re.search(
        rf"if\s*\(\s*schemaMismatch\(\s*{expected}\s*,[^()]*\bschema_version\b[^()]*\)\s*\)\s*\{{"
        r"\s*showBanner\([^;]*\);\s*return;\s*\}",
        main,
    )
    assert call, (
        "the gate is not `if (schemaMismatch(expected, <meta.schema_version>)) {banner; return}`"
    )
    # The expectation is never a literal in the script.
    assert not re.search(rf"\b{expected}\s*=\s*[\"'`]", js)


def test_the_banner_is_unhidden_by_script_and_not_defeated_by_css() -> None:
    banner = _by_id(PAGE, "schema-banner")
    assert banner.attrs.get("role") == "alert" and "hidden" in banner.attrs
    assert "banner.hidden = false" in _js_function("showBanner")
    # An author `display:` rule outranks the UA's [hidden] rule; this one must win back.
    assert re.search(r"\[hidden\]\s*\{\s*display:\s*none\s*!important", STYLES.read_text())


def test_a_failed_gate_withholds_bound_sentences_as_well_as_charts() -> None:
    body = _js_function("showBanner")
    assert '"main > section, [data-bind]"' in body and ".hidden = true" in body


# --- the page does not depend on what it was told not to --------------------------------


def test_page_does_not_read_the_row_count_field_under_review() -> None:
    for path in (INDEX, APP_JS):
        assert "n_rows_raw" not in path.read_text()


def test_every_key_the_page_reads_exists_in_summary_json() -> None:
    js = APP_JS.read_text()
    for view in set(re.findall(r"views\.(\w+)", js)):
        assert view in SUMMARY["views"], view
    for key in set(re.findall(r"\bmeta\.(\w+)", js)):
        assert key in SUMMARY["meta"], key
    bound = [n.attrs["data-bind"] for n in _walk(PAGE) if "data-bind" in n.attrs]
    assert bound, "no data-bind slots: the page would show no vintage numbers"
    for path in bound:
        _bound_value(path, {"data-trunc": "1"})


# --- order constants re-stated from Python ----------------------------------------------


def test_drug_class_order_matches_the_manuscript_figures() -> None:
    assert tuple(_js_frozen_array("DRUG_CLASS_ORDER")) == DRUG_CLASS_ORDER
    assert [r["drug_class"] for r in SUMMARY["views"]["prevalence_by_drug_class"]] == list(
        DRUG_CLASS_ORDER
    )


def test_level_order_matches_the_export_and_each_level_has_a_style() -> None:
    assert tuple(_js_frozen_array("LEVEL_ORDER")) == LEVEL_ORDER
    css = STYLES.read_text()
    for level in LEVEL_ORDER:
        assert re.search(rf"\.lvl-{level}\s*\{{\s*fill:", css), level
        assert f"--level-{level}:" in css, level


def test_genes_are_presented_in_gene_order_and_integrase_has_no_mutation_panel() -> None:
    assert _js_frozen_array("GENE_ORDER") == ["PR", "RT", "IN"]
    assert _js_frozen_array("SDRM_GENES") == ["PR", "RT"]
    assert "GENE_ORDER.forEach" in _js_function("renderCoverage")


def test_prevalence_keeps_the_exported_order_and_shows_each_denominator() -> None:
    body = _js_function("renderPrevalence")
    assert ".sort(" not in body and "toSorted(" not in body
    assert "rows.forEach" in body
    assert "n_sequences" in body and "n_resistant" in body


def test_class_to_gene_map_agrees_with_the_denominators_the_data_exports() -> None:
    # The page prints "N of M sequences scored for <gene>" from CLASS_GENE, a hand-typed copy
    # of which class lives on which gene. A wrong entry would mislabel a denominator and
    # nothing else would fail, so each class's exported denominator must equal the coverage
    # of the gene the map names. Exact on this vintage, where every sequence with a gene
    # score is scored for every class on it; a vintage where that stops being true should
    # fail here and be looked at.
    class_gene = _js_frozen_map("CLASS_GENE")
    assert set(class_gene) == set(DRUG_CLASS_ORDER)
    coverage = {r["gene"]: r["n_sequences"] for r in SUMMARY["views"]["gene_coverage"]}
    for row in SUMMARY["views"]["prevalence_by_drug_class"]:
        gene = class_gene[row["drug_class"]]
        assert row["n_sequences"] == coverage[gene], (row["drug_class"], gene)


def test_sdrm_panel_ranks_by_sequence_count_not_by_the_exported_name_order() -> None:
    body = _js_function("topSdrm")
    assert re.search(r"\.sort\(\(a, b\) => b\.n_sequences - a\.n_sequences", body)
    assert "SDRM_TOP_N" in body
    counts = [r["n_sequences"] for r in SUMMARY["views"]["sdrm_frequency"]]
    assert counts != sorted(counts, reverse=True), "export order is no longer by name"


# --- copy decisions, each pinned to the element that carries it -------------------------


def test_sections_appear_in_the_designed_order() -> None:
    ids = [
        n.attrs["id"]
        for n in _walk(PAGE)
        if n.attrs.get("id")
        in {"provenance", "prevalence", "landscape", "levels", "sdrm", "coverage", "footer"}
    ]
    assert ids == [
        "provenance",
        "prevalence",
        "landscape",
        "levels",
        "sdrm",
        "coverage",
        "footer",
    ]


_REGION_WORDS = re.compile(
    r"africa|kenya|uganda|tanzania|rwanda|burundi|ethiopia|somalia|south sudan|saharan", re.I
)


def test_only_the_landscape_section_names_a_country_and_nothing_claims_the_cohort_is_regional() -> (
    None
):
    # Country names belong to the landscape view, which shows them by construction. Everywhere
    # else, and in the script and stylesheet, they stay out: the cohort is described by its
    # retrieval process, and the page must not call it East African, because 3,000-odd of its
    # accessions state a country elsewhere.
    landscape = _by_id(PAGE, "landscape")
    elsewhere = _squash(_raw_text(PAGE)).replace(_squash(_raw_text(landscape)), "")
    assert not _REGION_WORDS.search(elsewhere.replace(_RETRIEVAL_PHRASE, "")), elsewhere
    for path in (STYLES, APP_JS):
        text = path.read_text().replace(_RETRIEVAL_PHRASE, "")
        assert not _REGION_WORDS.search(text), f"{path.name}: {_REGION_WORDS.search(text)}"
    assert "country-name" not in APP_JS.read_text()


def test_cohort_is_described_by_its_retrieval_process() -> None:
    cohort = _note("cohort")
    assert f"retrieved by an {_RETRIEVAL_PHRASE} LANL/GenBank query" in cohort
    assert "That names the query, not the sequences." in cohort


def test_cohort_line_states_scored_of_ingested_and_the_shortfall() -> None:
    meta = SUMMARY["meta"]
    scored, ingested = meta["n_accessions_scored"], meta["n_accessions_ingested"]
    cohort = _note("cohort")
    assert f"{scored:,} of {ingested:,} ingested accessions returned a Sierra result" in cohort
    assert f"({ingested - scored:,} did not)" in cohort
    assert f"giving {meta['n_calls']:,} resistance calls across" in cohort


def test_prevalence_is_caveated_as_not_a_pretreatment_resistance_estimate() -> None:
    prevalence = _by_id(PAGE, "prevalence")
    note = _by_id(PAGE, "treatment-history")
    chart = next(n for n in _walk(prevalence) if n.attrs.get("data-chart") == "prevalence")
    in_page_order = [n for n in _walk(prevalence) if n is note or n is chart]
    assert in_page_order == [note, chart], (
        "the caveat must sit in the section and before the bars"
    )
    text = _note("treatment-history")
    assert text.startswith("Treatment history is not recorded for these sequences.")
    assert "not pretreatment-resistance estimates" in text
    assert "surveillance, treatment-failure and trial sequences in unknown proportion" in text
    assert not re.search(r"\d", text)


def test_nothing_else_on_the_page_presents_the_figures_as_pretreatment_or_drug_naive() -> None:
    # The two caveats are the only places the word appears, and both use it to deny the
    # reading. Everything else is scanned for it.
    caveats = [_note("treatment-history"), _note("landscape-cohort-note")]
    elsewhere = _squash(_rendered_text(PAGE))
    for caveat in caveats:
        elsewhere = elsewhere.replace(caveat, "")
    assert not re.search(
        r"pre-?treatment|drug[- ]naive|treatment[- ]naive|\bPDR\b", elsewhere, re.I
    )


def test_unit_is_stated_as_accessions_with_the_measured_count_of_repeated_sequences() -> None:
    meta = SUMMARY["meta"]
    scored, distinct = meta["n_accessions_scored"], meta["n_distinct_sequences"]
    assert 0 < distinct <= scored
    note = _note("unit-note")
    assert note.startswith(
        "Every prevalence, level and coverage figure on this page counts accessions"
    )
    assert "“sequences” means accessions" in note
    expected = (
        f"The {scored:,} scored accessions hold {distinct:,} distinct sequences, "
        f"so {scored - distinct:,} repeat a sequence already present under another accession."
    )
    assert expected in note
    assert "Patient identity is not recorded, so these are not independent patients." in note
    assert not re.search(r"\d", _squash(_raw_text(_by_id(PAGE, "unit-note")))), (
        "a count is typed into the unit note instead of bound"
    )


def test_denominator_is_defined_in_words_and_is_not_confused_with_codon_coverage() -> None:
    note = _note("denominator-note")
    assert "“Sierra returned a score for this gene”" in note
    assert "not the same as “the surveillance codons were sequenced”" in note
    assert "re-calling Sierra" in note
    integrase = SUMMARY["views"]["gene_coverage"]
    in_n = next(r["n_sequences"] for r in integrase if r["gene"] == "IN")
    assert f"Only {in_n:,} of {SUMMARY['meta']['n_accessions_scored']:,} sequences" in note


def test_integrase_absence_is_explained_in_the_mutation_panel_note() -> None:
    note = _note("sdrm-scope").lower()
    assert note.startswith("integrase mutations are not shown")
    assert "dolutegravir" in note
    assert "bennett" in note and "no integrase component" in note
    assert "r263k" in note and "g118r" in note and "e138k" in note


def test_mixture_convention_and_its_measured_limit_are_both_stated() -> None:
    note = _note("mixture-convention")
    assert (
        "a mixture counts as carrying a surveillance mutation when any of its states" in note
    )
    assert "attributed to the alphabetically first matching state" in note
    assert "That affects" in note and "of observations" in note
    # The direction is part of the disclosure: the earlier state is right, the later one short.
    assert "the earlier state’s count is correct" in note
    assert "the later state’s count is understated" in note
    assert "inflated" not in note, (
        "a mixture does carry the earlier state, so it is not inflated"
    )


def test_mixture_percentage_is_computed_from_the_summary_counts() -> None:
    meta = SUMMARY["meta"]
    multi, total = meta["n_sdrm_observations_multistate"], meta["n_sdrm_observations"]
    assert 0 < multi < total, "the vintage has no multi-state observations to describe"
    # Expected text is built here from the raw counts, independently of the binding twin.
    expected = (
        f"That affects {multi / total * 100:.1f}% of observations ({multi:,} of {total:,})"
    )
    assert expected in _note("mixture-convention")


def test_no_mixture_figure_is_typed_into_the_page() -> None:
    raw = _squash(_raw_text(_by_id(PAGE, "mixture-convention")))
    assert not re.search(r"\d", raw), f"digits typed into the convention note: {raw}"
    js = _js_function("fmtShare")
    assert "(part / whole) * 100).toFixed(1)" in js
    assert "fmtShare(meta.n_sdrm_observations_multistate, meta.n_sdrm_observations)" in (
        _js_function("deriveFigures")
    )


def test_limitation_note_links_the_provenance_evaluation_that_exists() -> None:
    limitation = _by_id(PAGE, "limitation")
    hrefs = [n.attrs["href"] for n in _walk(limitation) if n.tag == "a"]
    target = "docs/evaluations/ingest-provenance.md"
    assert any(h.endswith(f"/blob/main/{target}") for h in hrefs), hrefs
    assert (REPO / target).is_file()


def test_limitation_says_geography_is_shown_and_no_longer_says_it_is_absent() -> None:
    note = _note("limitation")
    assert note.startswith("Geography is now shown")
    assert "rebuilt" in note and "record" in note
    assert "source kept per row" in note
    # The old note's claim, which is false now that the vintage is rebuilt.
    assert "deliberately absent" not in note
    assert "return after a provenance rebuild" not in note


def test_limitation_keeps_the_history_of_the_two_defects_as_separate_things() -> None:
    note = _note("limitation")
    assert "wrong source" not in note, "the two fields were not written from one wrong source"
    sentences = re.split(r"(?<=\.) ", note)
    (country,) = [sent for sent in sentences if "country field" in sent]
    assert "search query rather than from the record" in country
    assert "misparsed" in country and "every CRF was stored as subtype C" in country
    # A reader who follows the link must find what the page says: the evaluation's section on
    # the subtype defect is about CRFs being stored as subtype C.
    document = (REPO / "docs/evaluations/ingest-provenance.md").read_text()
    assert re.search(r"^## 2\. .*every CRF was stored as subtype C", document, re.MULTILINE), (
        "the evaluation no longer describes the subtype defect the page cites"
    )
    assert re.search(r"^## 3\. .*the country was the search query", document, re.MULTILINE)


def test_subtype_breakdowns_are_said_to_remain_absent_and_the_reason_is_the_rollup_not_the_field() -> (
    None
):
    note = _note("limitation")
    assert "subtype breakdowns are not" in note.lower()
    assert "remain absent although the subtype field is now correct" in note
    assert "sums across subtypes" in note
    intro = _note("landscape-intro")
    assert "sums across subtypes" in intro
    assert "Figure 2 stratifies by subtype and this page does not" in intro


def test_the_old_out_of_region_citation_is_gone_and_the_live_count_is_bound_not_typed() -> (
    None
):
    # 3,108 was the evaluation's count on the earlier vintage. On this one the count is
    # measured, so a typed copy of the old figure would contradict the bound one beside it.
    for path in (INDEX, STYLES, APP_JS):
        assert "3,108" not in path.read_text(), path.name
    raw = _squash(_raw_text(_by_id(PAGE, "landscape-exclusions-note")))
    assert not re.search(r"\d", raw), f"a count is typed into the exclusions note: {raw}"


def test_footer_carries_repo_citation_algorithm_vintage_and_commit() -> None:
    footer = _by_id(PAGE, "footer")
    hrefs = {n.attrs["href"] for n in _walk(footer) if n.tag == "a" and "href" in n.attrs}
    assert "https://github.com/Neoforge-Labs/SAWA" in hrefs
    assert "https://github.com/Neoforge-Labs/SAWA/blob/main/CITATION.cff" in hrefs
    meta = SUMMARY["meta"]
    text = _squash(_rendered_text(footer))
    assert f"Algorithm {meta['algorithm_version']}" in text
    assert f"data vintage {meta['called_at']}" in text
    assert f"generated at commit {meta['git_sha'][:7]}" in text
    commit_link = next(n for n in _walk(footer) if n.attrs.get("data-bind") == "meta.git_sha")
    assert commit_link.attrs["data-href"] == (
        "https://github.com/Neoforge-Labs/SAWA/commit/{value}"
    )
    assert "CITATION.cff" in text


def test_footer_labels_the_pointer_digest_and_the_file_digest_separately() -> None:
    meta = SUMMARY["meta"]
    for key in ("input_dvc_md5", "input_file_md5", "curated_dvc_md5", "curated_file_md5"):
        assert re.fullmatch(r"[0-9a-f]{32}", meta[key]), key
    text = _squash(_rendered_text(_by_id(PAGE, "footer")))
    assert f"DVC pointer md5 {meta['input_dvc_md5'][:8]}" in text
    assert f"input file md5 {meta['input_file_md5'][:8]}" in text
    assert f"curated table pointer md5 {meta['curated_dvc_md5'][:8]}" in text
    assert f"curated table file md5 {meta['curated_file_md5'][:8]}" in text
    assert "input md5" not in text, "an unlabelled digest is ambiguous about what it hashes"


def test_digest_mismatch_notice_is_hidden_until_the_digests_differ() -> None:
    note = _by_id(PAGE, "digest-note")
    assert "hidden" in note.attrs, "the notice would show even when the digests agree"
    said = _squash(_raw_text(note))
    assert "These digests differ." in said
    assert "does not hash to the md5 its DVC pointer records" in said
    assert any(n is note for n in _walk(_by_id(PAGE, "footer"))), (
        "notice sits outside the footer"
    )
    # The branch that is written but rarely run: the comparison, the unhide, and the call.
    assert "digestsDiffer(" in _js_function("syncDigestNote")
    sync = _js_function("syncDigestNote")
    assert 'getElementById("digest-note").hidden = !digestsDiffer(meta)' in sync
    main = _js_function("main")
    assert main.index("syncDigestNote(summary.meta)") > main.index("applyBindings(")


# --- null handling, run in node because the browser is not available -----------------------

_NODE = shutil.which("node")
_HARNESS_TAIL = """
globalThis.__page = {
  resolveBinding, boundText, deriveFigures, digestsDiffer, landscapeModel, columnLayout, yearLabels,
  rampColor, textTone, renderLandscape, describeLandscapeCell, RAMP_STOPS,
};
"""

# A stand-in for the DOM, enough for the SVG builders: elements that remember what they were given.
_FAKE_DOM = """
class FakeEl {
  constructor(tag) { this.tag = tag; this.attrs = {}; this.children = []; this.text = null; this.handlers = {}; }
  setAttribute(name, value) { this.attrs[name] = value; }
  append(...kids) { this.children.push(...kids); }
  set textContent(value) { this.text = value; }
  get textContent() { return this.text; }
  addEventListener(type, fn) { this.handlers[type] = fn; }
}
globalThis.document = { createElementNS: (_ns, tag) => new FakeEl(tag) };
const walk = (el, out = []) => { out.push(el); for (const k of el.children) if (typeof k === "object") walk(k, out); return out; };
const cellsOf = (svg) => svg.children.filter((el) => el.attrs.class === "hm-cell").map((g) => ({
  rect: g.children.find((k) => k.tag === "rect").attrs,
  slash: g.children.some((k) => k.tag === "line"),
  n: g.children.find((k) => k.tag === "text")?.text ?? null,
  title: g.children.find((k) => k.tag === "title")?.text ?? null,
  click: Boolean(g.handlers.click),
}));
"""


def _run_page_functions(script: str) -> Any:
    """Load app.js as a module with its `main()` call removed, run `script`, return its JSON."""
    source, replaced = re.subn(r"^main\(\);\s*$", "", APP_JS.read_text(), flags=re.MULTILINE)
    assert replaced == 1, "app.js no longer ends in a bare main() call; update this harness"
    harness = (
        source + _HARNESS_TAIL + _FAKE_DOM + f"\nconst page = globalThis.__page;\n{script}\n"
    )
    done = subprocess.run(
        [_NODE or "node", "--input-type=module"],
        input=harness,
        capture_output=True,
        text=True,
        check=False,
    )
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


needs_node = pytest.mark.skipif(_NODE is None, reason="node is not installed")


@needs_node
def test_null_in_a_required_meta_value_is_rejected_rather_than_rendered() -> None:
    required = sorted(
        {
            n.attrs["data-bind"]
            for n in _walk(PAGE)
            if n.attrs.get("data-bind", "").startswith("meta.")
        }
        - {"meta.algorithm_publish_date"}
    )
    assert {"meta.algorithm_version", "meta.called_at"} <= set(required)
    outcome = _run_page_functions(
        f"""
        const summary = {json.dumps(SUMMARY)};
        const out = {{}};
        for (const path of {json.dumps(required)}) {{
          const broken = structuredClone(summary);
          broken.meta[path.split(".")[1]] = null;
          try {{ out[path] = "rendered " + String(page.resolveBinding(path, broken, {{}})); }}
          catch (error) {{ out[path] = "threw"; }}
        }}
        console.log(JSON.stringify(out));
        """
    )
    assert set(outcome.values()) == {"threw"}, outcome


@needs_node
def test_a_null_publish_date_reads_as_a_dash_and_never_as_the_word_null() -> None:
    outcome = _run_page_functions(
        f"""
        const summary = {json.dumps(SUMMARY)};
        summary.meta.algorithm_publish_date = null;
        const value = page.resolveBinding("meta.algorithm_publish_date", summary, {{}});
        console.log(JSON.stringify({{ value, text: page.boundText(value, undefined) }}));
        """
    )
    assert outcome["value"] is None
    assert outcome["text"] == "\u2014"
    assert "null" not in outcome["text"].lower()


@needs_node
def test_a_value_absent_from_the_summary_is_still_an_error() -> None:
    outcome = _run_page_functions(
        f"""
        const summary = {json.dumps(SUMMARY)};
        delete summary.meta.called_at;
        try {{ page.resolveBinding("meta.called_at", summary, {{}}); console.log('"returned"'); }}
        catch (error) {{ console.log('"threw"'); }}
        """
    )
    assert outcome == "threw"


@needs_node
def test_bound_text_formats_numbers_truncates_strings_and_leaves_other_strings_whole() -> None:
    outcome = _run_page_functions(
        """
        console.log(JSON.stringify([
          page.boundText(35777, undefined), page.boundText(0, undefined),
          page.boundText("0123456789abcdef", "7"), page.boundText("HIVDB_10.2", undefined),
        ]));
        """
    )
    assert outcome == ["35,777", "0", "0123456", "HIVDB_10.2"]


# --- no external requests, no framework --------------------------------------------------


def test_page_is_one_html_one_stylesheet_one_module_with_no_external_loads() -> None:
    nodes = list(_walk(PAGE))
    scripts = [n for n in nodes if n.tag == "script"]
    assert [(s.attrs.get("type"), s.attrs.get("src")) for s in scripts] == [
        ("module", "app.js")
    ]
    sheets = [n for n in nodes if n.tag == "link" and n.attrs.get("rel") == "stylesheet"]
    assert [s.attrs["href"] for s in sheets] == ["styles.css"]
    loading = {"script": "src", "link": "href", "img": "src", "source": "src", "iframe": "src"}
    for node in nodes:
        attr = loading.get(node.tag)
        if attr and attr in node.attrs:
            assert not re.match(r"(https?:)?//", node.attrs[attr]), node.attrs[attr]


def test_script_and_styles_make_no_external_requests() -> None:
    js, css = APP_JS.read_text(), STYLES.read_text()
    assert set(re.findall(r"https?://[^\s\"'`)]+", js)) == {_SVG_NAMESPACE}
    assert not re.findall(r"https?://|@import|@font-face|url\(", css)
    assert re.findall(r"\bfetch\(([^)]*)\)", js) == ["SUMMARY_URL"]
    assert 'SUMMARY_URL = "./summary.json"' in js
    assert not re.search(r"\bimport\s*[({]|\bfrom\s+[\"']|\bXMLHttpRequest|\bWebSocket\b", js)


def test_no_work_markers_in_shipped_files() -> None:
    for path in (INDEX, STYLES, APP_JS):
        assert not re.search(r"TODO|FIXME|HACK", path.read_text()), path.name


# --- appearance --------------------------------------------------------------------------


def _root_block(css: str) -> str:
    match = re.search(r":root\s*\{(.*?)\n\}", css, re.S)
    assert match
    return match.group(1)


def test_palette_is_okabe_ito_from_the_manuscript_style_in_order() -> None:
    root = _root_block(STYLES.read_text())
    declared = re.findall(r"--okabe-[a-z]+:\s*(#[0-9A-Fa-f]{6})", root)
    assert [c.lower() for c in declared] == [c.lower() for c in COLORBLIND_PALETTE]


def test_class_and_panel_colours_come_from_the_palette_not_new_hex_values() -> None:
    root = _root_block(STYLES.read_text())
    for name in (*DRUG_CLASS_ORDER, "PR", "RT_NRTI", "RT_NNRTI"):
        assert re.search(rf"--(class|panel)-{name}:\s*var\(--okabe-[a-z]+\)", root), name
    for level in LEVEL_ORDER[1:]:
        assert re.search(rf"--level-{level}:\s*var\(--okabe-[a-z]+\)", root), level
    for name, colour in zip(
        DRUG_CLASS_ORDER, ("blue", "orange", "green", "purple"), strict=True
    ):
        assert f"--class-{name}: var(--okabe-{colour})" in root


def test_dark_mode_redefines_the_surface_and_body_paints_an_explicit_background() -> None:
    css = STYLES.read_text()
    dark = re.search(r"@media \(prefers-color-scheme: dark\)\s*\{(.*?)\n\}", css, re.S)
    assert dark
    for token in ("--bg:", "--text:", "--track:"):
        assert token in dark.group(1), token
    assert re.search(r"\nbody\s*\{[^}]*background:\s*var\(--bg\)", css)
    assert re.search(r"\nbody\s*\{[^}]*color:\s*var\(--text\)", css)


def test_phone_layout_has_a_viewport_meta_and_a_sixteen_pixel_gutter() -> None:
    viewport = next(
        n for n in _walk(PAGE) if n.tag == "meta" and n.attrs.get("name") == "viewport"
    )
    assert "width=device-width" in viewport.attrs["content"]
    css = STYLES.read_text()
    assert "--gutter: 16px" in css
    assert re.search(r"main, footer, \.banner\s*\{[^}]*padding-inline:\s*var\(--gutter\)", css)
    assert re.search(r"\.chart svg\s*\{[^}]*width:\s*100%", css)


# --- the landscape heatmap ------------------------------------------------------------------

_LANDSCAPE = SUMMARY["views"]["resistance_landscape"]
_FLOOR = SUMMARY["meta"]["landscape_min_cell_n"]


def _rgb(hex_colour: str) -> tuple[int, int, int]:
    return int(hex_colour[1:3], 16), int(hex_colour[3:5], 16), int(hex_colour[5:7], 16)


def _luminance(hex_colour: str) -> float:
    def channel(c: int) -> float:
        v = c / 255
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4

    r, g, b = (channel(c) for c in _rgb(hex_colour))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def test_the_committed_landscape_view_is_internally_consistent() -> None:
    # Structural properties that hold for any vintage, so a regenerated file passes and a
    # corrupted one does not. The page refuses a cell where these disagree; this says the
    # shipped file is not one.
    targeted = [c["code"] for c in SUMMARY["meta"]["landscape_countries"]]
    assert _LANDSCAPE, "the landscape view is empty"
    for row in _LANDSCAPE:
        assert row["country"] in targeted and row["drug_class"] in DRUG_CLASS_ORDER
        assert (row["n_sequences"] < _FLOOR) == (row["frac_resistant"] is None), row
        assert 0 <= row["n_resistant"] <= row["n_sequences"], row
        if row["frac_resistant"] is not None:
            assert row["frac_resistant"] == pytest.approx(
                row["n_resistant"] / row["n_sequences"], abs=5e-5
            )
    keys = [(r["drug_class"], r["country"], r["collection_year"]) for r in _LANDSCAPE]
    assert len(keys) == len(set(keys)), "a (class, country, year) cell appears twice"


def test_the_landscape_accounting_in_meta_sums_to_the_scored_cohort() -> None:
    meta = SUMMARY["meta"]
    excluded = sum(v for k, v in meta.items() if k.startswith("landscape_n_excluded_"))
    assert meta["landscape_n_accessions"] + excluded == meta["n_accessions_scored"]


def test_the_page_has_one_heatmap_per_class_in_manuscript_order_each_with_its_own_readout() -> (
    None
):
    figure = _by_id(PAGE, "landscape-figure")
    charts = [n for n in _walk(figure) if n.attrs.get("data-chart") == "landscape"]
    assert [n.attrs["data-class"] for n in charts] == list(DRUG_CLASS_ORDER)
    for node in charts:
        cls = node.attrs["data-class"]
        assert node.attrs["id"] == f"chart-landscape-{cls}"
        _by_id(PAGE, f"readout-landscape-{cls}")
        _by_id(PAGE, f"omitted-landscape-{cls}")


def test_the_cohort_composition_caveat_sits_before_the_heatmaps_and_the_caption_points_back_to_it() -> (
    None
):
    section = _by_id(PAGE, "landscape")
    note = _by_id(PAGE, "landscape-cohort-note")
    figure = _by_id(PAGE, "landscape-figure")
    first_chart = next(n for n in _walk(section) if n.attrs.get("data-chart") == "landscape")
    order = [n for n in _walk(section) if n in (note, first_chart)]
    assert order == [note, first_chart], "the caveat must precede the first heatmap"
    assert any(n is note for n in _walk(section)) and not any(n is note for n in _walk(figure))
    text = _note("landscape-cohort-note")
    assert text.startswith(
        "Cell-to-cell differences may reflect cohort composition rather than epidemiology."
    )
    assert "records no treatment history" in text
    assert "not as a trend over time or a comparison between countries" in text
    caption = _note("landscape-caption")
    assert "cohort caveat above applies to every cell" in caption
    assert any(
        n.tag == "figcaption" and n.attrs.get("id") == "landscape-caption"
        for n in _walk(figure)
    )


def test_the_cohort_caveat_quotes_two_real_cells_and_the_spread_between_them_is_real() -> None:
    low = _landscape_cell("NNRTI", "KE", "2010")
    high = _landscape_cell("NNRTI", "UG", "2016")
    for row in (low, high):
        assert row["frac_resistant"] is not None and row["n_sequences"] >= 100
    # The sentence says "a roughly N-fold spread between two large cells": it must stay true
    # of whatever vintage is shipped, or the example no longer illustrates the point.
    fold = high["frac_resistant"] / low["frac_resistant"]
    assert fold >= 10
    text = _note("landscape-cohort-note")
    expected = (
        f"Kenya in 2010 shows {low['frac_resistant'] * 100:.1f}% of {low['n_sequences']:,} "
        f"sequences resistant to an NNRTI; Uganda in 2016 shows "
        f"{high['frac_resistant'] * 100:.1f}% of {high['n_sequences']:,}. "
        f"That is a roughly {round(fold)}-fold spread"
    )
    assert expected in text
    raw = _squash(_raw_text(_by_id(PAGE, "landscape-cohort-note")))
    undated = re.sub(
        r"\b(?:19|20)\d\d\b", "", raw
    )  # the years name the cells; they are not vintage figures
    assert not re.search(r"\d", undated), "a figure is typed into the caveat instead of bound"


def test_a_binding_that_cites_a_cell_under_the_floor_is_an_error_not_a_blank() -> None:
    assert "frac_resistant === null ? null" in _js_function("cellBinding")
    assert "disagrees with the floor" in _js_function("landscapeModel")


def test_the_exclusions_note_states_every_measured_count_and_why_the_out_of_region_ones_exist() -> (
    None
):
    meta = SUMMARY["meta"]
    text = _note("landscape-exclusions-note")
    assert (
        f"Of {meta['n_accessions_scored']:,} scored accessions, "
        f"{meta['landscape_n_accessions']:,} are drawn."
    ) in text
    assert (
        f"{meta['landscape_n_excluded_out_of_region']:,} state a country outside the targeted "
        "region"
    ) in text
    assert f"{meta['landscape_n_excluded_no_country']:,} state no country." in text
    assert (
        f"{meta['landscape_n_excluded_no_year']:,} are from the targeted region but record no "
        "collection year, so they have no column"
    ) in text
    assert (
        "the search over-matched" in text and "these records state their own country" in text
    )
    assert "The prevalence figures above include every scored accession" in text
    names = ", ".join(c["name"] for c in meta["landscape_countries"][:-1])
    assert f"({names}, and {meta['landscape_countries'][-1]['name']})" in text
    links = [
        n.attrs["href"]
        for n in _walk(_by_id(PAGE, "landscape-exclusions-note"))
        if n.tag == "a"
    ]
    assert any(h.endswith("docs/evaluations/ingest-provenance.md") for h in links)


def test_the_undated_share_is_of_the_regions_accessions_not_of_everything_scored() -> None:
    meta = SUMMARY["meta"]
    undated = meta["landscape_n_excluded_no_year"]
    in_region = meta["landscape_n_accessions"] + undated
    assert undated > 0, "this vintage has no undated accessions; the sentence needs rewording"
    assert f"({_share(undated, in_region)} of the region" in _note("landscape-exclusions-note")


def test_the_sampling_note_states_the_range_and_median_of_cell_sizes() -> None:
    sizes = _drawn_cell_sizes()
    text = _note("landscape-sampling-note")
    assert text.startswith("Sampling effort is not uniform.")
    assert f"between {_FLOOR} and {max(sizes):,} sequences" in text
    assert f"half of them hold {_median(sizes):,.0f} or fewer" in text
    assert min(sizes) >= _FLOOR


def test_the_absent_countries_sentence_is_hidden_by_script_when_there_are_none() -> None:
    sentence = _by_id(PAGE, "landscape-absent")
    assert "hidden" in sentence.attrs
    sync = _js_function("syncLandscapeAbsent")
    assert 'getElementById("landscape-absent").hidden = absent.length === 0' in sync
    assert "syncLandscapeAbsent(summary)" in _js_function("main")


def test_the_heatmap_styles_keep_not_enough_data_apart_from_the_ramp() -> None:
    css = STYLES.read_text()
    assert re.search(r"\.hm-nodata\s*\{\s*fill:\s*var\(--nodata\)", css)
    assert re.search(r"\.hm-slash\s*\{[^}]*stroke:", css)
    root = _root_block(css)
    nodata = re.search(r"--nodata:\s*(#[0-9a-fA-F]{6})", root)
    assert nodata
    r, g, b = _rgb(nodata.group(1))
    assert max(r, g, b) - min(r, g, b) <= 16, "the not-enough-data fill must be a neutral grey"
    dark = re.search(r"@media \(prefers-color-scheme: dark\)\s*\{(.*?)\n\}", css, re.S)
    assert dark and "--nodata:" in dark.group(1)


def _js_ramp_stops() -> list[str]:
    return _js_frozen_array("RAMP_STOPS")


def test_the_ramp_is_built_from_the_manuscript_blues_and_darkens_towards_one() -> None:
    stops = _js_ramp_stops()
    palette = {c.lower() for c in COLORBLIND_PALETTE}
    assert {"#56b4e9", "#0072b2"} <= palette, (
        "the manuscript palette no longer has these blues"
    )
    assert {"#56b4e9", "#0072b2"} <= {c.lower() for c in stops}
    luminances = [_luminance(c) for c in stops]
    assert luminances == sorted(luminances, reverse=True), (
        "the ramp must get darker towards 100%, the convention resistance_cmap() documents"
    )
    # Every stop is blue-dominant: one hue, so the only thing that varies is lightness, which
    # is what survives colour-vision deficiency.
    for colour in stops:
        r, _g, b = _rgb(colour)
        assert b > r + 10, colour


@needs_node
def test_ramp_colours_get_strictly_darker_as_the_share_rises() -> None:
    colours = _run_page_functions(
        "console.log(JSON.stringify(Array.from({ length: 21 }, (_, i) => page.rampColor(i / 20))));"
    )
    luminances = [_luminance(c) for c in colours]
    assert all(a > b for a, b in pairwise(luminances)), luminances
    assert (
        colours[0] == _js_ramp_stops()[0].lower()
        and colours[-1] == _js_ramp_stops()[-1].lower()
    )


@needs_node
def test_cell_text_flips_to_white_on_dark_fills_and_stays_dark_on_pale_ones() -> None:
    tones = _run_page_functions(
        """
        const out = Array.from({ length: 21 }, (_, i) => page.textTone(page.rampColor(i / 20)));
        console.log(JSON.stringify(out));
        """
    )
    assert tones[0] == "dark" and tones[-1] == "light"
    assert tones == sorted(tones, key=lambda t: t == "light"), (
        "tone must flip once, not flicker"
    )


_MODEL_ROWS = [
    # PI: KE has a coloured cell in 2001 and a sparse one in 2002; UG's only 2002 cell is
    # sparse, then a coloured one in 2005; TZ has nothing at or above the floor.
    ("PI", "KE", 2001, 10, 1),
    ("PI", "KE", 2002, 3, 1),
    ("PI", "UG", 2002, 2, 0),
    ("PI", "UG", 2005, 7, 2),
    ("PI", "TZ", 2003, 1, 0),
    ("NNRTI", "TZ", 2009, 50, 5),
]


def _model_rows_json() -> str:
    rows = [
        {
            "drug_class": cls,
            "country": country,
            "collection_year": year,
            "n_sequences": n,
            "n_resistant": k,
            "frac_resistant": None if n < _FLOOR else k / n,
        }
        for cls, country, year, n, k in _MODEL_ROWS
    ]
    return json.dumps(rows)


@needs_node
def test_the_year_axis_keeps_only_years_with_a_coloured_cell_and_the_omissions_are_counted() -> (
    None
):
    model = _run_page_functions(
        f"""
        const m = page.landscapeModel({_model_rows_json()}, "PI", ["KE","UG","TZ"], {_FLOOR});
        console.log(JSON.stringify({{ years: m.years, countries: m.countries,
          cells: [...m.cells.keys()], omittedCells: m.omittedCells, omittedSequences: m.omittedSequences }}));
        """
    )
    assert model["years"] == [2001, 2005], "a year with only sparse cells must not be a column"
    assert model["countries"] == ["KE", "UG"], (
        "a country with only sparse cells must not be a row"
    )
    assert sorted(model["cells"]) == ["KE|2001", "UG|2005"]
    # KE 2002 (3), UG 2002 (2) and TZ 2003 (1) are dropped from the figure and said to be.
    assert (model["omittedCells"], model["omittedSequences"]) == (3, 6)


@needs_node
def test_rows_follow_the_exported_country_order_not_the_order_cells_arrive_in() -> None:
    order = _run_page_functions(
        f"""
        const rows = {_model_rows_json()}.reverse();
        console.log(JSON.stringify(page.landscapeModel(rows, "PI", ["UG","KE","TZ"], {_FLOOR}).countries));
        """
    )
    assert order == ["UG", "KE"]


@needs_node
def test_a_cell_whose_share_disagrees_with_the_floor_is_refused() -> None:
    outcome = _run_page_functions(
        f"""
        const row = (n, frac) => [{{ drug_class: "PI", country: "KE", collection_year: 2001,
          n_sequences: n, n_resistant: 0, frac_resistant: frac }}];
        const out = {{}};
        for (const [name, rows] of Object.entries({{
          coloured_below_floor: row({_FLOOR - 1}, 0.1), null_above_floor: row({_FLOOR}, null),
          fine_above: row({_FLOOR}, 0.2), fine_below: row({_FLOOR - 1}, null) }})) {{
          try {{ page.landscapeModel(rows, "PI", ["KE"], {_FLOOR}); out[name] = "ok"; }}
          catch (error) {{ out[name] = "threw"; }}
        }}
        console.log(JSON.stringify(out));
        """
    )
    assert outcome == {
        "coloured_below_floor": "threw",
        "null_above_floor": "threw",
        "fine_above": "ok",
        "fine_below": "ok",
    }


@needs_node
def test_non_adjacent_years_are_separated_by_a_break_and_adjacent_ones_are_not() -> None:
    layout = _run_page_functions(
        """
        const { cellWidth, xs } = page.columnLayout([2001, 2002, 2010], 300, 10);
        console.log(JSON.stringify({ cellWidth, xs: [...xs.entries()] }));
        """
    )
    xs = {int(year): x for year, x in layout["xs"]}
    width = layout["cellWidth"]
    assert xs[2002] - xs[2001] == pytest.approx(width)
    assert xs[2010] - xs[2002] == pytest.approx(width + 10)
    assert xs[2010] + width == pytest.approx(300), "the columns and the break fill the plot"


@needs_node
def test_year_labels_never_crowd_each_other_and_always_name_the_first_and_last_year() -> None:
    years = [1986, *range(2000, 2026)]
    labelled = _run_page_functions(
        f"""
        const years = {json.dumps(years)};
        const {{ cellWidth, xs }} = page.columnLayout(years, 313, 8);
        const labels = page.yearLabels(years, xs, cellWidth);
        console.log(JSON.stringify({{ labels, cellWidth, xs: years.map((y) => xs.get(y)) }}));
        """
    )
    labels = labelled["labels"]
    assert labels[0] == 1986 and labels[-1] == 2025
    centres = [labelled["xs"][years.index(y)] + labelled["cellWidth"] / 2 for y in labels]
    gaps = [b - a for a, b in pairwise(centres)]
    assert min(gaps) >= 34


def _render(rows_json: str, width: int, drug_class: str = "PI") -> Any:
    return _run_page_functions(
        f"""
        const ctx = {{ floor: {_FLOOR}, countries: [{{code:"KE",name:"Kenya"}},{{code:"UG",name:"Uganda"}},{{code:"TZ",name:"Tanzania"}}] }};
        const model = page.landscapeModel({rows_json}, "{drug_class}", ["KE","UG","TZ"], {_FLOOR});
        const clicked = [];
        const svg = page.renderLandscape(model, "{drug_class}", ctx, {{ clientWidth: {width} }}, (t) => clicked.push(t));
        const cells = cellsOf(svg);
        const first = svg.children.find((el) => el.attrs.class === "hm-cell");
        first.handlers.click();
        console.log(JSON.stringify({{ cells, clicked, viewBox: svg.attrs.viewBox }}));
        """
    )


def _wide_rows() -> str:
    # KE has coloured cells in 2001 and 2005; UG in 2005 only, so UG 2001 is an empty cell,
    # and KE 2005's sibling UG 2002 is sparse but in a year with no coloured cell elsewhere.
    rows = [
        ("PI", "KE", 2001, 1200, 120),
        ("PI", "KE", 2005, 7, 7),
        ("PI", "UG", 2005, 40, 0),
        ("PI", "UG", 2001, 4, 2),
    ]
    return json.dumps(
        [
            {
                "drug_class": c,
                "country": co,
                "collection_year": y,
                "n_sequences": n,
                "n_resistant": k,
                "frac_resistant": None if n < _FLOOR else k / n,
            }
            for c, co, y, n, k in rows
        ]
    )


@needs_node
def test_a_cell_under_the_floor_is_drawn_grey_with_a_slash_and_never_in_the_ramp() -> None:
    out = _render(_wide_rows(), 700)
    by_title = {c["title"].split(":")[0]: c for c in out["cells"]}
    sparse = by_title["Uganda 2001, PI"]
    assert sparse["rect"]["class"] == "hm-nodata" and sparse["slash"]
    assert "fill" not in sparse["rect"], "a sparse cell must not carry a ramp colour"
    coloured = by_title["Kenya 2001, PI"]
    assert coloured["rect"]["class"] == "hm-fill" and not coloured["slash"]
    assert coloured["rect"]["fill"] == _run_page_functions(
        "console.log(JSON.stringify(page.rampColor(120 / 1200)));"
    )
    # Even a sparse cell that is 50% resistant stays grey: "not enough data" is not a colour.
    assert sparse["title"].endswith("so no percentage is given.")


@needs_node
def test_every_drawn_cell_carries_its_counts_in_its_title_and_the_empty_ones_say_so() -> None:
    out = _render(_wide_rows(), 700)
    titles = {c["title"] for c in out["cells"]}
    assert (
        "Kenya 2001, PI: 120 of 1,200 sequences scored for PR are resistant (10.0%)." in titles
    )
    assert "Kenya 2005, PI: 7 of 7 sequences scored for PR are resistant (100.0%)." in titles
    assert "Uganda 2005, PI: 0 of 40 sequences scored for PR are resistant (0.0%)." in titles
    assert len(out["cells"]) == 4  # 2 countries x 2 years, none missing
    assert all(c["title"] for c in out["cells"])
    only_empty = [c for c in out["cells"] if c["title"].endswith("no sequences")]
    assert only_empty == [], "all four cells have rows in this fixture"


@needs_node
def test_a_country_year_with_no_sequences_is_an_empty_cell_with_a_title() -> None:
    rows = json.dumps(
        [
            {"drug_class": "PI", "country": "KE", "collection_year": 2001, "n_sequences": 9,
             "n_resistant": 1, "frac_resistant": 1 / 9},
            {"drug_class": "PI", "country": "UG", "collection_year": 2005, "n_sequences": 9,
             "n_resistant": 1, "frac_resistant": 1 / 9},
        ]
    )  # fmt: skip
    out = _render(rows, 700)
    empties = [c for c in out["cells"] if c["rect"]["class"] == "hm-empty"]
    assert sorted(c["title"] for c in empties) == [
        "Kenya 2005, PI: no sequences",
        "Uganda 2001, PI: no sequences",
    ]
    assert not any(c["click"] for c in empties), "an empty cell has nothing to read out"


@needs_node
def test_n_is_printed_in_the_cell_only_when_it_fits_and_large_counts_are_shortened() -> None:
    assert {c["n"] for c in _render(_wide_rows(), 700)["cells"]} == {"1.2k", "7", "40", "4"}
    many_years = json.dumps(
        [
            {
                "drug_class": "PI",
                "country": "KE",
                "collection_year": year,
                "n_sequences": 50,
                "n_resistant": 5,
                "frac_resistant": 0.1,
            }
            for year in range(2001, 2026)
        ]
    )
    wide, narrow = _render(many_years, 700), _render(many_years, 280)
    assert {c["n"] for c in wide["cells"]} == {"50"}
    assert {c["n"] for c in narrow["cells"]} == {None}, (
        "no room: the title and tap readout carry n"
    )
    assert len(narrow["cells"]) == len(wide["cells"]) == 25, "narrowing must not drop any cell"
    assert all(c["title"] for c in narrow["cells"])


@needs_node
def test_tapping_a_cell_hands_its_description_to_the_readout() -> None:
    out = _render(_wide_rows(), 700)
    assert out["clicked"] == [out["cells"][0]["title"]]


@needs_node
def test_a_class_with_no_cell_at_the_floor_says_so_instead_of_drawing_an_empty_grid() -> None:
    out = _run_page_functions(
        f"""
        const ctx = {{ floor: {_FLOOR}, countries: [{{code:"KE",name:"Kenya"}}] }};
        const rows = [{{ drug_class: "INSTI", country: "KE", collection_year: 2001,
          n_sequences: 2, n_resistant: 0, frac_resistant: null }}];
        const model = page.landscapeModel(rows, "INSTI", ["KE"], {_FLOOR});
        const svg = page.renderLandscape(model, "INSTI", ctx, {{ clientWidth: 400 }}, () => {{}});
        console.log(JSON.stringify(svg.children.filter((el) => el.tag === "text").map((el) => el.text)));
        """
    )
    assert out == [f"No INSTI cell reaches {_FLOOR} sequences."]


@needs_node
def test_bindings_that_cite_cells_resolve_and_a_missing_or_sparse_cell_throws() -> None:
    outcome = _run_page_functions(
        f"""
        const summary = {json.dumps(SUMMARY)};
        const out = {{}};
        const attempt = (path, s) => {{
          try {{ return page.boundText(page.resolveBinding(path, s, page.deriveFigures(s)), undefined); }}
          catch (error) {{ return "threw"; }}
        }};
        out.pct = attempt("cell.NNRTI.KE.2010.pct", summary);
        out.n = attempt("cell.NNRTI.KE.2010.n", summary);
        out.fold = attempt("ratio.NNRTI.UG.2016.KE.2010", summary);
        out.missing = attempt("cell.NNRTI.KE.1900.pct", summary);
        const sparse = structuredClone(summary);
        const target = sparse.views.resistance_landscape.find(
          (r) => r.drug_class === "NNRTI" && r.country === "KE" && r.collection_year === 2010);
        target.n_sequences = 2; target.frac_resistant = null;
        out.sparse_pct = attempt("cell.NNRTI.KE.2010.pct", sparse);
        out.sparse_n = attempt("cell.NNRTI.KE.2010.n", sparse);
        out.sparse_fold = attempt("ratio.NNRTI.UG.2016.KE.2010", sparse);
        console.log(JSON.stringify(out));
        """
    )
    low = _landscape_cell("NNRTI", "KE", "2010")
    assert outcome["pct"] == f"{low['frac_resistant'] * 100:.1f}%"
    assert outcome["n"] == f"{low['n_sequences']:,}"
    assert outcome["missing"] == "threw"
    assert outcome["sparse_pct"] == "threw" and outcome["sparse_fold"] == "threw"
    assert outcome["sparse_n"] == "2"


@needs_node
def test_the_digest_notice_shows_when_either_inputs_file_disagrees_with_its_pointer() -> None:
    outcome = _run_page_functions(
        """
        const agree = { input_dvc_md5: "a", input_file_md5: "a", curated_dvc_md5: "c", curated_file_md5: "c" };
        const out = {
          agree: page.digestsDiffer(agree),
          sierra: page.digestsDiffer({ ...agree, input_file_md5: "x" }),
          curated: page.digestsDiffer({ ...agree, curated_file_md5: "x" }),
        };
        console.log(JSON.stringify(out));
        """
    )
    assert outcome == {"agree": False, "sierra": True, "curated": True}
