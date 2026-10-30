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
from collections.abc import Iterator
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

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


def _bound_value(path: str, attrs: dict[str, str]) -> str:
    """Python twin of `resolveBinding` + `applyBindings` in app.js."""
    scope, key = path.split(".")
    views = SUMMARY["views"]
    if scope == "meta":
        value = SUMMARY["meta"][key]
    elif scope == "gene":
        value = next(r["n_sequences"] for r in views["gene_coverage"] if r["gene"] == key)
    elif scope == "derived":
        value = {
            "n_unscored": SUMMARY["meta"]["n_accessions_ingested"]
            - SUMMARY["meta"]["n_accessions_scored"],
            "n_drugs": len({r["drug"] for r in views["levels_by_drug"]}),
            "n_classes": len({r["drug_class"] for r in views["levels_by_drug"]}),
        }[key]
    else:
        raise AssertionError(f"unknown binding scope {scope!r}")
    if isinstance(value, int):
        return f"{value:,}"
    text = str(value)
    return text[: int(attrs["data-trunc"])] if "data-trunc" in attrs else text


def _rendered_text(node: _Node) -> str:
    """The element's text as a reader sees it once app.js has filled the bindings."""
    if node.tag in {"script", "style", "noscript"}:
        return ""
    if "data-bind" in node.attrs:
        return _bound_value(node.attrs["data-bind"], node.attrs)
    return "".join(_rendered_text(c) if isinstance(c, _Node) else c for c in node.children)


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


PAGE = _parse(INDEX.read_text())


# --- schema gate ----------------------------------------------------------------------


def test_body_declares_the_schema_version_summary_json_carries() -> None:
    body = next(n for n in _walk(PAGE) if n.tag == "body")
    assert body.attrs.get("data-summary-schema") == SUMMARY_SCHEMA_VERSION
    assert SUMMARY["meta"]["schema_version"] == SUMMARY_SCHEMA_VERSION


def test_schema_gate_runs_before_anything_is_drawn_and_shows_a_banner() -> None:
    main = _js_function("main")
    gate = main.index("schemaMismatch(")
    assert gate < main.index("applyBindings(") < main.index("drawAll(")
    mismatch_branch = main[gate : main.index("return;", gate)]
    assert "showBanner(" in mismatch_branch
    assert re.search(r'getAttribute\(\s*"data-summary-schema"\s*\)', main)
    assert re.search(r"summary\??\.meta\??\.schema_version", main)


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
        in {"provenance", "prevalence", "levels", "sdrm", "coverage", "footer"}
    ]
    assert ids == ["provenance", "prevalence", "levels", "sdrm", "coverage", "footer"]


def test_page_makes_no_region_claim_in_any_shipped_file() -> None:
    banned = re.compile(r"africa|kenya|uganda|tanzania|rwanda|burundi|ethiopia|saharan", re.I)
    for path in (INDEX, STYLES, APP_JS):
        text = path.read_text().replace(_RETRIEVAL_PHRASE, "")
        assert not banned.search(text), f"{path.name}: {banned.search(text)}"


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
    assert "attributed to a single one" in note
    assert "1.6% of observations (427 of 27,391)" in note


def test_limitation_note_links_the_provenance_evaluation_that_exists() -> None:
    limitation = _by_id(PAGE, "limitation")
    hrefs = [n.attrs["href"] for n in _walk(limitation) if n.tag == "a"]
    target = "docs/evaluations/ingest-provenance.md"
    assert any(h.endswith(f"/blob/main/{target}") for h in hrefs), hrefs
    assert (REPO / target).is_file()
    assert "deliberately absent" in _note("limitation")


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
