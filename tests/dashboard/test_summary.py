"""Tests for the deterministic `summary.json` export.

The summary is committed and reviewed in diffs, and a CI guard compares the
input hash recorded in it against the DVC pointer, so these tests pin byte
stability and provenance. Exact data values are pinned only for the synthetic
fixture, never for the real vintage, which is regenerated over time.

The synthetic fixture is shaped so that rounding is observable (RT has three
sequences, so fractions are thirds) and so that a duplicate-dropping or
integrase-leaking bug changes a number: PR and IN each have one sequence, one
raw call is stored three times, and one is stored twice.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

import polars as pl
import pytest

from hivdrwatch.dashboard.export import __main__ as cli
from hivdrwatch.dashboard.export.aggregate import prevalence_by_drug_class
from hivdrwatch.dashboard.export.dedup import NonRedundantDuplicateError, dedup_dr_results
from hivdrwatch.dashboard.export.drm import SdrmPanel, load_sdrm_panels
from hivdrwatch.dashboard.export.summary import (
    SUMMARY_SCHEMA_VERSION,
    MissingDvcHashError,
    NonSingleValuedProvenanceError,
    build_summary,
    input_dvc_md5,
    render_summary,
    write_summary,
)
from hivdrwatch.resistance.schema import DR_SCHEMA_VERSION

REPO = Path(__file__).resolve().parents[2]
POINTER = REPO / "data" / "interim" / "dr_results.parquet.dvc"
REFERENCE_PATH = REPO / "data" / "reference" / "drm_positions.json"
FIXTURE_DIR = REPO / "tests" / "fixtures" / "dashboard"
FIXTURE_PATH = FIXTURE_DIR / "_synthetic_dr_results.parquet"
GOLDEN_PATH = FIXTURE_DIR / "_synthetic_summary.json"

_GENE_OF_CLASS = {"PI": "PR", "NRTI": "RT", "NNRTI": "RT", "INSTI": "IN"}

_SCHEMA = {
    "accession": pl.String,
    "sequence_sha256": pl.String,
    "algorithm_name": pl.String,
    "algorithm_version": pl.String,
    "algorithm_publish_date": pl.String,
    "drug_class": pl.String,
    "drug": pl.String,
    "drug_display": pl.String,
    "score": pl.Float64,
    "level": pl.String,
    "gene": pl.String,
    "mutations": pl.List(pl.String),
    "called_at": pl.String,
    "schema_version": pl.String,
}


def _row(
    accession: str, drug_class: str, drug: str, level: str, mutations: list[str]
) -> dict[str, Any]:
    return {
        "accession": accession,
        "sequence_sha256": hashlib.sha256(accession.encode()).hexdigest(),
        "algorithm_name": "HIVDB",
        "algorithm_version": "HIVDB_10.2",
        "algorithm_publish_date": "2026-04-26",
        "drug_class": drug_class,
        "drug": drug,
        "drug_display": drug.lower(),
        "score": 0.0,
        "level": level,
        "gene": _GENE_OF_CLASS[drug_class],
        "mutations": mutations,
        "called_at": "2026-08-17",
        "schema_version": DR_SCHEMA_VERSION,
    }


def _synthetic() -> pl.DataFrame:
    """3 accessions, 7 distinct calls, 10 stored rows.

    SYN1 carries PR and RT; SYN2 carries RT; SYN3 carries RT and IN. Polymorphisms
    (R41K, V35T) and integrase mutations (E138K, G148H, R263K) are present so that a
    panel leak shows up as an extra mutation row.
    """
    pr = ["L90M", "M46I", "R41K"]
    rt1 = ["M184V", "K103N", "V35T"]
    rt2 = ["K103S", "V35T"]
    integrase = ["E138K", "G148H", "R263K"]
    rows = [
        _row("SYN000001.1", "PI", "ATV/r", "high", pr),
        _row("SYN000001.1", "NRTI", "3TC", "high", rt1),
        _row("SYN000001.1", "NNRTI", "EFV", "high", rt1),
        _row("SYN000002.1", "NRTI", "3TC", "susceptible", rt2),
        _row("SYN000002.1", "NNRTI", "EFV", "low", rt2),
        _row("SYN000003.1", "NRTI", "3TC", "susceptible", []),
        _row("SYN000003.1", "INSTI", "DTG", "intermediate", integrase),
    ]
    stored = [*rows, rows[0], rows[0], rows[6]]
    return pl.DataFrame(stored, schema=_SCHEMA)


def _panels() -> tuple[SdrmPanel, ...]:
    return load_sdrm_panels(REFERENCE_PATH)[0]


def _build(df: pl.DataFrame | None = None, **overrides: Any) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "panels": _panels(),
        "sdrm_version": "SDRM-2019",
        "input_md5": "abc",
        "git_sha": "deadbee",
        "n_accessions_ingested": 4,
    }
    kwargs.update(overrides)
    return build_summary(_synthetic() if df is None else df, **kwargs)


def _assert_golden(summary: dict[str, Any], path: Path, *, update: bool) -> None:
    rendered = render_summary(summary)
    if update:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rendered, encoding="utf-8")
    assert path.exists(), f"{path} is missing; regenerate it with --update-golden"
    assert rendered == path.read_text(encoding="utf-8")


# --- input_dvc_md5 -----------------------------------------------------------------


def test_input_dvc_md5_reads_the_pointer() -> None:
    assert input_dvc_md5(POINTER) == "e51cbb39d5ebf3bdcaff05d4df52b515"


def test_input_dvc_md5_ignores_the_hash_algorithm_line(tmp_path: Path) -> None:
    pointer = tmp_path / "x.dvc"
    pointer.write_text(
        "outs:\n- hash: md5\n  md5: 0123456789abcdef0123456789abcdef\n  size: 1\n"
    )
    assert input_dvc_md5(pointer) == "0123456789abcdef0123456789abcdef"


@pytest.mark.parametrize(
    "text",
    [
        "outs:\n- size: 1\n  hash: md5\n  path: x.parquet\n",
        "outs:\n- md5: not-a-digest\n  hash: md5\n",
        "outs:\n- md5: abc\n  hash: md5\n",
        "outs:\n- md5: 0123456789abcdef0123456789abcdef.dir\n  hash: md5\n",
        "",
    ],
)
def test_input_dvc_md5_raises_when_no_digest(tmp_path: Path, text: str) -> None:
    pointer = tmp_path / "x.dvc"
    pointer.write_text(text)
    with pytest.raises(MissingDvcHashError):
        input_dvc_md5(pointer)


# --- build_summary -----------------------------------------------------------------


def test_meta_records_provenance_without_asserting_data_values() -> None:
    meta = _build()["meta"]
    assert meta["schema_version"] == SUMMARY_SCHEMA_VERSION
    assert meta["dr_schema_version"] == DR_SCHEMA_VERSION
    assert meta["algorithm_version"] == "HIVDB_10.2"
    assert meta["algorithm_publish_date"] == "2026-04-26"
    assert meta["called_at"] == "2026-08-17"
    assert meta["sdrm_list_version"] == "SDRM-2019"
    assert meta["input_dvc_md5"] == "abc"
    assert meta["git_sha"] == "deadbee"
    assert meta["n_accessions_scored"] == 3
    assert meta["n_accessions_ingested"] == 4
    assert meta["n_calls"] == 7


def test_meta_counts_calls_after_deduplication() -> None:
    meta = _build()["meta"]
    assert meta["n_rows_raw"] == 10
    assert meta["n_calls"] == 7


def test_views_are_the_four_aggregations_in_row_form() -> None:
    views = _build()["views"]
    assert set(views) == {
        "prevalence_by_drug_class",
        "levels_by_drug",
        "gene_coverage",
        "sdrm_frequency",
    }
    assert [v["drug_class"] for v in views["prevalence_by_drug_class"]] == [
        "PI",
        "NRTI",
        "NNRTI",
        "INSTI",
    ]
    assert views["gene_coverage"] == [
        {"gene": "IN", "n_sequences": 1},
        {"gene": "PR", "n_sequences": 1},
        {"gene": "RT", "n_sequences": 3},
    ]


def test_non_redundant_duplicates_propagate_from_dedup() -> None:
    df = _synthetic()
    conflicting = pl.concat([df, df.head(1).with_columns(score=pl.lit(9.0))])
    with pytest.raises(NonRedundantDuplicateError):
        _build(conflicting)


@pytest.mark.parametrize(
    "column", ["algorithm_version", "called_at", "algorithm_publish_date"]
)
def test_provenance_that_is_not_single_valued_is_rejected(column: str) -> None:
    df = _synthetic()
    # Distinct values on a distinct accession keep dedup from absorbing the change.
    other = df.filter(pl.col("accession") == "SYN000002.1").with_columns(
        pl.lit("9999-other").alias(column)
    )
    mixed = pl.concat([df.filter(pl.col("accession") != "SYN000002.1"), other])
    with pytest.raises(NonSingleValuedProvenanceError) as excinfo:
        _build(mixed)
    message = str(excinfo.value)
    assert column in message
    assert "9999-other" in message


def test_summary_carries_no_region_claim() -> None:
    blob = json.dumps(_build()).lower()
    for banned in ("east africa", "kenya", "uganda", "tanzania", "country", "region"):
        assert banned not in blob


def test_summary_excludes_integrase_drms() -> None:
    sdrm = _build()["views"]["sdrm_frequency"]
    assert {v["panel"] for v in sdrm} == {"PR", "RT_NRTI", "RT_NNRTI"}
    assert {v["mutation"] for v in sdrm}.isdisjoint({"E138K", "G148H", "R263K"})


def test_summary_excludes_natural_polymorphisms() -> None:
    mutations = {v["mutation"] for v in _build()["views"]["sdrm_frequency"]}
    assert mutations == {"L90M", "M46I", "M184V", "K103N", "K103S"}


def test_fractions_are_rounded_to_four_places() -> None:
    unrounded = prevalence_by_drug_class(dedup_dr_results(_synthetic())[0])
    nrti = unrounded.filter(pl.col("drug_class") == "NRTI")["frac_resistant"][0]
    # Guards the fixture: a frac that is already 4-place exact would pass without rounding.
    assert nrti != round(nrti, 4)

    prevalence = _build()["views"]["prevalence_by_drug_class"]
    for v in prevalence:
        assert v["frac_resistant"] == round(v["frac_resistant"], 4)
    assert {v["drug_class"]: v["frac_resistant"] for v in prevalence}["NRTI"] == 0.3333


def test_counts_stay_integers() -> None:
    views = _build()["views"]
    for rows in views.values():
        for row in rows:
            for key, value in row.items():
                if key.startswith("n_"):
                    assert type(value) is int, (key, value)


# --- write_summary / render_summary -------------------------------------------------


def test_write_summary_is_byte_identical_across_runs(tmp_path: Path) -> None:
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    write_summary(_build(), a)
    write_summary(_build(), b)
    assert a.read_bytes() == b.read_bytes()


def test_output_does_not_depend_on_input_row_order() -> None:
    df = _synthetic()
    shuffled = df.sample(fraction=1.0, shuffle=True, seed=7)
    assert shuffled["accession"].to_list() != df["accession"].to_list()
    assert render_summary(_build(shuffled)) == render_summary(_build(df))
    assert render_summary(_build(df.reverse())) == render_summary(_build(df))


def test_write_summary_rounds_floats_it_is_given(tmp_path: Path) -> None:
    out = tmp_path / "s.json"
    write_summary({"meta": {}, "views": {"v": [{"frac": 1 / 3, "n": 3}]}}, out)
    assert json.loads(out.read_text())["views"]["v"][0] == {"frac": 0.3333, "n": 3}


def test_rendered_json_has_no_long_floats() -> None:
    assert not re.search(r"\d\.\d{5,}", render_summary(_build()))


def test_render_is_sorted_indented_unescaped_and_newline_terminated() -> None:
    text = render_summary({"b": 1, "a": "é"})
    assert text == '{\n  "a": "é",\n  "b": 1\n}\n'


def test_write_summary_creates_missing_parent_directories(tmp_path: Path) -> None:
    out = tmp_path / "a" / "b" / "summary.json"
    write_summary(_build(), out)
    assert json.loads(out.read_text())["meta"]["n_calls"] == 7


def test_write_summary_refuses_non_finite_floats(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Out of range float"):
        write_summary({"x": float("nan")}, tmp_path / "s.json")


# --- golden ------------------------------------------------------------------------


def test_committed_fixture_matches_its_builder() -> None:
    assert pl.read_parquet(FIXTURE_PATH).equals(_synthetic())


def test_golden_summary_matches(golden_update: bool) -> None:
    if golden_update:
        FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
        _synthetic().write_parquet(FIXTURE_PATH)
    summary = _build(pl.read_parquet(FIXTURE_PATH))
    _assert_golden(summary, GOLDEN_PATH, update=golden_update)


def test_golden_comparison_fails_on_a_changed_summary(tmp_path: Path) -> None:
    golden = tmp_path / "g.json"
    _assert_golden(_build(), golden, update=True)
    changed = _build(n_accessions_ingested=5)
    with pytest.raises(AssertionError):
        _assert_golden(changed, golden, update=False)


# --- CLI ---------------------------------------------------------------------------


def _cli_args(tmp_path: Path, out: Path) -> list[str]:
    pointer = tmp_path / "in.parquet.dvc"
    pointer.write_text("outs:\n- md5: 0123456789abcdef0123456789abcdef\n  hash: md5\n")
    return [
        "--input",
        str(FIXTURE_PATH),
        "--pointer",
        str(pointer),
        "--reference",
        str(REFERENCE_PATH),
        "--out",
        str(out),
        "--ingested-count",
        "4",
    ]


def test_cli_writes_only_summary_json(tmp_path: Path) -> None:
    out = tmp_path / "static" / "summary.json"
    assert cli.main(_cli_args(tmp_path, out)) == 0
    assert [p.name for p in out.parent.iterdir()] == ["summary.json"]
    meta = json.loads(out.read_text())["meta"]
    assert meta["input_dvc_md5"] == "0123456789abcdef0123456789abcdef"
    assert meta["n_accessions_ingested"] == 4
    assert meta["n_calls"] == 7
    assert meta["sdrm_list_version"] == "SDRM-2019"


def test_cli_output_is_stable_across_runs(tmp_path: Path) -> None:
    a, b = tmp_path / "a" / "summary.json", tmp_path / "b" / "summary.json"
    cli.main(_cli_args(tmp_path, a))
    cli.main(_cli_args(tmp_path, b))
    assert a.read_bytes() == b.read_bytes()


def test_cli_defaults_point_inside_the_package_and_repo() -> None:
    assert cli.DEFAULT_OUT.parts[-4:] == ("hivdrwatch", "dashboard", "static", "summary.json")
    assert cli.DEFAULT_INGESTED_COUNT == 36510
    args = cli.parse_args([])
    assert args.input == Path("data/interim/dr_results.parquet")
    assert args.pointer is None
    assert args.reference == Path("data/reference/drm_positions.json")
    assert args.out == cli.DEFAULT_OUT


def test_cli_pointer_defaults_to_input_plus_dvc(tmp_path: Path) -> None:
    args = cli.parse_args(["--input", str(tmp_path / "dr.parquet")])
    assert cli.pointer_path(args) == tmp_path / "dr.parquet.dvc"


def test_git_sha_falls_back_when_git_is_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(*_args: Any, **_kwargs: Any) -> None:
        raise FileNotFoundError("git")

    monkeypatch.setattr(cli.subprocess, "run", boom)
    assert cli.resolve_git_sha() == "unknown"


def test_git_sha_falls_back_on_nonzero_exit(monkeypatch: pytest.MonkeyPatch) -> None:
    def failed(*_args: Any, **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(["git"], 128, stdout="deadbeef\n", stderr="x")

    monkeypatch.setattr(cli.subprocess, "run", failed)
    assert cli.resolve_git_sha() == "unknown"


def test_git_sha_reads_head(monkeypatch: pytest.MonkeyPatch) -> None:
    def ok(*_args: Any, **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(["git"], 0, stdout="deadbeef\n", stderr="")

    monkeypatch.setattr(cli.subprocess, "run", ok)
    assert cli.resolve_git_sha() == "deadbeef"
