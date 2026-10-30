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
import logging
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
    IngestedCountBelowScoredError,
    MissingDvcHashError,
    NonSingleValuedProvenanceError,
    NullProvenanceError,
    build_summary,
    check_input_not_stale,
    file_md5,
    input_dvc_md5,
    render_summary,
    sdrm_observation_counts,
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
        "input_file_md5": "def",
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
    # Format only: the digest changes whenever the parquet is regenerated.
    assert re.fullmatch(r"[0-9a-f]{32}", input_dvc_md5(POINTER))


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


# --- check_input_not_stale ---------------------------------------------------------

_MD5_A = "0123456789abcdef0123456789abcdef"
_MD5_B = "fedcba9876543210fedcba9876543210"
_MD5_C = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"


def _drift_inputs(
    tmp_path: Path, *, summary_meta: dict[str, str], pointer_md5: str
) -> tuple[Path, Path]:
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({"meta": summary_meta, "views": {}}))
    pointer = tmp_path / "x.parquet.dvc"
    pointer.write_text(
        f"outs:\n- md5: {pointer_md5}\n  size: 1\n  hash: md5\n  path: x.parquet\n"
    )
    return summary, pointer


def test_drift_check_passes_when_hashes_agree(tmp_path: Path) -> None:
    summary, pointer = _drift_inputs(
        tmp_path, summary_meta={"input_dvc_md5": _MD5_A}, pointer_md5=_MD5_A
    )
    assert check_input_not_stale(summary, pointer) == 0


def test_drift_check_fails_when_the_pointer_moved(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    summary, pointer = _drift_inputs(
        tmp_path, summary_meta={"input_dvc_md5": _MD5_A}, pointer_md5=_MD5_B
    )
    with caplog.at_level(logging.ERROR):
        assert check_input_not_stale(summary, pointer) == 1
    errors = [r.getMessage() for r in caplog.records if r.levelno == logging.ERROR]
    assert len(errors) == 1
    assert _MD5_A in errors[0]
    assert _MD5_B in errors[0]


def test_drift_check_compares_the_pointer_digest_not_the_file_digest(tmp_path: Path) -> None:
    # The committed vintage's two digests differ (pointer vs hashed file), so passing here
    # is what lets the guard be green on it. Only the pointer-to-pointer claim is checked.
    summary, pointer = _drift_inputs(
        tmp_path,
        summary_meta={"input_dvc_md5": _MD5_A, "input_file_md5": _MD5_C},
        pointer_md5=_MD5_A,
    )
    assert check_input_not_stale(summary, pointer) == 0


def test_drift_check_is_not_satisfied_by_a_matching_file_digest(tmp_path: Path) -> None:
    summary, pointer = _drift_inputs(
        tmp_path,
        summary_meta={"input_dvc_md5": _MD5_A, "input_file_md5": _MD5_B},
        pointer_md5=_MD5_B,
    )
    assert check_input_not_stale(summary, pointer) == 1


def test_drift_check_fails_rather_than_passes_when_the_summary_records_no_digest(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    summary, pointer = _drift_inputs(tmp_path, summary_meta={}, pointer_md5=_MD5_A)
    with caplog.at_level(logging.ERROR):
        assert check_input_not_stale(summary, pointer) == 1
    assert any(r.levelno == logging.ERROR for r in caplog.records)


def test_drift_check_propagates_an_unreadable_pointer(tmp_path: Path) -> None:
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({"meta": {"input_dvc_md5": _MD5_A}}))
    pointer = tmp_path / "x.dvc"
    pointer.write_text("outs: []\n")
    with pytest.raises(MissingDvcHashError):
        check_input_not_stale(summary, pointer)


# --- file_md5 ----------------------------------------------------------------------


@pytest.mark.parametrize("chunk_bytes", [1, 7, 4096, 1 << 20])
def test_file_md5_matches_hashlib_whatever_the_chunk_size(
    tmp_path: Path, chunk_bytes: int
) -> None:
    path = tmp_path / "blob.bin"
    path.write_bytes(bytes(range(256)) * 41)
    expected = hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()
    assert file_md5(path, chunk_bytes=chunk_bytes) == expected


def test_file_md5_of_an_empty_file_is_the_empty_digest(tmp_path: Path) -> None:
    path = tmp_path / "empty"
    path.write_bytes(b"")
    assert file_md5(path) == "d41d8cd98f00b204e9800998ecf8427e"


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
    assert meta["input_file_md5"] == "def"
    assert meta["git_sha"] == "deadbee"
    assert meta["n_accessions_scored"] == 3
    assert meta["n_accessions_ingested"] == 4
    assert meta["n_calls"] == 7
    assert meta["n_distinct_sequences"] == 3


def test_distinct_sequences_are_counted_apart_from_accessions() -> None:
    # SYN1 and SYN2 are given one sequence between them: three accessions, two sequences.
    df = _synthetic().with_columns(
        pl.when(pl.col("accession") == "SYN000002.1")
        .then(pl.lit(hashlib.sha256(b"SYN000001.1").hexdigest()))
        .otherwise(pl.col("sequence_sha256"))
        .alias("sequence_sha256")
    )
    meta = _build(df)["meta"]
    assert meta["n_accessions_scored"] == 3
    assert meta["n_distinct_sequences"] == 2


def test_pointer_digest_and_file_digest_are_recorded_under_distinct_keys() -> None:
    meta = _build(input_md5="from-pointer", input_file_md5="from-file")["meta"]
    assert meta["input_dvc_md5"] == "from-pointer"
    assert meta["input_file_md5"] == "from-file"


def test_meta_records_sdrm_observation_counts() -> None:
    # Fixture observations: PR L90M, M46I; RT NRTI M184V (SYN1); RT NNRTI K103N (SYN1),
    # K103S (SYN2). No mixtures, so nothing is multi-state.
    meta = _build()["meta"]
    assert meta["n_sdrm_observations"] == 5
    assert meta["n_sdrm_observations_multistate"] == 0


def test_sdrm_observation_counts_separate_single_and_multistate_matches() -> None:
    both_rt = ["M184VI", "K103N", "V35T", "I220*"]
    other_rt = ["K103NS", "V35T"]
    pr = ["M46IL", "L90M", "R41K"]
    rows = [
        _row("A", "PI", "ATV/r", "high", pr),
        _row("A", "NRTI", "3TC", "high", both_rt),
        _row("A", "NNRTI", "EFV", "high", both_rt),
        _row("B", "NRTI", "3TC", "low", other_rt),
        _row("B", "NNRTI", "EFV", "low", other_rt),
        _row("B", "INSTI", "DTG", "low", ["E138K", "G148H"]),
        _row("C", "PI", "ATV/r", "low", ["L90M"]),
    ]
    frame = pl.DataFrame(rows, schema=_SCHEMA)
    # PR: M46IL (I and L both resistant) is multi-state; L90M is single and on A and C.
    # RT NRTI: A's M184VI (V and I) is multi-state. RT NNRTI: A's K103N single, B's K103NS multi.
    # A's two RT drug rows repeat the same strings and must not double-count; R41K, V35T,
    # the stop codon and the integrase strings match no panel.
    assert sdrm_observation_counts(frame, _panels()) == (6, 3)


def test_sdrm_observation_counts_are_zero_for_a_frame_without_surveillance_mutations() -> None:
    frame = pl.DataFrame([_row("A", "PI", "ATV/r", "susceptible", ["R41K"])], schema=_SCHEMA)
    assert sdrm_observation_counts(frame, _panels()) == (0, 0)


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
    "column", ["algorithm_version", "called_at", "algorithm_publish_date", "schema_version"]
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


def _all_null(column: str) -> pl.DataFrame:
    return _synthetic().with_columns(pl.lit(None, dtype=pl.String).alias(column))


@pytest.mark.parametrize("column", ["algorithm_version", "called_at", "schema_version"])
def test_provenance_the_page_must_name_is_rejected_when_null(column: str) -> None:
    with pytest.raises(NullProvenanceError, match=column):
        _build(_all_null(column))


def test_publish_date_is_nullable_by_design_and_recorded_as_null() -> None:
    meta = _build(_all_null("algorithm_publish_date"))["meta"]
    assert "algorithm_publish_date" in meta
    assert meta["algorithm_publish_date"] is None


def test_ingested_count_below_the_scored_count_is_rejected() -> None:
    with pytest.raises(IngestedCountBelowScoredError, match=r"2 .*3"):
        _build(n_accessions_ingested=2)


def test_ingested_count_equal_to_the_scored_count_is_accepted() -> None:
    meta = _build(n_accessions_ingested=3)["meta"]
    assert meta["n_accessions_ingested"] == meta["n_accessions_scored"] == 3


def test_dr_schema_version_is_the_one_the_data_carries_not_the_code_constant(
    caplog: pytest.LogCaptureFixture,
) -> None:
    older = _synthetic().with_columns(pl.lit("0.0.1").alias("schema_version"))
    assert older["schema_version"][0] != DR_SCHEMA_VERSION
    with caplog.at_level(logging.WARNING):
        meta = _build(older)["meta"]
    assert meta["dr_schema_version"] == "0.0.1"
    assert any(
        "0.0.1" in r.getMessage() for r in caplog.records if r.levelno == logging.WARNING
    )


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
    # The file's own digest, not the pointer's: the two are different claims.
    assert meta["input_file_md5"] == hashlib.md5(FIXTURE_PATH.read_bytes()).hexdigest()
    assert meta["input_file_md5"] != meta["input_dvc_md5"]
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
    args = cli.parse_args(["--ingested-count", "4"])
    assert args.input == Path("data/interim/dr_results.parquet")
    assert args.pointer is None
    assert args.reference == Path("data/reference/drm_positions.json")
    assert args.out == cli.DEFAULT_OUT


def test_cli_requires_the_ingested_count(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        cli.parse_args([])
    assert excinfo.value.code == 2
    assert "--ingested-count" in capsys.readouterr().err


def test_cli_pointer_defaults_to_input_plus_dvc(tmp_path: Path) -> None:
    args = cli.parse_args(["--input", str(tmp_path / "dr.parquet"), "--ingested-count", "4"])
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


def _fake_git(
    monkeypatch: pytest.MonkeyPatch, *, status_out: str, status_code: int = 0
) -> None:
    def run(args: list[str], *_a: Any, **_k: Any) -> subprocess.CompletedProcess[str]:
        if args[1] == "rev-parse":
            return subprocess.CompletedProcess(args, 0, stdout="deadbeef\n", stderr="")
        assert args[1:3] == ["status", "--porcelain"], args
        return subprocess.CompletedProcess(args, status_code, stdout=status_out, stderr="")

    monkeypatch.setattr(cli.subprocess, "run", run)


def test_git_sha_reads_head_when_the_tree_is_clean(monkeypatch: pytest.MonkeyPatch) -> None:
    _fake_git(monkeypatch, status_out="")
    assert cli.resolve_git_sha() == "deadbeef"


def test_git_sha_is_marked_dirty_when_status_lists_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _fake_git(monkeypatch, status_out=" M src/hivdrwatch/dashboard/export/aggregate.py\n")
    assert cli.resolve_git_sha() == "deadbeef-dirty"


def test_git_sha_is_marked_dirty_when_status_cannot_be_read(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Clean cannot be shown, so it is not claimed.
    _fake_git(monkeypatch, status_out="", status_code=128)
    assert cli.resolve_git_sha() == "deadbeef-dirty"


def _git_repo(path: Path) -> str:
    def git(*args: str) -> str:
        done = subprocess.run(
            ["git", "-c", "user.name=t", "-c", "user.email=t@example.org", *args],
            cwd=path,
            capture_output=True,
            text=True,
            check=True,
        )
        return done.stdout.strip()

    git("init", "-q")
    (path / "tracked.txt").write_text("one\n")
    git("add", "tracked.txt")
    git("commit", "-q", "-m", "init")
    return git("rev-parse", "HEAD")


def test_real_repository_is_reported_by_its_bare_sha_only_while_clean(tmp_path: Path) -> None:
    head = _git_repo(tmp_path)
    assert cli.resolve_git_sha(tmp_path) == head

    (tmp_path / "tracked.txt").write_text("two\n")
    assert cli.resolve_git_sha(tmp_path) == f"{head}-dirty"


def test_real_repository_with_only_an_untracked_file_is_dirty(tmp_path: Path) -> None:
    head = _git_repo(tmp_path)
    (tmp_path / "new_module.py").write_text("x = 1\n")
    assert cli.resolve_git_sha(tmp_path) == f"{head}-dirty"


def test_committed_summary_names_a_commit_from_a_clean_tree() -> None:
    sha = json.loads((REPO / "src/hivdrwatch/dashboard/static/summary.json").read_text())[
        "meta"
    ]["git_sha"]
    assert re.fullmatch(r"[0-9a-f]{40}", sha), (
        f"summary.json records git_sha={sha!r}: regenerate it from a clean checkout"
    )
