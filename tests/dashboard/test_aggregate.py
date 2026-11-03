"""Tests for the dashboard's four aggregated views.

Each function receives an already-deduplicated frame with one row per
(accession, drug), while `mutations` is per (accession, gene) and so repeats
on every drug row of a gene. Most of the traps these tests pin come from that
shape: counting a mutation once per drug, or losing mutation-free sequences
from a denominator.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl
import pytest

from hivdrwatch.dashboard.export.aggregate import (
    LEVEL_ORDER,
    UnknownDrugClassError,
    UnknownLevelError,
    gene_coverage,
    levels_by_drug,
    prevalence_by_drug_class,
    sdrm_frequency,
)
from hivdrwatch.dashboard.export.drm import SdrmPanel, load_sdrm_panels
from hivdrwatch.paper1.style import DRUG_CLASS_ORDER

REFERENCE_PATH = (
    Path(__file__).resolve().parents[2] / "data" / "reference" / "drm_positions.json"
)

_GENE_OF_CLASS = {"PI": "PR", "NRTI": "RT", "NNRTI": "RT", "INSTI": "IN"}


def _row(
    accession: str,
    drug_class: str,
    drug: str,
    level: str = "susceptible",
    *,
    mutations: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "accession": accession,
        "drug_class": drug_class,
        "drug": drug,
        "level": level,
        "gene": _GENE_OF_CLASS[drug_class],
        "mutations": mutations if mutations is not None else [],
    }


def _frame(rows: list[dict[str, Any]]) -> pl.DataFrame:
    return pl.DataFrame(
        rows,
        schema={
            "accession": pl.String,
            "drug_class": pl.String,
            "drug": pl.String,
            "level": pl.String,
            "gene": pl.String,
            "mutations": pl.List(pl.String),
        },
    )


def _panels() -> tuple[SdrmPanel, ...]:
    return load_sdrm_panels(REFERENCE_PATH)[0]


def _frame_two_nnrti_calls() -> pl.DataFrame:
    return _frame(
        [_row("A1", "NNRTI", "NVP", "high"), _row("A1", "NNRTI", "EFV", "susceptible")]
    )


def _frame_with_level(level: str) -> pl.DataFrame:
    return _frame([_row("A1", "NNRTI", "NVP", level), _row("A1", "NNRTI", "EFV", level)])


def _frame_without_insti() -> pl.DataFrame:
    return _frame([_row("A1", "PI", "LPV", "high"), _row("A1", "NRTI", "AZT", "low")])


def _frame_uneven_gene_coverage() -> pl.DataFrame:
    # Two accessions scored for PR, only one of them also for IN.
    return _frame(
        [
            _row("A1", "PI", "LPV", "high"),
            _row("A1", "PI", "DRV", "susceptible"),
            _row("A2", "PI", "LPV", "susceptible"),
            _row("A2", "PI", "DRV", "susceptible"),
            _row("A1", "INSTI", "DTG", "low"),
        ]
    )


def _frame_mixed_levels() -> pl.DataFrame:
    levels = ["susceptible", "potential_low", "low", "intermediate", "high"]
    return _frame([_row(f"A{i}", "NNRTI", "NVP", level) for i, level in enumerate(levels)])


def _frame_repeated_mutation() -> pl.DataFrame:
    return _frame(
        [
            _row("A1", "NRTI", "3TC", "high", mutations=["M184V"]),
            _row("A1", "NRTI", "FTC", "high", mutations=["M184V"]),
            _row("A1", "NRTI", "AZT", "low", mutations=["M184V"]),
        ]
    )


def _frame_one_drm_one_clean() -> pl.DataFrame:
    return _frame(
        [
            _row("A1", "NRTI", "3TC", "high", mutations=["M184V"]),
            _row("A1", "NRTI", "AZT", "low", mutations=["M184V"]),
            _row("A2", "NRTI", "3TC", "susceptible", mutations=[]),
            _row("A2", "NRTI", "AZT", "susceptible", mutations=[]),
        ]
    )


def test_prevalence_counts_sequence_resistant_if_any_drug_in_class_is():
    # One accession: high to NVP, susceptible to EFV -> NNRTI-resistant once.
    out = prevalence_by_drug_class(_frame_two_nnrti_calls())
    row = out.filter(pl.col("drug_class") == "NNRTI").row(0, named=True)
    assert row["n_sequences"] == 1 and row["n_resistant"] == 1
    assert row["frac_resistant"] == pytest.approx(1.0)


def test_prevalence_uses_resistant_levels_from_fig02():
    from hivdrwatch.paper1.fig02_aggregate import RESISTANT_LEVELS

    assert frozenset({"high", "intermediate"}) == RESISTANT_LEVELS
    out = prevalence_by_drug_class(_frame_with_level("low"))
    assert out.filter(pl.col("drug_class") == "NNRTI").row(0, named=True)["n_resistant"] == 0
    out = prevalence_by_drug_class(_frame_with_level("intermediate"))
    assert out.filter(pl.col("drug_class") == "NNRTI").row(0, named=True)["n_resistant"] == 1


def test_prevalence_omits_a_drug_class_absent_from_the_vintage():
    out = prevalence_by_drug_class(_frame_without_insti())
    assert "INSTI" not in out["drug_class"].to_list()
    assert out["frac_resistant"].is_not_nan().all()


def test_denominators_are_per_class_not_global():
    # PR scored for 2 accessions, INSTI for 1.
    out = prevalence_by_drug_class(_frame_uneven_gene_coverage())
    by = {r["drug_class"]: r["n_sequences"] for r in out.iter_rows(named=True)}
    assert by["PI"] == 2 and by["INSTI"] == 1


def test_level_order_knows_every_resistant_level():
    from hivdrwatch.paper1.fig02_aggregate import RESISTANT_LEVELS

    assert set(LEVEL_ORDER) >= RESISTANT_LEVELS


def test_levels_by_drug_covers_all_five_levels_per_drug():
    out = levels_by_drug(_frame_mixed_levels())
    assert set(out.filter(pl.col("drug") == "NVP")["level"]) == {
        "susceptible",
        "potential_low",
        "low",
        "intermediate",
        "high",
    }


def test_levels_by_drug_zero_fills_levels_a_drug_never_received():
    out = levels_by_drug(_frame_with_level("high")).filter(pl.col("drug") == "NVP")
    by = {r["level"]: r["n_sequences"] for r in out.iter_rows(named=True)}
    assert by == {
        "susceptible": 0,
        "potential_low": 0,
        "low": 0,
        "intermediate": 0,
        "high": 1,
    }


def test_levels_by_drug_fractions_sum_to_one_per_drug():
    out = levels_by_drug(_frame_mixed_levels())
    assert out["frac"].sum() == pytest.approx(1.0)
    assert out["frac"].to_list() == pytest.approx([0.2] * 5)


def test_levels_by_drug_orders_levels_by_severity_not_alphabetically():
    out = levels_by_drug(_frame_mixed_levels())
    assert out["level"].to_list() == [
        "susceptible",
        "potential_low",
        "low",
        "intermediate",
        "high",
    ]


def test_levels_by_drug_rejects_an_unrecognised_level():
    frame = _frame([_row("A1", "NNRTI", "NVP", "very_high")])
    with pytest.raises(UnknownLevelError, match="very_high"):
        levels_by_drug(frame)


def test_gene_coverage_counts_distinct_accessions_per_gene():
    assert dict(gene_coverage(_frame_uneven_gene_coverage()).iter_rows()) == {"PR": 2, "IN": 1}


def test_sdrm_frequency_counts_each_accession_once_per_mutation():
    # Same accession, M184V listed on both its RT rows.
    out = sdrm_frequency(_frame_repeated_mutation(), _panels())
    row = out.filter(pl.col("mutation") == "M184V").row(0, named=True)
    assert row["n_sequences"] == 1
    # One accession over three drug rows: the denominator is where a missing reduction shows.
    assert row["n_panel_denominator"] == 1


def test_sdrm_denominator_includes_sequences_with_no_mutations():
    # Two RT accessions: one carries M184V, one has an empty mutation list.
    out = sdrm_frequency(_frame_one_drm_one_clean(), _panels())
    row = out.filter(pl.col("mutation") == "M184V").row(0, named=True)
    assert row["n_sequences"] == 1
    assert row["n_panel_denominator"] == 2
    assert row["frac"] == pytest.approx(0.5)


def test_sdrm_frequency_drops_polymorphisms_but_keeps_them_in_the_denominator():
    frame = _frame(
        [
            _row("A1", "NRTI", "3TC", mutations=["V35T", "M184V"]),
            _row("A2", "NRTI", "3TC", mutations=["V35T"]),
        ]
    )
    out = sdrm_frequency(frame, _panels())
    assert out["mutation"].to_list() == ["M184V"]
    assert out["n_panel_denominator"].to_list() == [2]


def test_sdrm_panel_only_sees_rows_of_its_own_gene():
    # K103N is an RT mutation; a PR row listing it must not feed the NNRTI panel,
    # and the PR sequence must not enter the RT denominator.
    frame = _frame(
        [
            _row("A1", "NNRTI", "NVP", mutations=["K103N"]),
            _row("A2", "PI", "LPV", mutations=["K103N"]),
        ]
    )
    out = sdrm_frequency(frame, _panels())
    row = out.filter((pl.col("panel") == "RT_NNRTI") & (pl.col("mutation") == "K103N")).row(
        0, named=True
    )
    assert row["n_sequences"] == 1 and row["n_panel_denominator"] == 1
    assert out.filter(pl.col("panel") == "PR").height == 0


def test_both_rt_panels_share_the_rt_denominator():
    frame = _frame(
        [
            _row("A1", "NRTI", "3TC", mutations=["M184V"]),
            _row("A1", "NNRTI", "NVP", mutations=["M184V"]),
            _row("A2", "NRTI", "3TC", mutations=[]),
            _row("A3", "NNRTI", "NVP", mutations=["K103N"]),
        ]
    )
    out = sdrm_frequency(frame, _panels())
    assert set(out["n_panel_denominator"]) == {3}
    assert set(out["panel"]) == {"RT_NRTI", "RT_NNRTI"}


def test_sdrm_frequency_with_no_matches_keeps_the_output_schema():
    frame = _frame([_row("A1", "NRTI", "3TC", mutations=["V35T"])])
    out = sdrm_frequency(frame, _panels())
    assert out.height == 0
    assert out.columns == ["panel", "mutation", "n_sequences", "n_panel_denominator", "frac"]


def test_sdrm_frequency_with_a_null_mutation_list_counts_the_sequence_once():
    frame = _frame([_row("A1", "NRTI", "3TC", mutations=["M184V"])]).vstack(
        pl.DataFrame(
            {
                "accession": ["A2"],
                "drug_class": ["NRTI"],
                "drug": ["3TC"],
                "level": ["susceptible"],
                "gene": ["RT"],
                "mutations": [None],
            },
            schema=_frame([]).schema,
        )
    )
    row = sdrm_frequency(frame, _panels()).row(0, named=True)
    assert row["n_panel_denominator"] == 2


@pytest.mark.parametrize(
    "build",
    [
        lambda f, p: prevalence_by_drug_class(f),
        lambda f, p: levels_by_drug(f),
        lambda f, p: gene_coverage(f),
        lambda f, p: sdrm_frequency(f, p),
    ],
    ids=["prevalence", "levels", "coverage", "sdrm"],
)
def test_output_does_not_depend_on_input_row_order(build):
    rows = []
    for i in range(30):
        for cls, drug in (("PI", "LPV"), ("NRTI", "3TC"), ("NNRTI", "NVP"), ("INSTI", "DTG")):
            level = ("susceptible", "low", "intermediate", "high")[(i + len(drug)) % 4]
            muts = [("M184V", "K103N", "V35T", "L90M")[(i + j) % 4] for j in range(i % 3)]
            rows.append(_row(f"A{i:02d}", cls, drug, level, mutations=muts))
    frame = _frame(rows)
    reference = build(frame, _panels())
    for seed in range(5):
        shuffled = frame.sample(fraction=1.0, shuffle=True, seed=seed)
        assert build(shuffled, _panels()).equals(reference)


def test_drug_classes_follow_the_manuscript_order_not_alphabetical():
    rows = [
        _row("A1", "INSTI", "DTG", "low"),
        _row("A1", "NNRTI", "NVP", "low"),
        _row("A1", "NRTI", "3TC", "low"),
        _row("A1", "PI", "LPV", "low"),
    ]
    for build in (prevalence_by_drug_class, levels_by_drug):
        out = build(_frame(rows))
        seen = list(dict.fromkeys(out["drug_class"].to_list()))
        assert seen == list(DRUG_CLASS_ORDER)


def test_an_unknown_drug_class_is_rejected_rather_than_sorted_to_the_end():
    frame = _frame([_row("A1", "NRTI", "3TC")]).with_columns(drug_class=pl.lit("NCINI"))
    with pytest.raises(UnknownDrugClassError, match="NCINI"):
        prevalence_by_drug_class(frame)


def test_missing_columns_are_named():
    frame = _frame_two_nnrti_calls().drop("level")
    with pytest.raises(ValueError, match="level"):
        prevalence_by_drug_class(frame)
    with pytest.raises(ValueError, match="gene"):
        gene_coverage(frame.drop("gene"))
    with pytest.raises(ValueError, match="mutations"):
        sdrm_frequency(frame.drop("mutations"), _panels())
