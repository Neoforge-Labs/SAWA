"""Unit tests for causal.data_prep — the Paper1→Paper2 adapter."""

from __future__ import annotations

import polars as pl
import pytest

from causal.data_prep import (
    RESISTANT_LEVELS,
    _parse_mutation,
    build_dr_causal_frames,
)

# --- Helpers ---


def _dr(rows: list[dict]) -> pl.DataFrame:
    """Build a Sierra-shaped dr_results frame from dict rows."""
    return pl.from_dicts(rows)


def _curated(rows: list[dict]) -> pl.DataFrame:
    return pl.from_dicts(rows)


DEFAULT_DR_ROW = {
    "accession": "A1",
    "gene": "RT",
    "mutations": ["M184V"],
    "drug_class": "NRTI",
    "drug": "3TC",
    "level": "high",
}


# --- Constants ---


class TestConstants:
    def test_resistant_levels_matches_paper1(self) -> None:
        assert frozenset({"high", "intermediate"}) == RESISTANT_LEVELS


# --- Mutation parser ---


class TestParseMutation:
    def test_basic(self) -> None:
        assert _parse_mutation("RT", "M184V") == "RT_184_V"

    def test_all_three_genes(self) -> None:
        assert _parse_mutation("PR", "L90M") == "PR_90_M"
        assert _parse_mutation("RT", "K103N") == "RT_103_N"
        assert _parse_mutation("IN", "Q148H") == "IN_148_H"

    def test_gene_case_insensitive(self) -> None:
        assert _parse_mutation("rt", "M184V") == "RT_184_V"

    def test_stop_codon_variant(self) -> None:
        # Sierra emits early stop codons as '*'
        assert _parse_mutation("RT", "K65*") == "RT_65_*"

    def test_none_on_missing_inputs(self) -> None:
        assert _parse_mutation(None, "M184V") is None
        assert _parse_mutation("RT", None) is None
        assert _parse_mutation("", "M184V") is None

    def test_none_on_unknown_gene(self) -> None:
        assert _parse_mutation("VIF", "K22R") is None

    def test_none_on_insertions_and_mixed_variants(self) -> None:
        # Sierra encodes insertions as e.g. 'T69S_SS' — not a single variant
        assert _parse_mutation("RT", "T69S_SS") is None
        # Mixed calls like 'M184MV' → skip
        assert _parse_mutation("RT", "M184MV") is None
        # Deletions with '-'
        assert _parse_mutation("RT", "M184-") is None


# --- Frame builder — shape ---


class TestBuildDrCausalFrames:
    def test_returns_one_frame_per_drug_class_by_default(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR_ROW, "accession": "A1", "drug_class": "NRTI", "drug": "3TC"},
                {
                    **DEFAULT_DR_ROW,
                    "accession": "A1",
                    "drug_class": "NNRTI",
                    "drug": "EFV",
                    "mutations": ["K103N"],
                },
                {
                    **DEFAULT_DR_ROW,
                    "accession": "A2",
                    "drug_class": "NRTI",
                    "drug": "3TC",
                    "level": "susceptible",
                    "mutations": [],
                },
            ]
        )
        curated = _curated(
            [
                {"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "A2", "country": "KE", "collection_year": 2020, "subtype": "A1"},
            ]
        )
        frames, counts = build_dr_causal_frames(dr, curated)
        assert set(frames) == {"NRTI", "NNRTI"}
        assert counts.scope_kind == "drug_class"
        assert counts.n_scopes == 2

    def test_per_drug_mode(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR_ROW, "drug_class": "NRTI", "drug": "3TC"},
                {**DEFAULT_DR_ROW, "drug_class": "NRTI", "drug": "FTC", "mutations": []},
            ]
        )
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        frames, counts = build_dr_causal_frames(dr, curated, per_drug=True)
        assert set(frames) == {"3TC", "FTC"}
        assert counts.scope_kind == "drug"

    def test_frame_has_all_dag_mutation_columns(self) -> None:
        from causal.dag import DRDag

        dr = _dr([DEFAULT_DR_ROW])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        for mut in DRDag().mutation_features:
            assert mut in frames["NRTI"].columns

    def test_frame_has_covariates_and_label(self) -> None:
        dr = _dr([DEFAULT_DR_ROW])
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        cols = set(frames["NRTI"].columns)
        assert {"accession", "label_binary", "subtype", "country", "collection_year"} <= cols


# --- Frame builder — semantics ---


class TestSemantics:
    def test_high_and_intermediate_count_as_resistant(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR_ROW, "accession": "A1", "level": "high"},
                {**DEFAULT_DR_ROW, "accession": "A2", "level": "intermediate"},
                {**DEFAULT_DR_ROW, "accession": "A3", "level": "susceptible", "mutations": []},
                {**DEFAULT_DR_ROW, "accession": "A4", "level": "low", "mutations": []},
            ]
        )
        curated = _curated(
            [
                {
                    "accession": f"A{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "A1",
                }
                for i in range(1, 5)
            ]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        labels = dict(
            zip(
                frames["NRTI"]["accession"].to_list(),
                frames["NRTI"]["label_binary"].to_list(),
                strict=True,
            )
        )
        assert labels == {"A1": 1, "A2": 1, "A3": 0, "A4": 0}

    def test_any_drug_in_class_resistant_flags_accession(self) -> None:
        # A1 is high to 3TC but susceptible to AZT → still counts as NRTI-resistant.
        dr = _dr(
            [
                {**DEFAULT_DR_ROW, "accession": "A1", "drug": "3TC", "level": "high"},
                {
                    **DEFAULT_DR_ROW,
                    "accession": "A1",
                    "drug": "AZT",
                    "level": "susceptible",
                    "mutations": [],
                },
            ]
        )
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        assert frames["NRTI"]["label_binary"].to_list() == [1]

    def test_mutation_present_when_sierra_lists_it(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR_ROW, "accession": "A1", "mutations": ["M184V", "K103N"]},
            ]
        )
        # K103N is on RT, but this row's gene is RT so both parse to RT_*.
        curated = _curated(
            [{"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"}]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        row = frames["NRTI"].to_dicts()[0]
        assert row["RT_184_V"] == 1
        assert row["RT_103_N"] == 1
        # Unlisted mutations stay zero.
        assert row["PR_46_L"] == 0
        assert row["IN_148_H"] == 0

    def test_missing_covariates_dropped(self) -> None:
        dr = _dr(
            [
                {**DEFAULT_DR_ROW, "accession": "A1"},
                {**DEFAULT_DR_ROW, "accession": "A2"},
            ]
        )
        curated = _curated(
            [
                {"accession": "A1", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "A2", "country": None, "collection_year": 2020, "subtype": "A1"},
            ]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        assert frames["NRTI"]["accession"].to_list() == ["A1"]


# --- Validation ---


class TestValidation:
    def test_missing_dr_columns_raises(self) -> None:
        with pytest.raises(ValueError, match="dr_results missing"):
            build_dr_causal_frames(
                pl.DataFrame({"accession": ["A"]}),
                pl.DataFrame(
                    {
                        "accession": ["A"],
                        "country": ["KE"],
                        "collection_year": [2020],
                        "subtype": ["A1"],
                    }
                ),
            )

    def test_missing_curated_columns_raises(self) -> None:
        with pytest.raises(ValueError, match="curated missing"):
            build_dr_causal_frames(
                _dr([DEFAULT_DR_ROW]),
                pl.DataFrame({"accession": ["A1"]}),
            )


# --- End-to-end contract with the estimator ---


class TestEstimatorContract:
    def test_produced_frame_slots_into_estimate_dr_effect(self) -> None:
        """The single most important test — the frame's shape must be
        exactly what estimate_dr_effect() expects, so a downstream user
        can pipe adapter → estimator without gymnastics."""
        dr = _dr(
            [
                {
                    **DEFAULT_DR_ROW,
                    "accession": f"A{i}",
                    "mutations": ["M184V"] if i % 2 == 0 else [],
                }
                for i in range(1, 41)
            ]
        )
        curated = _curated(
            [
                {
                    "accession": f"A{i}",
                    "country": "KE",
                    "collection_year": 2015 + (i % 5),
                    "subtype": ["A1", "C"][i % 2],
                }
                for i in range(1, 41)
            ]
        )
        frames, _ = build_dr_causal_frames(dr, curated)
        nrti = frames["NRTI"]

        # The estimator requires all of these columns; missing any is a
        # broken contract.
        for col in ("RT_184_V", "label_binary", "subtype", "country", "collection_year"):
            assert col in nrti.columns, f"missing {col}"

        # Types must be usable by dowhy — mutation + label as ints,
        # subtype + country as categorical/string.
        assert nrti["label_binary"].dtype == pl.Int8
        assert nrti["RT_184_V"].dtype == pl.Int8
