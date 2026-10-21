"""Tests for hivdrwatch.cluster.features.build_outbreak_features."""

from __future__ import annotations

import math

import polars as pl
import pytest

from hivdrwatch.cluster.features import (
    REPORTED_SUBTYPES,
    build_outbreak_features,
)


def _seq_meta(rows: list[dict]) -> pl.DataFrame:
    return (
        pl.from_dicts(rows)
        if rows
        else pl.DataFrame(
            schema={
                "accession": pl.Utf8,
                "country": pl.Utf8,
                "collection_year": pl.Int64,
                "subtype": pl.Utf8,
            }
        )
    )


def _cluster(rows: list[dict]) -> pl.DataFrame:
    return (
        pl.from_dicts(rows)
        if rows
        else pl.DataFrame(schema={"accession": pl.Utf8, "cluster_id": pl.Int64})
    )


def _growth(rows: list[dict]) -> pl.DataFrame:
    return (
        pl.from_dicts(rows)
        if rows
        else pl.DataFrame(schema={"cluster_id": pl.Int64, "growth_rate": pl.Float64})
    )


def _subtype_mix(rows: list[dict]) -> pl.DataFrame:
    return (
        pl.from_dicts(rows)
        if rows
        else pl.DataFrame(
            schema={"cluster_id": pl.Int64, "subtype": pl.Utf8, "fraction": pl.Float64}
        )
    )


class TestSchema:
    def test_output_has_all_outbreak_dag_columns(self) -> None:
        expected = {
            "country",
            "year",
            "n_sequences_deposited",
            "n_new_cluster_members",
            "n_active_clusters",
            "top_cluster_size",
            "mean_growth_rate",
            *(f"subtype_{s.lower()}_fraction" for s in REPORTED_SUBTYPES),
            "mean_subtype_entropy",
            "deposit_rate",
        }
        seq = _seq_meta(
            [
                {"accession": "A", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "B", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "C", "country": "UG", "collection_year": 2021, "subtype": "D"},
            ]
        )
        clusters = _cluster([{"accession": "A", "cluster_id": 1}])
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=clusters,
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        )
        assert set(out.columns) == expected


class TestCounts:
    def test_deposits_and_new_members(self) -> None:
        seq = _seq_meta(
            [
                {"accession": "A", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "B", "country": "KE", "collection_year": 2020, "subtype": "A1"},
                {"accession": "C", "country": "KE", "collection_year": 2020, "subtype": "C"},
            ]
        )
        clusters = _cluster(
            [
                {"accession": "A", "cluster_id": 1},
                {"accession": "B", "cluster_id": 1},
                # C stays a singleton (no cluster_id)
            ]
        )
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=clusters,
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        )
        row = out.filter((pl.col("country") == "KE") & (pl.col("year") == 2020)).row(
            0, named=True
        )
        assert row["n_sequences_deposited"] == 3
        assert row["n_new_cluster_members"] == 2
        assert row["n_active_clusters"] == 1
        assert row["top_cluster_size"] == 2


class TestSubtypeFractions:
    def test_named_fractions_sum_to_typed_ratio(self) -> None:
        # 4 A1, 2 D, 2 C in Uganda 2020
        seq = _seq_meta(
            [
                {
                    "accession": f"A{i}",
                    "country": "UG",
                    "collection_year": 2020,
                    "subtype": "A1",
                }
                for i in range(4)
            ]
            + [
                {
                    "accession": f"D{i}",
                    "country": "UG",
                    "collection_year": 2020,
                    "subtype": "D",
                }
                for i in range(2)
            ]
            + [
                {
                    "accession": f"C{i}",
                    "country": "UG",
                    "collection_year": 2020,
                    "subtype": "C",
                }
                for i in range(2)
            ]
        )
        clusters = _cluster([])
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=clusters,
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        )
        row = out.row(0, named=True)
        assert row["subtype_a1_fraction"] == pytest.approx(0.5)
        assert row["subtype_c_fraction"] == pytest.approx(0.25)
        assert row["subtype_d_fraction"] == pytest.approx(0.25)


class TestShannonEntropy:
    def test_zero_when_only_one_subtype(self) -> None:
        seq = _seq_meta(
            [
                {
                    "accession": f"A{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "A1",
                }
                for i in range(10)
            ]
        )
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=_cluster([]),
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        )
        assert out.row(0, named=True)["mean_subtype_entropy"] == pytest.approx(0.0)

    def test_ln2_when_two_equal_subtypes(self) -> None:
        seq = _seq_meta(
            [
                {
                    "accession": f"A{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "A1",
                }
                for i in range(5)
            ]
            + [
                {
                    "accession": f"D{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "D",
                }
                for i in range(5)
            ]
        )
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=_cluster([]),
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        )
        assert out.row(0, named=True)["mean_subtype_entropy"] == pytest.approx(math.log(2))


class TestNullHandling:
    def test_missing_country_or_year_rows_dropped(self) -> None:
        seq = _seq_meta(
            [
                {"accession": "A", "country": None, "collection_year": 2020, "subtype": "A1"},
                {"accession": "B", "country": "KE", "collection_year": None, "subtype": "A1"},
                {"accession": "C", "country": "KE", "collection_year": 2020, "subtype": "A1"},
            ]
        )
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=_cluster([]),
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        )
        assert out.height == 1
        assert out.row(0, named=True)["n_sequences_deposited"] == 1


class TestDepositRate:
    def test_normalised_to_country_max(self) -> None:
        seq = _seq_meta(
            [
                {
                    "accession": f"A{i}",
                    "country": "KE",
                    "collection_year": 2020,
                    "subtype": "A1",
                }
                for i in range(10)
            ]
            + [
                {
                    "accession": f"B{i}",
                    "country": "KE",
                    "collection_year": 2021,
                    "subtype": "A1",
                }
                for i in range(5)
            ]
        )
        out = build_outbreak_features(
            sequence_metadata=seq,
            cluster_assignments=_cluster([]),
            cluster_growth=_growth([]),
            cluster_subtypes=_subtype_mix([]),
        ).sort("year")
        assert out.row(0, named=True)["deposit_rate"] == pytest.approx(1.0)
        assert out.row(1, named=True)["deposit_rate"] == pytest.approx(0.5)


class TestValidation:
    def test_missing_required_columns_raise(self) -> None:
        with pytest.raises(ValueError, match="sequence_metadata missing"):
            build_outbreak_features(
                sequence_metadata=pl.DataFrame({"x": [1]}),
                cluster_assignments=_cluster([]),
                cluster_growth=_growth([]),
                cluster_subtypes=_subtype_mix([]),
            )
        with pytest.raises(ValueError, match="cluster_assignments missing"):
            build_outbreak_features(
                sequence_metadata=_seq_meta(
                    [
                        {
                            "accession": "A",
                            "country": "KE",
                            "collection_year": 2020,
                            "subtype": "A1",
                        }
                    ]
                ),
                cluster_assignments=pl.DataFrame({"x": [1]}),
                cluster_growth=_growth([]),
                cluster_subtypes=_subtype_mix([]),
            )
