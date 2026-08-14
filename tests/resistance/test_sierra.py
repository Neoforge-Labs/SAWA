"""Unit tests for the Sierra client — cache, parsing, batching, offline mode.

Network access is mocked via pytest-httpx. If Sierra's GraphQL shape ever
changes, the parser tests will fail loudly against the fixture payloads.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import httpx
import pytest

from hivdrwatch.resistance.cache import SierraCache, sequence_sha256
from hivdrwatch.resistance.schema import DrugClass, ResistanceLevel
from hivdrwatch.resistance.sierra import (
    AlgorithmVersion,
    SierraClient,
    SierraError,
    _flatten_mutations,
    _maybe_date,
    _parse_algorithm,
    _parse_dr_results,
    _parse_subtype,
)


# ---------- fixtures ----------


ONE_SEQ_PAYLOAD = {
    "currentVersion": {"text": "9.6", "publishDate": "2024-08-15"},
    "analysis": {
        "inputSequence": {"header": "TESTSEQ"},
        "bestMatchingSubtype": {"display": "A1"},
        "drugResistance": [
            {
                "gene": {"name": "RT"},
                "drugScores": [
                    {
                        "drugClass": {"name": "NRTI"},
                        "drug": {
                            "name": "LMV",
                            "displayAbbr": "3TC",
                            "fullName": "lamivudine",
                        },
                        "score": 60.0,
                        "text": "High-Level Resistance",
                    },
                    {
                        "drugClass": {"name": "NNRTI"},
                        "drug": {
                            "name": "EFV",
                            "displayAbbr": "EFV",
                            "fullName": "efavirenz",
                        },
                        "score": 0.0,
                        "text": "Susceptible",
                    },
                ],
                "mutationsByTypes": [
                    {
                        "mutationType": "NRTI",
                        "mutations": [{"text": "M184V"}],
                    },
                    {
                        "mutationType": "NNRTI",
                        "mutations": [{"text": "K103N"}],
                    },
                ],
            }
        ],
    },
}


BATCH_RESPONSE = {
    "data": {
        "currentVersion": {"text": "9.6", "publishDate": "2024-08-15"},
        "viewer": {
            "sequenceAnalysis": [
                ONE_SEQ_PAYLOAD["analysis"],
                ONE_SEQ_PAYLOAD["analysis"],
            ]
        },
    }
}


# ---------- cache ----------


class TestSierraCache:
    def test_put_get_roundtrip(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        sha = sequence_sha256("ACGT")
        cache.put(sha, {"hello": "world"})
        assert cache.has(sha)
        assert cache.get(sha) == {"hello": "world"}

    def test_sharded_layout(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        sha = sequence_sha256("ACGT")
        cache.put(sha, {"x": 1})
        assert (tmp_path / sha[:2] / f"{sha}.json").exists()

    def test_missing_lists_gaps(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        present = sequence_sha256("A")
        cache.put(present, {"x": 1})
        absent = sequence_sha256("B")
        assert cache.missing([present, absent]) == [absent]

    def test_algorithm_version_roundtrip(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        cache.write_algorithm("9.6", "2024-08-15")
        assert cache.read_algorithm() == {
            "version": "9.6",
            "publish_date": "2024-08-15",
        }

    def test_corrupt_entry_ignored(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        sha = sequence_sha256("A")
        p = cache.path_for(sha)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("{not json")
        assert cache.get(sha) is None


class TestSequenceSha256:
    def test_normalization(self) -> None:
        # U→T, case insensitive, whitespace-stripped → same hash
        assert sequence_sha256("acgt") == sequence_sha256("ACGT")
        assert sequence_sha256("ACGU") == sequence_sha256("ACGT")
        assert sequence_sha256("  ACGT  ") == sequence_sha256("ACGT")

    def test_different_sequences_hash_differently(self) -> None:
        assert sequence_sha256("ACGT") != sequence_sha256("ACGA")


# ---------- pure parsers ----------


class TestParseAlgorithm:
    def test_extracts_version_and_date(self) -> None:
        algo = _parse_algorithm(ONE_SEQ_PAYLOAD)
        assert algo.version == "9.6"
        assert algo.publish_date == date(2024, 8, 15)

    def test_missing_current_version(self) -> None:
        algo = _parse_algorithm({})
        assert algo.version == "unknown"
        assert algo.publish_date is None


class TestParseDrResults:
    def test_yields_row_per_drug(self) -> None:
        algo = AlgorithmVersion("9.6", date(2024, 8, 15))
        rows = _parse_dr_results(
            ONE_SEQ_PAYLOAD,
            accession="AF004885",
            sha="a" * 64,
            algo=algo,
            called_at=date(2026, 8, 11),
        )
        assert len(rows) == 2
        by_drug = {r.drug: r for r in rows}
        assert set(by_drug) == {"3TC", "EFV"}
        assert by_drug["3TC"].drug_class == DrugClass.NRTI
        assert by_drug["3TC"].level == ResistanceLevel.HIGH
        assert by_drug["3TC"].score == 60.0
        assert "M184V" in by_drug["3TC"].mutations
        assert by_drug["EFV"].level == ResistanceLevel.SUSCEPTIBLE

    def test_empty_analysis_yields_nothing(self) -> None:
        algo = AlgorithmVersion("9.6", None)
        assert (
            _parse_dr_results(
                {"analysis": {"drugResistance": []}},
                accession="MK123456",
                sha="a" * 64,
                algo=algo,
                called_at=date(2026, 8, 11),
            )
            == []
        )

    def test_skips_unknown_drug_class(self) -> None:
        algo = AlgorithmVersion("9.6", None)
        payload = {
            "analysis": {
                "drugResistance": [
                    {
                        "gene": {"name": "RT"},
                        "drugScores": [
                            {
                                "drugClass": {"name": "FUTUREBLOCKER"},
                                "drug": {"displayAbbr": "XYZ"},
                                "score": 10,
                                "text": "Susceptible",
                            }
                        ],
                        "mutationsByTypes": [],
                    }
                ]
            }
        }
        assert (
            _parse_dr_results(
                payload,
                accession="MK123456",
                sha="a" * 64,
                algo=algo,
                called_at=date(2026, 8, 11),
            )
            == []
        )


class TestParseSubtype:
    def test_yields_subtype(self) -> None:
        algo = AlgorithmVersion("9.6", None)
        st = _parse_subtype(
            ONE_SEQ_PAYLOAD, "AF004885", "a" * 64, algo, date(2026, 8, 11)
        )
        assert st is not None
        assert st.subtype == "A1"

    def test_none_when_missing(self) -> None:
        algo = AlgorithmVersion("9.6", None)
        assert (
            _parse_subtype(
                {"analysis": {}}, "AF004885", "a" * 64, algo, date(2026, 8, 11)
            )
            is None
        )


class TestFlattenMutations:
    def test_dedupes_and_sorts(self) -> None:
        m = _flatten_mutations(
            [
                {"mutationType": "NRTI", "mutations": [{"text": "M184V"}, {"text": "K65R"}]},
                {"mutationType": "NNRTI", "mutations": [{"text": "M184V"}]},
            ]
        )
        assert m == ["K65R", "M184V"]

    def test_empty(self) -> None:
        assert _flatten_mutations([]) == []


class TestMaybeDate:
    def test_iso_date(self) -> None:
        assert _maybe_date("2024-08-15") == date(2024, 8, 15)

    def test_iso_datetime(self) -> None:
        assert _maybe_date("2024-08-15T00:00:00Z") == date(2024, 8, 15)

    def test_none_on_garbage(self) -> None:
        assert _maybe_date("not a date") is None
        assert _maybe_date(None) is None


# ---------- integration: client with mocked HTTP ----------


class TestSierraClientOfflineMode:
    def test_offline_cache_hit(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        seq = "ACGT" * 60
        sha = sequence_sha256(seq)
        cache.put(sha, ONE_SEQ_PAYLOAD)
        cache.write_algorithm("9.6", "2024-08-15")

        client = SierraClient(cache=cache, offline=True)
        dr, subtypes = client.analyze_sequences([("AF004885", seq)])
        assert len(dr) == 2
        assert len(subtypes) == 1
        assert dr[0].algorithm_version == "9.6"

    def test_offline_cache_miss_raises(self, tmp_path: Path) -> None:
        cache = SierraCache(tmp_path)
        client = SierraClient(cache=cache, offline=True)
        with pytest.raises(SierraError, match="offline"):
            client.analyze_sequences([("AF004885", "ACGT" * 60)])


class TestSierraClientHTTP:
    def test_fetch_and_cache(self, tmp_path: Path, httpx_mock) -> None:
        httpx_mock.add_response(
            url="https://hivdb.stanford.edu/graphql",
            json=BATCH_RESPONSE,
        )
        cache = SierraCache(tmp_path)
        with httpx.Client() as http:
            client = SierraClient(cache=cache, client=http)
            dr, subtypes = client.analyze_sequences(
                [("AF004885", "ACGT" * 60), ("MK123456", "GGGG" * 60)]
            )
        assert len(dr) == 4  # 2 drugs × 2 sequences
        assert len(subtypes) == 2
        # Both entries cached — a second call would be pure cache.
        assert cache.has(sequence_sha256("ACGT" * 60))
        assert cache.has(sequence_sha256("GGGG" * 60))

    def test_dedups_same_sequence_across_accessions(
        self, tmp_path: Path, httpx_mock
    ) -> None:
        """Two accessions with the same sequence should only hit Sierra once."""
        single_response = {
            "data": {
                "currentVersion": BATCH_RESPONSE["data"]["currentVersion"],
                "viewer": {"sequenceAnalysis": [ONE_SEQ_PAYLOAD["analysis"]]},
            }
        }
        httpx_mock.add_response(
            url="https://hivdb.stanford.edu/graphql", json=single_response
        )
        cache = SierraCache(tmp_path)
        with httpx.Client() as http:
            client = SierraClient(cache=cache, client=http)
            dr, subtypes = client.analyze_sequences(
                [("AF004885", "ACGT" * 60), ("MK123456", "ACGT" * 60)]
            )
        # 2 drugs × 2 accessions expanded from 1 cached payload
        assert len(dr) == 4
        assert {r.accession for r in dr} == {"AF004885", "MK123456"}
        # httpx_mock will raise on unmatched requests, so only 1 POST happened.
