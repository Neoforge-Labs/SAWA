"""Unit tests for the HIV-TRACE wrapper.

We don't hit Docker in tests — the wrapper's job is to shell out and
parse a JSON. We test the parser directly with fixtures shaped like
hivnetworkcsv's `-J` output, plus the availability guard.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hivdrwatch.cluster.hivtrace import (
    DEFAULT_IMAGE,
    HIVTraceNotAvailable,
    parse_network_json,
    run_hivtrace,
)


def _write_network_legacy(path: Path, nodes: list[dict]) -> Path:
    """Legacy list-of-dicts shape."""
    path.write_text(json.dumps({"Nodes": nodes, "Edges": []}))
    return path


def _write_network_compact(
    path: Path, ids: list[str], keys_map: dict, values: list[int]
) -> Path:
    """Compact `-J` shape emitted by hivclustering 1.9.x."""
    payload = {
        "Nodes": {
            "id": ids,
            "cluster": {"keys": keys_map, "values": values},
        },
        "Edges": [],
    }
    path.write_text(json.dumps(payload))
    return path


class TestParseNetworkJsonLegacy:
    """Coverage for the old list-of-dicts `-j` (pretty-print) shape."""

    def test_basic_two_clusters(self, tmp_path: Path) -> None:
        p = _write_network_legacy(
            tmp_path / "net.json",
            [
                {"id": "S0", "cluster": 1},
                {"id": "S1", "cluster": 1},
                {"id": "S2", "cluster": 2},
                {"id": "S3", "cluster": 2},
            ],
        )
        mem = parse_network_json(p, ["S0", "S1", "S2", "S3"])
        assert mem["S0"] == mem["S1"]
        assert mem["S2"] == mem["S3"]
        assert mem["S0"] != mem["S2"]

    def test_singletons_added_back(self, tmp_path: Path) -> None:
        p = _write_network_legacy(
            tmp_path / "net.json",
            [
                {"id": "S0", "cluster": 1},
                {"id": "S1", "cluster": 1},
            ],
        )
        mem = parse_network_json(p, ["S0", "S1", "S2", "S3"])
        assert set(mem.keys()) == {"S0", "S1", "S2", "S3"}
        assert mem["S2"].startswith("SINGLETON_")
        assert mem["S3"].startswith("SINGLETON_")
        assert mem["S2"] != mem["S3"]

    def test_ids_are_stable_stringified(self, tmp_path: Path) -> None:
        p = _write_network_legacy(
            tmp_path / "net.json",
            [{"id": "S0", "cluster": 0}, {"id": "S1", "cluster": 0}],
        )
        mem = parse_network_json(p, ["S0", "S1"])
        assert mem["S0"] == "HT-0000"
        assert mem["S1"] == "HT-0000"

    def test_missing_fields_are_skipped(self, tmp_path: Path) -> None:
        p = _write_network_legacy(
            tmp_path / "net.json",
            [
                {"id": "S0"},              # no cluster
                {"cluster": 1},              # no id
                {"id": "S1", "cluster": 1},
            ],
        )
        mem = parse_network_json(p, ["S0", "S1", "S2"])
        assert mem["S1"] == "HT-0001"
        assert mem["S0"].startswith("SINGLETON_")
        assert mem["S2"].startswith("SINGLETON_")


class TestParseNetworkJsonCompact:
    """Coverage for the compact `-J` shape emitted by hivclustering 1.9.x."""

    def test_basic_two_clusters(self, tmp_path: Path) -> None:
        # Cluster 1 has S0, S1; cluster 2 has S2, S3. Compact idx 0 → 1,
        # compact idx 1 → 2.
        p = _write_network_compact(
            tmp_path / "net.json",
            ids=["S0", "S1", "S2", "S3"],
            keys_map={"0": 1, "1": 2},
            values=[0, 0, 1, 1],
        )
        mem = parse_network_json(p, ["S0", "S1", "S2", "S3"])
        assert mem["S0"] == "HT-0001"
        assert mem["S1"] == "HT-0001"
        assert mem["S2"] == "HT-0002"
        assert mem["S3"] == "HT-0002"

    def test_singletons_added_back(self, tmp_path: Path) -> None:
        p = _write_network_compact(
            tmp_path / "net.json",
            ids=["S0", "S1"],           # only these appear in Nodes
            keys_map={"0": 1},
            values=[0, 0],
        )
        mem = parse_network_json(p, ["S0", "S1", "S2", "S3"])
        assert mem["S0"] == mem["S1"] == "HT-0001"
        assert mem["S2"].startswith("SINGLETON_")
        assert mem["S3"].startswith("SINGLETON_")
        assert mem["S2"] != mem["S3"]

    def test_falls_back_to_compact_idx_when_keys_missing(self, tmp_path: Path) -> None:
        # No `keys` mapping — parser should use the compact idx directly.
        p = tmp_path / "net.json"
        p.write_text(
            json.dumps(
                {
                    "Nodes": {
                        "id": ["S0", "S1"],
                        "cluster": {"values": [3, 3]},   # no `keys`
                    },
                    "Edges": [],
                }
            )
        )
        mem = parse_network_json(p, ["S0", "S1"])
        assert mem["S0"] == mem["S1"] == "HT-0003"

    def test_cluster_as_plain_list(self, tmp_path: Path) -> None:
        # Some builds emit `cluster` as a plain list (no keys/values dict).
        p = tmp_path / "net.json"
        p.write_text(
            json.dumps(
                {
                    "Nodes": {"id": ["S0", "S1"], "cluster": [7, 7]},
                    "Edges": [],
                }
            )
        )
        mem = parse_network_json(p, ["S0", "S1"])
        assert mem["S0"] == mem["S1"] == "HT-0007"


class TestAvailabilityGuard:
    def test_run_raises_when_image_missing(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Force is_available → False by monkeypatching the module symbol.
        import hivdrwatch.cluster.hivtrace as ht

        monkeypatch.setattr(ht, "is_available", lambda image=DEFAULT_IMAGE: False)

        distances = tmp_path / "distances.csv"
        distances.write_text("id1,id2,distance\nS0,S1,0.01\n")
        network = tmp_path / "network.json"

        with pytest.raises(HIVTraceNotAvailable):
            run_hivtrace(distances, network, ["S0", "S1"])

    def test_run_raises_on_missing_csv(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import hivdrwatch.cluster.hivtrace as ht

        monkeypatch.setattr(ht, "is_available", lambda image=DEFAULT_IMAGE: True)

        with pytest.raises(FileNotFoundError):
            run_hivtrace(tmp_path / "nope.csv", tmp_path / "out.json", ["S0"])
