"""Route-level tests — read against a seeded SQLite."""

from __future__ import annotations

from fastapi.testclient import TestClient


class TestSequences:
    def test_list_all(self, client: TestClient) -> None:
        r = client.get("/sequences")
        assert r.status_code == 200
        body = r.json()
        assert body["total"] == 5
        assert body["limit"] == 100
        assert body["offset"] == 0
        assert len(body["items"]) == 5

    def test_filter_by_country(self, client: TestClient) -> None:
        r = client.get("/sequences", params={"country": "KE"})
        assert r.status_code == 200
        assert r.json()["total"] == 3

    def test_filter_by_subtype(self, client: TestClient) -> None:
        r = client.get("/sequences", params={"subtype": "C"})
        assert r.status_code == 200
        assert r.json()["total"] == 2

    def test_filter_by_year(self, client: TestClient) -> None:
        r = client.get("/sequences", params={"collection_year": 2022})
        assert r.status_code == 200
        assert r.json()["total"] == 1

    def test_pagination(self, client: TestClient) -> None:
        r = client.get("/sequences", params={"limit": 2, "offset": 0})
        assert r.status_code == 200
        body = r.json()
        assert body["total"] == 5
        assert len(body["items"]) == 2
        assert body["items"][0]["accession"] == "AB000000"

        r = client.get("/sequences", params={"limit": 2, "offset": 2})
        assert r.json()["items"][0]["accession"] == "AB000002"

    def test_invalid_country_length_rejected(self, client: TestClient) -> None:
        r = client.get("/sequences", params={"country": "KEN"})
        assert r.status_code == 422


class TestDRCalls:
    def test_list_all(self, client: TestClient) -> None:
        r = client.get("/dr-calls")
        assert r.status_code == 200
        assert r.json()["total"] == 5

    def test_filter_by_level(self, client: TestClient) -> None:
        r = client.get("/dr-calls", params={"level": "high"})
        assert r.status_code == 200
        # 3 sequences at index 0, 2, 4 → level=high
        assert r.json()["total"] == 3

    def test_filter_by_drug_class(self, client: TestClient) -> None:
        r = client.get("/dr-calls", params={"drug_class": "PI"})
        assert r.status_code == 200
        assert r.json()["total"] == 0  # seed only has NNRTI

    def test_filter_by_accession_uppercased(self, client: TestClient) -> None:
        r = client.get("/dr-calls", params={"accession": "ab000000"})
        assert r.status_code == 200
        body = r.json()
        assert body["total"] == 1
        assert body["items"][0]["accession"] == "AB000000"


class TestClusters:
    def test_list_all(self, client: TestClient) -> None:
        r = client.get("/clusters")
        assert r.status_code == 200
        assert r.json()["total"] == 4

    def test_filter_by_cluster_id(self, client: TestClient) -> None:
        r = client.get("/clusters", params={"cluster_id": "C-0001"})
        assert r.status_code == 200
        assert r.json()["total"] == 2

    def test_min_cluster_size(self, client: TestClient) -> None:
        r = client.get("/clusters", params={"min_cluster_size": 3})
        assert r.status_code == 200
        # All seed clusters have size 2 — 3+ filter yields none.
        assert r.json()["total"] == 0
