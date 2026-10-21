"""Health, version, docs endpoints."""

from __future__ import annotations

from fastapi.testclient import TestClient

from hivdrwatch.cluster.schema import CLUSTER_SCHEMA_VERSION
from hivdrwatch.ingest.schema import SCHEMA_VERSION as INGEST_SCHEMA_VERSION
from hivdrwatch.storage.schema import STORAGE_SCHEMA_VERSION


class TestHealth:
    def test_returns_ok(self, client: TestClient) -> None:
        r = client.get("/health")
        assert r.status_code == 200
        body = r.json()
        assert body["status"] == "ok"
        assert "checked_at" in body

    def test_request_id_header_echoed(self, client: TestClient) -> None:
        r = client.get("/health", headers={"X-Request-Id": "test-rid-abc"})
        assert r.status_code == 200
        assert r.headers["X-Request-Id"] == "test-rid-abc"

    def test_generated_request_id_when_absent(self, client: TestClient) -> None:
        r = client.get("/health")
        assert r.status_code == 200
        assert r.headers.get("X-Request-Id")
        assert len(r.headers["X-Request-Id"]) >= 8


class TestVersion:
    def test_returns_versions(self, client: TestClient) -> None:
        r = client.get("/version")
        assert r.status_code == 200
        body = r.json()
        # In a full 3.12 env with every schema module importable these
        # report their module's own SCHEMA_VERSION. In restricted envs
        # (missing bio deps, older Python) they degrade to "unknown" per
        # the /version endpoint's graceful-degradation contract — but the
        # response must still be 200 and every key must be present.
        expected = {
            "ingest_schema": INGEST_SCHEMA_VERSION,
            "cluster_schema": CLUSTER_SCHEMA_VERSION,
            "storage_schema": STORAGE_SCHEMA_VERSION,
        }
        for key, version in expected.items():
            assert body[key] in {version, "unknown"}, f"{key}={body[key]!r}"
        # storage schema is pure-Python and always importable — this one
        # is the canary for "did the version endpoint even work?"
        assert body["storage_schema"] == STORAGE_SCHEMA_VERSION
        assert body["api_schema"] == "0.1.0"


class TestOpenAPI:
    def test_openapi_json_available(self, client: TestClient) -> None:
        r = client.get("/openapi.json")
        assert r.status_code == 200
        spec = r.json()
        assert spec["openapi"].startswith("3.")
        # Every route we registered should appear.
        paths = set(spec["paths"])
        assert "/health" in paths
        assert "/version" in paths
        assert "/sequences" in paths
        assert "/dr-calls" in paths
        assert "/clusters" in paths

    def test_swagger_docs_available(self, client: TestClient) -> None:
        r = client.get("/docs")
        assert r.status_code == 200
        assert "swagger" in r.text.lower()
