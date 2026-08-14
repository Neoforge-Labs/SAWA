"""DuckDB analytics wrapper tests."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from hivdrwatch.storage.duckdb import DuckDBAnalytics, open_analytics


def _write_parquet(path: Path, data: dict) -> Path:
    pl.DataFrame(data).write_parquet(path)
    return path


class TestParquetView:
    def test_open_and_query(self, tmp_path: Path) -> None:
        p = _write_parquet(
            tmp_path / "dr.parquet",
            {"accession": ["A", "B", "C"], "level": ["high", "susceptible", "high"]},
        )
        with open_analytics() as db:
            db.open_parquet_view(p, "dr")
            result = db.sql("SELECT COUNT(*) AS n FROM dr WHERE level = 'high'")
            assert result["n"][0] == 2

    def test_missing_parquet_raises(self, tmp_path: Path) -> None:
        with open_analytics() as db, pytest.raises(FileNotFoundError):
            db.open_parquet_view(tmp_path / "nope.parquet", "x")

    def test_multiple_views_join(self, tmp_path: Path) -> None:
        _write_parquet(tmp_path / "a.parquet", {"k": [1, 2, 3], "va": [10, 20, 30]})
        _write_parquet(tmp_path / "b.parquet", {"k": [2, 3, 4], "vb": [200, 300, 400]})
        with open_analytics() as db:
            db.open_parquet_view(tmp_path / "a.parquet", "a")
            db.open_parquet_view(tmp_path / "b.parquet", "b")
            result = db.sql("SELECT SUM(va + vb) AS s FROM a JOIN b USING (k)")
            assert result["s"][0] == (20 + 200) + (30 + 300)


class TestIdentifierSafety:
    def test_bad_ident_rejected(self, tmp_path: Path) -> None:
        p = _write_parquet(tmp_path / "a.parquet", {"k": [1]})
        with open_analytics() as db, pytest.raises(ValueError, match="unsafe identifier"):
            db.open_parquet_view(p, "bad.name")
        with open_analytics() as db, pytest.raises(ValueError):
            db.open_parquet_view(p, 'bad"name')


class TestContextManager:
    def test_closes_on_exit(self, tmp_path: Path) -> None:
        p = _write_parquet(tmp_path / "a.parquet", {"k": [1]})
        db = DuckDBAnalytics()
        db.open_parquet_view(p, "a")
        db.close()
        # Second close should not crash.
        db.close()
