"""Read-side analytics wrapper — DuckDB against Parquet or Postgres.

DuckDB is our SQL analytics engine: it reads Parquet directly (no
loader), joins Parquets to each other cheaply, and can attach a
Postgres schema via the `postgres_scan` extension when the storage
layer is populated.

Design:
  * `open_parquet_view(path, name)` — register a Parquet as a view.
  * `attach_postgres(url, alias='pg')` — attach Postgres tables.
  * `sql(query)` — raw SQL passthrough returning a Polars DataFrame.

We deliberately don't hide SQL — analysts already know SQL; a thin
wrapper is more valuable than a domain-specific query builder.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class DuckDBAnalytics:
    """Thin wrapper around a DuckDB connection.

    Usage::

        db = DuckDBAnalytics()
        db.open_parquet_view(Path("data/interim/dr_results.parquet"), "dr")
        db.open_parquet_view(Path("data/interim/curated/aligned.parquet"), "curated")
        df = db.sql(
            \"\"\"
            SELECT country, subtype, COUNT(*) AS n
            FROM dr JOIN curated USING (accession)
            WHERE level IN ('high', 'intermediate')
            GROUP BY 1, 2 ORDER BY n DESC
            \"\"\"
        )
    """

    def __init__(self, database: str = ":memory:") -> None:
        try:
            import duckdb  # noqa: PLC0415
        except ImportError as e:  # pragma: no cover
            raise ImportError(
                "duckdb not installed; add `duckdb>=1.0` to dependencies"
            ) from e
        self._duckdb = duckdb
        self._con = duckdb.connect(database=database)

    def close(self) -> None:
        self._con.close()

    def __enter__(self) -> DuckDBAnalytics:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def open_parquet_view(self, path: Path, name: str) -> None:
        """Register `path` as SQL view `name`."""
        path = Path(path).resolve()
        if not path.exists():
            raise FileNotFoundError(f"parquet not found: {path}")
        # `?` bind params in DuckDB DDL don't inline into read_parquet('...')
        # cleanly across versions — use a properly quoted literal instead.
        # `self._quote` gates identifiers; path is quoted with SQL-safe
        # single-quote escaping.
        safe_path = str(path).replace("'", "''")
        self._con.execute(
            f"CREATE OR REPLACE VIEW {self._quote(name)} AS "
            f"SELECT * FROM read_parquet('{safe_path}')"
        )
        logger.info("view %s ← %s", name, path.name)

    def attach_postgres(self, url: str, *, alias: str = "pg", read_only: bool = True) -> None:
        """Attach a Postgres database via the postgres_scan extension.

        The extension autoloads on modern DuckDB; older builds may need
        an explicit `INSTALL postgres; LOAD postgres;`.
        """
        self._con.execute("INSTALL postgres")
        self._con.execute("LOAD postgres")
        mode = "READ_ONLY" if read_only else ""
        self._con.execute(
            f"ATTACH ? AS {self._quote(alias)} (TYPE postgres{', ' + mode if mode else ''})",
            [url],
        )
        logger.info("attached postgres as %s (read_only=%s)", alias, read_only)

    def sql(self, query: str):
        """Run `query`, return the result as a Polars DataFrame."""
        try:
            import polars as pl  # noqa: PLC0415
        except ImportError as e:  # pragma: no cover
            raise ImportError("polars required for sql()") from e
        arrow = self._con.execute(query).arrow()
        return pl.from_arrow(arrow)

    def table_exists(self, name: str) -> bool:
        row = self._con.execute(
            "SELECT COUNT(*) FROM duckdb_tables() WHERE table_name = ?", [name]
        ).fetchone()
        return bool(row and row[0])

    @staticmethod
    def _quote(ident: str) -> str:
        """Quote an SQL identifier — no dots allowed for safety."""
        if "." in ident or '"' in ident:
            raise ValueError(f"unsafe identifier: {ident}")
        return f'"{ident}"'


def open_analytics(database: Optional[str] = None) -> DuckDBAnalytics:
    """Convenience factory — returns an in-memory DuckDB by default."""
    return DuckDBAnalytics(database or ":memory:")


__all__ = ["DuckDBAnalytics", "open_analytics"]
