"""Standalone CLI for Sierra — thin wrapper around the Prefect flow.

For orchestrated runs, prefer `python -m hivdrwatch.pipeline.dr_ingest` which
runs inside a Prefect flow context (retries, caching, UI). This module
exists so the pyproject `hivdrwatch-sierra` entrypoint has something to
point at.
"""

from __future__ import annotations

from ..pipeline.dr_ingest import main

if __name__ == "__main__":
    raise SystemExit(main())
