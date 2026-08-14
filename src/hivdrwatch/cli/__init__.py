"""Unified `hivdrwatch` CLI.

Every pipeline stage has its own `python -m hivdrwatch.<module>` entry
point; this Typer app wraps them behind a single `hivdrwatch
<subcommand>` so a fresh clone (or a JOSS reviewer) hits ONE command
instead of hunting through modules.

Subcommands mirror the pipeline stages 1:1 — ingest, sierra, curate,
cluster, phylo — plus operational helpers (validate, bench, fig,
version). Each subcommand shells out to the underlying module's `main`
so behavior is identical whether you invoke via CLI or module path.
"""

from .app import app

__all__ = ["app"]
