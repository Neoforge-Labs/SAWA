"""The `hivdrwatch` Typer application — dispatch to underlying modules.

We keep this layer *thin* on purpose: subcommands delegate to the same
`main(argv)` functions that `python -m hivdrwatch.<module>` uses. That
means:

  * No new behavior lives here; anything you can do via CLI you can
    also do via module path, and vice versa.
  * The module's argparse remains the single source of truth for
    flags — Typer forwards `--` args through with `context_settings=
    {"allow_extra_args": True, "ignore_unknown_options": True}`.

The result is a top-level UX that looks like `docker`, `git`, or `dvc`:
one entry point + discoverable subcommands.
"""

from __future__ import annotations

import importlib
from typing import Callable

import typer

app = typer.Typer(
    name="hivdrwatch",
    help=(
        "SAWA — Sub-Saharan Africa Watch for Antiretroviral resistance. "
        "Open real-time HIV-1 drug-resistance surveillance. Every "
        "subcommand mirrors a pipeline stage. "
        "(CLI is `hivdrwatch`, brand is SAWA.)"
    ),
    no_args_is_help=True,
    add_completion=False,
)


# Map subcommand → (module path, help string). Kept as a table so
# adding a new stage is a one-line change.
_SUBCOMMANDS: dict[str, tuple[str, str]] = {
    "ingest": (
        "hivdrwatch.ingest.__main__",
        "Merge LANL Table exports + GenBank Entrez pulls → canonical bronze Parquet.",
    ),
    "sierra": (
        "hivdrwatch.pipeline.dr_ingest",
        "Call Stanford HIVdb Sierra resistance for every ingested sequence.",
    ),
    "curate": (
        "hivdrwatch.pipeline.curate",
        "Align to HXB2, trim, QC, dedup → silver Parquet + FASTA.",
    ),
    "cluster": (
        "hivdrwatch.pipeline.cluster",
        "tn93 pairwise distances + connected-components clustering → gold cluster set.",
    ),
    "validate-clusters": (
        "hivdrwatch.pipeline.cluster_validate",
        "Score our clusterer (and optionally HIV-TRACE) against simulated ground truth.",
    ),
    "phylo": (
        "hivdrwatch.pipeline.phylo",
        "Per-subtype IQ-TREE 2 + TreeTime → dated per-subtype trees.",
    ),
    "bench-ingest": (
        "benchmarks.scale_ingest",
        "Scale benchmark for the ingest stage — synthetic LANL TSV parser throughput.",
    ),
    "bench-align": (
        "benchmarks.scale_align",
        "Scale benchmark for the align stage — mappy pairwise-to-HXB2 throughput.",
    ),
    "bench-cluster": (
        "benchmarks.scale_cluster",
        "Scale benchmark for the cluster stage — tn93 + assembly wall time and peak RSS.",
    ),
}


def _make_forwarder(module_path: str, help_text: str) -> Callable[..., None]:
    """Wrap `module.main(argv)` in a Typer command that forwards extra args.

    The underlying module returns an int exit code; we bubble it up
    via `typer.Exit` so the shell sees the same code as
    `python -m <module>` would produce.
    """

    def _cmd(
        ctx: typer.Context = typer.Context,  # type: ignore[assignment]
    ) -> None:
        module = importlib.import_module(module_path)
        if not hasattr(module, "main"):
            raise typer.BadParameter(
                f"underlying module {module_path!r} has no `main(argv)`; "
                "please open an issue."
            )
        code = module.main(list(ctx.args))
        raise typer.Exit(code=int(code or 0))

    _cmd.__doc__ = help_text
    return _cmd


for name, (module_path, help_text) in _SUBCOMMANDS.items():
    app.command(
        name=name,
        help=help_text,
        context_settings={
            "allow_extra_args": True,
            "ignore_unknown_options": True,
        },
    )(_make_forwarder(module_path, help_text))


def _safe_import_version(module_path: str, attr: str) -> str:
    """Import `module_path` and return `attr`; return 'unknown' on any failure.

    Version-dump commands should never crash because a dependency isn't
    installed — degrade gracefully to `unknown` per component instead.
    """
    try:
        module = importlib.import_module(module_path)
        return str(getattr(module, attr))
    except Exception:  # noqa: BLE001 - defensive at report boundary
        return "unknown"


@app.command(help="Print the installed hivdrwatch version and schema versions, then exit.")
def version() -> None:
    """Version and schema-versions dump — the first thing a reviewer types."""
    import importlib.metadata

    try:
        pkg_version = importlib.metadata.version("hivdrwatch")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover - editable dev only
        pkg_version = "0.0.0+dev"

    ingest_version = _safe_import_version("hivdrwatch.ingest.schema", "SCHEMA_VERSION")
    cluster_version = _safe_import_version("hivdrwatch.cluster.schema", "CLUSTER_SCHEMA_VERSION")
    storage_version = _safe_import_version("hivdrwatch.storage.schema", "STORAGE_SCHEMA_VERSION")
    bench_version = _safe_import_version("benchmarks", "BENCHMARK_SCHEMA_VERSION")

    lines = [
        f"hivdrwatch          {pkg_version}",
        f"ingest schema       {ingest_version}",
        f"cluster schema      {cluster_version}",
        f"storage schema      {storage_version}",
        f"benchmark schema    {bench_version}",
    ]
    typer.echo("\n".join(lines))


@app.command(help="List every registered subcommand with a one-line summary.")
def commands() -> None:
    """Nice-to-have alternative to `--help` for a wide subcommand table."""
    width = max(len(name) for name in _SUBCOMMANDS) + 2
    for name, (_, help_text) in _SUBCOMMANDS.items():
        typer.echo(f"  {name.ljust(width)}{help_text}")
    typer.echo(f"  {'version'.ljust(width)}Print version and schema info.")


__all__ = ["app"]
