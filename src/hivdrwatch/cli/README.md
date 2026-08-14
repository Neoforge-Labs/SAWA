# `cli` — the unified `hivdrwatch` entry point

Every pipeline stage still has its own `python -m hivdrwatch.<module>` entry point; this Typer app wraps them behind one command:

```bash
hivdrwatch --help
hivdrwatch ingest --lanl-dir data/raw/lanl --out data/interim/east_africa_pol.parquet
hivdrwatch sierra --input data/interim/east_africa_pol.parquet
hivdrwatch curate --input data/interim/east_africa_pol.parquet --lanl-dir data/raw/lanl
hivdrwatch cluster --aligned-fasta data/interim/curate/aligned.fasta ...
hivdrwatch validate-clusters --with-hivtrace
hivdrwatch phylo --curated data/interim/curated/aligned.parquet
hivdrwatch bench-ingest / bench-align / bench-cluster
hivdrwatch version
hivdrwatch commands
```

## Design

Thin dispatch layer. Every subcommand forwards all extra args to the underlying module's `main(argv)`:

```
hivdrwatch ingest --skip-genbank ...
    ⇓
hivdrwatch.ingest.__main__.main(["--skip-genbank", ...])
```

That means:

- **One source of truth** — argparse flags stay in the module.
- **Discoverability** — one command reviewers grep for, `hivdrwatch --help` lists everything.
- **Backwards compatibility** — `python -m hivdrwatch.pipeline.dr_ingest ...` still works.

## Tests

```bash
pytest tests/cli -v
```
