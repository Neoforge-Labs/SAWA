"""Prefect flow: LANL + GenBank → merged bronze parquet.

Wraps the existing `hivdrwatch.ingest.__main__` logic in a Prefect flow
with per-country GenBank checkpointing (via watermark.py) and MLflow
tracking.

Incremental behavior:
- Reads `<genbank_raw_dir>/.watermarks.json`.
- For each country, passes `since=<watermark>` to GenBank so only newly
  deposited sequences come back.
- On success, bumps the watermark to today.

If a country fetch fails mid-way, its watermark isn't bumped; the next
run re-fetches from the same starting point.

Run locally:
    prefect server start &                      # if not already running
    python -m hivdrwatch.pipeline.ingest \\
        --email you@example.com \\
        --countries KE UG TZ \\
        --lanl-dir data/raw/lanl \\
        --genbank-dir data/raw/genbank \\
        --out data/interim/east_africa_pol.parquet
"""

from __future__ import annotations

import argparse
import logging
import os
from datetime import date, datetime
from pathlib import Path

import polars as pl
from prefect import flow, get_run_logger, task

from ..ingest.genbank import EAST_AFRICA_COUNTRIES, GenBankIngester
from ..ingest.lanl import LANLIngester
from ..ingest.merge import merge_sources
from ..ingest.watermark import (
    bump_watermark,
    get_since,
    read_watermarks,
    write_watermarks,
)
from ..tracking import tracked_run

FIRST_EVER_PULL_FLOOR = date(1980, 1, 1)


# ------------------------------------------------------------------
# Tasks
# ------------------------------------------------------------------


@task
def load_lanl(lanl_dir: Path, min_length: int) -> pl.DataFrame:
    log = get_run_logger()
    df = LANLIngester(lanl_dir).load_directory(min_length=min_length)
    log.info("LANL: %d rows from %s", df.height, lanl_dir)
    return df


@task(retries=2, retry_delay_seconds=60)
def fetch_genbank_country(
    country: str,
    email: str,
    genbank_dir: Path,
    since: date | None,
    until: date | None,
    min_length: int,
) -> list[dict]:
    """Fetch one country's GenBank records incrementally.

    Retried by Prefect on transient failures (E-utilities does hiccup).
    The `raw_dir` cache means retries don't re-download batches that
    already landed successfully in the previous attempt.
    """
    log = get_run_logger()
    api_key = os.environ.get("NCBI_API_KEY")
    ingester = GenBankIngester(email=email, raw_dir=genbank_dir, api_key=api_key)
    records = []
    for rec in ingester.fetch_country(
        country_code=country,
        since=since,
        until=until,
        min_length=min_length,
    ):
        records.append(rec.model_dump(mode="json"))
    log.info(
        "GenBank %s: %d validated records (since=%s until=%s)",
        country,
        len(records),
        since,
        until,
    )
    return records


@task
def merge_and_write(
    genbank_records: list[dict],
    lanl_df: pl.DataFrame,
    out_path: Path,
) -> Path:
    log = get_run_logger()
    if genbank_records:
        genbank_df = pl.from_dicts(genbank_records)
    else:
        genbank_df = pl.DataFrame({"accession": []})
    merged = merge_sources(genbank_df, lanl_df, out_path)
    log.info("merged: %d rows → %s", merged.height, out_path)
    return out_path


@task
def bump_country_watermarks(
    per_country_records: dict[str, list[dict]],
    genbank_dir: Path,
    watermarks: dict[str, date],
) -> Path:
    """Update watermark for every country that succeeded this run."""
    log = get_run_logger()
    updated = watermarks
    today = date.today()
    for country, records in per_country_records.items():
        # Only bump on success (records may be an empty list — that's fine,
        # it means GenBank had nothing new since the last watermark).
        if records is not None:
            updated = bump_watermark(updated, country, today)
    path = write_watermarks(genbank_dir, updated)
    log.info("wrote watermarks → %s", path)
    return path


# ------------------------------------------------------------------
# Flow
# ------------------------------------------------------------------


@flow(name="ingest")
def ingest_flow(
    email: str | None,
    countries: list[str],
    lanl_dir: Path,
    genbank_dir: Path,
    out_path: Path,
    *,
    since: date | None = None,
    until: date | None = None,
    min_length: int = 800,
    skip_genbank: bool = False,
    to_postgres: bool | None = None,
) -> Path:
    """Top-level ingest — LANL, optional GenBank, merge, write parquet."""
    log = get_run_logger()

    lanl_df = load_lanl(lanl_dir, min_length)

    per_country_records: dict[str, list[dict]] = {}
    genbank_records: list[dict] = []
    watermarks: dict[str, date] = {}

    if skip_genbank:
        log.info("skipping GenBank (skip_genbank=True)")
    else:
        if not email:
            raise ValueError("email is required for GenBank ingestion")
        watermarks = read_watermarks(genbank_dir)
        for country in countries:
            eff_since = since or get_since(watermarks, country, fallback=FIRST_EVER_PULL_FLOOR)
            records = fetch_genbank_country(
                country=country,
                email=email,
                genbank_dir=genbank_dir,
                since=eff_since,
                until=until,
                min_length=min_length,
            )
            per_country_records[country] = records
            genbank_records.extend(records)

    out = merge_and_write(genbank_records, lanl_df, out_path)

    if per_country_records:
        bump_country_watermarks(per_country_records, genbank_dir, watermarks)

    # MLflow: one run per ingest with per-country counts.
    with tracked_run(
        experiment="ingest",
        run_name=f"ingest-{date.today().isoformat()}",
        tags={
            "flow": "ingest",
            "skip_genbank": str(skip_genbank),
        },
    ) as run:
        run.log_params(
            {
                "countries": ",".join(countries),
                "min_length": min_length,
                "since": str(since or ""),
                "until": str(until or ""),
            }
        )
        run.log_metrics(
            {
                "n_lanl": float(lanl_df.height),
                "n_genbank_total": float(len(genbank_records)),
            }
        )
        for country, recs in per_country_records.items():
            run.log_metrics({f"n_genbank_{country}": float(len(recs))})
        if out.exists():
            run.log_artifact(out)

    # Storage tail — load bronze Parquet into Postgres if enabled.
    # Never crashes the flow; Parquet on disk remains source of truth.
    from ._storage_tail import storage_tail

    storage_tail("ingest", out, enabled=to_postgres)

    return out


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------


def _parse_date(s: str) -> date:
    return datetime.strptime(s, "%Y-%m-%d").date()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.pipeline.ingest")
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL"),
        help="contact email for NCBI E-utilities (required unless --skip-genbank; "
        "reads NCBI_EMAIL env var by default)",
    )
    parser.add_argument(
        "--countries",
        nargs="+",
        default=list(EAST_AFRICA_COUNTRIES.keys()),
        help="ISO 3166-1 alpha-2 country codes",
    )
    parser.add_argument("--since", type=_parse_date, default=None)
    parser.add_argument("--until", type=_parse_date, default=None)
    parser.add_argument("--lanl-dir", type=Path, default=Path("data/raw/lanl"))
    parser.add_argument("--genbank-dir", type=Path, default=Path("data/raw/genbank"))
    parser.add_argument(
        "--out", type=Path, default=Path("data/interim/east_africa_pol.parquet")
    )
    parser.add_argument("--min-length", type=int, default=800)
    parser.add_argument("--skip-genbank", action="store_true")
    parser.add_argument(
        "--to-postgres",
        action="store_true",
        help="After writing the Parquet, also load into Postgres "
        "(also honors HIVDRWATCH_LOAD_TO_POSTGRES=1 env var).",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    if not args.skip_genbank and not args.email:
        parser.error(
            "--email is required unless --skip-genbank is set (also honors NCBI_EMAIL env var)"
        )

    ingest_flow(
        email=args.email,
        countries=args.countries,
        lanl_dir=args.lanl_dir,
        genbank_dir=args.genbank_dir,
        out_path=args.out,
        since=args.since,
        until=args.until,
        min_length=args.min_length,
        skip_genbank=args.skip_genbank,
        to_postgres=args.to_postgres or None,  # None → fall back to env var
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
