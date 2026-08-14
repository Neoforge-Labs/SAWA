"""End-to-end ingestion CLI.

Usage (GenBank + LANL join):
    python -m hivdrwatch.ingest \\
        --email you@neoforge.io \\
        --countries KE UG TZ \\
        --since 2020-01-01 \\
        --lanl-dir data/raw/lanl \\
        --genbank-dir data/raw/genbank \\
        --out data/interim/east_africa_pol.parquet

Usage (LANL only — no NCBI network calls):
    python -m hivdrwatch.ingest \\
        --skip-genbank \\
        --lanl-dir data/raw/lanl \\
        --out data/interim/east_africa_pol.parquet

Expected LANL layout under --lanl-dir:
    <lanl_dir>/<cc>/results.tsv       (or results.txt — LANL's default)
    <lanl_dir>/<cc>/hiv-db.fasta      (optional alignment FASTA)

Environment:
    NCBI_API_KEY   optional; raises rate limit from 3 → 10 req/s.
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import date, datetime
from pathlib import Path

import polars as pl

from .genbank import EAST_AFRICA_COUNTRIES, GenBankIngester
from .lanl import LANLIngester
from .merge import merge_sources


def _configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def _parse_date(s: str) -> date:
    return datetime.strptime(s, "%Y-%m-%d").date()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="hivdrwatch.ingest")
    parser.add_argument(
        "--email",
        default=None,
        help="contact email for NCBI E-utilities (required unless --skip-genbank)",
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
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    _configure_logging(args.verbose)
    log = logging.getLogger("hivdrwatch.ingest")

    if not args.skip_genbank and not args.email:
        parser.error("--email is required unless --skip-genbank is set")

    ingester = LANLIngester(args.lanl_dir)

    # --- GenBank ---
    if args.skip_genbank:
        log.warning("--skip-genbank set; building canonical output from LANL alone.")
        # SequenceRecord path attaches the FASTA fallback and drops rows
        # without any sequence.
        lanl_records = [
            r.model_dump(mode="json")
            for r in ingester.to_sequence_records(min_length=args.min_length)
        ]
        if not lanl_records:
            log.error("LANL yielded no usable SequenceRecords; nothing to write.")
            return 1
        # Merge branch detects an empty GenBank frame and promotes LANL directly.
        genbank_df = pl.DataFrame({"accession": []})
        lanl_df = pl.from_dicts(lanl_records)
    else:
        # --- LANL first (fast, no network) ---
        log.info("loading LANL curation from %s", args.lanl_dir)
        lanl_df = ingester.load_directory(min_length=args.min_length)
        gb = GenBankIngester(email=args.email, raw_dir=args.genbank_dir)
        records: list[dict] = []
        for country in args.countries:
            log.info("fetching GenBank for %s", country)
            for rec in gb.fetch_country(
                country_code=country,
                since=args.since,
                until=args.until,
                min_length=args.min_length,
            ):
                records.append(rec.model_dump(mode="json"))
        genbank_df = pl.from_dicts(records) if records else pl.DataFrame({"accession": []})
        log.info("GenBank total: %d records", genbank_df.height)

    # --- Merge & write ---
    merge_sources(genbank_df, lanl_df, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
