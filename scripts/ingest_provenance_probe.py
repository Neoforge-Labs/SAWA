"""Measure what the three ingest provenance bugs did to the shipped corpus.

Exploratory, and deliberately not part of the package. Its numbers are the
ones quoted in `docs/evaluations/ingest-provenance.md`.

Subtype and country are re-derived from the local GenBank flat-file cache
rather than estimated, and the pre-fix implementations are reproduced
verbatim below so "the shipped value is the buggy one" is a measurement
and not an assumption.

    python scripts/ingest_provenance_probe.py scan-cache \\
        --genbank-dir data/raw/genbank --out /tmp/gb_meta.jsonl
    python scripts/ingest_provenance_probe.py report \\
        --cache /tmp/gb_meta.jsonl \\
        --curated data/interim/curated/aligned.parquet \\
        --fm-dir data/gold/fm --out /tmp/provenance.json
"""

from __future__ import annotations

import argparse
import collections
import json
import logging
import sys
from pathlib import Path
from typing import Any

import polars as pl

from hivdrwatch.ingest.genbank import (
    EAST_AFRICA_COUNTRIES,
    _normalize_country,
    _normalize_subtype,
)
from hivdrwatch.ingest.schema import Subtype

logger = logging.getLogger("ingest_provenance_probe")

# Source-feature qualifiers worth carrying out of the flat files.
_WANTED_QUALIFIERS = ("country", "geo_loc_name", "note", "collection_date")


# ------------------------------------------------------------------
# the pre-fix implementations, reproduced so the delta is measurable
# ------------------------------------------------------------------


def legacy_normalize_subtype(note: str) -> str:
    """`genbank._normalize_subtype` as it stood at 99c312b.

    `f"SUBTYPE: {token}" in upper` is a substring test, and the Subtype
    enum declares C before the CRF members, so every CRF matched C.
    """
    if not note:
        return Subtype.UNKNOWN.value
    upper = note.upper()
    for st in Subtype:
        if st in (Subtype.UNKNOWN, Subtype.RECOMBINANT, Subtype.A, Subtype.F):
            continue  # A and F did not exist at 99c312b
        token = st.value.upper()
        if f"SUBTYPE: {token}" in upper or f"SUBTYPE {token}" in upper:
            return st.value
    if "RECOMBINANT" in upper or "CRF" in upper:
        return Subtype.RECOMBINANT.value
    return Subtype.UNKNOWN.value


def legacy_normalize_country(value: str | None) -> str | None:
    """`genbank._normalize_country` as it stood at 99c312b: East Africa only."""
    if not value:
        return None
    name = value.split(":")[0].strip()
    reverse = {v.lower(): k for k, v in EAST_AFRICA_COUNTRIES.items()}
    return reverse.get(name.lower())


# ------------------------------------------------------------------
# step 1 — recover the source-feature metadata the parser actually saw
# ------------------------------------------------------------------


def _parse_flat_file(path: Path) -> list[dict[str, Any]]:
    """Source-feature qualifiers per record, first value of each.

    A hand-rolled scan rather than `Bio.SeqIO`: we need four qualifiers
    out of ~600 MB of flat files, and SeqIO builds a full record per
    entry. Mirrors `_extract_source_features`, which takes the source
    feature only and the first value of each qualifier.
    """
    out: list[dict[str, Any]] = []
    version: str | None = None
    locus_length: int | None = None
    in_features = False
    feature_type: str | None = None
    source: dict[str, str] = {}
    qualifier: str | None = None
    buffer: list[str] = []

    def flush() -> None:
        nonlocal qualifier, buffer
        if qualifier is None:
            return
        value = " ".join(buffer).strip().strip('"')
        if feature_type == "source" and qualifier not in source:
            source[qualifier] = value
        qualifier, buffer = None, []

    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("//"):
                flush()
                if version is not None:
                    out.append(
                        {
                            "accession": version,
                            "locus_length": locus_length,
                            "batch": path.name,
                            **{q: source.get(q) for q in _WANTED_QUALIFIERS},
                        }
                    )
                version, locus_length, in_features = None, None, False
                feature_type, source = None, {}
                continue
            if line.startswith("ORIGIN"):
                flush()
                in_features = False
                continue
            if line.startswith("LOCUS"):
                parts = line.split()
                locus_length = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else None
                continue
            if line.startswith("VERSION"):
                version = line.split()[1]
                continue
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue
            if line.startswith(" " * 21):
                body = line[21:].rstrip("\n")
                if body.startswith("/"):
                    flush()
                    name, _, rest = body[1:].partition("=")
                    qualifier, buffer = name, [rest]
                elif qualifier is not None:
                    buffer.append(body.strip())
            elif line.startswith("     ") and line[5] != " ":
                flush()
                feature_type = line.split()[0]
    flush()
    return out


def scan_cache(genbank_dir: Path, out: Path) -> None:
    files = sorted(genbank_dir.glob("*.gb"))
    if not files:
        raise SystemExit(f"no *.gb files under {genbank_dir}")
    total = 0
    with out.open("w") as fh:
        for ix, path in enumerate(files, start=1):
            for row in _parse_flat_file(path):
                fh.write(json.dumps(row) + "\n")
                total += 1
            if ix % 300 == 0:
                logger.info("%d/%d files, %d records", ix, len(files), total)
    logger.info("wrote %d cache records from %d files → %s", total, len(files), out)


# ------------------------------------------------------------------
# step 2 — compare the shipped corpus against the re-derived truth
# ------------------------------------------------------------------


def _load_truth(cache: Path) -> pl.DataFrame:
    """One row per accession: shipped-era and fixed subtype + country."""
    first: dict[str, dict[str, Any]] = {}
    batches: dict[str, set[str]] = collections.defaultdict(set)
    for line in cache.read_text().splitlines():
        row = json.loads(line)
        first.setdefault(row["accession"], row)
        batches[row["accession"]].add(row["batch"])

    records = []
    for accession, row in first.items():
        note = row.get("note") or ""
        stated = row.get("geo_loc_name") or row.get("country")
        records.append(
            {
                "accession": accession,
                "locus_length": row.get("locus_length"),
                "n_cache_batches": len(batches[accession]),
                "has_country_qualifier": row.get("country") is not None,
                "stated_country_raw": stated,
                "fixed_country": _normalize_country(stated),
                "legacy_country": legacy_normalize_country(row.get("country")),
                "note": note,
                "fixed_subtype": _normalize_subtype(note).value,
                "legacy_subtype": legacy_normalize_subtype(note),
            }
        )
    return pl.from_dicts(records)


def _counter(df: pl.DataFrame, *cols: str) -> list[dict[str, Any]]:
    """Descending counts, ties broken on the group columns so repeat runs
    of this script produce byte-identical JSON."""
    return (
        df.group_by(list(cols))
        .len()
        .sort(["len", *cols], descending=[True, *[False] * len(cols)], nulls_last=True)
        .to_dicts()
    )


def report(cache: Path, curated: Path, fm_dir: Path | None, out: Path) -> None:
    truth = _load_truth(cache)
    shipped = pl.read_parquet(curated)

    # Curated keys accessions with the version suffix stripped in some
    # vintages; join on both spellings to be safe.
    truth = truth.with_columns(accession_base=pl.col("accession").str.replace(r"\.\d+$", ""))
    shipped = shipped.with_columns(
        accession_base=pl.col("accession").str.replace(r"\.\d+$", "")
    )
    joined = shipped.join(
        truth.drop("accession"), on="accession_base", how="left", validate="m:1"
    )

    result: dict[str, Any] = {
        "inputs": {
            "curated": str(curated),
            "cache_accessions": truth.height,
            "curated_rows": shipped.height,
            "curated_accessions": shipped["accession"].n_unique(),
            "curated_sha256": shipped["sequence_sha256"].n_unique(),
            "joined_to_cache": int(joined["fixed_subtype"].is_not_null().sum()),
        }
    }

    # ---- bug 1: subtype ------------------------------------------------
    # A shipped row is provably a bug-1 casualty when its label is the one
    # the buggy parser produced from this record's own note and the fixed
    # parser disagrees. Where the shipped label differs from the legacy
    # output, LANL curation (or another vintage) set it, not bug 1.
    subtype = joined.filter(pl.col("fixed_subtype").is_not_null())
    legacy_set = subtype.filter(pl.col("subtype") == pl.col("legacy_subtype"))
    relabelled = legacy_set.filter(pl.col("legacy_subtype") != pl.col("fixed_subtype"))
    result["bug1_subtype"] = {
        "rows_with_a_note_derived_label": legacy_set.height,
        "rows_relabelled": relabelled.height,
        "accessions_relabelled": relabelled["accession"].n_unique(),
        "transitions": _counter(relabelled, "subtype", "fixed_subtype"),
        "shipped_distribution": _counter(joined, "subtype"),
        "fixed_distribution_note_only": _counter(subtype, "fixed_subtype"),
        "crf_rows_shipped_as_c": legacy_set.filter(
            (pl.col("subtype") == "C") & pl.col("fixed_subtype").str.starts_with("CRF")
        ).height,
    }

    # ---- bug 2: country ------------------------------------------------
    stated = joined.filter(pl.col("stated_country_raw").is_not_null())
    wrong = stated.filter(
        pl.col("fixed_country").is_not_null() & (pl.col("country") != pl.col("fixed_country"))
    )
    no_statement = joined.filter(
        pl.col("fixed_subtype").is_not_null() & pl.col("stated_country_raw").is_null()
    )
    result["bug2_country"] = {
        "accessions_with_legacy_country_qualifier": int(truth["has_country_qualifier"].sum()),
        "accessions_with_geo_loc_name": int(truth["stated_country_raw"].is_not_null().sum()),
        "rows_with_a_stated_country": stated.height,
        "rows_where_shipped_disagrees": wrong.height,
        "accessions_where_shipped_disagrees": wrong["accession"].n_unique(),
        "rows_stating_a_non_east_african_country": stated.filter(
            ~pl.col("fixed_country").is_in(list(EAST_AFRICA_COUNTRIES))
            & pl.col("fixed_country").is_not_null()
        ).height,
        "rows_with_no_stated_country": no_statement.height,
        "shipped_vs_stated": _counter(wrong, "country", "fixed_country"),
        "shipped_distribution": _counter(joined, "country"),
        "unmapped_stated_names": _counter(
            stated.filter(pl.col("fixed_country").is_null()), "stated_country_raw"
        ),
    }

    # ---- bug 3: duplicate rows ----------------------------------------
    per_accession = shipped.group_by("accession").agg(
        pl.len().alias("rows"), pl.col("country").n_unique().alias("countries")
    )
    result["bug3_duplicates"] = {
        "surplus_rows": shipped.height - shipped["accession"].n_unique(),
        "accessions_duplicated": int((per_accession["rows"] > 1).sum()),
        "rows_per_accession": _counter(per_accession, "rows"),
        "patient_groups": shipped["patient_group_id"].n_unique(),
        "rows_flagged_earliest": int(shipped["is_earliest_per_patient"].sum()),
        "duplicate_accessions_sharing_a_patient_group": int(
            shipped.filter(
                pl.col("accession").is_in(
                    per_accession.filter(pl.col("rows") > 1)["accession"]
                )
            )
            .group_by("accession")
            .agg(pl.col("patient_group_id").n_unique().alias("groups"))
            .filter(pl.col("groups") == 1)
            .height
        ),
    }

    # ---- figure 2: country x year x subtype cells ----------------------
    def cells(country: str, subtype: str) -> set[tuple]:
        return set(
            joined.filter(pl.col("fixed_subtype").is_not_null())
            .select(country, "collection_year", subtype)
            .unique()
            .iter_rows()
        )

    shipped_cells = cells("country", "subtype")
    fixed_cells = cells("fixed_country", "fixed_subtype")
    result["figure2_cells"] = {
        "shipped": len(shipped_cells),
        "fixed": len(fixed_cells),
        "unchanged": len(shipped_cells & fixed_cells),
        "only_shipped": len(shipped_cells - fixed_cells),
        "only_fixed": len(fixed_cells - shipped_cells),
    }

    # ---- gene_region --------------------------------------------------
    lengths = joined.filter(pl.col("locus_length").is_not_null())
    result["gene_region"] = {
        "rows_over_4000nt_labelled_pol": lengths.filter(pl.col("locus_length") > 4000).height,
        "accessions_over_4000nt": int((truth["locus_length"] > 4000).sum()),
        "accessions_over_8000nt": int((truth["locus_length"] > 8000).sum()),
    }

    # ---- Paper 3 splits ------------------------------------------------
    if fm_dir is not None and fm_dir.exists():
        splits: dict[str, Any] = {}
        for path in sorted(fm_dir.glob("*/subtype_holdout/*.parquet")):
            df = pl.read_parquet(path)
            if "subtype" not in df.columns:
                continue
            key = f"{path.parent.parent.name}/{path.stem}"
            splits[key] = {
                "rows": df.height,
                "by_subtype": _counter(df, "subtype"),
            }
        result["paper3_subtype_holdout"] = splits

    out.write_text(json.dumps(result, indent=2, default=str))
    logger.info("wrote %s", out)
    print(json.dumps(result, indent=2, default=str))


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    scan = sub.add_parser("scan-cache")
    scan.add_argument("--genbank-dir", type=Path, required=True)
    scan.add_argument("--out", type=Path, required=True)

    rep = sub.add_parser("report")
    rep.add_argument("--cache", type=Path, required=True)
    rep.add_argument("--curated", type=Path, required=True)
    rep.add_argument("--fm-dir", type=Path, default=None)
    rep.add_argument("--out", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "scan-cache":
        scan_cache(args.genbank_dir, args.out)
    else:
        report(args.cache, args.curated, args.fm_dir, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
