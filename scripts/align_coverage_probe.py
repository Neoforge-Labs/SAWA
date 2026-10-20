"""Measure how much pol the curated alignment places, and under which preset.

Exploratory, and deliberately not part of the package. It answers
`docs/evaluations/ligase-hivdb.md` §9 item 1 — why `curated/aligned.parquet`
covers less of pol than Sierra interprets from the same raw sequence — and
its numbers are the ones quoted in `docs/evaluations/align-coverage.md`.

The raw sequences come back out of the local GenBank flat-file cache rather
than from a fresh NCBI pull, and every one is checked against the curated
row by `sequence_sha256`, so what the aligner was handed is known rather
than assumed.

    python scripts/align_coverage_probe.py scan-cache --out /tmp/raw_all.parquet
    python scripts/align_coverage_probe.py compare \
        --raw /tmp/raw_all.parquet --out /tmp/coverage.parquet
    python scripts/align_coverage_probe.py report --compared /tmp/coverage.parquet
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import mappy as mp
import polars as pl

from hivdrwatch.align.align_mappy import _build_keeplength, _usable_hits
from hivdrwatch.align.reference import POL_REGIONS
from hivdrwatch.resistance.cache import sequence_sha256

logger = logging.getLogger("align_coverage_probe")

POL_LEN = 3012
GENE_START = {gene: start for gene, (start, _) in POL_REGIONS.items()}


# ------------------------------------------------------------------
# step 1 — recover the raw sequences the aligner actually saw
# ------------------------------------------------------------------


def scan_cache(genbank_dir: Path, out: Path) -> None:
    """Pull accession → sequence out of every cached GenBank flat file.

    A hand-rolled scan rather than `Bio.SeqIO`: we need two fields from
    ~600 MB of flat files, and SeqIO builds a full record per entry.
    """
    accs: list[str] = []
    seqs: list[str] = []
    country: list[str | None] = []
    note: list[str | None] = []
    defs: list[str] = []

    for path in sorted(genbank_dir.glob("*.gb")):
        ver: str | None = None
        defline: list[str] = []
        in_def = False
        ctry: str | None = None
        nt: str | None = None
        buf: list[str] = []
        in_origin = False
        with path.open(errors="replace") as fh:
            for line in fh:
                if line.startswith("ORIGIN"):
                    in_origin, buf = True, []
                    continue
                if line.startswith("//"):
                    if ver is not None and buf:
                        accs.append(ver)
                        seqs.append("".join(buf))
                        country.append(ctry)
                        note.append(nt)
                        defs.append(" ".join(defline)[:200])
                    ver, defline, ctry, nt, buf = None, [], None, None, []
                    in_origin = in_def = False
                    continue
                if in_origin:
                    buf.append("".join(line.split()[1:]).upper())
                    continue
                if line.startswith("VERSION"):
                    ver, in_def = line.split()[1], False
                elif line.startswith("DEFINITION"):
                    defline, in_def = [line[12:].strip()], True
                elif in_def and line.startswith(" " * 12):
                    defline.append(line.strip())
                else:
                    in_def = False
                    stripped = line.strip()
                    if stripped.startswith(('/country="', '/geo_loc_name="')):
                        ctry = stripped.split('="', 1)[1].rstrip('"')
                    elif stripped.startswith('/note="subtype:'):
                        nt = stripped.split('="', 1)[1].rstrip('"')

    df = pl.DataFrame(
        {
            "accession": accs,
            "sequence": seqs,
            "length": [len(s) for s in seqs],
            "definition": defs,
            "gb_country": country,
            "gb_note": note,
        }
    ).unique(subset=["accession"], keep="first")
    df.write_parquet(out)
    logger.info("recovered %d records from %s → %s", df.height, genbank_dir, out)


# ------------------------------------------------------------------
# step 2 — align each sequence under each preset and score the span
# ------------------------------------------------------------------


def _span(aligned: str) -> tuple[int, int, int]:
    idx = [i for i, c in enumerate(aligned) if c != "-"]
    return (idx[0] + 1, idx[-1] + 1, len(idx)) if idx else (0, 0, 0)


def _sierra_positions(cache_dir: Path, sha: str) -> tuple[set[int], set[int]] | None:
    """Sierra mutation positions as pol-relative nt: (all, scoring-only)."""
    path = cache_dir / sha[:2] / f"{sha}.json"
    if not path.exists():
        return None
    analysis = json.loads(path.read_text()).get("analysis") or {}
    every: set[int] = set()
    scoring: set[int] = set()
    for gene in analysis.get("drugResistance", []):
        name = gene.get("gene", {}).get("name")
        if name not in GENE_START:
            continue
        for group in gene.get("mutationsByTypes", []):
            is_scoring = group.get("mutationType") != "Other"
            for mut in group.get("mutations", []):
                digits = "".join(c for c in mut.get("text", "") if c.isdigit())
                if not digits:
                    continue
                pos = GENE_START[name] + 3 * (int(digits) - 1)
                every.add(pos)
                if is_scoring:
                    scoring.add(pos)
    return every, scoring


def _score_span(
    label: str, start: int, end: int, sierra: tuple[set[int], set[int]] | None
) -> dict[str, int]:
    """How much of what Sierra read falls outside the span [start, end]."""
    if sierra is None:
        return {}
    every, scoring = sierra
    outside = [p for p in every if not start <= p <= end]
    return {
        f"{label}_outside": len(outside),
        f"{label}_outside_scoring": sum(1 for p in scoring if not start <= p <= end),
        f"{label}_max_outside": max(
            (min(abs(p - start), abs(p - end)) for p in outside), default=0
        ),
    }


def _legacy_single_block(aligner: mp.Aligner, query: str) -> str:
    """The pre-fix projection: one hit, whichever has the most matches.

    Kept here, not imported, because `align_mappy` no longer does this.
    It is the control: running it at `asm10` must reproduce the shipped
    `aligned.parquet` span exactly, which is what licenses reading any
    other column in this table as a like-for-like comparison.
    """
    hits = list(aligner.map(query))
    if not hits:
        return "-" * POL_LEN
    return _build_keeplength(query, [max(hits, key=lambda h: h.mlen)], POL_LEN)


def compare(
    raw_parquet: Path, curated: Path, reference: Path, cache_dir: Path, out: Path
) -> None:
    raw = pl.read_parquet(raw_parquet)
    cur = pl.read_parquet(
        curated, columns=["accession", "sequence_sha256", "aligned_pol", "subtype"]
    ).unique(subset=["accession"], keep="first")

    # `asm10_legacy` is the shipped pipeline exactly; `asm20` is the fixed
    # path. `asm10` isolates the preset change from the block-merge change.
    aligners = {
        preset: mp.Aligner(str(reference), preset=preset) for preset in ("asm10", "asm20")
    }
    legacy = mp.Aligner(str(reference), preset="asm10")
    seqs = dict(zip(raw["accession"].to_list(), raw["sequence"].to_list(), strict=True))

    rows: list[dict] = []
    for row in cur.iter_rows(named=True):
        acc = row["accession"]
        raw_seq = seqs.get(acc)
        if raw_seq is None:
            continue
        query = raw_seq.upper().replace("U", "T").replace("-", "").strip()
        shipped_start, shipped_end, shipped_n = _span(row["aligned_pol"])
        out_row: dict = {
            "accession": acc,
            "subtype": row["subtype"],
            "raw_len": len(query),
            # The join that makes everything else meaningful.
            "sha_matches": sequence_sha256(raw_seq) == row["sequence_sha256"],
            "shipped_start": shipped_start,
            "shipped_end": shipped_end,
            "shipped_n": shipped_n,
        }
        sierra = _sierra_positions(cache_dir, row["sequence_sha256"])

        lstart, lend, ln = _span(_legacy_single_block(legacy, query))
        out_row |= {"legacy_start": lstart, "legacy_end": lend, "legacy_n": ln}
        out_row |= _score_span("legacy", lstart, lend, sierra)

        for preset, aligner in aligners.items():
            hits = _usable_hits(aligner, query)
            aligned = _build_keeplength(query, hits, POL_LEN) if hits else "-" * POL_LEN
            start, end, n = _span(aligned)
            out_row |= {
                f"{preset}_start": start,
                f"{preset}_end": end,
                f"{preset}_n": n,
                f"{preset}_hits": len(hits),
            }
            out_row |= _score_span(preset, start, end, sierra)
        out_row["has_sierra"] = sierra is not None
        rows.append(out_row)

    df = pl.DataFrame(rows)
    df.write_parquet(out)
    logger.info("compared %d sequences → %s", df.height, out)


# ------------------------------------------------------------------
# step 3 — the tables quoted in the write-up
# ------------------------------------------------------------------


def report(compared: Path) -> None:
    df = pl.read_parquet(compared)
    n = df.height
    print(f"sequences: {n}")
    print(f"sha256 matches the curated row: {df['sha_matches'].sum()} / {n}")
    exact = df.filter(
        (pl.col("shipped_start") == pl.col("legacy_start"))
        & (pl.col("shipped_end") == pl.col("legacy_end"))
    ).height
    print(f"legacy path reproduces the shipped span: {exact} / {n}\n")

    def pct(k: int) -> str:
        return f"{k} ({100 * k / n:.1f}%)"

    print(f"{'coverage':<30}{'shipped':>18}{'asm10+merge':>18}{'asm20+merge':>18}")
    for label, expr in [
        ("no alignment at all", lambda p: pl.col(f"{p}_n") == 0),
        ("coverage >= 800 nt", lambda p: pl.col(f"{p}_n") >= 800),
        ("coverage >= 2900 nt", lambda p: pl.col(f"{p}_n") >= 2900),
    ]:
        print(
            f"{label:<30}"
            f"{pct(df.filter(expr('legacy')).height):>18}"
            f"{pct(df.filter(expr('asm10')).height):>18}"
            f"{pct(df.filter(expr('asm20')).height):>18}"
        )
    for label, col in [("mean coverage (nt)", "mean"), ("median coverage (nt)", "median")]:
        vals = [getattr(df[f"{p}_n"], col)() for p in ("legacy", "asm10", "asm20")]
        print(f"{label:<30}" + "".join(f"{v:>18.0f}" for v in vals))
    for preset in ("asm10", "asm20"):
        k = df.filter(pl.col(f"{preset}_hits") > 1).height
        print(f"  {preset}: >1 alignment block: {pct(k)}")

    s = df.filter(pl.col("has_sierra"))
    if not s.height:
        return
    print(f"\nagainst the Sierra cache ({s.height} sequences)")
    cols = ("legacy", "asm10", "asm20")
    print(f"{'':<30}{'shipped':>18}{'asm10+merge':>18}{'asm20+merge':>18}")
    for label, col in [
        (">=1 mutation outside span", "outside"),
        ("  ...a scoring mutation", "outside_scoring"),
        ("outside by >150 nt", "max_outside"),
    ]:
        threshold = 150 if col == "max_outside" else 0
        counts = [s.filter(pl.col(f"{p}_{col}") > threshold).height for p in cols]
        print(
            f"{label:<30}"
            + "".join(f"{f'{k} ({100 * k / s.height:.1f}%)':>18}" for k in counts)
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    scan = sub.add_parser("scan-cache", help="recover raw sequences from data/raw/genbank")
    scan.add_argument("--genbank-dir", type=Path, default=Path("data/raw/genbank"))
    scan.add_argument("--out", type=Path, required=True)

    cmp_ = sub.add_parser("compare", help="align under each preset and score the span")
    cmp_.add_argument("--raw", type=Path, required=True)
    cmp_.add_argument(
        "--curated", type=Path, default=Path("data/interim/curated/aligned.parquet")
    )
    cmp_.add_argument("--reference", type=Path, default=Path("data/reference/hxb2_pol.fasta"))
    cmp_.add_argument("--cache-dir", type=Path, default=Path("data/raw/hivdb"))
    cmp_.add_argument("--out", type=Path, required=True)

    rep = sub.add_parser("report", help="print the tables used in the write-up")
    rep.add_argument("--compared", type=Path, required=True)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.cmd == "scan-cache":
        scan_cache(args.genbank_dir, args.out)
    elif args.cmd == "compare":
        compare(args.raw, args.curated, args.reference, args.cache_dir, args.out)
    else:
        report(args.compared)
    return 0


if __name__ == "__main__":
    sys.exit(main())
