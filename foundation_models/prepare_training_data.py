"""Prepare foundation-model training data from SAWA Parquets.

Reads the curated + dr_results Parquets, joins them, translates pol
nucleotide to protein (for ESM-2 and other protein LMs), splits into
train/val/test with either stratified or subtype-holdout strategy, and
writes one Parquet trio per drug class under
``data/gold/fm/<drug_class>/<strategy>/``.

The output frames carry the canonical columns every downstream trainer
expects:

    accession, sequence_nt, sequence_aa, subtype, country,
    collection_year, drug, drug_class, label_binary, split_tag

``split_tag`` is redundant with the file location (train / val / test
each written to a separate Parquet) but is preserved inside every file
so a downstream loader can concatenate + filter without needing to
know the directory layout.

Determinism:

* Splits are RNG-seeded (default seed=42).
* Sequences are gap-stripped when the source is aligned_pol.
* Translation uses the standard genetic code; premature stops are
  preserved in the AA sequence and reported in the summary.

Usage:

    python -m foundation_models.prepare_training_data \\
        --curated data/interim/curated/aligned.parquet \\
        --dr-results data/interim/dr_results.parquet \\
        --out-dir data/gold/fm \\
        --drug-classes INSTI NRTI NNRTI PI \\
        --strategy stratified \\
        --seed 42
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from .data import RESISTANT_LEVELS, build_fm_dataset

logger = logging.getLogger(__name__)

# Standard genetic code (nuc triplet → 1-letter AA).
_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_UNKNOWN_AA = "X"


def translate_pol_nt(nt: str) -> str:
    """Translate a nucleotide pol sequence to single-letter protein.

    Unknown / ambiguous codons (any character outside ACGT) map to
    ``X``. Premature stops (``*``) are preserved — downstream QC can
    decide whether to drop them; here we don't silently mask real data.
    Trailing 1-2 nt that don't complete a codon are dropped.
    """
    if not nt:
        return ""
    nt = nt.upper()
    n_full = len(nt) - (len(nt) % 3)
    aa_chars: list[str] = []
    for i in range(0, n_full, 3):
        codon = nt[i : i + 3]
        aa_chars.append(_CODON_TABLE.get(codon, _UNKNOWN_AA))
    return "".join(aa_chars)


@dataclass(frozen=True)
class ClassRunSummary:
    """Per-drug-class run record for the manifest."""

    drug_class: str
    strategy: str
    n_train: int
    n_val: int
    n_test: int
    n_unique_accessions: int
    n_subtypes: int
    positive_rate_train: float
    positive_rate_test: float
    sequence_column_used: str
    strip_gaps: bool
    aa_translation_applied: bool
    aa_stop_codon_rate: float
    n_dropped_high_ambiguity: int
    max_aa_x_fraction: float
    seed: int
    generated_at: str
    train_parquet: str
    val_parquet: str
    test_parquet: str


# Which drug classes need which pol subgene. Sierra scores PIs against
# PR only, NRTI/NNRTI against RT only, INSTI against IN only. Training
# on the wrong subgene is data leakage in reverse: the model has to
# learn to ignore irrelevant genes. We concatenate all three by default
# so a downstream trainer can subset if it wants, but the per-class
# **primary** subgene is recorded in the summary.
_PRIMARY_GENE: dict[str, str] = {
    "PI": "pr_nt",
    "NRTI": "rt_nt",
    "NNRTI": "rt_nt",
    "INSTI": "in_nt",
}


def _in_frame_pol_from_genes(row: dict) -> str:
    """Concatenate pr_nt + rt_nt + in_nt into a single in-frame pol NT.

    Missing sub-genes are silently omitted (their absence is already
    surfaced via curated.parquet's QC flags). The concatenation order
    matches the natural HIV-1 pol layout (PR → RT → IN), so
    frame-1 translation yields the canonical pol protein.
    """
    parts: list[str] = []
    for col in ("pr_nt", "rt_nt", "in_nt"):
        seq = row.get(col)
        if seq:
            parts.append(seq)
    return "".join(parts)


def prepare_drug_class(
    dr_results: pl.DataFrame,
    curated: pl.DataFrame,
    drug_class: str,
    out_dir: Path,
    *,
    strategy: str = "stratified",
    holdout_train_subtypes: tuple[str, ...] = ("A1", "C", "D"),
    seed: int = 42,
    translate_to_protein: bool = True,
    use_primary_gene_only: bool = True,
    max_aa_x_fraction: float = 1.0,
) -> ClassRunSummary:
    """Produce the train/val/test Parquets for one drug class.

    ``use_primary_gene_only`` (default True): emit only the pol subgene
    Sierra actually scores for this drug class:

      * PI    → PR  (99 AA)
      * NRTI  → RT  (~440 AA)
      * NNRTI → RT  (~440 AA)
      * INSTI → IN  (~289 AA)

    This is scientifically correct (matches how Sierra assigns levels —
    training on all three genes lets the model learn spurious
    cross-gene correlations) AND dramatically faster because
    transformer attention is O(n²) on sequence length. For INSTI the
    savings are 830→289 = 2.9× shorter → ~8× faster attention.

    Set to False to fall back to the concatenated PR+RT+IN sequence.
    """
    dc_dr = dr_results.filter(pl.col("drug_class") == drug_class)
    if dc_dr.is_empty():
        raise ValueError(f"no dr_results rows for drug_class={drug_class!r}")

    per_gene_cols = {"pr_nt", "rt_nt", "in_nt"}
    primary_gene = _PRIMARY_GENE.get(drug_class)

    if use_primary_gene_only and primary_gene and primary_gene in curated.columns:
        # Scientifically-aligned path: only the Sierra-scored subgene.
        curated_prepped = curated.with_columns(
            pl.col(primary_gene).alias("sequence")
        )
        sequence_source = f"{primary_gene} (primary gene only)"
    elif per_gene_cols.issubset(curated.columns):
        # Legacy path: concatenated pol. Kept for A/B comparison but no
        # longer default — see docstring rationale.
        primary_gene = "concatenated pol (PR+RT+IN)"
        curated_prepped = curated.with_columns(
            pl.struct(["pr_nt", "rt_nt", "in_nt"])
            .map_elements(
                lambda s: _in_frame_pol_from_genes(s),
                return_dtype=pl.Utf8,
            )
            .alias("sequence")
        )
        sequence_source = "pr_nt+rt_nt+in_nt (in-frame concat)"
    else:
        curated_prepped = curated
        primary_gene = "aligned_pol"
        sequence_source = "aligned_pol (gap-stripped)"

    ds = build_fm_dataset(
        dc_dr,
        curated_prepped,
        strategy=strategy,
        holdout_train_subtypes=holdout_train_subtypes,
        seed=seed,
    )

    logger.info(
        "%s: sequence source = %s, primary gene = %s",
        drug_class,
        sequence_source,
        primary_gene,
    )

    # Translate NT → AA in each split. Keep both columns so the trainer
    # can pick — protein LMs (ESM-2) use `sequence_aa`, nucleotide LMs
    # (DNABERT-2, Nucleotide Transformer) use `sequence_nt`.
    n_dropped_ambig = 0

    def _finalize(frame: pl.DataFrame, tag: str) -> pl.DataFrame:
        nonlocal n_dropped_ambig
        if translate_to_protein:
            aa = [translate_pol_nt(s) for s in frame["sequence"].to_list()]
            frame = frame.with_columns(
                pl.col("sequence").alias("sequence_nt"),
                pl.Series("sequence_aa", aa),
            ).drop("sequence")

            # Drop rows whose AA sequence is mostly ambiguous (X). These
            # are usually sequences with heavy NT ambiguity throughout —
            # they teach the model nothing useful.
            if max_aa_x_fraction < 1.0:
                before = frame.height
                # Vectorized: count 'X' chars per row / total length.
                frame = frame.with_columns(
                    (
                        pl.col("sequence_aa").str.count_matches("X")
                        / pl.col("sequence_aa").str.len_chars().cast(pl.Float64)
                    ).alias("_x_frac")
                )
                frame = frame.filter(pl.col("_x_frac") <= max_aa_x_fraction).drop(
                    "_x_frac"
                )
                dropped = before - frame.height
                if dropped > 0:
                    logger.info(
                        "%s: dropped %d/%d rows with AA X-fraction > %.2f",
                        tag,
                        dropped,
                        before,
                        max_aa_x_fraction,
                    )
                    n_dropped_ambig += dropped
        else:
            frame = frame.rename({"sequence": "sequence_nt"}).with_columns(
                pl.lit(None, dtype=pl.Utf8).alias("sequence_aa")
            )
        return frame.with_columns(pl.lit(tag).alias("split_tag"))

    train = _finalize(ds.train, "train")
    val = _finalize(ds.val, "val")
    test = _finalize(ds.test, "test")

    dc_dir = out_dir / drug_class.lower() / strategy
    dc_dir.mkdir(parents=True, exist_ok=True)
    train_path = dc_dir / "train.parquet"
    val_path = dc_dir / "val.parquet"
    test_path = dc_dir / "test.parquet"
    train.write_parquet(train_path, compression="zstd", compression_level=6)
    val.write_parquet(val_path, compression="zstd", compression_level=6)
    test.write_parquet(test_path, compression="zstd", compression_level=6)

    # Diagnostic: fraction of AA sequences with a premature stop.
    aa_stop_rate = 0.0
    if translate_to_protein and train.height > 0:
        stop_hits = sum(1 for s in train["sequence_aa"].to_list() if "*" in s)
        aa_stop_rate = stop_hits / train.height

    summary = ClassRunSummary(
        drug_class=drug_class,
        strategy=strategy,
        n_train=train.height,
        n_val=val.height,
        n_test=test.height,
        n_unique_accessions=ds.counts.n_unique_accessions,
        n_subtypes=ds.counts.n_subtypes,
        positive_rate_train=float(train["label_binary"].mean() or 0.0),
        positive_rate_test=float(test["label_binary"].mean() or 0.0),
        sequence_column_used=sequence_source,
        strip_gaps=True,
        aa_translation_applied=translate_to_protein,
        aa_stop_codon_rate=aa_stop_rate,
        n_dropped_high_ambiguity=n_dropped_ambig,
        max_aa_x_fraction=max_aa_x_fraction,
        seed=seed,
        generated_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        train_parquet=str(train_path),
        val_parquet=str(val_path),
        test_parquet=str(test_path),
    )
    logger.info(
        "%s/%s: train=%d val=%d test=%d, +rate train=%.2f%% test=%.2f%%, "
        "aa_stop_rate=%.2f%%",
        drug_class,
        strategy,
        summary.n_train,
        summary.n_val,
        summary.n_test,
        summary.positive_rate_train * 100,
        summary.positive_rate_test * 100,
        summary.aa_stop_codon_rate * 100,
    )
    return summary


def write_manifest(summaries: list[ClassRunSummary], out_path: Path) -> Path:
    """Emit a JSON manifest listing every produced Parquet + its shape."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps({"runs": [asdict(s) for s in summaries]}, indent=2),
        encoding="utf-8",
    )
    logger.info("FM training manifest → %s (%d runs)", out_path, len(summaries))
    return out_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="foundation_models.prepare_training_data")
    parser.add_argument("--curated", type=Path, required=True)
    parser.add_argument("--dr-results", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("data/gold/fm"))
    parser.add_argument(
        "--drug-classes",
        nargs="+",
        default=["INSTI", "NRTI", "NNRTI", "PI"],
        help="Which drug classes to prepare. INSTI first because it's smallest.",
    )
    parser.add_argument(
        "--strategy",
        choices=("stratified", "subtype_holdout"),
        default="stratified",
    )
    parser.add_argument(
        "--holdout-train-subtypes",
        nargs="+",
        default=["A1", "C", "D"],
        help="For subtype_holdout: subtypes to train on (test = all others).",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--no-translate",
        action="store_true",
        help="Skip nt→aa translation; keep only sequence_nt.",
    )
    parser.add_argument(
        "--concat-pol",
        action="store_true",
        help=(
            "Legacy: emit PR+RT+IN concatenated instead of the "
            "per-drug-class primary gene. Slower AND lets the model "
            "learn spurious cross-gene correlations. Kept only for A/B "
            "comparison against the default (primary-gene-only)."
        ),
    )
    parser.add_argument(
        "--max-aa-x-fraction",
        type=float,
        default=1.0,
        help=(
            "Drop AA sequences with more than this fraction of X "
            "(unknown / ambiguous residues). Default 1.0 = disabled. "
            "Opt in with e.g. 0.50 to drop majority-ambiguous sequences, "
            "or 0.20 for a stricter filter. NB: on the current East "
            "African cohort, 0.50 drops ~57%% of rows and 0.20 drops "
            "~60%% — usable if you also verify the positive rate holds."
        ),
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    dr = pl.read_parquet(args.dr_results)
    cur = pl.read_parquet(args.curated)
    if "passes_qc" in cur.columns:
        n_before = cur.height
        cur = cur.filter(pl.col("passes_qc"))
        logger.info(
            "curated QC filter: %d → %d rows (%d dropped)",
            n_before,
            cur.height,
            n_before - cur.height,
        )

    summaries: list[ClassRunSummary] = []
    for dc in args.drug_classes:
        try:
            s = prepare_drug_class(
                dr,
                cur,
                drug_class=dc,
                out_dir=args.out_dir,
                strategy=args.strategy,
                holdout_train_subtypes=tuple(args.holdout_train_subtypes),
                seed=args.seed,
                translate_to_protein=not args.no_translate,
                use_primary_gene_only=not args.concat_pol,
                max_aa_x_fraction=args.max_aa_x_fraction,
            )
            summaries.append(s)
        except ValueError as e:
            logger.warning("skipping %s: %s", dc, e)

    if summaries:
        write_manifest(summaries, args.out_dir / f"manifest_{args.strategy}.json")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
