"""Head-to-head: ESM-2 fine-tune vs kmer-LR / one-hot-RF baselines.

Loads the same train/val/test Parquets that ``prepare_training_data``
emitted, runs the two non-FM baselines, loads the ESM-2 checkpoint
saved by ``esm.py``, and produces a comparison table:

    data/gold/fm/comparison_<class>_<strategy>.parquet    (per-drug)
    data/gold/fm/comparison_<class>_<strategy>.json       (headline)

The Parquet has one row per (drug, model) with AUROC, AUPRC, and
n_positive. The JSON is the headline table an analyst can drop
straight into a manuscript figure caption.

This is the file that turns "ESM-2 got AUROC 0.71" into a defensible
Paper 3 result — without it the ESM-2 number is unmoored.

Usage:

    python -m foundation_models.compare_baselines \\
        --manifest data/gold/fm/manifest_stratified.json \\
        --drug-class INSTI --strategy stratified

The comparison uses the same *test* split as the ESM-2 run, so the
numbers are directly comparable. kmer-LR and one-hot-RF are fit on
the same train split.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import polars as pl

from .baselines import KmerLogisticRegression, OneHotRandomForest
from .evaluate import compute_metrics

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ModelRow:
    """One (drug, model) result row in the comparison Parquet."""

    drug_class: str
    strategy: str
    drug: str
    model: str  # "esm2" | "kmer_lr" | "onehot_rf"
    auroc: float | None
    auprc: float | None
    n_test: int
    n_positive_test: int


def _load_splits(manifest_path: Path, drug_class: str, strategy: str):
    manifest = json.loads(manifest_path.read_text())
    runs = manifest.get("runs", [])
    match = next(
        (r for r in runs if r["drug_class"] == drug_class and r["strategy"] == strategy),
        None,
    )
    if match is None:
        raise ValueError(
            f"manifest {manifest_path} has no entry for "
            f"drug_class={drug_class!r} strategy={strategy!r}"
        )
    train = pl.read_parquet(match["train_parquet"])
    test = pl.read_parquet(match["test_parquet"])
    return train, test, match


def _score_kmer_lr(train: pl.DataFrame, test: pl.DataFrame) -> np.ndarray:
    """Fit kmer-LR on nucleotide sequences; return per-test-row P(class=1)."""
    lr = KmerLogisticRegression(k=3, seed=42).fit(
        train["sequence_nt"].to_list(),
        train["label_binary"].to_list(),
    )
    return lr.predict_proba(test["sequence_nt"].to_list())


def _score_onehot_rf(train: pl.DataFrame, test: pl.DataFrame) -> np.ndarray:
    """Fit one-hot-RF on protein sequences; return per-test-row P(class=1)."""
    rf = OneHotRandomForest(n_estimators=100, seed=42).fit(
        train["sequence_aa"].to_list(),
        train["label_binary"].to_list(),
    )
    return rf.predict_proba(test["sequence_aa"].to_list())


def _score_esm2(
    ckpt_path: Path,
    test: pl.DataFrame,
    model_name: str,
    batch_size: int = 16,
    max_length: int = 512,
) -> np.ndarray:
    """Load ESM-2 checkpoint + score the test set. Returns P(class=1)."""
    import torch
    from torch.utils.data import DataLoader
    from transformers import (
        AutoModelForSequenceClassification,
        AutoTokenizer,
        DataCollatorWithPadding,
    )

    if not ckpt_path.exists():
        raise FileNotFoundError(
            f"ESM-2 checkpoint missing: {ckpt_path}. Run `make fm-esm-train` first."
        )
    device = (
        "cuda"
        if torch.cuda.is_available()
        else "mps"
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
        else "cpu"
    )
    logger.info("loading ESM-2 checkpoint on %s: %s", device, ckpt_path)

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    try:
        model = AutoModelForSequenceClassification.from_pretrained(
            model_name, num_labels=2, attn_implementation="sdpa"
        )
    except (ValueError, TypeError):
        model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    model.load_state_dict(ckpt["model_state_dict"])
    model = model.to(device).eval()

    seqs = test["sequence_aa"].to_list()
    collator = DataCollatorWithPadding(tokenizer=tokenizer)

    class _Ds:
        def __init__(self, xs):
            self.xs = xs

        def __len__(self):
            return len(self.xs)

        def __getitem__(self, i):
            return tokenizer(self.xs[i], truncation=True, max_length=max_length, padding=False)

    loader = DataLoader(_Ds(seqs), batch_size=batch_size, collate_fn=collator)
    probs: list[float] = []
    with torch.no_grad():
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            logits = model(**batch).logits
            probs.extend(torch.softmax(logits, dim=-1)[:, 1].float().cpu().tolist())
    return np.asarray(probs)


def compare(
    manifest_path: Path,
    drug_class: str,
    strategy: str,
    out_dir: Path,
    *,
    esm_ckpt_dir: Path | None = None,
    esm_model_name: str = "facebook/esm2_t6_8M_UR50D",
    skip_esm: bool = False,
) -> list[ModelRow]:
    """Run all three models on the same splits; return per-drug rows."""
    train, test, _meta = _load_splits(manifest_path, drug_class, strategy)
    logger.info(
        "loaded splits: train=%d test=%d for %s/%s",
        train.height,
        test.height,
        drug_class,
        strategy,
    )

    # kmer-LR is fit once per drug for a fair per-drug comparison
    # (baselines don't have a shared multi-drug head, unlike ESM-2).
    all_rows: list[ModelRow] = []

    esm_probs: np.ndarray | None = None
    if not skip_esm:
        ckpt_dir = esm_ckpt_dir or Path("data/gold/fm/checkpoints")
        ckpt = ckpt_dir / drug_class.lower() / strategy / "best.pt"
        try:
            esm_probs = _score_esm2(ckpt, test, esm_model_name)
        except FileNotFoundError as e:
            logger.warning("skipping ESM-2 comparison: %s", e)

    for drug in sorted(test["drug"].unique().to_list()):
        drug_train = train.filter(pl.col("drug") == drug)
        drug_test = test.filter(pl.col("drug") == drug)
        drug_test_idx = (
            test.with_row_index("__idx").filter(pl.col("drug") == drug)["__idx"].to_list()
        )
        y_test = drug_test["label_binary"].to_numpy()

        if drug_train.height == 0 or drug_test.height == 0:
            logger.info("skipping %s: no train or test rows", drug)
            continue

        # kmer-LR
        try:
            lr_probs = _score_kmer_lr(drug_train, drug_test)
            m = compute_metrics(y_test, lr_probs)
            all_rows.append(
                ModelRow(
                    drug_class=drug_class,
                    strategy=strategy,
                    drug=drug,
                    model="kmer_lr",
                    auroc=m.auroc,
                    auprc=m.pr_auc,
                    n_test=drug_test.height,
                    n_positive_test=int(y_test.sum()),
                )
            )
        except (ValueError, RuntimeError) as e:
            logger.warning("kmer_lr failed on %s: %s", drug, e)

        # one-hot RF (needs equal-length AA; per-drug per-gene sequences
        # are uniform length by construction — the primary-gene path
        # yields IN=289 AA for every INSTI row)
        try:
            rf_probs = _score_onehot_rf(drug_train, drug_test)
            m = compute_metrics(y_test, rf_probs)
            all_rows.append(
                ModelRow(
                    drug_class=drug_class,
                    strategy=strategy,
                    drug=drug,
                    model="onehot_rf",
                    auroc=m.auroc,
                    auprc=m.pr_auc,
                    n_test=drug_test.height,
                    n_positive_test=int(y_test.sum()),
                )
            )
        except (ValueError, RuntimeError) as e:
            logger.warning("onehot_rf failed on %s: %s", drug, e)

        # ESM-2 (already scored the full test set; index into it)
        if esm_probs is not None:
            m = compute_metrics(y_test, esm_probs[drug_test_idx])
            all_rows.append(
                ModelRow(
                    drug_class=drug_class,
                    strategy=strategy,
                    drug=drug,
                    model="esm2",
                    auroc=m.auroc,
                    auprc=m.pr_auc,
                    n_test=drug_test.height,
                    n_positive_test=int(y_test.sum()),
                )
            )
    return all_rows


def write_outputs(rows: list[ModelRow], out_dir: Path, drug_class: str, strategy: str):
    """Emit the per-drug Parquet + a headline JSON summary."""
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = out_dir / f"comparison_{drug_class.lower()}_{strategy}.parquet"
    json_path = out_dir / f"comparison_{drug_class.lower()}_{strategy}.json"

    df = pl.from_dicts([asdict(r) for r in rows])
    df.write_parquet(parquet_path, compression="zstd", compression_level=6)

    # Headline: mean AUROC per model across drugs (only cells with ≥5 positives).
    scored = df.filter(pl.col("n_positive_test") >= 5)
    headline = (
        scored.group_by("model")
        .agg(
            pl.col("auroc").mean().alias("mean_auroc"),
            pl.col("auprc").mean().alias("mean_auprc"),
            pl.len().alias("n_drugs_scored"),
        )
        .sort("mean_auroc", descending=True)
        .to_dicts()
    )

    summary = {
        "drug_class": drug_class,
        "strategy": strategy,
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "headline": headline,
        "per_drug": df.to_dicts(),
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    logger.info("wrote %s and %s", parquet_path.name, json_path.name)
    logger.info("headline:\n%s", json.dumps(headline, indent=2))
    return parquet_path, json_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="foundation_models.compare_baselines")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=(
            "Path to the FM training manifest JSON. If omitted, auto-"
            "derived as data/gold/fm/manifest_<strategy>.json — the "
            "default emitted by prepare_training_data. Override only for "
            "custom out-dirs."
        ),
    )
    parser.add_argument("--drug-class", required=True)
    parser.add_argument(
        "--strategy", choices=("stratified", "subtype_holdout"), default="stratified"
    )
    parser.add_argument("--out-dir", type=Path, default=Path("data/gold/fm"))
    parser.add_argument("--esm-ckpt-dir", type=Path, default=Path("data/gold/fm/checkpoints"))
    parser.add_argument("--esm-model-name", default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument(
        "--skip-esm",
        action="store_true",
        help="Skip ESM-2 scoring — useful when only comparing baselines "
        "or when the checkpoint doesn't exist yet.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    manifest_path = args.manifest or (args.out_dir / f"manifest_{args.strategy}.json")
    if not manifest_path.exists():
        parser.error(
            f"manifest not found: {manifest_path}. Run "
            f"`make fm-prepare` (stratified) or `make fm-esm-holdout` "
            f"(subtype_holdout) first."
        )

    rows = compare(
        manifest_path,
        args.drug_class,
        args.strategy,
        args.out_dir,
        esm_ckpt_dir=args.esm_ckpt_dir,
        esm_model_name=args.esm_model_name,
        skip_esm=args.skip_esm,
    )
    if not rows:
        logger.error("no comparison rows produced")
        return 1
    write_outputs(rows, args.out_dir, args.drug_class, args.strategy)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
