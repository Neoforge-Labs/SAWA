"""Synthetic-ground-truth baseline benchmark for Paper 3.

Generates a controlled sequence-and-label corpus with:

  * A **subtype-shared true motif** — inserted at position 20 in
    every class-1 sequence regardless of subtype; this is the
    signal both baselines must recover.
  * A **subtype-specific decoy motif** — inserted at position 60
    only in subtype-C sequences (regardless of class); this is
    what tests whether a model conflates subtype composition with
    the label.

Runs both baselines (`KmerLogisticRegression`,
`OneHotRandomForest`) under both split strategies
(`stratified`, `subtype_holdout`) at three sample sizes
(`n ∈ {2000, 8000, 20000}`). Writes the per-drug × per-subtype
metrics frame to `benchmarks/results/paper3_baselines_synth.json`
so Paper 3's Results section cites a real number.

The generator itself is preserved under `tests/fixtures/` with the
`_synthetic_` prefix per CLAUDE.md so the fixture can never be
mistaken for real data.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from foundation_models.baselines import KmerLogisticRegression, OneHotRandomForest
from foundation_models.data import build_fm_dataset
from foundation_models.evaluate import per_subtype_metrics

from .common import HostInfo, today_iso

logger = logging.getLogger(__name__)

DEFAULT_SIZES: tuple[int, ...] = (2000, 8000, 20000)
DEFAULT_SEED: int = 42
SEQ_LEN: int = 200
ALPHABET: tuple[str, ...] = ("A", "C", "G", "T")

# Fixed motifs — never edit without bumping the benchmark's schema
# version and regenerating the CI baseline.
TRUE_MOTIF: str = "AAAAAAAA"  # subtype-shared → causal for label
TRUE_MOTIF_POS: int = 20
DECOY_MOTIF: str = "GGGGGGGG"  # subtype-C only → confounder
DECOY_MOTIF_POS: int = 60
DECOY_SUBTYPE: str = "C"

# Class-1 base rate before decoy injection
POSITIVE_RATE: float = 0.35


@dataclass(frozen=True)
class SyntheticCorpus:
    dr: pl.DataFrame
    curated: pl.DataFrame


def generate_corpus(n: int, seed: int) -> SyntheticCorpus:
    """Build the synthetic (dr_results, curated) Parquet-shaped frames.

    Balanced across three subtypes (B, A1, C) and six countries. One
    drug (3TC, NRTI class) — the benchmark is about the baseline
    behaviour, not about drug diversity.
    """
    rng = np.random.default_rng(seed)
    subtypes = np.array(["B", "A1", "C"])
    countries = np.array(["KE", "UG", "TZ", "RW", "BI", "ET"])

    acc = [f"S{i:07d}" for i in range(n)]
    subtype = rng.choice(subtypes, size=n, p=[0.34, 0.33, 0.33])
    country = rng.choice(countries, size=n)
    labels = rng.random(n) < POSITIVE_RATE

    seqs: list[str] = []
    for i in range(n):
        base = rng.choice(ALPHABET, size=SEQ_LEN)
        if labels[i]:
            base[TRUE_MOTIF_POS : TRUE_MOTIF_POS + len(TRUE_MOTIF)] = list(TRUE_MOTIF)
        if subtype[i] == DECOY_SUBTYPE:
            base[DECOY_MOTIF_POS : DECOY_MOTIF_POS + len(DECOY_MOTIF)] = list(DECOY_MOTIF)
        seqs.append("".join(base.tolist()))

    curated = pl.DataFrame(
        {
            "accession": acc,
            "sequence": seqs,
            "subtype": subtype.tolist(),
            "country": country.tolist(),
            "collection_year": [2020] * n,
        }
    )
    dr = pl.DataFrame(
        {
            "accession": acc,
            "drug": ["3TC"] * n,
            "drug_class": ["NRTI"] * n,
            "level": ["high" if lab else "susceptible" for lab in labels],
        }
    )
    return SyntheticCorpus(dr=dr, curated=curated)


def evaluate_pair(corpus: SyntheticCorpus, *, strategy: str, seed: int) -> tuple[dict, dict]:
    """Fit both baselines on the split, return per-subtype metrics."""
    ds = build_fm_dataset(corpus.dr, corpus.curated, strategy=strategy, seed=seed)
    train_seqs = ds.train["sequence"].to_list()
    train_y = ds.train["label_binary"].to_list()
    test_seqs = ds.test["sequence"].to_list()

    # Baseline 1 — k-mer LR
    lr = KmerLogisticRegression(k=3, seed=seed).fit(train_seqs, train_y)
    lr_prob = lr.predict_proba(test_seqs)

    # Baseline 2 — one-hot RF (needs equal-length; synth data satisfies this)
    rf = OneHotRandomForest(n_estimators=100, seed=seed).fit(train_seqs, train_y)
    rf_prob = rf.predict_proba(test_seqs)

    test = ds.test.with_columns(
        pl.Series("prob_lr", lr_prob),
        pl.Series("prob_rf", rf_prob),
    )

    lr_metrics = per_subtype_metrics(test, prob_col="prob_lr", min_positives=1).to_dicts()
    rf_metrics = per_subtype_metrics(test, prob_col="prob_rf", min_positives=1).to_dicts()

    lr_summary = {
        "auroc_median": float(
            np.median([m["auroc"] for m in lr_metrics if m["auroc"] is not None])
        )
        if lr_metrics
        else None,
        "per_subtype": lr_metrics,
    }
    rf_summary = {
        "auroc_median": float(
            np.median([m["auroc"] for m in rf_metrics if m["auroc"] is not None])
        )
        if rf_metrics
        else None,
        "per_subtype": rf_metrics,
    }
    return lr_summary, rf_summary


def run(sizes: tuple[int, ...], seed: int) -> dict:
    rows: list[dict] = []
    for n in sizes:
        corpus = generate_corpus(n, seed=seed)
        for strategy in ("stratified", "subtype_holdout"):
            lr, rf = evaluate_pair(corpus, strategy=strategy, seed=seed)
            rows.append(
                {
                    "n": n,
                    "strategy": strategy,
                    "kmer_lr_auroc_median": lr["auroc_median"],
                    "kmer_lr_per_subtype": lr["per_subtype"],
                    "one_hot_rf_auroc_median": rf["auroc_median"],
                    "one_hot_rf_per_subtype": rf["per_subtype"],
                }
            )
            logger.info(
                "n=%d strategy=%s → lr_med=%.3f rf_med=%.3f",
                n,
                strategy,
                lr["auroc_median"] or float("nan"),
                rf["auroc_median"] or float("nan"),
            )

    return {
        "stage": "paper3_baselines_synth",
        "seed": seed,
        "generated_at": today_iso(),
        "host": {
            "platform": HostInfo.collect().platform,
            "machine": HostInfo.collect().machine,
            "python": HostInfo.collect().python,
            "cpu_count": HostInfo.collect().cpu_count,
        },
        "schema_version": "0.1.0",
        "notes": (
            "Baseline benchmark for Paper 3 v1 Results section. "
            "True motif AAAAAAAA at pos 20 (subtype-shared, causal); "
            "decoy motif GGGGGGGG at pos 60 in subtype C (non-causal)."
        ),
        "measurements": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=list(DEFAULT_SIZES))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "results" / "paper3_baselines_synth.json",
    )
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level, format="%(levelname)s %(name)s: %(message)s")

    result = run(sizes=tuple(sorted(set(args.sizes))), seed=args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=False) + "\n")
    logger.info("wrote → %s", args.out)


if __name__ == "__main__":
    main()
