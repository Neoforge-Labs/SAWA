"""Synthetic-ground-truth recovery benchmark for Paper 2 (CausalPhylo).

Simulates a corpus with **known** causal mutations plus **known**
decoy mutations that co-vary with subtype but do not drive the
label. Runs the full identify → estimate pipeline and reports
whether the estimator recovers the true causal edges and rejects
the decoy edges.

Because the causal pipeline in `causal/flows/causal.py` imports
DoWhy and Python 3.11+ StrEnum, this script is a *specification*
plus a lightweight, DoWhy-free estimator sanity check on the same
synthetic corpus. The full DoWhy identify → estimate → refute run
happens via `causal/flows/causal.py` when Kelyn runs it locally
against the same generator seed.

The generator lives under `tests/fixtures/_synthetic_paper2/` per
CLAUDE.md so the fixture can never be mistaken for real data.

Outputs `benchmarks/results/paper2_synth_recovery.json` for the
manuscript build to cite.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.linear_model import LogisticRegression

from .common import HostInfo, today_iso

logger = logging.getLogger(__name__)

DEFAULT_SIZES: tuple[int, ...] = (500, 2000, 8000)
DEFAULT_SEED: int = 42

# True causal mutations — the effect estimator should recover these.
# Each is (mutation_key, drug_class, effect_on_positive_rate).
TRUE_CAUSAL: tuple[tuple[str, str, float], ...] = (
    ("RT_184_V", "NRTI", 0.35),
    ("RT_103_N", "NNRTI", 0.30),
    ("IN_148_H", "INSTI", 0.25),
)

# Decoy mutations — subtype-associated but non-causal. Estimator
# should not flag these as causal after adjustment.
DECOY: tuple[tuple[str, str, float], ...] = (
    ("RT_215_Y", "C", 0.15),
    ("PR_10_I", "A1", 0.25),
    ("IN_140_S", "C", 0.20),
)

SUBTYPE_MIX: tuple[tuple[str, float], ...] = (
    ("B", 0.40),
    ("A1", 0.30),
    ("C", 0.30),
)


@dataclass(frozen=True)
class SyntheticFrames:
    mutations: pl.DataFrame  # per-accession binary mutation columns
    labels: pl.DataFrame  # per-(accession, drug_class) label
    covariates: pl.DataFrame  # accession-level


def generate(n: int, seed: int) -> SyntheticFrames:
    """Generate a synthetic corpus with known causal + decoy structure."""
    rng = np.random.default_rng(seed)
    subtypes, subtype_probs = zip(*SUBTYPE_MIX, strict=True)
    subtype = rng.choice(subtypes, size=n, p=subtype_probs)
    countries = rng.choice(["KE", "UG", "TZ", "RW", "BI", "ET"], size=n)
    years = rng.integers(2015, 2024, size=n)
    treated = rng.random(n) < 0.55  # latent treatment history

    # Base mutation prevalence (all subtypes)
    base_mut_rate = 0.10
    mut_cols = {
        m: (rng.random(n) < base_mut_rate).astype(np.int8)
        for m, _drug_class, _eff in TRUE_CAUSAL
    }
    # Decoy mutations: elevated in the target subtype
    for m, target_subtype, elevated_rate in DECOY:
        base = rng.random(n) < 0.04
        boosted = rng.random(n) < elevated_rate
        mut_cols[m] = np.where(subtype == target_subtype, boosted, base).astype(np.int8)

    # Treated patients acquire causal mutations at a boosted rate
    for m, _drug_class, effect in TRUE_CAUSAL:
        boost = rng.random(n) < effect
        mut_cols[m] = np.maximum(mut_cols[m], (treated & boost).astype(np.int8))

    accessions = [f"A{i:07d}" for i in range(n)]
    mutations = pl.DataFrame(
        {"accession": accessions, **{m: v.tolist() for m, v in mut_cols.items()}}
    )

    # Labels: per drug class, positive iff any causal mutation for that class
    # is present. This is the ground-truth causal structure.
    label_rows: list[dict] = []
    for drug_class in {dc for _m, dc, _e in TRUE_CAUSAL}:
        cs = [m for m, dc2, _e in TRUE_CAUSAL if dc2 == drug_class]
        pos = np.zeros(n, dtype=np.int8)
        for m in cs:
            pos = np.maximum(pos, mut_cols[m])
        label_rows.extend(
            {
                "accession": acc,
                "drug_class": drug_class,
                "label_binary": int(p),
            }
            for acc, p in zip(accessions, pos, strict=True)
        )
    labels = pl.DataFrame(label_rows)

    covariates = pl.DataFrame(
        {
            "accession": accessions,
            "subtype": subtype.tolist(),
            "country": countries.tolist(),
            "collection_year": years.tolist(),
        }
    )
    return SyntheticFrames(mutations=mutations, labels=labels, covariates=covariates)


def _fit_adjusted_effect(
    df: pl.DataFrame, *, mut: str, adj_cols: tuple[str, ...]
) -> tuple[float, float, float]:
    """Return (log_odds, ci_low, ci_high) for `mut → label` with
    logistic-regression adjustment on `adj_cols`, plus a
    500-bootstrap 95 % percentile CI. Sanity check only — the
    real analysis uses DoWhy through `causal/flows/causal.py`.
    """
    x = df.select([mut, *adj_cols]).to_pandas()
    # One-hot categorical covariates
    x = _one_hot(x, [c for c in adj_cols if not np.issubdtype(x[c].dtype, np.number)])
    y = df["label_binary"].to_numpy()
    lr = LogisticRegression(C=1.0, solver="liblinear", max_iter=500)
    lr.fit(x.values, y)
    idx = list(x.columns).index(mut)
    point = float(lr.coef_[0][idx])

    rng = np.random.default_rng(0)
    boots = np.empty(500)
    n = len(y)
    for i in range(500):
        ix = rng.integers(0, n, size=n)
        lr_b = LogisticRegression(C=1.0, solver="liblinear", max_iter=500)
        try:
            lr_b.fit(x.values[ix], y[ix])
            boots[i] = lr_b.coef_[0][idx]
        except ValueError:
            boots[i] = np.nan
    boots = boots[~np.isnan(boots)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def _one_hot(df, cols):
    import pandas as pd

    if not cols:
        return df
    return pd.get_dummies(df, columns=cols, drop_first=True, dtype=np.int8)


def run(sizes: tuple[int, ...], seed: int) -> dict:
    rows: list[dict] = []
    for n in sizes:
        f = generate(n, seed=seed)
        joined = f.mutations.join(f.covariates, on="accession").join(f.labels, on="accession")
        by_class = {
            dc: joined.filter(pl.col("drug_class") == dc) for dc in {"NRTI", "NNRTI", "INSTI"}
        }
        for m, dc, _expected in TRUE_CAUSAL:
            df = by_class[dc]
            unadj, ul, uh = _fit_adjusted_effect(df, mut=m, adj_cols=())
            adj, al, ah = _fit_adjusted_effect(
                df, mut=m, adj_cols=("subtype", "country", "collection_year")
            )
            rows.append(
                {
                    "n": n,
                    "mutation": m,
                    "drug_class": dc,
                    "kind": "true_causal",
                    "expected_effect_direction": "positive",
                    "unadjusted_log_odds": unadj,
                    "unadjusted_ci": [ul, uh],
                    "adjusted_log_odds": adj,
                    "adjusted_ci": [al, ah],
                    "recovers": bool(al > 0),
                }
            )
        for m, subtype, _ in DECOY:
            # Decoy is scored against every drug class (should be null)
            for dc in {"NRTI", "NNRTI", "INSTI"}:
                df = by_class[dc]
                if m not in df.columns:
                    df = df.with_columns(pl.lit(0).alias(m))
                unadj, ul, uh = _fit_adjusted_effect(df, mut=m, adj_cols=())
                adj, al, ah = _fit_adjusted_effect(
                    df, mut=m, adj_cols=("subtype", "country", "collection_year")
                )
                rows.append(
                    {
                        "n": n,
                        "mutation": m,
                        "drug_class": dc,
                        "kind": "decoy_subtype_confounded",
                        "confounded_by_subtype": subtype,
                        "unadjusted_log_odds": unadj,
                        "unadjusted_ci": [ul, uh],
                        "adjusted_log_odds": adj,
                        "adjusted_ci": [al, ah],
                        # Decoy passes if adjusted CI includes 0
                        "adjusted_ci_includes_zero": bool(al <= 0 <= ah),
                    }
                )

    return {
        "stage": "paper2_synth_recovery",
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
            "Sanity-check estimator on synthetic ground truth for Paper 2 "
            "v1 Results. The full DoWhy identify → estimate → refute run "
            "lives in causal/flows/causal.py; this benchmark checks that "
            "backdoor adjustment on {subtype, country, collection_year} "
            "correctly deconfounds the subtype-associated decoys."
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
        default=Path(__file__).parent / "results" / "paper2_synth_recovery.json",
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
