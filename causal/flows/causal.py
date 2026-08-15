"""Prefect flow: gold DR + outbreak parquets → causal effect estimates + HTML report.

Runs the committed domain DAG against both gold marts, producing:
    data/gold/causal_effects.parquet   — one row per (treatment, outcome, scope)
    data/gold/causal_report.html       — human-readable summary

Optionally logs to MLflow as an experiment (`causal-effect-estimation`)
with the parquet + HTML as artifacts. See tracking.py for the no-op path.

Run locally:
    python -m causal.flows.causal \\
        --dr-gold data/gold/dr_features.parquet \\
        --outbreak-gold data/gold/outbreak_features.parquet \\
        --out-dir data/gold
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import polars as pl
from prefect import flow, get_run_logger, task
from prefect.tasks import task_input_hash

from ..causal import (
    DAG_VERSION,
    build_dr_causal_frames,
    derive_prune_config,
    effects_to_frame,
    estimate_all_dr_effects,
    estimate_all_outbreak_effects,
    estimate_dr_effect,
    estimate_dr_effects_per_subtype,
    write_html_report,
)
from ..dag import DRDag
from ..tracking import tracked_run


@task(cache_key_fn=task_input_hash)
def load_dr_gold(path: Path) -> pl.DataFrame:
    log = get_run_logger()
    if not path.exists():
        log.warning("%s missing — DR causal analysis will be skipped", path)
        return pl.DataFrame()
    df = pl.read_parquet(path)
    log.info("loaded %d DR gold rows", df.height)
    return df


@task(cache_key_fn=task_input_hash)
def load_outbreak_gold(path: Path) -> pl.DataFrame:
    log = get_run_logger()
    if not path.exists():
        log.warning("%s missing — outbreak causal analysis will be skipped", path)
        return pl.DataFrame()
    df = pl.read_parquet(path)
    log.info("loaded %d outbreak gold rows", df.height)
    return df


@task
def run_dr_effects(
    df: pl.DataFrame,
    drugs: list[str] | None,
    bootstrap_samples: int,
    per_subtype: bool,
) -> list:
    log = get_run_logger()
    if df.is_empty():
        return []
    if drugs is None and "drug" in df.columns:
        drugs = df["drug"].unique().to_list()

    effects = []
    if drugs:
        for drug in drugs:
            effects.extend(
                estimate_all_dr_effects(
                    df,
                    drug=drug,
                    bootstrap_samples=bootstrap_samples,
                    per_subtype=per_subtype,
                )
            )
    else:
        effects.extend(
            estimate_all_dr_effects(
                df,
                bootstrap_samples=bootstrap_samples,
                per_subtype=per_subtype,
            )
        )
    log.info(
        "DR effects estimated: %d (bootstrap=%d per_subtype=%s)",
        len(effects),
        bootstrap_samples,
        per_subtype,
    )
    return effects


@task
def run_dr_effects_from_paper1(
    dr_parquet: Path,
    curated_parquet: Path,
    bootstrap_samples: int,
    per_subtype: bool,
    per_drug: bool,
    subtype_min_samples: int,
) -> list:
    """New path (Sprint 1 for Paper 2): build causal frames on the fly
    from Paper 1's dr_results + curated Parquet, then estimate per scope.

    This is the recommended entrypoint — the old `dr_gold` frame was
    produced by a module deleted in Phase 0. The adapter (data_prep.py)
    replaces it and stays in lockstep with the SAWA pipeline schema.
    """
    log = get_run_logger()
    if not dr_parquet.exists() or not curated_parquet.exists():
        log.warning(
            "missing dr_results (%s) or curated (%s) — DR causal analysis skipped",
            dr_parquet,
            curated_parquet,
        )
        return []

    dr = pl.read_parquet(dr_parquet)
    curated = pl.read_parquet(curated_parquet)
    frames, counts = build_dr_causal_frames(dr, curated, per_drug=per_drug)
    log.info(
        "adapter built %d %s frames (%d accessions × %d mutations)",
        counts.n_scopes,
        counts.scope_kind,
        counts.n_accessions,
        counts.n_mutation_features,
    )

    dag_mutations = DRDag().mutation_features
    all_effects: list = []
    for scope, frame in frames.items():
        scope_tag = f"DR/{counts.scope_kind}={scope}"
        if frame.height < subtype_min_samples:
            log.info("skipping %s (n=%d < %d)", scope_tag, frame.height, subtype_min_samples)
            continue
        for mutation in dag_mutations:
            all_effects.append(
                estimate_dr_effect(
                    frame,
                    mutation,
                    data_scope=scope_tag,
                    bootstrap_samples=bootstrap_samples,
                )
            )
            if per_subtype:
                per_st = estimate_dr_effects_per_subtype(
                    frame,
                    mutation,
                    min_samples=subtype_min_samples,
                    bootstrap_samples=bootstrap_samples,
                )
                # Prefix scope so cross-scope comparison stays disambiguated.
                for eff in per_st:
                    if not eff.data_scope.startswith(scope_tag):
                        eff = eff.model_copy(
                            update={
                                "data_scope": eff.data_scope.replace("DR/", f"{scope_tag}/"),
                            },
                        )
                    all_effects.append(eff)

    log.info(
        "Paper1→causal produced %d effects (%d %ss × %d mutations, "
        "bootstrap=%d, per_subtype=%s)",
        len(all_effects),
        counts.n_scopes,
        counts.scope_kind,
        len(dag_mutations),
        bootstrap_samples,
        per_subtype,
    )
    return all_effects


@task
def run_outbreak_effects(df: pl.DataFrame, bootstrap_samples: int) -> list:
    log = get_run_logger()
    if df.is_empty():
        return []
    effects = estimate_all_outbreak_effects(df, bootstrap_samples=bootstrap_samples)
    log.info("outbreak effects estimated: %d (bootstrap=%d)", len(effects), bootstrap_samples)
    return effects


@task
def write_outputs(
    dr_effects: list, outbreak_effects: list, out_dir: Path
) -> tuple[Path, Path]:
    log = get_run_logger()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_effects = list(dr_effects) + list(outbreak_effects)
    parquet_path = out_dir / "causal_effects.parquet"
    html_path = out_dir / "causal_report.html"

    dr_prune_path = out_dir / "prune_dr.json"
    outbreak_prune_path = out_dir / "prune_outbreak.json"
    if all_effects:
        frame = effects_to_frame(all_effects)
        frame.write_parquet(parquet_path, compression="zstd", compression_level=6)
        write_html_report(all_effects, html_path)
        # Emit per-domain prune configs — trainers can consume these.
        derive_prune_config(frame, dag_version=DAG_VERSION, scope_prefix="DR").write(
            dr_prune_path
        )
        derive_prune_config(frame, dag_version=DAG_VERSION, scope_prefix="outbreak").write(
            outbreak_prune_path
        )
    else:
        pl.DataFrame().write_parquet(parquet_path)
        html_path.write_text("<html><body><h1>no effects</h1></body></html>")

    log.info(
        "wrote %d causal effects → %s + %s",
        len(all_effects),
        parquet_path,
        html_path,
    )

    with tracked_run(
        experiment="causal-effect-estimation",
        run_name=f"dag-{DAG_VERSION}",
        tags={"flow": "causal", "dag_version": DAG_VERSION},
    ) as run:
        run.log_params(
            {
                "n_dr_effects": len(dr_effects),
                "n_outbreak_effects": len(outbreak_effects),
                "dag_version": DAG_VERSION,
            }
        )
        # High-level metrics: identified / passed-refutation counts.
        identified = sum(1 for e in all_effects if e.identified)
        passing = sum(
            1 for e in all_effects for r in e.refutations if r.verdict.value == "passes"
        )
        run.log_metrics(
            {
                "n_total": float(len(all_effects)),
                "n_identified": float(identified),
                "n_refutations_pass": float(passing),
            }
        )
        if parquet_path.exists():
            run.log_artifact(parquet_path)
        if html_path.exists():
            run.log_artifact(html_path)
        if dr_prune_path.exists():
            run.log_artifact(dr_prune_path)
        if outbreak_prune_path.exists():
            run.log_artifact(outbreak_prune_path)

    return parquet_path, html_path


@flow(name="causal")
def causal_flow(
    outbreak_gold: Path,
    out_dir: Path,
    *,
    dr_parquet: Path | None = None,
    curated_parquet: Path | None = None,
    dr_gold: Path | None = None,
    drugs: list[str] | None = None,
    bootstrap_samples: int = 0,
    per_subtype: bool = False,
    per_drug: bool = False,
    subtype_min_samples: int = 100,
) -> tuple[Path, Path]:
    """Causal effect estimation flow.

    Two DR input modes (mutually exclusive):
      * `dr_parquet` + `curated_parquet` (RECOMMENDED, Sprint 1):
        build frames on the fly via the data_prep adapter.
      * `dr_gold`: legacy pre-built DR gold parquet (kept for BC).
    """
    if dr_parquet and curated_parquet:
        dr_effects = run_dr_effects_from_paper1(
            dr_parquet,
            curated_parquet,
            bootstrap_samples,
            per_subtype,
            per_drug,
            subtype_min_samples,
        )
    elif dr_gold:
        dr_df = load_dr_gold(dr_gold)
        dr_effects = run_dr_effects(dr_df, drugs, bootstrap_samples, per_subtype)
    else:
        dr_effects = []

    outbreak_df = load_outbreak_gold(outbreak_gold)
    outbreak_effects = run_outbreak_effects(outbreak_df, bootstrap_samples)
    return write_outputs(dr_effects, outbreak_effects, out_dir)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="causal.flows.causal")
    parser.add_argument(
        "--dr-parquet",
        type=Path,
        default=Path("data/interim/dr_results.parquet"),
        help="Sierra DR results parquet (Paper 1 output). "
        "Used with --curated-parquet by the data_prep adapter.",
    )
    parser.add_argument(
        "--curated-parquet",
        type=Path,
        default=Path("data/interim/curated/aligned.parquet"),
        help="Silver curated parquet (Paper 1 output).",
    )
    parser.add_argument(
        "--dr-gold",
        type=Path,
        default=None,
        help="Legacy: pre-built DR gold parquet. Overrides --dr-parquet / "
        "--curated-parquet when supplied. Prefer the adapter path.",
    )
    parser.add_argument(
        "--outbreak-gold",
        type=Path,
        default=Path("data/gold/outbreak_features.parquet"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("data/gold"))
    parser.add_argument(
        "--drugs",
        nargs="+",
        default=None,
        help="Legacy path only — restrict DR analysis to these drug short names.",
    )
    parser.add_argument(
        "--bootstrap-samples",
        type=int,
        default=0,
        help="percentile-bootstrap resamples for CIs (0=off, ~200 recommended)",
    )
    parser.add_argument(
        "--per-subtype",
        action="store_true",
        help="also emit per-subtype (interaction proxy) DR effect rows",
    )
    parser.add_argument(
        "--per-drug",
        action="store_true",
        help="adapter path: emit one frame per drug (not per drug class). "
        "Adapter-only; ignored when --dr-gold is supplied.",
    )
    parser.add_argument(
        "--subtype-min-samples",
        type=int,
        default=100,
        help="per-subtype effect rows require ≥ this many rows in the subtype.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    causal_flow(
        outbreak_gold=args.outbreak_gold,
        out_dir=args.out_dir,
        dr_parquet=args.dr_parquet if args.dr_gold is None else None,
        curated_parquet=args.curated_parquet if args.dr_gold is None else None,
        dr_gold=args.dr_gold,
        drugs=args.drugs,
        bootstrap_samples=args.bootstrap_samples,
        per_subtype=args.per_subtype,
        per_drug=args.per_drug,
        subtype_min_samples=args.subtype_min_samples,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
