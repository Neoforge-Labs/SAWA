"""Paper 1 Figure 4 — DTG-era INSTI resistance case study.

Preregistered analysis: `paper1_tools/preregistration_fig04.md`.
Do NOT change the aggregation semantics without amending the preregistration
first.

Renders:
  * Left: pre vs post INSTI prevalence per (country × subtype),
    with 95% Wilson score intervals.
  * Right: delta prevalence (post − pre) with 95% Newcombe hybrid CI;
    strata with a signal (CI excludes zero) are highlighted.

Reads:
    data/interim/dr_results.parquet
    data/interim/curated/aligned.parquet

Writes:
    paper1_tools/figures/fig04_dtg_case_study.pdf
    paper1_tools/figures/fig04_dtg_case_study.png

Launch:
    marimo edit notebooks/fig04_dtg_case_study.py
    make fig4
"""

from __future__ import annotations

import marimo

__generated_with = "0.23.16"
app = marimo.App(width="full")


@app.cell
def _():
    from pathlib import Path

    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from hivdrwatch.paper1 import COLORBLIND_PALETTE, SUBTYPE_ORDER, apply_paper_style
    from hivdrwatch.paper1.fig04_dtg import (
        POST_YEARS,
        PRE_YEARS,
        ROLLOUT_YEAR,
        aggregate_dtg_landscape,
    )

    apply_paper_style()

    REPO_ROOT = Path(__file__).resolve().parent.parent
    DR_PATH = REPO_ROOT / "data" / "interim" / "dr_results.parquet"
    CURATED_PATH = REPO_ROOT / "data" / "interim" / "curated" / "aligned.parquet"
    OUT_DIR = REPO_ROOT / "paper1_tools" / "figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return (
        COLORBLIND_PALETTE,
        CURATED_PATH,
        DR_PATH,
        OUT_DIR,
        POST_YEARS,
        PRE_YEARS,
        ROLLOUT_YEAR,
        SUBTYPE_ORDER,
        aggregate_dtg_landscape,
        np,
        pl,
        plt,
    )


@app.cell
def _(CURATED_PATH, DR_PATH, aggregate_dtg_landscape, pl):
    dr = pl.read_parquet(DR_PATH)
    curated = pl.read_parquet(CURATED_PATH)
    per_epoch, delta, summary = aggregate_dtg_landscape(dr, curated)
    per_epoch, delta, summary
    return curated, delta, dr, per_epoch, summary


@app.cell
def _(
    COLORBLIND_PALETTE,
    OUT_DIR,
    POST_YEARS,
    PRE_YEARS,
    SUBTYPE_ORDER,
    delta,
    np,
    per_epoch,
    pl,
    plt,
    summary,
):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), constrained_layout=True)
    ax_pair, ax_delta = axes

    # ----- Left: pre vs post pairs -----
    strata = delta.sort(["country", "subtype"]).to_dicts()
    y = np.arange(len(strata))
    labels = [f"{s['country']}·{s['subtype']}" for s in strata]

    for row in per_epoch.to_dicts():
        key = (row["country"], row["subtype"])
        matching = [i for i, s in enumerate(strata) if (s["country"], s["subtype"]) == key]
        if not matching:
            continue
        yi = matching[0]
        color = COLORBLIND_PALETTE[0 if row["epoch"] == "pre" else 1]
        ax_pair.errorbar(
            row["prev"],
            yi,
            xerr=[[row["prev"] - row["prev_ci_low"]], [row["prev_ci_high"] - row["prev"]]],
            fmt="o",
            color=color,
            markersize=5,
            capsize=3,
            label="pre"
            if row["epoch"] == "pre" and yi == 0
            else "post"
            if row["epoch"] == "post" and yi == 0
            else None,
        )

    ax_pair.set_yticks(y)
    ax_pair.set_yticklabels(labels)
    ax_pair.set_xlabel("INSTI major-mutation prevalence (95% CI)")
    ax_pair.set_title(
        f"Pre ({PRE_YEARS[0]}–{PRE_YEARS[1]}) vs post ({POST_YEARS[0]}–) DTG rollout"
    )
    ax_pair.grid(axis="x", linestyle="--", alpha=0.3)

    handles = [
        plt.Line2D(
            [0], [0], marker="o", color=COLORBLIND_PALETTE[0], linestyle="", label="pre"
        ),
        plt.Line2D(
            [0], [0], marker="o", color=COLORBLIND_PALETTE[1], linestyle="", label="post"
        ),
    ]
    ax_pair.legend(handles=handles, loc="lower right")

    # ----- Right: delta with signal highlight -----
    for i, s in enumerate(strata):
        color = COLORBLIND_PALETTE[2] if s["signal"] else "#888888"
        ax_delta.errorbar(
            s["delta_prev"],
            i,
            xerr=[
                [s["delta_prev"] - s["delta_ci_low"]],
                [s["delta_ci_high"] - s["delta_prev"]],
            ],
            fmt="o",
            color=color,
            markersize=5,
            capsize=3,
        )
    ax_delta.axvline(0, color="black", linewidth=0.7, linestyle="--")
    ax_delta.set_yticks(y)
    ax_delta.set_yticklabels([])
    ax_delta.set_xlabel("Δ prevalence (post − pre, 95% CI)")
    ax_delta.set_title(
        f"{summary.n_headline_after_fdr}/{summary.n_strata_primary} primary strata BH-significant ({summary.n_signal_secondary}/{summary.n_strata_secondary} secondary exploratory)"
    )
    ax_delta.grid(axis="x", linestyle="--", alpha=0.3)

    fig.suptitle(
        "Emerging INSTI resistance in East Africa following the 2019 WHO DTG rollout",
        fontsize=11,
    )

    fig_path_pdf = OUT_DIR / "fig04_dtg_case_study.pdf"
    fig_path_png = OUT_DIR / "fig04_dtg_case_study.png"
    fig.savefig(fig_path_pdf)
    fig.savefig(fig_path_png)
    fig_path_pdf, fig_path_png
    return fig, fig_path_pdf, fig_path_png


@app.cell
def _(fig_path_pdf, fig_path_png, summary):
    print("DTG case study:")
    print(f"  N pre  = {summary.n_pre_seqs:,}")
    print(f"  N post = {summary.n_post_seqs:,}")
    print(f"  primary strata: {summary.n_strata_primary} (signal: {summary.n_signal_primary}, BH-headline: {summary.n_headline_after_fdr})")
    print(f"  secondary (exploratory) strata: {summary.n_strata_secondary} (signal: {summary.n_signal_secondary})")
    print()
    print("figures:")
    print(f"  PDF → {fig_path_pdf}")
    print(f"  PNG → {fig_path_png}")
    return


if __name__ == "__main__":
    app.run()
