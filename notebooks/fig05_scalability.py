"""Paper 1 Figure 5 — pipeline scalability (wall time vs sequence count).

Reads `benchmarks/results/scale_*.json` (produced by `make benchmarks`)
and plots wall time per stage on a log-log scale. Also annotates each
stage with the fitted power-law exponent so a reviewer can see whether
the stage is roughly linear (b≈1) or super-linear (tn93 is ~O(N²)).

Reads:
    benchmarks/results/scale_ingest.json
    benchmarks/results/scale_align.json
    benchmarks/results/scale_cluster.json

Writes:
    paper1_tools/figures/fig05_scalability.pdf
    paper1_tools/figures/fig05_scalability.png

Launch:
    marimo edit notebooks/fig05_scalability.py
    make fig5
"""

from __future__ import annotations

import marimo

__generated_with = "0.23.16"
app = marimo.App(width="full")


@app.cell
def _():
    from pathlib import Path

    import matplotlib.pyplot as plt
    import polars as pl

    from hivdrwatch.paper1 import COLORBLIND_PALETTE, apply_paper_style
    from hivdrwatch.paper1.fig05_scalability import (
        STAGE_ORDER,
        fit_power_law,
        load_scalability_frame,
    )

    apply_paper_style()

    REPO_ROOT = Path(__file__).resolve().parent.parent
    RESULTS_DIR = REPO_ROOT / "benchmarks" / "results"
    OUT_DIR = REPO_ROOT / "paper1_tools" / "figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return (
        COLORBLIND_PALETTE,
        OUT_DIR,
        RESULTS_DIR,
        STAGE_ORDER,
        fit_power_law,
        load_scalability_frame,
        pl,
        plt,
    )


@app.cell
def _(RESULTS_DIR, load_scalability_frame):
    frame, summary = load_scalability_frame(RESULTS_DIR)
    frame, summary
    return frame, summary


@app.cell
def _(COLORBLIND_PALETTE, OUT_DIR, STAGE_ORDER, fit_power_law, frame, pl, plt):
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)

    stages_present = [s for s in STAGE_ORDER if s in frame["stage"].unique().to_list()]

    for i, stage in enumerate(stages_present):
        sub = frame.filter(pl.col("stage") == stage).sort("n_sequences")
        color = COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)]
        ax.plot(
            sub["n_sequences"].to_list(),
            sub["wall_seconds"].to_list(),
            marker="o",
            label=stage,
            color=color,
            linewidth=1.6,
        )
        fit = fit_power_law(frame, stage)
        if fit is not None:
            _a, b = fit
            ax.annotate(
                f"{stage} — b≈{b:.2f}",
                xy=(sub["n_sequences"][-1], sub["wall_seconds"][-1]),
                xytext=(6, 0),
                textcoords="offset points",
                fontsize=8,
                color=color,
                va="center",
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("input sequences (N)")
    ax.set_ylabel("wall time (seconds)")
    ax.set_title("HIV-DR-Watch pipeline scalability")
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    ax.legend(loc="upper left")

    fig_path_pdf = OUT_DIR / "fig05_scalability.pdf"
    fig_path_png = OUT_DIR / "fig05_scalability.png"
    fig.savefig(fig_path_pdf)
    fig.savefig(fig_path_png)
    fig_path_pdf, fig_path_png
    return fig, fig_path_pdf, fig_path_png


@app.cell
def _(fig_path_pdf, fig_path_png, summary):
    print("scalability:")
    print(f"  {summary.n_stages} stages, {summary.n_points} points, N up to {summary.max_n:,}")
    print(f"  wall range: {summary.min_wall_seconds:.3f}s .. {summary.max_wall_seconds:.3f}s")
    print()
    print("figures:")
    print(f"  PDF → {fig_path_pdf}")
    print(f"  PNG → {fig_path_png}")
    return


if __name__ == "__main__":
    app.run()
