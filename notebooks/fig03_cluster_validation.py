"""Paper 1 Figure 3 — cluster-detection parity vs HIV-TRACE.

Reads `data/gold/cluster_validation.json` (produced by
`make cluster-validate WITH_HIVTRACE=1`) and renders:

  Left panel:  grouped bar chart — precision/recall/F1/ARI for our
               clusterer vs HIV-TRACE against simulated ground truth.
  Right panel: pair-agreement fraction between the two clusterers,
               as a single annotated bar.

If HIV-TRACE wasn't run (e.g. the Docker image isn't built yet), the
figure gracefully renders only our clusterer and prints a warning.

Launch:
    marimo edit notebooks/fig03_cluster_validation.py
    make fig3
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

    from hivdrwatch.paper1 import COLORBLIND_PALETTE, apply_paper_style
    from hivdrwatch.paper1.fig03_parity import (
        CLUSTERER_ORDER,
        FIG3_METRICS,
        load_parity_payload,
        metrics_frame,
        summarize,
    )

    apply_paper_style()

    REPO_ROOT = Path(__file__).resolve().parent.parent
    RESULTS_JSON = REPO_ROOT / "data" / "gold" / "cluster_validation.json"
    OUT_DIR = REPO_ROOT / "paper1_tools" / "figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return (
        CLUSTERER_ORDER,
        COLORBLIND_PALETTE,
        FIG3_METRICS,
        OUT_DIR,
        RESULTS_JSON,
        load_parity_payload,
        metrics_frame,
        np,
        pl,
        plt,
        summarize,
    )


@app.cell
def _(RESULTS_JSON, load_parity_payload, metrics_frame, summarize):
    payload = load_parity_payload(RESULTS_JSON)
    m_frame = metrics_frame(payload)
    summary = summarize(payload)
    m_frame, summary
    return m_frame, payload, summary


@app.cell
def _(
    CLUSTERER_ORDER,
    COLORBLIND_PALETTE,
    FIG3_METRICS,
    OUT_DIR,
    m_frame,
    np,
    pl,
    plt,
    summary,
):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.0, 3.6),
        gridspec_kw={"width_ratios": [3, 1]},
        constrained_layout=True,
    )
    ax_metrics, ax_parity = axes

    clusterers_present = [
        c for c in CLUSTERER_ORDER if c in m_frame["clusterer"].unique().to_list()
    ]
    metrics_present = [m for m in FIG3_METRICS if m in m_frame["metric"].unique().to_list()]

    x = np.arange(len(metrics_present))
    bar_width = 0.8 / max(1, len(clusterers_present))

    for i, clust in enumerate(clusterers_present):
        vals = []
        for metric in metrics_present:
            row = m_frame.filter((pl.col("clusterer") == clust) & (pl.col("metric") == metric))
            vals.append(row["value"][0] if row.height else 0.0)
        offset = (i - (len(clusterers_present) - 1) / 2) * bar_width
        ax_metrics.bar(
            x + offset,
            vals,
            width=bar_width,
            color=COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)],
            label="ours" if clust.startswith("tn93") else "HIV-TRACE",
            edgecolor="black",
            linewidth=0.4,
        )
        for xi, v in zip(x + offset, vals, strict=True):
            ax_metrics.text(
                xi,
                v + 0.02,
                f"{v:.2f}",
                ha="center",
                va="bottom",
                fontsize=7,
            )

    ax_metrics.set_xticks(x)
    ax_metrics.set_xticklabels(metrics_present, rotation=0)
    ax_metrics.set_ylim(0, 1.1)
    ax_metrics.set_ylabel("score")
    ax_metrics.set_title(
        f"Cluster-detection parity  ·  N={summary.n_sequences} sequences, "
        f"{summary.n_truth_clusters} truth clusters"
    )
    ax_metrics.legend(loc="lower right")
    ax_metrics.grid(axis="y", linestyle="--", alpha=0.3)

    # Parity panel
    if summary.pair_agreement is not None:
        ax_parity.bar(
            [0],
            [summary.pair_agreement],
            width=0.6,
            color=COLORBLIND_PALETTE[2],
            edgecolor="black",
            linewidth=0.4,
        )
        ax_parity.text(
            0,
            summary.pair_agreement + 0.02,
            f"{summary.pair_agreement:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
        ax_parity.set_ylim(0, 1.1)
        ax_parity.set_xticks([0])
        ax_parity.set_xticklabels(["pair agreement"])
        ax_parity.set_title("ours vs HIV-TRACE")
        ax_parity.grid(axis="y", linestyle="--", alpha=0.3)
    else:
        ax_parity.text(
            0.5,
            0.5,
            "HIV-TRACE not run\n(run `make cluster-validate WITH_HIVTRACE=1`)",
            ha="center",
            va="center",
            fontsize=8,
            wrap=True,
        )
        ax_parity.set_axis_off()

    fig_path_pdf = OUT_DIR / "fig03_cluster_validation.pdf"
    fig_path_png = OUT_DIR / "fig03_cluster_validation.png"
    fig.savefig(fig_path_pdf)
    fig.savefig(fig_path_png)
    fig_path_pdf, fig_path_png
    return fig, fig_path_pdf, fig_path_png


@app.cell
def _(fig_path_pdf, fig_path_png, summary):
    print("parity:")
    print(f"  clusterers: {summary.n_clusterers}")
    print(
        f"  simulated: N={summary.n_sequences} across {summary.n_truth_clusters} truth clusters"
    )
    print(f"  ours F1:      {summary.ours_f1:.3f}")
    if summary.hivtrace_f1 is not None:
        print(f"  HIV-TRACE F1: {summary.hivtrace_f1:.3f}")
        print(f"  pair agreement: {summary.pair_agreement:.3f}")
    print()
    print("figures:")
    print(f"  PDF → {fig_path_pdf}")
    print(f"  PNG → {fig_path_png}")
    return


if __name__ == "__main__":
    app.run()
