"""Paper 1 Figure 1 — HIV-DR-Watch pipeline architecture.

Publication-quality swimlane diagram:
  * horizontal phase bands as soft background
  * generously spaced rounded boxes with subtle drop shadow
  * curved Bezier arrows with tier labels
  * legend at the bottom, phase labels inside each band

Reads:  nothing (spec lives in `hivdrwatch.paper1.fig01_pipeline`)
Writes: paper1_tools/figures/fig01_pipeline.{pdf,png}

Launch:
    marimo edit notebooks/fig01_pipeline_architecture.py
    make fig1
"""

from __future__ import annotations

import marimo

__generated_with = "0.23.16"
app = marimo.App(width="full")


@app.cell
def _():
    from pathlib import Path

    import matplotlib.patches as mpatches
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt

    from hivdrwatch.paper1 import apply_paper_style
    from hivdrwatch.paper1.fig01_pipeline import (
        PHASE_COLORS,
        pipeline_diagram,
        stage_summary,
    )

    apply_paper_style()

    REPO_ROOT = Path(__file__).resolve().parent.parent
    OUT_DIR = REPO_ROOT / "paper1_tools" / "figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return (
        OUT_DIR,
        PHASE_COLORS,
        mpatches,
        pe,
        pipeline_diagram,
        plt,
        stage_summary,
    )


@app.cell
def _(pipeline_diagram, stage_summary):
    diagram = pipeline_diagram()
    summary = stage_summary()
    diagram, summary
    return diagram, summary


@app.cell
def _(OUT_DIR, PHASE_COLORS, diagram, mpatches, pe, plt):
    # --- Layout constants — tuned generously for a two-column journal figure ---
    BOX_WIDTH = 2.4
    BOX_HEIGHT = 1.05
    ROW_SPACING = 2.4      # vertical distance between row centers
    COL_SPACING = 3.2      # horizontal distance per col unit
    BAND_PAD_X = 0.6       # extra background band width beyond outermost boxes
    BAND_PAD_Y = 0.35      # band half-height beyond box half-height

    n_rows = max(b.row for b in diagram.blocks) + 1
    all_cols = [b.col for b in diagram.blocks]
    min_col, max_col = min(all_cols), max(all_cols)

    fig_w = (max_col - min_col + 2) * COL_SPACING * 0.85
    fig_h = n_rows * ROW_SPACING * 0.95 + 1.1

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    ax.set_axis_off()
    ax.set_aspect("equal")

    # --- Position lookup for every block ---
    positions: dict[str, tuple[float, float]] = {}
    for b in diagram.blocks:
        x = b.col * COL_SPACING
        y = -b.row * ROW_SPACING
        positions[b.key] = (x, y)

    # --- Soft phase bands as horizontal backgrounds ---
    x_left = min_col * COL_SPACING - COL_SPACING * BAND_PAD_X - BOX_WIDTH / 2
    x_right = max_col * COL_SPACING + COL_SPACING * BAND_PAD_X + BOX_WIDTH / 2
    band_width = x_right - x_left

    for row, label in diagram.row_labels:
        # Pick this row's phase from the first block on the row (rows
        # are single-phase in this spec).
        row_phase = next(b.phase for b in diagram.blocks if b.row == row)
        color = PHASE_COLORS[row_phase]
        y_center = -row * ROW_SPACING
        band = mpatches.FancyBboxPatch(
            (x_left, y_center - BOX_HEIGHT / 2 - BAND_PAD_Y),
            band_width,
            BOX_HEIGHT + 2 * BAND_PAD_Y,
            boxstyle="round,pad=0.0,rounding_size=0.25",
            facecolor=color,
            edgecolor="none",
            alpha=0.10,
            zorder=0,
        )
        ax.add_patch(band)

        # Row label sits inside the band on the far left, uppercase and muted.
        ax.text(
            x_left + 0.25,
            y_center,
            label.upper(),
            fontsize=8,
            fontweight="bold",
            color="#555555",
            ha="left",
            va="center",
            family="sans-serif",
            zorder=1,
        )

    # --- Blocks ---
    for b in diagram.blocks:
        x, y = positions[b.key]
        color = PHASE_COLORS[b.phase]

        # Subtle drop shadow via path effects on a same-shape shadow patch.
        shadow = mpatches.FancyBboxPatch(
            (x - BOX_WIDTH / 2 + 0.06, y - BOX_HEIGHT / 2 - 0.08),
            BOX_WIDTH, BOX_HEIGHT,
            boxstyle="round,pad=0.04,rounding_size=0.15",
            facecolor="#000000",
            edgecolor="none",
            alpha=0.10,
            zorder=2,
        )
        ax.add_patch(shadow)

        rect = mpatches.FancyBboxPatch(
            (x - BOX_WIDTH / 2, y - BOX_HEIGHT / 2),
            BOX_WIDTH, BOX_HEIGHT,
            boxstyle="round,pad=0.04,rounding_size=0.15",
            facecolor="white",
            edgecolor=color,
            linewidth=1.6,
            zorder=3,
        )
        ax.add_patch(rect)

        # A thin colored bar across the top of each box — subtle way to
        # signal the phase without overwhelming the label.
        bar_h = 0.14
        top_bar = mpatches.FancyBboxPatch(
            (x - BOX_WIDTH / 2, y + BOX_HEIGHT / 2 - bar_h),
            BOX_WIDTH, bar_h,
            boxstyle="round,pad=0.0,rounding_size=0.15",
            facecolor=color,
            edgecolor="none",
            alpha=0.55,
            zorder=4,
        )
        ax.add_patch(top_bar)

        # Split label — first line bold (name), rest small (subtitle).
        parts = b.label.split("\n", 1)
        title = parts[0]
        subtitle = parts[1] if len(parts) > 1 else ""

        ax.text(
            x, y + 0.10,
            title,
            ha="center", va="center",
            fontsize=9.5,
            fontweight="bold",
            color="#1a1a1a",
            zorder=5,
        )
        if subtitle:
            ax.text(
                x, y - 0.22,
                subtitle,
                ha="center", va="center",
                fontsize=7.5,
                color="#555555",
                zorder=5,
            )

    # --- Edges (curved bezier) ---
    for e in diagram.edges:
        (x0, y0) = positions[e.src]
        (x1, y1) = positions[e.dst]

        # Always attach at box edges — top/bottom for vertical/diagonal
        # arrows, left/right for horizontal siblings.
        row_gap = abs(y1 - y0) > 0.1
        if row_gap:
            sy = y0 - BOX_HEIGHT / 2
            ey = y1 + BOX_HEIGHT / 2
            sx, ex = x0, x1
        else:
            if x1 > x0:
                sx, ex = x0 + BOX_WIDTH / 2, x1 - BOX_WIDTH / 2
            else:
                sx, ex = x0 - BOX_WIDTH / 2, x1 + BOX_WIDTH / 2
            sy, ey = y0, y1

        # Curvature: straight-down is straight; small arcs for adjacent
        # cols; more pronounced curve for anything spanning ≥2 cols so
        # the arrow arcs around any box that might sit between them.
        dx_units = abs(x1 - x0) / COL_SPACING
        if dx_units < 0.1:
            rad = 0.0
        elif dx_units < 1.1:
            rad = 0.18
        else:
            rad = 0.32
        # Alternate arc direction so parallel arrows don't overlap.
        if x1 < x0:
            rad = -rad
        connection = f"arc3,rad={rad}"

        arrow = mpatches.FancyArrowPatch(
            (sx, sy), (ex, ey),
            arrowstyle="-|>,head_length=0.55,head_width=0.42",
            mutation_scale=18,          # bigger overall arrow, incl. head
            linewidth=1.1,
            color="#3a3a3a",
            connectionstyle=connection,
            # Leave a visible gap between the arrow tip and the box so
            # the head isn't hidden by the destination box's border /
            # top-bar. shrinkA/B are in pixel units.
            shrinkA=2,
            shrinkB=7,
            # Above the box bodies so heads always render on top; still
            # below tier badges (zorder=6).
            zorder=5,
        )
        ax.add_patch(arrow)

        if e.tier:
            # Place the tier badge near the SOURCE end of the arrow
            # (1/3 of the way along), not the midpoint — midpoints of
            # curved arrows drift into empty space between rows and
            # look disconnected from either endpoint.
            t = 0.35
            mx = sx + t * (ex - sx)
            my = sy + t * (ey - sy)
            ax.text(
                mx, my,
                e.tier,
                fontsize=6.5,
                color="#333333",
                ha="center", va="center",
                bbox={
                    "boxstyle": "round,pad=0.28",
                    "facecolor": "white",
                    "edgecolor": "#bbbbbb",
                    "linewidth": 0.5,
                },
                zorder=6,
            )

    # --- Legend (phase colors) at the bottom ---
    handles = [
        mpatches.Patch(color=color, label=phase, alpha=0.85)
        for phase, color in PHASE_COLORS.items()
    ]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=len(PHASE_COLORS),
        frameon=False,
        fontsize=8,
        columnspacing=1.6,
        handletextpad=0.4,
    )

    # --- Bounds — pad so nothing clips ---
    ax.set_xlim(x_left - 0.4, x_right + 0.4)
    ax.set_ylim(-(n_rows - 1) * ROW_SPACING - 1.5, 1.4)

    fig.suptitle(
        "HIV-DR-Watch — pipeline data flow from public sources to figures and API",
        fontsize=11,
        fontweight="bold",
        color="#1a1a1a",
    )

    fig_path_pdf = OUT_DIR / "fig01_pipeline.pdf"
    fig_path_png = OUT_DIR / "fig01_pipeline.png"
    fig.savefig(fig_path_pdf)
    fig.savefig(fig_path_png)
    fig_path_pdf, fig_path_png
    return fig, fig_path_pdf, fig_path_png


@app.cell
def _(diagram, fig_path_pdf, fig_path_png, summary):
    print("pipeline:")
    print(f"  {len(diagram.blocks)} blocks · {len(diagram.edges)} edges · {len(diagram.row_labels)} rows")
    print("  per-phase counts:")
    for phase, count in summary.items():
        print(f"    {phase:8s} {count}")
    print()
    print("figures:")
    print(f"  PDF → {fig_path_pdf}")
    print(f"  PNG → {fig_path_png}")
    return


if __name__ == "__main__":
    app.run()
