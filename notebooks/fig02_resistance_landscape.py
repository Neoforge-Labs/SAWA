"""Paper 1 Figure 2 — HIV-1 drug-resistance landscape across East Africa.

Rows: (country × year). Columns: drug class × subtype panel.
Cell value: fraction of sequences with Sierra resistance ≥ intermediate.
Cell annotation: N (dropped when < min_cell_n).

Reads:
    data/interim/dr_results.parquet    (Sierra output)
    data/interim/curated/aligned.parquet (country + subtype + year)

Writes:
    paper1_tools/figures/fig02_resistance_landscape.pdf
    paper1_tools/figures/fig02_resistance_landscape.png

Launch:
    marimo edit notebooks/fig02_resistance_landscape.py       # interactive
    python -m marimo run notebooks/fig02_resistance_landscape.py  # headless render
    make fig2                                                    # via Makefile
"""

from __future__ import annotations

import marimo

__generated_with = "0.23.16"
app = marimo.App(width="full")


@app.cell
def _():
    from pathlib import Path
    import polars as pl
    import matplotlib.pyplot as plt
    import numpy as np

    from hivdrwatch.paper1 import (
        DRUG_CLASS_ORDER,
        SUBTYPE_ORDER,
        apply_paper_style,
        resistance_cmap,
    )
    from hivdrwatch.paper1.fig02_aggregate import (
        MIN_CELL_N,
        aggregate_resistance_landscape,
    )

    apply_paper_style()

    REPO_ROOT = Path(__file__).resolve().parent.parent
    DR_PATH = REPO_ROOT / "data" / "interim" / "dr_results.parquet"
    CURATED_PATH = REPO_ROOT / "data" / "interim" / "curated" / "aligned.parquet"
    OUT_DIR = REPO_ROOT / "paper1_tools" / "figures"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return (
        CURATED_PATH,
        DRUG_CLASS_ORDER,
        DR_PATH,
        MIN_CELL_N,
        OUT_DIR,
        Path,
        SUBTYPE_ORDER,
        aggregate_resistance_landscape,
        np,
        pl,
        plt,
        resistance_cmap,
    )


@app.cell
def _(CURATED_PATH, DR_PATH, aggregate_resistance_landscape, pl):
    dr = pl.read_parquet(DR_PATH)
    curated = pl.read_parquet(CURATED_PATH)
    landscape, counts = aggregate_resistance_landscape(dr, curated)
    landscape, counts
    return counts, curated, dr, landscape


@app.cell
def _(
    DRUG_CLASS_ORDER,
    OUT_DIR,
    SUBTYPE_ORDER,
    landscape,
    np,
    pl,
    plt,
    resistance_cmap,
):
    # Country × year rows we'll plot — sort by (country, year) and keep
    # only combinations that actually appear in the data.
    country_years = (
        landscape.select(["country", "collection_year"])
        .unique()
        .sort(["country", "collection_year"])
    )
    countries = country_years["country"].unique(maintain_order=True).to_list()
    years = country_years["collection_year"].unique(maintain_order=True).to_list()

    row_labels = [f"{c} {y}" for c in countries for y in years]
    row_index = {(c, y): i for i, (c, y) in enumerate((c, y) for c in countries for y in years)}
    n_rows = len(row_labels)

    # Columns: (drug_class × subtype) panels. Restrict to subtypes we
    # care about + drop empty combos to keep the figure compact.
    present_subtypes = [
        s for s in SUBTYPE_ORDER
        if landscape.filter(pl.col("subtype") == s).height > 0
    ]

    fig, axes = plt.subplots(
        1, len(DRUG_CLASS_ORDER),
        figsize=(3.2 * len(DRUG_CLASS_ORDER), 0.35 * n_rows + 1.5),
        sharey=True,
        constrained_layout=True,
    )
    if len(DRUG_CLASS_ORDER) == 1:
        axes = [axes]

    cmap = resistance_cmap()
    vmin, vmax = 0.0, 1.0

    for ax, drug_class in zip(axes, DRUG_CLASS_ORDER, strict=True):
        # Grid: rows = country_years, cols = present_subtypes
        grid = np.full((n_rows, len(present_subtypes)), np.nan)
        annot = [[""] * len(present_subtypes) for _ in range(n_rows)]

        class_df = landscape.filter(pl.col("drug_class") == drug_class)
        for row in class_df.iter_rows(named=True):
            key = (row["country"], row["collection_year"])
            if key not in row_index:
                continue
            r = row_index[key]
            s = row["subtype"]
            if s not in present_subtypes:
                continue
            c = present_subtypes.index(s)
            if row["frac_resistant"] is not None:
                grid[r, c] = row["frac_resistant"]
            annot[r][c] = str(row["n_sequences"])

        im = ax.imshow(
            grid,
            aspect="auto",
            cmap=cmap,
            vmin=vmin, vmax=vmax,
            interpolation="nearest",
        )
        # Grey out NaN (below-min-N cells)
        ax.set_facecolor("#e6e6e6")

        # Cell annotations — N (small text)
        for r in range(n_rows):
            for c in range(len(present_subtypes)):
                if annot[r][c]:
                    val = grid[r, c]
                    txt_color = "white" if (not np.isnan(val) and val > 0.5) else "black"
                    ax.text(
                        c, r, annot[r][c],
                        ha="center", va="center",
                        fontsize=6, color=txt_color,
                    )

        ax.set_xticks(range(len(present_subtypes)))
        ax.set_xticklabels(present_subtypes, rotation=0)
        ax.set_title(drug_class)
        ax.set_xlabel("subtype")

    axes[0].set_yticks(range(n_rows))
    axes[0].set_yticklabels(row_labels)
    axes[0].set_ylabel("country · year")

    cbar = fig.colorbar(
        im, ax=axes, shrink=0.8, aspect=40, pad=0.02,
        label="fraction with Sierra ≥ intermediate resistance",
    )
    cbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])

    fig.suptitle(
        "HIV-1 pol drug-resistance landscape, East Africa"
        "  ·  cells with N < min are shown grey; N annotated in each cell",
        fontsize=10,
    )

    fig_path_pdf = OUT_DIR / "fig02_resistance_landscape.pdf"
    fig_path_png = OUT_DIR / "fig02_resistance_landscape.png"
    fig.savefig(fig_path_pdf)
    fig.savefig(fig_path_png)
    fig_path_pdf, fig_path_png
    return fig, fig_path_pdf, fig_path_png


@app.cell
def _(counts, fig_path_pdf, fig_path_png):
    print(f"aggregation:")
    print(f"  {counts.n_dr_rows:,} DR rows  ×  {counts.n_curated_rows:,} curated rows")
    print(f"  → {counts.n_joined_rows:,} per-sequence-per-class rows")
    print(f"  → {counts.n_output_cells:,} landscape cells")
    print(f"     of which {counts.n_cells_above_min:,} have N ≥ min_cell_n")
    print()
    print(f"figures:")
    print(f"  PDF → {fig_path_pdf}")
    print(f"  PNG → {fig_path_png}")
    return


if __name__ == "__main__":
    app.run()
