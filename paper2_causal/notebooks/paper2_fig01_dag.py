"""Render Paper 2 Figure 1: the DR domain DAG.

Reads the DAG from `causal/dag.py::DRDag.build_graph()` — no
duplication of edge lists here. Writes both PDF (for the manuscript
build) and PNG (for HTML preview) under `paper2_causal/figures/`.

Solid boxes for observed nodes, dashed for latent. Node colour by
kind (mutation, covariate, outcome, latent, mediator). Layout is
graphviz's `dot` when available, networkx spring layout otherwise.

Run: `python paper2_causal/notebooks/paper2_fig01_dag.py`
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

from causal.dag import DRDag

logger = logging.getLogger(__name__)

OUT_DIR = Path(__file__).resolve().parent.parent / "figures"

# Node styling by kind
_KIND_COLOUR: dict[str, str] = {
    "outcome": "#ffd7a0",  # warm amber
    "mutation": "#c6dbef",  # muted blue
    "covariate": "#c7e9c0",  # muted green
    "latent": "#f4d0d0",  # pale red
    "mediator": "#e1d5f0",  # pale purple
}
_LATENT_KINDS: frozenset[str] = frozenset({"latent", "mediator"})


def _kind_of(g: nx.DiGraph, node: str) -> str:
    return g.nodes[node].get("kind", "covariate")


def render(out_dir: Path = OUT_DIR) -> tuple[Path, Path]:
    """Render the DR DAG. Returns (pdf_path, png_path)."""
    g = DRDag().build_graph()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prefer graphviz dot layout — it makes causal graphs readable.
    # Fall back to spring layout if pygraphviz isn't installed.
    try:
        from networkx.drawing.nx_agraph import graphviz_layout

        pos = graphviz_layout(g, prog="dot")
    except ImportError:
        pos = nx.spring_layout(g, seed=42, k=1.5)

    fig, ax = plt.subplots(figsize=(9, 6))
    edge_widths = [1.2] * g.number_of_edges()

    nx.draw_networkx_edges(
        g,
        pos,
        ax=ax,
        edge_color="#555555",
        arrows=True,
        width=edge_widths,
        arrowsize=15,
        connectionstyle="arc3,rad=0.05",
    )

    # Solid boxes for observed, dashed for latent
    for node in g.nodes:
        kind = _kind_of(g, node)
        style = "dashed" if kind in _LATENT_KINDS else "solid"
        nx.draw_networkx_nodes(
            g,
            pos,
            nodelist=[node],
            ax=ax,
            node_color=_KIND_COLOUR.get(kind, "#eeeeee"),
            node_shape="s",
            node_size=2200,
            edgecolors="black",
            linewidths=1.2,
        )
        # Overlay dashed outline for latents by drawing a rectangle
        # patch directly — draw_networkx_nodes doesn't accept linestyle
        # on all matplotlib/networkx versions.
        if style == "dashed":
            from matplotlib.patches import FancyBboxPatch

            x, y = pos[node]
            # Rectangle sized to visually surround the ~2400-pt node
            w, h = 0.11, 0.06
            xlim = ax.get_xlim() or (x - 1, x + 1)
            ylim = ax.get_ylim() or (y - 1, y + 1)
            xspan = max(xlim[1] - xlim[0], 1.0)
            yspan = max(ylim[1] - ylim[0], 1.0)
            rect = FancyBboxPatch(
                (x - w * xspan / 2, y - h * yspan / 2),
                w * xspan,
                h * yspan,
                boxstyle="round,pad=0.005",
                linewidth=1.5,
                linestyle="dashed",
                edgecolor="black",
                facecolor="none",
                zorder=3,
            )
            ax.add_patch(rect)

    nx.draw_networkx_labels(g, pos, ax=ax, font_size=8, font_family="sans-serif")

    # Legend
    from matplotlib.patches import Patch  # local import to keep module top light

    legend_items = [
        Patch(
            facecolor=_KIND_COLOUR["outcome"],
            edgecolor="black",
            label="Outcome (Sierra level)",
        ),
        Patch(
            facecolor=_KIND_COLOUR["mutation"],
            edgecolor="black",
            label="Mutation feature (observed)",
        ),
        Patch(
            facecolor=_KIND_COLOUR["covariate"],
            edgecolor="black",
            label="Covariate (observed)",
        ),
        Patch(
            facecolor=_KIND_COLOUR["mediator"],
            edgecolor="black",
            label="Mediator (latent)",
            linestyle="dashed",
        ),
        Patch(
            facecolor=_KIND_COLOUR["latent"],
            edgecolor="black",
            label="Latent confounder",
            linestyle="dashed",
        ),
    ]
    ax.legend(handles=legend_items, loc="upper left", fontsize=7, frameon=True)

    ax.set_axis_off()
    ax.set_title(
        "SAWA drug-resistance domain DAG\n"
        "(latent nodes → dashed outline; adjustment set = {subtype, country, collection_year})",
        fontsize=10,
    )
    fig.tight_layout()

    pdf_path = out_dir / "fig01_dr_dag.pdf"
    png_path = out_dir / "fig01_dr_dag.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("wrote %s and %s", pdf_path, png_path)
    return pdf_path, png_path


if __name__ == "__main__":
    logging.basicConfig(level="INFO", format="%(levelname)s %(name)s: %(message)s")
    render()
