"""Pipeline-architecture spec for Paper 1 Figure 1.

Kept out of the notebook so the layout is unit-testable and so we can
re-use the block definitions in the paper's methods-section text.

The diagram is stage-by-stage top-to-bottom. Each row is a phase, each
box is a concrete module in `src/hivdrwatch/`. Arrows show data flow;
group boxes show what maps to which colour band.

If you add a stage, add it here — the notebook draws whatever this
module returns, no separate list to keep in sync.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

# --- Phase colour bands (Okabe-Ito palette, matches paper1.style) ---
PHASE_COLORS: dict[str, str] = {
    "source": "#F0E442",  # yellow  — external inputs
    "ingest": "#0072B2",  # blue    — bronze Parquet
    "call": "#009E73",  # green   — Sierra resistance calls
    "curate": "#56B4E9",  # sky     — align + trim + QC + dedup → silver
    "analyze": "#E69F00",  # orange  — phylo + cluster (gold)
    "serve": "#CC79A7",  # purple  — storage + API + dashboard
}

Phase = Literal["source", "ingest", "call", "curate", "analyze", "serve"]


@dataclass(frozen=True)
class Block:
    """One rendered box in the diagram."""

    key: str  # unique id used in edges
    label: str  # visible text
    module: str  # e.g. "hivdrwatch/ingest"
    phase: Phase
    row: int  # 0 = top
    col: float  # horizontal placement (float ok for centering)


@dataclass(frozen=True)
class Edge:
    """One arrow from one block to another."""

    src: str
    dst: str
    # Optional label ("bronze", "silver", "gold") for the Parquet
    # tier the arrow carries. Rendered as small text on the arrow.
    tier: str | None = None


@dataclass(frozen=True)
class PipelineDiagram:
    """The whole diagram spec — blocks + edges."""

    blocks: tuple[Block, ...]
    edges: tuple[Edge, ...]
    row_labels: tuple[tuple[int, str], ...] = field(default_factory=tuple)


def pipeline_diagram() -> PipelineDiagram:
    """Return the canonical Paper 1 architecture.

    Layout: 5-col grid (0..4). Column placement chosen so most edges
    are straight-down or short diagonals; only one arrow spans more
    than one column laterally.

      col:      1              2              3              4
      row 0:    LANL                          GenBank
      row 1:                   ingest
      row 2:    align                         sierra
      row 3:    cluster                       phylo             (spread wide)
      row 4:                   figures 2-5    storage          api

    The Products row is deliberately shifted right so `storage → api`
    is a short adjacent edge (they're the serving pair), and `sierra
    → storage` and `phylo → storage` become clean straight-down edges.
    Figs sits directly under align, so `align → cluster → figs` reads
    as a straight vertical spine.

    We drop the `phylo → figures` edge — no Paper 1 figure directly
    consumes phylogenetic trees. Phylo output reaches storage (and
    eventually the API dashboard) but not the paper's figures.
    """

    blocks: tuple[Block, ...] = (
        # Row 0 — external inputs
        Block("lanl", "LANL Table\nTSV + FASTA", "public", "source", 0, 1.0),
        Block("gb", "GenBank\nNCBI Entrez", "public", "source", 0, 3.0),
        # Row 1 — ingest merge (centered)
        Block(
            "ingest",
            "ingest\nmerge · dedup · watermark",
            "src/hivdrwatch/ingest",
            "ingest",
            1,
            2.0,
        ),
        # Row 2 — align (left, feeds cluster+phylo+figs), sierra (right, feeds storage+figs)
        Block(
            "align",
            "align\nmappy → HXB2 · trim · QC",
            "src/hivdrwatch/align",
            "curate",
            2,
            1.0,
        ),
        Block(
            "sierra",
            "resistance\nStanford HIVdb Sierra",
            "src/hivdrwatch/resistance",
            "call",
            2,
            3.0,
        ),
        # Row 3 — gold analytics: cluster (left, downstream of align), phylo (right, downstream of align)
        Block(
            "cluster",
            "cluster\ntn93 · connected components",
            "src/hivdrwatch/cluster",
            "analyze",
            3,
            1.0,
        ),
        Block(
            "phylo", "phylo\nIQ-TREE 2 · TreeTime", "src/hivdrwatch/phylo", "analyze", 3, 3.0
        ),
        # Row 4 — serve: figs under the analytical spine, storage+api on the right
        Block("figs", "figures 2-5\nMarimo notebooks", "notebooks/", "serve", 4, 2.0),
        Block(
            "storage",
            "storage\nPostgres · DuckDB · DVC",
            "src/hivdrwatch/storage",
            "serve",
            4,
            3.0,
        ),
        Block("api", "api\nFastAPI · OpenAPI", "src/hivdrwatch/api", "serve", 4, 4.0),
    )

    edges: tuple[Edge, ...] = (
        # Sources → ingest
        Edge("lanl", "ingest"),
        Edge("gb", "ingest"),
        # Ingest → parallel analytics (bronze)
        Edge("ingest", "align", tier="bronze"),
        Edge("ingest", "sierra", tier="bronze"),
        # Align → downstream analytics (silver)
        Edge("align", "cluster", tier="silver"),
        Edge("align", "phylo", tier="silver"),
        # Gold → products
        Edge("cluster", "figs"),  # cluster_validation.json → Fig 3
        Edge("sierra", "figs"),  # dr_results.parquet → Fig 2, Fig 4
        Edge("sierra", "storage", tier="gold"),
        Edge("cluster", "storage", tier="gold"),
        Edge("phylo", "storage", tier="gold"),
        # Serve
        Edge("storage", "api"),
    )

    row_labels: tuple[tuple[int, str], ...] = (
        (0, "Sources"),
        (1, "Bronze"),
        (2, "Silver"),
        (3, "Gold"),
        (4, "Products"),
    )

    return PipelineDiagram(blocks=blocks, edges=edges, row_labels=row_labels)


def stage_summary() -> dict[str, int]:
    """Count blocks per phase — used by the notebook's caption line."""
    diagram = pipeline_diagram()
    counts: dict[str, int] = {phase: 0 for phase in PHASE_COLORS}
    for b in diagram.blocks:
        counts[b.phase] += 1
    return counts


__all__ = [
    "PHASE_COLORS",
    "Block",
    "Edge",
    "PipelineDiagram",
    "pipeline_diagram",
    "stage_summary",
]
