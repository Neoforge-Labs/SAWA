"""Paper 1 style constants — colorblind palette, ordering, matplotlib rc.

Every Paper 1 figure imports `apply_paper_style()` at the top of its
notebook so the fonts, sizes, and color choices are consistent across
figures (a reviewer complaint otherwise).

Palette source: Okabe & Ito 2008 — the standard colorblind-safe
qualitative palette used in Nature journals + nf-core. Ordered so the
first 4 colors are also monochrome-print distinguishable.
"""

from __future__ import annotations

from typing import Any

# Okabe-Ito colorblind-safe qualitative palette. Reordered so the
# first four are also visually distinct in grayscale (blue → orange →
# green → purple).
COLORBLIND_PALETTE: tuple[str, ...] = (
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#F0E442",  # yellow
    "#56B4E9",  # sky blue
    "#D55E00",  # vermilion
    "#000000",  # black
)

# Canonical drug-class ordering used across figures. PI, then reverse
# transcriptase classes (NRTI, NNRTI), then integrase. Matches Sierra's
# `drugClass.name` enum values.
DRUG_CLASS_ORDER: tuple[str, ...] = ("PI", "NRTI", "NNRTI", "INSTI")

# Canonical subtype ordering (top East African subtypes first, then
# recombinants, then everything else). Kept short so faceted figures
# stay readable.
SUBTYPE_ORDER: tuple[str, ...] = ("A1", "C", "D", "recombinant")


def apply_paper_style() -> None:
    """Apply matplotlib rcParams for Paper 1 figures.

    Call once at the top of every figure notebook. Idempotent.
    """
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            # Fonts sized for a two-column journal spread; downstream can
            # scale with `plt.rcParams["figure.dpi"] = 150` etc.
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.titlesize": 11,
            # Publication-quality
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            # Cleaner spines
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.6,
            # No gridlines by default (heatmaps get their own; scatter
            # plots turn them on explicitly)
            "axes.grid": False,
            # PDF: keep text as text, not paths — reviewers can select it
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def resistance_cmap() -> Any:
    """Sequential colormap for resistance prevalence (0-1).

    We use viridis reversed so darker = higher prevalence, matching the
    convention in most HIV-DR surveillance papers. Colorblind-safe by
    construction.
    """
    import matplotlib as mpl

    return mpl.colormaps["viridis"]
