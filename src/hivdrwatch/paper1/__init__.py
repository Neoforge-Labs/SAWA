"""Paper 1 (tools paper) — shared helpers used by figure notebooks.

The figures themselves live in `notebooks/fig0N_*.py` per CLAUDE.md
("no figure lives outside notebooks/"). This package holds the pure-
Python building blocks those notebooks import: color palettes, data
aggregators, canonical column names. Keeping them here makes them
unit-testable and reusable across figures.
"""

from .style import (
    COLORBLIND_PALETTE,
    DRUG_CLASS_ORDER,
    SUBTYPE_ORDER,
    apply_paper_style,
    resistance_cmap,
)

__all__ = [
    "COLORBLIND_PALETTE",
    "DRUG_CLASS_ORDER",
    "SUBTYPE_ORDER",
    "apply_paper_style",
    "resistance_cmap",
]
