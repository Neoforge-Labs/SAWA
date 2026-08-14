"""Canonical schemas for phylogenetic tree artifacts.

The `phylo/` module produces two things per subtype:
  1. A raw Newick file for the ML tree from IQ-TREE 2
  2. A time-scaled Newick + node-metadata JSON from TreeTime

We keep the raw files as artifacts on disk (`data/interim/trees/<subtype>/`)
and record a `TreeMetadata` row per subtype in a manifest parquet so
downstream consumers (cluster validation, DTG case study, Auspice
export) can look up "which tree, what vintage, built how" without
re-parsing Newick.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, Field

PHYLO_SCHEMA_VERSION = "0.1.0"


class TreeKind(StrEnum):
    """The two flavors of tree we build per subtype."""

    ML = "ml"  # maximum-likelihood, IQ-TREE 2
    TIME_SCALED = "time_scaled"  # dated, TreeTime


class SubstitutionModel(StrEnum):
    """Substitution models we allow — pinned per MODULES.md.

    HIV pol is well-served by GTR+F+I+G4:
      * GTR: general time-reversible substitution
      * +F: empirical base frequencies
      * +I: proportion of invariant sites
      * +G4: 4-category discrete gamma rate heterogeneity
    """

    GTR_F_I_G4 = "GTR+F+I+G4"


class TreeMetadata(BaseModel):
    """One row per (subtype, kind) tree artifact — indexed manifest."""

    subtype: str = Field(..., min_length=1, max_length=20)
    kind: TreeKind
    n_taxa: int = Field(..., ge=1)

    # --- artifact paths (relative to repo root) ---
    newick_path: str
    log_path: str | None = None
    dates_csv_path: str | None = None  # only for TIME_SCALED

    # --- IQ-TREE specifics ---
    substitution_model: SubstitutionModel | None = None
    log_likelihood: float | None = None
    bootstrap_replicates: int | None = None

    # --- TreeTime specifics ---
    clock_rate: float | None = None  # substitutions per site per year
    clock_rate_stddev: float | None = None
    time_span_years: float | None = None  # tip-date range

    # --- reproducibility ---
    random_seed: int
    tool_version: str  # from `iqtree2 --version` etc.
    wall_seconds: float | None = None

    # --- provenance ---
    built_at: date
    schema_version: str = PHYLO_SCHEMA_VERSION


class SubtypePartitionRow(BaseModel):
    """Records which sequences went into each subtype's tree build."""

    subtype: str
    accession: str
    collection_year: int | None = None
    country: str | None = None


def default_tree_dir(root: Path, subtype: str) -> Path:
    """Standard output directory for per-subtype tree artifacts."""
    return root / "interim" / "trees" / subtype
