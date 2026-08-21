"""Nextstrain / Auspice JSON export for SAWA per-subtype trees.

CONTEXT.md commits SAWA to *embedding* Auspice (Nextstrain) rather
than reinventing tree visualisation. That means every per-subtype
time-scaled tree emitted by `phylo/treetime.py` must be exportable
to the Auspice v2 JSON schema
(https://nextstrain.org/docs/bioinformatics/data-formats#version-2)
so the dashboard can `<auspice-client>` the tree directly.

## What this module exports

For each subtype we produce **one** Auspice-v2 dataset JSON:

```
data/gold/auspice/hivdrwatch_<subtype>.json
```

The dataset includes:

* `tree` — the nested-JSON tree with branch_length, num_date, and
  node_attrs {country, subtype, resistance_high (bool), sample_year}
* `meta` — dataset title, colour scales for country and resistance,
  build metadata (SAWA version, git SHA, algorithm version)
* `version` — always 2 (Auspice v2 schema)

The exporter does **not** try to be a general phylodynamics
visualiser — it is a thin adapter from SAWA's per-subtype tree +
DR + curated Parquets to the shape Auspice expects. Colour scales
follow SAWA's `paper1.style` palette so the dashboard matches
the manuscript figures.

## Contract

```python
from pathlib import Path
from hivdrwatch.phylo.auspice import export_subtype_dataset

dataset_path = export_subtype_dataset(
    subtype="A1",
    tree_newick=Path("data/interim/trees/A1/time_scaled.nwk"),
    curated_parquet=Path("data/interim/curated/aligned.parquet"),
    dr_parquet=Path("data/interim/dr_results.parquet"),
    out_dir=Path("data/gold/auspice"),
    sawa_version="v0.1.0",
    git_sha="abc123",
)
# → data/gold/auspice/hivdrwatch_A1.json
```

The result is directly consumable by Auspice as either a hosted
dataset (`nextstrain.org/community/user/repo/<subtype>`) or via
`auspice view --datasetDir data/gold/auspice/`.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

logger = logging.getLogger(__name__)

AUSPICE_JSON_VERSION: int = 2
DATASET_TITLE_TEMPLATE: str = "SAWA HIV-1 subtype {subtype} — East African cohort"
DEFAULT_MAINTAINERS: tuple[dict[str, str], ...] = (
    {"name": "Kelyn Njeri", "url": "https://github.com/Neoforge-Labs/SAWA"},
)
DEFAULT_BUILD_URL: str = "https://github.com/Neoforge-Labs/SAWA"

# Sierra levels that count as "resistant" in the boolean node attr.
# Matches the RESISTANT_LEVELS convention used throughout SAWA.
RESISTANT_LEVELS: frozenset[str] = frozenset({"high", "intermediate"})

# Country colour palette — chosen for accessibility (Okabe-Ito subset).
COUNTRY_COLOUR: dict[str, str] = {
    "KE": "#0072B2",  # blue
    "UG": "#D55E00",  # vermillion
    "TZ": "#009E73",  # green
    "RW": "#CC79A7",  # magenta
    "BI": "#F0E442",  # yellow
    "ET": "#56B4E9",  # sky blue
    "unknown": "#999999",
}


@dataclass(frozen=True)
class AuspiceExportCounts:
    """Reproducibility summary for a single subtype export."""

    subtype: str
    n_tips_in_tree: int
    n_tips_with_metadata: int
    n_tips_with_dr: int
    n_resistant: int
    output_path: str


def export_subtype_dataset(
    *,
    subtype: str,
    tree_newick: Path,
    curated_parquet: Path,
    dr_parquet: Path,
    out_dir: Path,
    sawa_version: str,
    git_sha: str,
    resistance_drug_class: str | None = None,
) -> AuspiceExportCounts:
    """Export one subtype's time-scaled tree + metadata as Auspice-v2 JSON.

    `resistance_drug_class` (optional) restricts the `resistance_high`
    boolean node attr to the given class (e.g. "INSTI"); when None,
    a tip is marked resistant if any drug's Sierra level is in
    RESISTANT_LEVELS.
    """
    try:
        from Bio import Phylo
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "Biopython required for Auspice export; install with `uv add biopython`"
        ) from e

    tree = Phylo.read(tree_newick, "newick")
    tips = [t for t in tree.get_terminals()]
    tip_accessions = [_strip_accession(t.name) for t in tips]

    curated = pl.read_parquet(curated_parquet).select(
        ["accession", "country", "subtype", "collection_year"]
    )
    dr = pl.read_parquet(dr_parquet)
    if resistance_drug_class:
        dr = dr.filter(pl.col("drug_class") == resistance_drug_class)

    resistant_accessions = set(
        dr.filter(pl.col("level").is_in(list(RESISTANT_LEVELS)))["accession"]
        .unique()
        .to_list()
    )
    metadata_lookup = {
        r["accession"]: r
        for r in curated.filter(pl.col("subtype") == subtype).iter_rows(named=True)
    }
    # Also permit metadata lookup for accessions whose Sierra subtype
    # differs from the LANL curated subtype — surveillance data is
    # noisy, don't drop those tips silently.
    all_metadata = {r["accession"]: r for r in curated.iter_rows(named=True)}

    n_with_metadata = 0
    n_with_dr = 0
    n_resistant = 0
    for tip, acc in zip(tips, tip_accessions, strict=True):
        meta = metadata_lookup.get(acc, all_metadata.get(acc, {}))
        country = meta.get("country", "unknown") or "unknown"
        year = meta.get("collection_year")
        is_resistant = acc in resistant_accessions
        _annotate_tip(
            tip, country=country, year=year, is_resistant=is_resistant, subtype=subtype
        )
        if meta:
            n_with_metadata += 1
        if (
            any(r["accession"] == acc for r in dr.head(0).iter_rows(named=True))
            or acc in resistant_accessions
        ):
            n_with_dr += 1
        if is_resistant:
            n_resistant += 1

    _propagate_internal_dates(tree)

    dataset = {
        "version": AUSPICE_JSON_VERSION,
        "meta": _build_meta(
            subtype=subtype,
            sawa_version=sawa_version,
            git_sha=git_sha,
            resistance_drug_class=resistance_drug_class,
        ),
        "tree": _biopython_to_auspice(tree.clade),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"hivdrwatch_{_slug(subtype)}.json"
    out_path.write_text(
        json.dumps(dataset, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    logger.info(
        "wrote Auspice dataset → %s (tips=%d, resistant=%d)",
        out_path,
        len(tips),
        n_resistant,
    )
    return AuspiceExportCounts(
        subtype=subtype,
        n_tips_in_tree=len(tips),
        n_tips_with_metadata=n_with_metadata,
        n_tips_with_dr=n_with_dr,
        n_resistant=n_resistant,
        output_path=str(out_path),
    )


# ------------------------------------------------------------------
# internals
# ------------------------------------------------------------------


def _strip_accession(name: str) -> str:
    """Recover the GenBank accession from a tree tip name.

    Tip names in SAWA trees typically carry the LANL FASTA header
    convention, which is dot-delimited with fields:

        Subtype.Country.Year.Name.PatientID.CDID.SEID.Accession.Length.Start.End.HIV-1

    So we split on '.' and take position [7]. If the name has fewer
    fields (e.g. bare accession, or a GenBank-only tip name), fall
    back to splitting on '|' (older format) or returning the whole
    string.
    """
    if not name:
        return ""
    parts = name.split(".")
    if len(parts) >= 8:
        return parts[7]
    if "|" in name:
        return name.split("|", 1)[0]
    return name


def _annotate_tip(
    tip: Any, *, country: str, year: int | None, is_resistant: bool, subtype: str
) -> None:
    """Attach Auspice-compatible node_attrs to a Biopython terminal clade."""
    if not hasattr(tip, "node_attrs"):
        tip.node_attrs = {}
    tip.node_attrs = {
        "country": {"value": country},
        "subtype": {"value": subtype},
        "resistance_high": {"value": is_resistant},
    }
    if year is not None:
        tip.node_attrs["num_date"] = {"value": float(year)}


def _propagate_internal_dates(tree: Any) -> None:
    """Assign `num_date` to internal nodes as the min of descendant tip dates.

    TreeTime normally sets these on the time-scaled tree; if Biopython
    lost them during round-trip, this ensures Auspice still renders.
    """
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            continue
        child_dates = []
        for child in clade.clades:
            if hasattr(child, "node_attrs") and "num_date" in child.node_attrs:
                child_dates.append(child.node_attrs["num_date"]["value"])
        if child_dates:
            if not hasattr(clade, "node_attrs"):
                clade.node_attrs = {}
            clade.node_attrs["num_date"] = {"value": min(child_dates)}


def _biopython_to_auspice(clade: Any) -> dict:
    """Convert a Biopython Clade (subtree root) to Auspice-v2 tree JSON."""
    node: dict = {
        "name": clade.name or _auto_internal_name(clade),
        "node_attrs": getattr(clade, "node_attrs", {}),
    }
    if clade.branch_length is not None:
        node["branch_attrs"] = {"length": float(clade.branch_length)}
    children = list(clade.clades)
    if children:
        node["children"] = [_biopython_to_auspice(c) for c in children]
    return node


_INTERNAL_COUNTER = {"n": 0}


def _auto_internal_name(clade: Any) -> str:
    _INTERNAL_COUNTER["n"] += 1
    return f"NODE_{_INTERNAL_COUNTER['n']:05d}"


def _build_meta(
    *,
    subtype: str,
    sawa_version: str,
    git_sha: str,
    resistance_drug_class: str | None,
) -> dict:
    """Auspice meta block — title, colourings, filters, build info."""
    resistance_label = (
        f"{resistance_drug_class} resistance (Sierra intermediate or high)"
        if resistance_drug_class
        else "Any-class resistance (Sierra intermediate or high)"
    )
    return {
        "title": DATASET_TITLE_TEMPLATE.format(subtype=subtype),
        "updated": date.today().isoformat(),
        "genome_annotations": {},
        "colorings": [
            {
                "key": "country",
                "title": "Country",
                "type": "categorical",
                "scale": [[k, v] for k, v in COUNTRY_COLOUR.items()],
            },
            {
                "key": "resistance_high",
                "title": resistance_label,
                "type": "categorical",
                "scale": [[True, "#D55E00"], [False, "#0072B2"]],
            },
            {
                "key": "num_date",
                "title": "Sample year",
                "type": "continuous",
            },
        ],
        "display_defaults": {
            "color_by": "resistance_high",
            "branch_label": "none",
            "map_triplicate": False,
        },
        "filters": ["country", "resistance_high"],
        "panels": ["tree", "map"],
        "maintainers": list(DEFAULT_MAINTAINERS),
        "build_url": DEFAULT_BUILD_URL,
        "extensions": {
            "sawa": {
                "version": sawa_version,
                "git_sha": git_sha,
                "resistance_drug_class": resistance_drug_class,
            }
        },
    }


def _slug(subtype: str) -> str:
    """Filesystem-safe subtype identifier for Auspice URL routing."""
    return "".join(c if c.isalnum() else "_" for c in subtype).strip("_") or "unknown"


__all__ = [
    "AUSPICE_JSON_VERSION",
    "COUNTRY_COLOUR",
    "RESISTANT_LEVELS",
    "AuspiceExportCounts",
    "export_subtype_dataset",
]
