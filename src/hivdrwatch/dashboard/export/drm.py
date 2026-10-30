"""Surveillance drug-resistance mutation (SDRM) matching for Sierra mutation strings.

Sierra reports every variant relative to consensus, so the most frequent
mutation strings are natural polymorphisms (PR R41K, RT V35T, ...), not
resistance. A "top mutations" view built on the raw column would show
background variation, so everything the dashboard publishes about mutations
goes through `match_sdrm`, which keeps only the Stanford surveillance list.

Integrase is deliberately not a panel. The surveillance list is cited to
Bennett et al. 2009, which has no integrase component, and the reference
file's `IN` block is the raltegravir/elvitegravir-era list, missing
dolutegravir-era mutations such as R263K, G118R and E138K.

The reference file stores positions and resistant states but no consensus
residues, so the reference residue in a label is the leading letter of the
mutation string itself.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

SDRM_PANELS: tuple[str, ...] = ("PR", "RT_NRTI", "RT_NNRTI")

# The parquet `gene` column names the enzyme, while the reference file splits RT by drug class.
_PANEL_GENE: dict[str, str] = {"PR": "PR", "RT_NRTI": "RT", "RT_NNRTI": "RT"}

_INSERTION_STATE = "ins"

_SUBSTITUTION = re.compile(r"^([A-Z])(\d+)([A-Z]+)$")
_INSERTION = re.compile(r"^([A-Z])(\d+)[A-Z]_[A-Z]+$")


class MissingPanelError(ValueError):
    """The reference file has no block for a panel the dashboard requires."""


@dataclass(frozen=True)
class SdrmPanel:
    """One surveillance panel: resistant amino-acid states per position."""

    name: str
    gene: str
    positions: dict[int, frozenset[str]]


def load_sdrm_panels(path: Path) -> tuple[tuple[SdrmPanel, ...], str]:
    """Load the `SDRM_PANELS` blocks of the reference file and its pinned list version."""
    reference = json.loads(path.read_text())
    genes: dict[str, dict[str, dict[str, list[str]]]] = reference["genes"]

    absent = [name for name in SDRM_PANELS if name not in genes]
    if absent:
        raise MissingPanelError(
            f"{path} has no block for panel(s) {absent}; found {sorted(genes)}"
        )

    panels = tuple(
        SdrmPanel(
            name=name,
            gene=_PANEL_GENE[name],
            positions={
                int(pos): frozenset(states) for pos, states in genes[name]["positions"].items()
            },
        )
        for name in SDRM_PANELS
    )
    version: str = reference["algorithm_version_pinned"]
    logger.info(
        "loaded SDRM panels %s (%s)",
        {p.name: len(p.positions) for p in panels},
        version,
    )
    return panels, version


def _parse(raw: str) -> tuple[str, int, frozenset[str]] | None:
    # Stop codons and deletions never match either pattern: `*` and lowercase `del`
    # fall outside [A-Z].
    insertion = _INSERTION.match(raw)
    if insertion:
        return insertion[1], int(insertion[2]), frozenset({_INSERTION_STATE})
    substitution = _SUBSTITUTION.match(raw)
    if substitution:
        return substitution[1], int(substitution[2]), frozenset(substitution[3])
    return None


def parse_mutation(raw: str) -> tuple[int, frozenset[str]] | None:
    """Position and observed states of a Sierra mutation string, or None if unparseable."""
    parsed = _parse(raw)
    if parsed is None:
        return None
    _, position, states = parsed
    return position, states


def match_sdrm(raw: str, panel: SdrmPanel) -> str | None:
    """Return the canonical `<ref><pos><state>` label if `raw` carries an SDRM in `panel`.

    A mixture carries an SDRM when any of its states is a resistant state at that
    position (the surveillance convention). When several match, the first in the
    panel's sorted state order is used so the label is identical on every run.
    """
    parsed = _parse(raw)
    if parsed is None:
        return None
    reference, position, observed = parsed

    resistant = panel.positions.get(position)
    if resistant is None:
        return None
    matched = sorted(resistant & observed)
    if not matched:
        return None
    return f"{reference}{position}{matched[0]}"
