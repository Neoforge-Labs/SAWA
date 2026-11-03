"""Ingestion module — sources HIV pol sequences from public databases.

The source ingesters are imported on first use. `genbank` pulls in Biopython, so importing
the schema (as the dashboard does for its country list) must not pay for it.
"""

from __future__ import annotations

import importlib
from typing import Any

from .schema import EAST_AFRICA_COUNTRIES, SCHEMA_VERSION, SequenceRecord, Source, Subtype

__all__ = [
    "EAST_AFRICA_COUNTRIES",
    "SCHEMA_VERSION",
    "GenBankIngester",
    "LANLIngester",
    "SequenceRecord",
    "Source",
    "Subtype",
    "merge_sources",
]

_LAZY_MODULES: dict[str, str] = {
    "GenBankIngester": ".genbank",
    "LANLIngester": ".lanl",
    "merge_sources": ".merge",
}


def __getattr__(name: str) -> Any:
    module = _LAZY_MODULES.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module, __name__), name)
    globals()[name] = value
    return value
