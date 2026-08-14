"""Per-country GenBank ingest watermarks.

E-utilities `since` filters by *modification date* (MDAT), which is the
right semantics for incremental pulls: re-running with `since=<last>`
returns anything newly deposited since we last checked, plus any older
records that were re-uploaded / annotated.

Layout:
    <genbank_raw_dir>/.watermarks.json
        {
          "KE": "2026-08-11",
          "UG": "2026-08-11",
          "TZ": "2026-07-30"
        }

The file is atomic-written (tmp + rename) so a crash mid-write can't
leave a truncated JSON that would trip the next run.
"""

from __future__ import annotations

import json
import logging
from datetime import date
from pathlib import Path

logger = logging.getLogger(__name__)

WATERMARK_FILENAME = ".watermarks.json"


def watermark_path(genbank_raw_dir: Path) -> Path:
    return genbank_raw_dir / WATERMARK_FILENAME


def read_watermarks(genbank_raw_dir: Path) -> dict[str, date]:
    """Return `{country_code: last_ingest_date}`. Empty dict if the file
    doesn't exist yet (first-ever run)."""
    p = watermark_path(genbank_raw_dir)
    if not p.exists():
        return {}
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        logger.warning("corrupt watermark file at %s; treating as empty", p)
        return {}
    out: dict[str, date] = {}
    for cc, iso in raw.items():
        try:
            out[cc.upper()] = date.fromisoformat(iso)
        except ValueError:
            logger.warning("bad watermark date for %s: %r; skipping", cc, iso)
    return out


def write_watermarks(
    genbank_raw_dir: Path, watermarks: dict[str, date]
) -> Path:
    p = watermark_path(genbank_raw_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = {cc.upper(): d.isoformat() for cc, d in watermarks.items()}
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(p)
    return p


def get_since(
    watermarks: dict[str, date], country: str, *, fallback: date | None = None
) -> date | None:
    """Return the appropriate `since` date for a country's incremental pull.

    Falls back to `fallback` (e.g. an all-time floor like 1980-01-01)
    when there's no watermark yet — that's the "first ever pull" behavior.
    """
    cc = country.upper()
    if cc in watermarks:
        return watermarks[cc]
    return fallback


def bump_watermark(
    watermarks: dict[str, date], country: str, new_date: date | None = None
) -> dict[str, date]:
    """Return a new dict with `country`'s watermark set to `new_date`
    (defaults to today). Doesn't mutate the input; caller re-writes."""
    updated = dict(watermarks)
    updated[country.upper()] = new_date or date.today()
    return updated
