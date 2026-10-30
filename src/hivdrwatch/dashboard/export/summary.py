"""Deterministic `summary.json`: the single file the static dashboard page draws.

The summary is committed and read in diffs, and CI compares the input hash
recorded in `meta` against the DVC pointer, so two runs over the same input
must agree byte for byte. That rests on three things kept here: every
aggregation already ends in an explicit sort, keys are serialised sorted, and
floats are rounded to a fixed number of places before they reach the encoder.

`meta` records the measured counts and provenance of the input vintage but
never asserts them. The parquet is regenerated over time, so only properties
that hold for any vintage fail the export: lossless deduplication (raised by
`dedup_dr_results`) and a single algorithm version, publish date and call date
across the frame. The summary makes no geographic claim about the sequences.

Summary schema 1.2.0: `meta` (provenance and counts) and `views` (four lists of
row dicts: prevalence_by_drug_class, levels_by_drug, gene_coverage, sdrm_frequency).

Migration 1.0.0 -> 1.1.0 (minor: new required fields). `meta` gained
`n_sdrm_observations` and `n_sdrm_observations_multistate`, which let the page state
the measured size of the single-state attribution rule instead of a typed-in figure,
and `input_file_md5`, the digest of the file actually read, kept apart from
`input_dvc_md5`, which is only what the DVC pointer claims. The static page requires
all three, so a 1.0.0 file is refused by its schema gate rather than half-drawn.

Migration 1.1.0 -> 1.2.0 (minor: new required field). `meta` gained `n_distinct_sequences`,
the number of distinct `sequence_sha256` values among the scored accessions. Every count
in the summary is of accessions, and several accessions can carry one sequence, so the
page needs this number to say how far the unit differs from distinct sequences. The static
page requires it, so a 1.1.0 file is refused by its schema gate rather than half-drawn.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Any

import polars as pl

from hivdrwatch.dashboard.export.aggregate import (
    gene_coverage,
    levels_by_drug,
    prevalence_by_drug_class,
    sdrm_frequency,
)
from hivdrwatch.dashboard.export.dedup import dedup_dr_results
from hivdrwatch.dashboard.export.drm import SdrmPanel, parse_mutation
from hivdrwatch.resistance.schema import DR_SCHEMA_VERSION

logger = logging.getLogger(__name__)

SUMMARY_SCHEMA_VERSION: str = "1.2.0"

FLOAT_PLACES: int = 4

_HASH_CHUNK_BYTES: int = 1 << 20

_PROVENANCE_COLUMNS: tuple[str, ...] = (
    "algorithm_version",
    "algorithm_publish_date",
    "called_at",
)

# A DVC file-output pointer stores a 32-hex md5. A directory output ends in `.dir`, which
# this deliberately rejects: the dashboard input is a single parquet file.
_MD5_LINE = re.compile(r"^\s*(?:-\s+)?md5:\s*([0-9a-f]{32})\s*$", re.MULTILINE)


class MissingDvcHashError(ValueError):
    """A DVC pointer has no `md5:` line holding a 32-character digest."""


class NonSingleValuedProvenanceError(ValueError):
    """A provenance column takes more than one value, so `meta` cannot name one."""


def input_dvc_md5(pointer: Path) -> str:
    """Return the md5 recorded in a `.dvc` pointer file."""
    match = _MD5_LINE.search(pointer.read_text())
    if match is None:
        raise MissingDvcHashError(f"{pointer} has no 'md5:' line with a 32-character digest")
    return match.group(1)


def check_input_not_stale(summary: Path, pointer: Path) -> int:
    """Return 0 if `meta.input_dvc_md5` equals the pointer's md5, else log an ERROR and return 1.

    Narrow by construction. It compares one pointer value against another, so it detects
    that `dvc add` moved the pointer after the summary was exported. It does not detect a
    parquet that disagrees with its own pointer, and it does not re-derive any aggregate:
    the only DVC remote is a localhost MinIO, so CI never holds the data. `input_file_md5`
    is the digest of the file actually read at export time; it cannot be checked here for
    the same reason, and is deliberately not consulted.
    """
    meta = json.loads(summary.read_text()).get("meta", {})
    recorded = meta.get("input_dvc_md5")
    current = input_dvc_md5(pointer)
    if recorded == current:
        return 0
    logger.error(
        "%s records input_dvc_md5=%s but %s now says md5=%s; re-run the export against the "
        "current input",
        summary,
        recorded,
        pointer,
        current,
    )
    return 1


def file_md5(path: Path, *, chunk_bytes: int = _HASH_CHUNK_BYTES) -> str:
    """Return the md5 of the file's bytes, read in chunks so the whole file is never in memory."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def _single_value(df: pl.DataFrame, column: str) -> str | None:
    if column not in df.columns:
        raise ValueError(f"dr_results missing required column: {column}")
    values = df[column].unique().sort().to_list()
    if len(values) != 1:
        raise NonSingleValuedProvenanceError(
            f"{column} must take exactly one value across the frame, found {values}"
        )
    value = values[0]
    return None if value is None else str(value)


def sdrm_observation_counts(
    df: pl.DataFrame, panels: tuple[SdrmPanel, ...]
) -> tuple[int, int]:
    """Count surveillance-mutation observations, and how many matched more than one state.

    An observation is one mutation string on one accession in one panel that matches
    at least one of the panel's resistant states at its position; it is multi-state when
    it matches two or more. `match_sdrm` attributes a multi-state observation to a single
    label, so the second count is the size of that attribution rule's effect. Rows are
    reduced to one per accession and gene first, because `mutations` repeats on every
    drug row of a gene.
    """
    total = 0
    multistate = 0
    for panel in panels:
        per_accession = df.filter(pl.col("gene") == panel.gene).unique(
            subset="accession", keep="first", maintain_order=True
        )
        raw_counts = (
            per_accession["mutations"]
            .explode(empty_as_null=False)
            .drop_nulls()
            .value_counts()
            .iter_rows()
        )
        for raw, count in raw_counts:
            parsed = parse_mutation(raw)
            if parsed is None:
                continue
            position, observed = parsed
            matched = len(panel.positions.get(position, frozenset()) & observed)
            if matched:
                total += count
            if matched > 1:
                multistate += count
    return total, multistate


def _round_floats(value: Any) -> Any:
    if isinstance(value, float):
        return round(value, FLOAT_PLACES)
    if isinstance(value, dict):
        return {key: _round_floats(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_round_floats(item) for item in value]
    return value


def build_summary(
    df: pl.DataFrame,
    *,
    panels: tuple[SdrmPanel, ...],
    sdrm_version: str,
    input_md5: str,
    input_file_md5: str,
    git_sha: str,
    n_accessions_ingested: int,
) -> dict[str, Any]:
    """Deduplicate `df`, aggregate it, and return `{"meta": ..., "views": ...}`.

    Exceptions from `dedup_dr_results` propagate unchanged. Floats in the result
    are already rounded, so what a caller inspects is what `write_summary` writes.
    """
    deduped, report = dedup_dr_results(df)

    provenance = {column: _single_value(deduped, column) for column in _PROVENANCE_COLUMNS}
    n_observations, n_multistate = sdrm_observation_counts(deduped, panels)

    meta: dict[str, Any] = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "dr_schema_version": DR_SCHEMA_VERSION,
        **provenance,
        "sdrm_list_version": sdrm_version,
        "input_dvc_md5": input_md5,
        "input_file_md5": input_file_md5,
        "git_sha": git_sha,
        "n_rows_raw": report.n_rows_in,
        "n_calls": report.n_rows_out,
        "n_accessions_scored": deduped["accession"].n_unique(),
        "n_distinct_sequences": report.n_distinct_sequences,
        "n_accessions_ingested": n_accessions_ingested,
        "n_sdrm_observations": n_observations,
        "n_sdrm_observations_multistate": n_multistate,
    }
    views: dict[str, list[dict[str, Any]]] = {
        "prevalence_by_drug_class": prevalence_by_drug_class(deduped).to_dicts(),
        "levels_by_drug": levels_by_drug(deduped).to_dicts(),
        "gene_coverage": gene_coverage(deduped).to_dicts(),
        "sdrm_frequency": sdrm_frequency(deduped, panels).to_dicts(),
    }
    logger.info(
        "summary: %d calls over %d accessions scored (%d ingested)",
        meta["n_calls"],
        meta["n_accessions_scored"],
        n_accessions_ingested,
    )
    result: dict[str, Any] = _round_floats({"meta": meta, "views": views})
    return result


def render_summary(summary: dict[str, Any]) -> str:
    """Serialise to the exact text `write_summary` stores."""
    # NaN has no JSON spelling; refuse it rather than emit a file other parsers reject.
    text = json.dumps(
        _round_floats(summary), sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False
    )
    return text + "\n"


def write_summary(summary: dict[str, Any], out: Path) -> None:
    """Write `summary` to `out`, creating parent directories."""
    out.parent.mkdir(parents=True, exist_ok=True)
    # Bytes, not text mode, so the newline convention does not vary by platform.
    out.write_bytes(render_summary(summary).encode("utf-8"))
    logger.info("wrote %s", out)
