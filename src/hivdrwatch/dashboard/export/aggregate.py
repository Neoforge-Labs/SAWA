"""Aggregated views of deduplicated Sierra results for the resistance dashboard.

Each function takes a frame already passed through `dedup_dr_results` (one row
per (accession, drug)) and returns a small summary frame. Row order is fixed by
an explicit trailing sort on every function, because Polars' `group_by` and
`unique` return rows in a different order on every call and the frames are
serialised into a committed `summary.json` that must not churn between runs.
Drug classes sort by the manuscript's `DRUG_CLASS_ORDER`, not alphabetically.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable

import polars as pl

from hivdrwatch.dashboard.export.drm import SdrmPanel, match_sdrm
from hivdrwatch.paper1.fig02_aggregate import RESISTANT_LEVELS
from hivdrwatch.paper1.style import DRUG_CLASS_ORDER

logger = logging.getLogger(__name__)

# Least to most severe, so stacked bars and legends read in a fixed order.
LEVEL_ORDER: tuple[str, ...] = ("susceptible", "potential_low", "low", "intermediate", "high")

_CLASS_RANK: dict[str, int] = {name: i for i, name in enumerate(DRUG_CLASS_ORDER)}
_LEVEL_RANK: dict[str, int] = {name: i for i, name in enumerate(LEVEL_ORDER)}

_SDRM_SCHEMA = pl.Schema(
    {
        "panel": pl.String,
        "mutation": pl.String,
        "n_sequences": pl.UInt32,
        "n_panel_denominator": pl.UInt32,
        "frac": pl.Float64,
    }
)


class UnknownDrugClassError(ValueError):
    """A drug class outside `DRUG_CLASS_ORDER`, which has no canonical position."""


class UnknownLevelError(ValueError):
    """A Sierra level outside `LEVEL_ORDER`, which would silently vanish from the breakdown."""


def _require(df: pl.DataFrame, columns: Iterable[str], *, view: str) -> None:
    missing = set(columns) - set(df.columns)
    if missing:
        raise ValueError(f"{view}: dr_results missing required columns: {sorted(missing)}")


def _check_drug_classes(df: pl.DataFrame, *, view: str) -> None:
    unknown = sorted(set(df["drug_class"].drop_nulls().unique()) - set(DRUG_CLASS_ORDER))
    if unknown:
        raise UnknownDrugClassError(
            f"{view}: drug classes {unknown} are not in the canonical order {DRUG_CLASS_ORDER}"
        )


def _class_rank() -> pl.Expr:
    return pl.col("drug_class").replace_strict(_CLASS_RANK, return_dtype=pl.UInt8)


def prevalence_by_drug_class(df: pl.DataFrame) -> pl.DataFrame:
    """Sequences scored per drug class and how many are resistant to any drug in it.

    Output: [drug_class, n_sequences, n_resistant, frac_resistant]. The
    denominator is per class because most deposited sequences cover PR and RT
    but not integrase. A class with no scored sequence has no row.
    """
    _require(df, ("accession", "drug_class", "level"), view="prevalence")
    _check_drug_classes(df, view="prevalence")

    per_sequence = df.group_by(["accession", "drug_class"]).agg(
        pl.col("level").is_in(list(RESISTANT_LEVELS)).any().alias("is_resistant")
    )
    out = (
        per_sequence.group_by("drug_class")
        .agg(
            n_sequences=pl.len(),
            n_resistant=pl.col("is_resistant").sum(),
        )
        .with_columns(frac_resistant=pl.col("n_resistant") / pl.col("n_sequences"))
        .sort(_class_rank())
        .select("drug_class", "n_sequences", "n_resistant", "frac_resistant")
    )
    logger.info(
        "prevalence: %d drug classes over %d accessions",
        out.height,
        df["accession"].n_unique(),
    )
    return out


def levels_by_drug(df: pl.DataFrame) -> pl.DataFrame:
    """Per-drug breakdown across every Sierra level, zero-filled.

    Output: [drug_class, drug, level, n_sequences, frac]. Every drug carries a
    row for each of `LEVEL_ORDER` so a renderer can draw stacked bars without
    guessing which levels a drug never received.
    """
    _require(df, ("accession", "drug_class", "drug", "level"), view="levels")
    _check_drug_classes(df, view="levels")
    unknown = sorted(set(df["level"].drop_nulls().unique()) - set(LEVEL_ORDER))
    if unknown or df["level"].null_count():
        raise UnknownLevelError(
            f"levels: found levels {unknown} and {df['level'].null_count()} null levels; "
            f"expected only {LEVEL_ORDER}"
        )

    drugs = df.select("drug_class", "drug").unique()
    grid = drugs.join(pl.DataFrame({"level": list(LEVEL_ORDER)}), how="cross")
    counts = df.group_by(["drug_class", "drug", "level"]).agg(
        n_sequences=pl.col("accession").n_unique()
    )
    out = (
        grid.join(counts, on=["drug_class", "drug", "level"], how="left")
        .with_columns(pl.col("n_sequences").fill_null(0).cast(pl.UInt32))
        .with_columns(
            frac=pl.col("n_sequences")
            / pl.col("n_sequences").sum().over(["drug_class", "drug"])
        )
        .sort(
            _class_rank(),
            "drug",
            pl.col("level").replace_strict(_LEVEL_RANK, return_dtype=pl.UInt8),
        )
        .select("drug_class", "drug", "level", "n_sequences", "frac")
    )
    logger.info("levels: %d drugs, %d rows", drugs.height, out.height)
    return out


def gene_coverage(df: pl.DataFrame) -> pl.DataFrame:
    """Distinct accessions scored per gene, as [gene, n_sequences], sorted by gene name."""
    _require(df, ("accession", "gene"), view="coverage")
    out = (
        df.group_by("gene")
        .agg(n_sequences=pl.col("accession").n_unique())
        .sort("gene")
        .select("gene", "n_sequences")
    )
    logger.info("coverage: %s", dict(out.iter_rows()))
    return out


def _panel_frequency(df: pl.DataFrame, panel: SdrmPanel) -> pl.DataFrame:
    # `mutations` repeats on every drug row of a gene, so reduce to one row per
    # accession before exploding or each mutation is counted once per drug.
    per_accession = df.filter(pl.col("gene") == panel.gene).unique(
        subset="accession", keep="first", maintain_order=True
    )
    denominator = per_accession.height
    if denominator == 0:
        logger.warning("panel %s: no sequences scored for gene %s", panel.name, panel.gene)
        return pl.DataFrame(schema=_SDRM_SCHEMA)

    exploded = (
        per_accession.select("accession", pl.col("mutations").alias("raw"))
        .explode("raw", empty_as_null=False)
        .drop_nulls("raw")
    )
    labels = pl.DataFrame(
        {"raw": (raws := sorted(exploded["raw"].unique()))},
        schema={"raw": pl.String},
    ).with_columns(
        mutation=pl.Series([match_sdrm(raw, panel) for raw in raws], dtype=pl.String)
    )
    return (
        exploded.join(labels, on="raw")
        .drop_nulls("mutation")
        .group_by("mutation")
        .agg(n_sequences=pl.col("accession").n_unique())
        .with_columns(
            panel=pl.lit(panel.name),
            n_panel_denominator=pl.lit(denominator, dtype=pl.UInt32),
            frac=pl.col("n_sequences") / denominator,
        )
        .select(_SDRM_SCHEMA.names())
        .cast(_SDRM_SCHEMA)
    )


def sdrm_frequency(df: pl.DataFrame, panels: tuple[SdrmPanel, ...]) -> pl.DataFrame:
    """Surveillance-DRM frequency per panel.

    Output: [panel, mutation, n_sequences, n_panel_denominator, frac]. A panel
    sees only rows of its own gene, and its denominator is every accession
    scored for that gene, including those carrying no mutations. Both RT panels
    therefore share one denominator. Output is sorted by panel (in the order
    given) and then mutation label.
    """
    _require(df, ("accession", "gene", "mutations"), view="sdrm")
    frames = [_panel_frequency(df, panel) for panel in panels]
    rank = {panel.name: i for i, panel in enumerate(panels)}
    out = (
        pl.concat(frames)
        .sort(pl.col("panel").replace_strict(rank, return_dtype=pl.UInt8), "mutation")
        .select(_SDRM_SCHEMA.names())
    )
    logger.info("sdrm: %d (panel, mutation) rows across %d panels", out.height, len(panels))
    return out
