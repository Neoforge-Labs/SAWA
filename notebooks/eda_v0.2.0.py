"""SAWA v0.2.0 — Comprehensive exploratory data-analysis report.

All-in-one Marimo notebook covering the SAWA v0.2.0 data vintage
(6 East African countries, ~42k HIV-1 pol sequences post 6-country
expansion). Reviewer-grade: every panel is sourced from a live
Parquet / JSON artifact and captions its source path + `n`.

Sections:

  A. Data provenance + filter waterfall
  B. Cohort composition (country, subtype, year, is_drug_naive)
  C. Sequence quality (length, ambiguity, APOBEC, per-gene coverage, QC)
  D. Resistance landscape — per-drug-class trends, country×year heatmap,
     INSTI/DTG spotlight card, top-10 mutations per drug class
  E. Phylogeny — manifest, per-subtype partition audit trail
     (`data/gold/trees/partition_report.parquet`), stale-manifest
     detector, per-tree quality panel (|logL|/tip, clock rate, wall time)
  F. Transmission clusters — size distribution, singleton banner if
     applicable, cross-border %, subtype composition, HIV-TRACE parity +
     synthetic-truth validation metrics
  G. Outbreak features — country × year feature matrix that feeds
     Paper 2's causal DAG (deposit_rate, subtype entropy)
  H. Pipeline performance — wall time + peak RSS from committed
     benchmarks
  I. Reproducibility — env, per-artifact schema_version + mtime table,
     tool versions from tree manifest, DVC pointer

Fig 4 DTG headline (Fisher exact + Clopper-Pearson + BH-FDR) lives in
its own notebook at `notebooks/dtg_deep_dive.py` — kept separate so it
can be rendered independently as the Fig 4 supplement.

Run:
    marimo edit notebooks/eda_v0.2.0.py
    marimo export html notebooks/eda_v0.2.0.py -o \
        paper1_tools/figures/eda_v0.2.0.html
"""

from __future__ import annotations

import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


# ------------------------------------------------------------------
# §0 imports + constants
# ------------------------------------------------------------------


@app.cell
def _():
    import json
    import logging
    import platform
    import subprocess
    import sys
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    logging.basicConfig(level="INFO", format="%(levelname)s %(name)s: %(message)s")

    REPO_ROOT = Path(__file__).resolve().parent.parent

    # SAWA convention across the codebase
    RESISTANT_LEVELS = frozenset({"high", "intermediate"})
    TARGET_COUNTRIES = ("KE", "UG", "TZ", "RW", "BI", "ET")

    # Okabe–Ito palette for country
    COUNTRY_COLOUR = {
        "KE": "#0072B2",
        "UG": "#D55E00",
        "TZ": "#009E73",
        "RW": "#CC79A7",
        "BI": "#F0E442",
        "ET": "#56B4E9",
    }
    return (
        COUNTRY_COLOUR,
        Path,
        REPO_ROOT,
        RESISTANT_LEVELS,
        TARGET_COUNTRIES,
        json,
        mo,
        np,
        platform,
        plt,
        pl,
        subprocess,
        sys,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
        # SAWA v0.2.0 — exploratory data-analysis report

        Companion to the Paper 1 manuscript. Reviewers, clinical users, and
        the maintainer can inspect every stage of the pipeline visually here.
        All numbers regenerate deterministically from the DVC-frozen data
        vintage at git tag `v0.2.0-sawa`.

        **Every figure caption names the source Parquet and the n it was
        computed against.** No claim in this notebook lives outside a
        committed data path.
        """
    )
    return


# ------------------------------------------------------------------
# §A Data provenance + filter waterfall
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §A Data provenance

        * Git SHA + tag
        * LANL / GenBank pull dates
        * Sierra algorithm version
        * Per-stage record counts (raw → curated → dr_results)
        * Filter waterfall with per-step drop counts and reasons
        """
    )
    return


@app.cell
def _(REPO_ROOT, mo, subprocess, sys):
    try:
        _sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
        _tag = subprocess.check_output(
            ["git", "describe", "--tags", "--always"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        _sha, _tag = "unknown", "unknown"

    mo.md(
        f"""
        | Field | Value |
        |---|---|
        | git SHA | `{_sha}` |
        | git tag | `{_tag}` |
        | Python | {sys.version.split()[0]} |
        """
    )
    return


@app.cell
def _(REPO_ROOT, mo, pl):
    _raw_dir = REPO_ROOT / "data" / "raw" / "lanl"
    _rows = []
    _total_raw = 0
    for _cd in sorted([p for p in _raw_dir.iterdir() if p.is_dir()]):
        _fasta = _cd / "hiv-db.fasta"
        if _fasta.exists():
            _n = sum(1 for line in _fasta.open() if line.startswith(">"))
            _rows.append({"country_code": _cd.name.upper(), "raw_fasta_records": _n})
            _total_raw += _n
    raw_counts = pl.DataFrame(_rows) if _rows else pl.DataFrame()

    mo.md(
        f"""
        ### Raw ingest — per-country LANL FASTA record counts

        Total raw records across all six countries: **{_total_raw:,}**
        """
    )
    raw_counts
    return (raw_counts,)


@app.cell
def _(REPO_ROOT, mo, pl):
    _interim = REPO_ROOT / "data" / "interim"
    ea = pl.read_parquet(_interim / "east_africa_pol.parquet")
    curated = pl.read_parquet(_interim / "curated" / "aligned.parquet")
    dr = pl.read_parquet(_interim / "dr_results.parquet")

    mo.md(
        f"""
        ### Interim stages

        | Stage | Rows | Unique accessions |
        |---|---|---|
        | `east_africa_pol.parquet` | {ea.height:,} | {ea["accession"].n_unique():,} |
        | `curated/aligned.parquet` | {curated.height:,} | {curated["accession"].n_unique():,} |
        | `dr_results.parquet` | {dr.height:,} | {dr["accession"].n_unique():,} |

        Sierra algorithm versions:
        {sorted(dr["algorithm_version"].unique().to_list())}
        """
    )
    return curated, dr, ea


@app.cell
def _(curated, dr, ea, mo, raw_counts):
    _raw_total = int(raw_counts["raw_fasta_records"].sum()) if not raw_counts.is_empty() else 0
    _ingest_n = ea["accession"].n_unique()
    _curated_n = curated["accession"].n_unique()
    _dr_n = dr["accession"].n_unique()

    _mkline = lambda a, b, label: (
        f"| {label} | {b:,} | −{max(a - b, 0):,} |"
        if a >= b
        else f"| {label} | {b:,} | +{b - a:,} |"
    )
    mo.md(
        f"""
        ### Filter waterfall

        | Stage | Accessions | Δ from previous |
        |---|---|---|
        | Raw FASTA (LANL) | {_raw_total:,} | — |
        {_mkline(_raw_total, _ingest_n, "Post-ingest (dedup, parse OK)")}
        {_mkline(_ingest_n, _curated_n, "Post-QC (align, ambiguity, APOBEC, stop codons)")}
        {_mkline(_curated_n, _dr_n, "Sierra-scored")}

        Common drop reasons at each step are logged during the pipeline
        run to `pipeline.log`; the totals above are the observable
        outcomes.
        """
    )
    return


# ------------------------------------------------------------------
# §B Cohort composition
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §B Cohort composition

        Countries, subtypes, years, and `is_drug_naive` completeness.
        The stacked subtype-per-country bar is the visual "non-B-native"
        story.
        """
    )
    return


@app.cell
def _(COUNTRY_COLOUR, curated, mo, plt):
    _fig, _ax = plt.subplots(figsize=(7, 3.5))
    _by_country = (
        curated.group_by("country").agg([pl.len().alias("n")]).sort("n", descending=True)
    )
    _ax.bar(
        _by_country["country"].to_list(),
        _by_country["n"].to_list(),
        color=[COUNTRY_COLOUR.get(c, "#999") for c in _by_country["country"].to_list()],
    )
    _ax.set_ylabel("Accessions (post-QC)")
    _ax.set_title(f"Records per country — curated/aligned.parquet (n={curated.height:,})")
    for _i, _v in enumerate(_by_country["n"].to_list()):
        _ax.text(_i, _v, f"{_v:,}", ha="center", va="bottom", fontsize=9)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(curated, plt):
    # Sampling-year distribution per country (small multiples)
    _countries = sorted(c for c in curated["country"].unique().to_list() if c is not None)
    _n_c = len(_countries)
    _fig, _axes = plt.subplots(1, _n_c, figsize=(2.2 * _n_c, 3), sharey=True)
    if _n_c == 1:
        _axes = [_axes]
    for _ax, _c in zip(_axes, _countries, strict=True):
        _years = (
            curated.filter(pl.col("country") == _c)
            .filter(pl.col("collection_year").is_not_null())["collection_year"]
            .to_list()
        )
        if _years:
            _ax.hist(
                _years, bins=range(min(_years), max(_years) + 2), color="#0072B2", alpha=0.85
            )
        else:
            _ax.text(0.5, 0.5, "no year data", ha="center", va="center", fontsize=8)
        _ax.set_title(f"{_c}  (n={len(_years):,})", fontsize=9)
        _ax.tick_params(labelsize=7)
    _axes[0].set_ylabel("Accessions")
    _fig.suptitle("Sampling-year distribution per country", fontsize=10)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(curated, plt):
    # Subtype composition per country (stacked bar, top 8 subtypes)
    _sub_col = "subtype" if "subtype" in curated.columns else None
    if _sub_col:
        _top = (
            curated.group_by(_sub_col)
            .agg(pl.len().alias("n"))
            .sort("n", descending=True)
            .head(8)[_sub_col]
            .to_list()
        )
        _cts = sorted(c for c in curated["country"].unique().to_list() if c is not None)
        _fig, _ax = plt.subplots(figsize=(7, 3.5))
        _bottom = [0.0] * len(_cts)
        _palette = [
            "#0072B2",
            "#D55E00",
            "#009E73",
            "#CC79A7",
            "#F0E442",
            "#56B4E9",
            "#E69F00",
            "#999999",
        ]
        for _s, _col in zip(_top, _palette, strict=False):
            _vals = []
            for _c in _cts:
                _n = curated.filter(
                    (pl.col("country") == _c) & (pl.col(_sub_col) == _s)
                ).height
                _tot = curated.filter(pl.col("country") == _c).height
                _vals.append(_n / _tot if _tot else 0.0)
            _ax.bar(_cts, _vals, bottom=_bottom, label=_s, color=_col)
            _bottom = [b + v for b, v in zip(_bottom, _vals, strict=True)]
        _ax.set_ylabel("Fraction of accessions")
        _ax.set_title("Top-8 subtype composition per country")
        _ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)
        _ax.set_ylim(0, 1.02)
        _fig.tight_layout()
    else:
        _fig = None
    _fig
    return


# ------------------------------------------------------------------
# §C Sequence quality
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §C Sequence quality

        Per-country distributions of length, ambiguity fraction, APOBEC
        signature ratio, per-gene coverage, and QC pass rate. Reviewers
        of a tools paper look at these before anything else.
        """
    )
    return


@app.cell
def _(curated, plt):
    _num_cols = [
        c
        for c in ("aligned_length", "ambiguity_fraction", "apobec_signature_ratio")
        if c in curated.columns
    ]
    if not _num_cols:
        _fig = None
    else:
        _fig, _axes = plt.subplots(1, len(_num_cols), figsize=(3.5 * len(_num_cols), 3))
        if len(_num_cols) == 1:
            _axes = [_axes]
        for _ax, _c in zip(_axes, _num_cols, strict=True):
            _vals = curated.filter(pl.col(_c).is_not_null())[_c].to_numpy()
            # Strip inf / NaN so matplotlib's hist doesn't reject the range —
            # apobec_signature_ratio in particular can produce inf on the
            # small handful of sequences with zero-count reference denominators.
            _vals = _vals[np.isfinite(_vals)]
            if _vals.size == 0:
                _ax.text(0.5, 0.5, "no finite values", ha="center", va="center")
            else:
                # Clip to the 1-99 percentile to keep the histogram legible
                # when a few outliers dominate the range.
                _lo, _hi = np.percentile(_vals, [1, 99])
                if _hi <= _lo:
                    _hi = _lo + 1e-9
                _ax.hist(_vals, bins=40, range=(_lo, _hi), color="#0072B2", alpha=0.85)
                _ax.set_xlim(_lo, _hi)
            _ax.set_title(_c, fontsize=9)
            _ax.tick_params(labelsize=7)
        _fig.suptitle("Sequence-quality distributions (post-QC)", fontsize=10)
        _fig.tight_layout()
    _fig
    return


@app.cell
def _(curated, mo, pl):
    _qc_cols = [c for c in curated.columns if c.startswith("has_") or c == "passes_qc"]
    _rows = []
    for _c in _qc_cols:
        _true = curated.filter(pl.col(_c) == True).height  # noqa: E712
        _rows.append(
            {"flag": _c, "n_true": _true, "pct": _true / max(curated.height, 1) * 100}
        )
    _qc_summary = pl.DataFrame(_rows)
    mo.md(f"QC flag prevalence across **{curated.height:,}** curated accessions:")
    _qc_summary
    return


# ------------------------------------------------------------------
# §D Resistance landscape
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §D Resistance landscape

        Deeper cuts than the manuscript Fig 2 heatmap:

        * Per-drug prevalence over time, faceted by country
        * Per-drug-class prevalence heatmap (country × year)
        * Top mutations per drug class (when the mutations column is present)
        """
    )
    return


@app.cell
def _(RESISTANT_LEVELS, curated, dr, plt):
    _drug_class_col = "drug_class"
    if _drug_class_col not in dr.columns:
        _fig = None
    else:
        _joined = dr.join(
            curated.select(["accession", "country", "collection_year"]),
            on="accession",
            how="inner",
        ).filter(pl.col("country").is_not_null() & pl.col("collection_year").is_not_null())
        _prev = (
            _joined.with_columns(
                pl.col("level").is_in(list(RESISTANT_LEVELS)).cast(pl.Int8).alias("R")
            )
            .group_by(["drug_class", "country", "collection_year"])
            .agg([pl.col("R").mean().alias("prev"), pl.len().alias("n")])
            .filter(pl.col("n") >= 10)
        )
        _classes = sorted(c for c in _prev["drug_class"].unique().to_list() if c is not None)
        _fig, _axes = plt.subplots(
            1, len(_classes), figsize=(3.5 * len(_classes), 3), sharey=True
        )
        if len(_classes) == 1:
            _axes = [_axes]
        _country_palette = {
            "KE": "#0072B2",
            "UG": "#D55E00",
            "TZ": "#009E73",
            "RW": "#CC79A7",
            "BI": "#F0E442",
            "ET": "#56B4E9",
        }
        for _ax, _cls in zip(_axes, _classes, strict=True):
            _sub = _prev.filter(pl.col("drug_class") == _cls).sort("collection_year")
            for _country in sorted(
                c for c in _sub["country"].unique().to_list() if c is not None
            ):
                _c_sub = _sub.filter(pl.col("country") == _country).sort("collection_year")
                _ax.plot(
                    _c_sub["collection_year"].to_list(),
                    _c_sub["prev"].to_list(),
                    marker="o",
                    color=_country_palette.get(_country, "#999"),
                    label=_country,
                    markersize=3,
                    linewidth=1,
                )
            _ax.set_title(_cls, fontsize=10)
            _ax.set_xlabel("Year")
            _ax.grid(alpha=0.3)
        _axes[0].set_ylabel("P(intermediate or high)")
        _axes[-1].legend(
            loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False
        )
        _fig.suptitle(
            "Per-drug-class resistance prevalence over time, by country", fontsize=11
        )
        _fig.tight_layout()
    _fig
    return


@app.cell
def _(RESISTANT_LEVELS, curated, dr, mo, np, pl, plt):
    if "drug_class" not in dr.columns:
        _fig = None
    else:
        _joined = dr.join(
            curated.select(["accession", "country", "collection_year"]),
            on="accession",
            how="inner",
        ).filter(
            # Drop rows with missing country or year before aggregation —
            # otherwise `sorted(unique)` fails on None comparisons.
            pl.col("country").is_not_null() & pl.col("collection_year").is_not_null()
        )
        _grid = (
            _joined.with_columns(
                pl.col("level").is_in(list(RESISTANT_LEVELS)).cast(pl.Int8).alias("R")
            )
            .group_by(["country", "collection_year"])
            .agg([pl.col("R").mean().alias("prev"), pl.len().alias("n")])
            .filter(pl.col("n") >= 20)
        )
        _cts = sorted(c for c in _grid["country"].unique().to_list() if c is not None)
        _years = sorted(
            y for y in _grid["collection_year"].unique().to_list() if y is not None
        )
        _mat = np.full((len(_cts), len(_years)), np.nan)
        for _r in _grid.iter_rows(named=True):
            _i = _cts.index(_r["country"])
            _j = _years.index(_r["collection_year"])
            _mat[_i, _j] = _r["prev"]
        _fig, _ax = plt.subplots(figsize=(9, 3))
        _im = _ax.imshow(_mat, cmap="YlOrRd", aspect="auto", vmin=0, vmax=0.3)
        _ax.set_yticks(range(len(_cts)))
        _ax.set_yticklabels(_cts)
        _ax.set_xticks(range(0, len(_years), max(len(_years) // 12, 1)))
        _ax.set_xticklabels(
            [str(_years[_i]) for _i in range(0, len(_years), max(len(_years) // 12, 1))],
            rotation=45,
            fontsize=8,
        )
        _ax.set_title(
            "Any-class resistance prevalence, country × year (cells n≥20)", fontsize=10
        )
        _cbar = _fig.colorbar(_im, ax=_ax, shrink=0.9)
        _cbar.set_label("P(intermediate or high)")
        _fig.tight_layout()

    mo.md("Grey cells: stratum had fewer than 20 accessions; not scored.")
    _fig
    return


@app.cell
def _(RESISTANT_LEVELS, curated, dr, mo, pl):
    # ---------- INSTI spotlight (feeds Paper 1 Fig 4 headline) ----------
    if "drug_class" not in dr.columns:
        _insti_summary = mo.md("> `drug_class` missing from dr_results.")
    else:
        _insti = dr.filter(pl.col("drug_class") == "INSTI").join(
            curated.select(["accession", "country", "collection_year"]),
            on="accession",
            how="inner",
        ).filter(
            pl.col("country").is_not_null() & pl.col("collection_year").is_not_null()
        )
        if _insti.height == 0:
            _insti_summary = mo.md("> No INSTI rows after join.")
        else:
            _insti = _insti.with_columns(
                pl.when(pl.col("collection_year") <= 2018)
                .then(pl.lit("pre (2010-2018)"))
                .when(pl.col("collection_year") >= 2020)
                .then(pl.lit("post (2020+)"))
                .otherwise(pl.lit("rollout (2019)"))
                .alias("epoch"),
                pl.col("level")
                .is_in(list(RESISTANT_LEVELS))
                .cast(pl.Int8)
                .alias("is_resistant"),
            )
            # Per-drug × epoch prevalence (excludes 2019 rollout buffer).
            _rollup = (
                _insti.filter(pl.col("epoch") != "rollout (2019)")
                .group_by(["drug", "epoch"])
                .agg(
                    pl.len().alias("n_calls"),
                    pl.col("is_resistant").sum().alias("n_resistant"),
                    (pl.col("is_resistant").mean() * 100).round(2).alias("prev_pct"),
                )
                .sort(["drug", "epoch"])
            )
            _insti_summary = mo.vstack([
                mo.md(
                    "### INSTI spotlight — pre vs post DTG rollout\n\n"
                    "Feeds Paper 1 Fig 4 headline. Excludes rollout year 2019 "
                    "as an implementation buffer. `n_calls` = number of "
                    "(sequence, drug) pairs Sierra scored; `n_resistant` = "
                    "calls at Sierra level ∈ {intermediate, high}."
                ),
                _rollup,
            ])
    _insti_summary
    return


@app.cell
def _(RESISTANT_LEVELS, dr, mo, pl):
    # ---------- Top mutations per drug class (co-occurring with R calls) ----------
    if "mutations" not in dr.columns:
        _top_muts_out = mo.md("> `mutations` column missing from dr_results.")
    else:
        _resistant = dr.filter(pl.col("level").is_in(list(RESISTANT_LEVELS)))
        if _resistant.height == 0:
            _top_muts_out = mo.md("> No resistant calls in dr_results.")
        else:
            # Deduplicate (accession, drug_class, mutation) so a single
            # sequence tested against 6 INSTI drugs doesn't count its
            # mutations 6× — we want mutation prevalence at the sequence
            # level, not the call level.
            _exploded = (
                _resistant.select(["accession", "drug_class", "mutations"])
                .explode("mutations")
                .drop_nulls("mutations")
                .unique(subset=["accession", "drug_class", "mutations"])
            )
            _top = (
                _exploded.group_by(["drug_class", "mutations"])
                .agg(pl.len().alias("n_sequences"))
                .sort(["drug_class", "n_sequences"], descending=[False, True])
            )
            _top_by_class = (
                _top.with_columns(
                    pl.col("n_sequences")
                    .rank(method="ordinal", descending=True)
                    .over("drug_class")
                    .alias("rank_in_class")
                )
                .filter(pl.col("rank_in_class") <= 10)
                .sort(["drug_class", "rank_in_class"])
                .rename({"mutations": "mutation"})
            )
            _top_muts_out = mo.vstack([
                mo.md(
                    "### Top-10 mutations per drug class (in resistant-calling sequences)\n\n"
                    "Ranked by number of distinct sequences carrying the mutation "
                    "among sequences whose Sierra call was intermediate or high "
                    "for at least one drug in that class. Includes accessory + "
                    "known DRAMs — Sierra's `mutations` column is the full "
                    "mutation profile, not a filtered DRAM list."
                ),
                _top_by_class,
            ])
    _top_muts_out
    return


# ------------------------------------------------------------------
# §E Phylogeny sanity checks
# ------------------------------------------------------------------


@app.cell
def _(REPO_ROOT, mo, pl):
    # Manifest lives under `data/gold/` per phylo/build.py convention;
    # fall back to `data/interim/` if a legacy vintage is present.
    _candidates = [
        REPO_ROOT / "data" / "gold" / "trees" / "manifest.parquet",
        REPO_ROOT / "data" / "interim" / "trees" / "manifest.parquet",
    ]
    _manifest_path = next((p for p in _candidates if p.exists()), None)
    _report_path = REPO_ROOT / "data" / "gold" / "trees" / "partition_report.parquet"

    _parts: list = []
    _parts.append(mo.md("## §E Phylogeny sanity checks"))

    if _manifest_path is None:
        _parts.append(
            mo.md(
                "> No tree manifest found. Expected at "
                "`data/gold/trees/manifest.parquet`. Run `make phylo`."
            )
        )
    else:
        _tree_manifest = pl.read_parquet(_manifest_path)
        _parts.append(mo.md("### Tree manifest (built artifacts)"))
        _parts.append(_tree_manifest)

    # Partition report — full per-subtype accounting including subtypes
    # that were NOT built (excluded, below thresholds). Emitted by
    # `hivdrwatch.phylo.build_all_trees` when `partition_report_path`
    # is passed; the Prefect flow does this by default.
    if _report_path.exists():
        _partition = pl.read_parquet(_report_path)
        _parts.append(
            mo.md(
                "### Partition audit trail\n\n"
                "One row per candidate subtype label seen in curated data — "
                "including labels that did NOT produce a tree. This is the "
                "reproducibility record: an auditor should be able to explain "
                "manifest membership from this table alone, without reading "
                "IQ-TREE logs."
            )
        )
        _parts.append(_partition)

        # Stale-tree detector: cross-check manifest n_taxa against what the
        # partition report says would be built from the current curated
        # data. If they diverge, the on-disk trees were built from an
        # earlier curated vintage — the manifest is provenance-mismatched.
        if _manifest_path is not None:
            _ml_only = _tree_manifest.filter(pl.col("kind") == "ml").select(
                ["subtype", "n_taxa"]
            )
            _expected = _partition.filter(pl.col("decision") == "built").select(
                pl.col("subtype"),
                pl.col("n_after_coverage").alias("expected_n_taxa"),
            )
            _joined = _ml_only.join(_expected, on="subtype", how="full").select(
                pl.coalesce([pl.col("subtype"), pl.col("subtype_right")]).alias(
                    "subtype"
                ),
                pl.col("n_taxa").alias("manifest_n_taxa"),
                pl.col("expected_n_taxa"),
            )
            _stale = _joined.filter(
                (pl.col("manifest_n_taxa").is_null())
                | (pl.col("expected_n_taxa").is_null())
                | (pl.col("manifest_n_taxa") != pl.col("expected_n_taxa"))
            )
            if _stale.height > 0:
                _parts.append(
                    mo.md(
                        "> ⚠️ **Stale manifest detected.** The on-disk trees "
                        "do not match what the current curated data would "
                        "produce. Re-run `make phylo` to rebuild. Discrepant "
                        "subtypes:"
                    )
                )
                _parts.append(_stale)
            else:
                _parts.append(
                    mo.md(
                        "> ✅ Manifest matches the current curated data — "
                        "trees are provenance-consistent."
                    )
                )
    else:
        _parts.append(
            mo.md(
                "> No partition report at "
                "`data/gold/trees/partition_report.parquet`. Re-run "
                "`make phylo` to emit the audit trail alongside the "
                "manifest — required for reviewer-grade reproducibility."
            )
        )

    mo.vstack(_parts)
    return


@app.cell
def _(REPO_ROOT, mo, pl, plt):
    # ---------- Per-tree quality panel (ML + time-scaled together) ----------
    _mp = REPO_ROOT / "data" / "gold" / "trees" / "manifest.parquet"
    if not _mp.exists():
        _fig = None
        _out = mo.md("> No manifest — skipping tree quality panel.")
    else:
        _m = pl.read_parquet(_mp)
        if _m.height == 0:
            _fig = None
            _out = mo.md("> Manifest empty — no tree quality to plot.")
        else:
            _ml = _m.filter(pl.col("kind") == "ml").sort("subtype")
            _tt = _m.filter(pl.col("kind") == "time_scaled").sort("subtype")
            _subtypes = _ml["subtype"].to_list()

            _fig, _axes = plt.subplots(1, 3, figsize=(11, 3.2))

            # Panel 1 — log-likelihood normalized by n_taxa (larger n has
            # larger |logL| by construction; per-tip logL is comparable).
            _per_tip_ll = [
                (r["log_likelihood"] / r["n_taxa"]) if r["log_likelihood"] else None
                for r in _ml.iter_rows(named=True)
            ]
            _axes[0].bar(
                _subtypes,
                [abs(v) if v is not None else 0 for v in _per_tip_ll],
                color="#0072B2",
                alpha=0.85,
            )
            _axes[0].set_title("|log-likelihood| per tip (ML)", fontsize=10)
            _axes[0].set_ylabel("|logL| / n_taxa")
            _axes[0].grid(axis="y", alpha=0.3)

            # Panel 2 — TreeTime clock rate ± stddev.
            _rates = [r["clock_rate"] for r in _tt.iter_rows(named=True)]
            _sd = [r["clock_rate_stddev"] or 0.0 for r in _tt.iter_rows(named=True)]
            _sub2 = _tt["subtype"].to_list()
            _axes[1].bar(
                _sub2,
                [r or 0.0 for r in _rates],
                yerr=_sd,
                color="#D55E00",
                alpha=0.85,
                capsize=4,
            )
            _axes[1].set_title("TreeTime clock rate ± SD", fontsize=10)
            _axes[1].set_ylabel("substitutions / site / year")
            _axes[1].grid(axis="y", alpha=0.3)

            # Panel 3 — wall seconds (ML vs time-scaled stacked).
            _ml_walls = [r["wall_seconds"] or 0.0 for r in _ml.iter_rows(named=True)]
            _tt_walls = [r["wall_seconds"] or 0.0 for r in _tt.iter_rows(named=True)]
            _x = list(range(len(_subtypes)))
            _axes[2].bar(_x, _ml_walls, color="#0072B2", label="IQ-TREE (ML)")
            _axes[2].bar(
                _x, _tt_walls, bottom=_ml_walls, color="#009E73", label="TreeTime"
            )
            _axes[2].set_xticks(_x)
            _axes[2].set_xticklabels(_subtypes)
            _axes[2].set_title("Wall-time per subtype", fontsize=10)
            _axes[2].set_ylabel("seconds")
            _axes[2].legend(fontsize=8, frameon=False)
            _axes[2].grid(axis="y", alpha=0.3)

            _fig.suptitle(
                "Per-tree quality — internal statistics from the manifest",
                fontsize=11,
            )
            _fig.tight_layout()
            _out = mo.md(
                "*|logL|/tip* is comparable across subtypes; large values indicate "
                "either sequence divergence or model misfit. Clock rates for HIV-1 "
                "pol typically fall in **1e-3 to 3e-3 substitutions/site/year** — "
                "flag anything outside that band."
            )
    _out
    _fig
    return


# ------------------------------------------------------------------
# §F Transmission clusters
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §F Transmission clusters

        Cluster size distribution, single- vs cross-border composition,
        subtype composition, year span. Populated when `data/gold/clusters.parquet`
        (or equivalent) exists — skipped cleanly otherwise.
        """
    )
    return


@app.cell
def _(REPO_ROOT, curated, mo, pl, plt):
    # ---------- Cluster size distribution + honest singleton banner ----------
    _cluster_paths = [
        REPO_ROOT / "data" / "gold" / "clusters.parquet",
        REPO_ROOT / "data" / "gold" / "cluster_assignments.parquet",
        REPO_ROOT / "data" / "interim" / "cluster_assignments.parquet",
    ]
    _cp = next((p for p in _cluster_paths if p.exists()), None)
    if _cp is None:
        _size_fig = None
        _size_out = mo.md("> No cluster assignments Parquet found — skipping §F.")
        _clusters = pl.DataFrame()
    else:
        _clusters = pl.read_parquet(_cp)
        _sizes = (
            _clusters.group_by("cluster_id")
            .agg(pl.len().alias("size"))
            .sort("size", descending=True)
        )
        _n_singletons = _sizes.filter(pl.col("size") == 1).height
        _n_multi = _sizes.filter(pl.col("size") > 1).height
        _threshold = (
            _clusters["tn93_threshold"].unique().to_list()[0]
            if "tn93_threshold" in _clusters.columns
            else None
        )
        _size_fig, _ax = plt.subplots(figsize=(7, 3))
        _ax.hist(_sizes["size"].to_list(), bins=40, color="#0072B2", alpha=0.85)
        _ax.set_yscale("log")
        _ax.set_xlabel("Cluster size")
        _ax.set_ylabel("Count (log)")
        _ax.set_title(
            f"Cluster size distribution — n_clusters={_sizes.height:,}, "
            f"n_assigned={_clusters.height:,}",
            fontsize=10,
        )
        _size_fig.tight_layout()

        if _n_multi == 0:
            _size_out = mo.md(
                f"> ⚠️ **All {_n_singletons:,} clusters are singletons at "
                f"TN93 threshold = {_threshold}.** No transmission links "
                "formed on the current cohort. Likely causes: (1) the "
                "full-cohort TN93 pair file at `data/interim/cluster/"
                "distances.csv` is empty — the tn93 stage may not have "
                "run over the full cohort; (2) 0.015 is too tight for the "
                "genetic diversity in the current pol alignment (predominantly "
                "single-gene submissions after coverage filtering). "
                "**Recommended:** verify the pipeline reached the tn93 stage "
                "on the current curated vintage, and consider a threshold "
                "sensitivity analysis at 0.015 / 0.02 / 0.03 before publishing "
                "any cluster-level claim."
            )
        else:
            _size_out = mo.md(
                f"n_singletons = {_n_singletons:,} · "
                f"n_multi_member_clusters = {_n_multi:,} · "
                f"tn93_threshold = {_threshold}"
            )
    _size_out
    _size_fig
    return


@app.cell
def _(REPO_ROOT, curated, mo, pl, plt):
    # ---------- Cluster subtype + country composition (multi-member only) ----------
    _cp = REPO_ROOT / "data" / "gold" / "cluster_assignments.parquet"
    if not _cp.exists():
        _comp_out = mo.md("> No cluster assignments — skipping composition panel.")
    else:
        _clusters = pl.read_parquet(_cp)
        _multi = _clusters.filter(pl.col("cluster_size") > 1)
        if _multi.height == 0:
            _comp_out = mo.md(
                "> No multi-member clusters — subtype / cross-border composition "
                "panel skipped. See warning above."
            )
        else:
            _joined = _multi.join(
                curated.select(["accession", "country", "subtype"]),
                on="accession",
                how="left",
            )
            # Cross-border: cluster spans ≥ 2 distinct countries.
            _cluster_countries = (
                _joined.group_by("cluster_id")
                .agg(pl.col("country").n_unique().alias("n_countries"))
            )
            _n_cross = _cluster_countries.filter(pl.col("n_countries") > 1).height
            _n_total = _cluster_countries.height
            # Subtype mix per cluster (dominant subtype).
            _subtype_mix = (
                _joined.group_by(["cluster_id", "subtype"])
                .agg(pl.len().alias("n"))
                .sort(["cluster_id", "n"], descending=[False, True])
                .group_by("cluster_id")
                .agg(pl.col("subtype").first().alias("dominant_subtype"))
            )
            _mix_summary = _subtype_mix.group_by("dominant_subtype").agg(
                pl.len().alias("n_clusters")
            ).sort("n_clusters", descending=True)
            _comp_out = mo.vstack([
                mo.md(
                    f"### Cluster geographic + subtype composition\n\n"
                    f"* multi-member clusters: **{_n_total:,}**\n"
                    f"* cross-border (≥ 2 countries): **{_n_cross:,}** "
                    f"({100 * _n_cross / _n_total:.1f}%)"
                ),
                _mix_summary,
            ])
    _comp_out
    return


@app.cell
def _(REPO_ROOT, json, mo):
    # ---------- HIV-TRACE parity + validation-harness metrics ----------
    _vp = REPO_ROOT / "data" / "gold" / "cluster_validation.json"
    if not _vp.exists():
        _val_out = mo.md("> No cluster_validation.json.")
    else:
        _v = json.loads(_vp.read_text())
        _clusterer = _v.get("clusterer", "unknown")
        _threshold = _v.get("threshold", "unknown")
        _sim = _v.get("sim", {})
        _parity = _v.get("parity", {})
        _metrics = _v.get("metrics", {})
        _md_lines = [
            "### Cluster validation — synthetic ground truth",
            "",
            f"* clusterer: `{_clusterer}` at threshold `{_threshold}`",
            f"* synthetic corpus: **n_sequences = {_sim.get('n_sequences', '?')}**, "
            f"true clusters = {_sim.get('n_clusters', '?')}, "
            f"singletons = {_sim.get('n_singletons', '?')} "
            f"(seed = {_sim.get('seed', '?')})",
            "",
            "**Primary clusterer metrics vs synthetic truth:**",
            "",
            f"* precision = **{_metrics.get('precision', 0):.3f}** · "
            f"recall = **{_metrics.get('recall', 0):.3f}** · "
            f"F1 = **{_metrics.get('f1', 0):.3f}** · "
            f"ARI = **{_metrics.get('adjusted_rand_index', 0):.3f}**",
            f"* n_predicted_clusters = {_metrics.get('n_predicted_clusters', '?')} · "
            f"n_truth_clusters = {_metrics.get('n_truth_clusters', '?')}",
            "",
            "**HIV-TRACE parity (pair-level agreement between our clusterer and "
            "HIV-TRACE's hivnetworkcsv on the same synthetic corpus):**",
            "",
            f"* n_pairs = {_parity.get('n_pairs', '?'):,} · "
            f"agreement = **{100 * _parity.get('fraction_agree', 0):.2f}%** "
            f"({_parity.get('n_agree', 0):,} / {_parity.get('n_pairs', 0):,})",
        ]
        _val_out = mo.md("\n".join(_md_lines))
    _val_out
    return


# ------------------------------------------------------------------
# §G Outbreak features (Paper 2 causal DAG input)
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §G Outbreak features — causal-DAG input

        The `outbreak_features` table is the country × year feature matrix
        that feeds Paper 2's causal DAG. Every row summarises one
        (country, year) cell with cluster-level deposits, subtype fractions,
        Shannon entropy, and normalized deposit rate. Reviewers reading
        Paper 2 will trace claims back to these numbers, so they belong in
        the EDA next to the raw cohort data.

        > **Note on DTG-era deep-dive:** Fig 4 (§G in the paper) has its
        > own dedicated notebook at `notebooks/dtg_deep_dive.py` — it renders
        > the Fisher / Clopper-Pearson / BH-FDR headline with per-stratum
        > forest plots.
        """
    )
    return


@app.cell
def _(COUNTRY_COLOUR, REPO_ROOT, mo, pl, plt):
    _of = REPO_ROOT / "data" / "gold" / "outbreak_features.parquet"
    if not _of.exists():
        _feat_out = mo.md(
            "> No `data/gold/outbreak_features.parquet` — run "
            "`make cluster` (or the cluster stage of the Prefect flow) "
            "to generate."
        )
        _feat_fig = None
        _peek = None
    else:
        _feat = pl.read_parquet(_of).filter(
            pl.col("country").is_not_null() & pl.col("year").is_not_null()
        )
        _peek = _feat.sort(["country", "year"]).head(20)

        # Two-panel plot: deposit_rate over time (per country) + mean subtype
        # entropy over time (per country).
        _feat_fig, _axes = plt.subplots(1, 2, figsize=(11, 3.2))
        for _country in sorted(
            c for c in _feat["country"].unique().to_list() if c is not None
        ):
            _sub = _feat.filter(pl.col("country") == _country).sort("year")
            _c = COUNTRY_COLOUR.get(_country, "#999")
            _axes[0].plot(
                _sub["year"].to_list(),
                _sub["deposit_rate"].to_list(),
                "o-",
                color=_c,
                label=_country,
                markersize=3,
                linewidth=1,
            )
            _axes[1].plot(
                _sub["year"].to_list(),
                _sub["mean_subtype_entropy"].to_list(),
                "o-",
                color=_c,
                label=_country,
                markersize=3,
                linewidth=1,
            )
        _axes[0].set_title("Deposit rate (normalized to country max)", fontsize=10)
        _axes[0].set_ylabel("deposit_rate")
        _axes[0].set_xlabel("Year")
        _axes[0].grid(alpha=0.3)
        _axes[1].set_title("Mean subtype Shannon entropy", fontsize=10)
        _axes[1].set_ylabel("H")
        _axes[1].set_xlabel("Year")
        _axes[1].grid(alpha=0.3)
        _axes[1].legend(
            loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False
        )
        _feat_fig.suptitle(
            "Outbreak features over time — inputs to the Paper 2 causal DAG",
            fontsize=11,
        )
        _feat_fig.tight_layout()
        _feat_out = mo.md(
            f"### First 20 (country × year) rows\n\n"
            f"Full table: **{_feat.height} rows** covering "
            f"{_feat['country'].n_unique()} countries × "
            f"{_feat['year'].n_unique()} distinct years."
        )
    _feat_out
    _peek
    _feat_fig
    return


# ------------------------------------------------------------------
# §H Pipeline performance
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §H Pipeline performance

        Wall time, throughput, and peak RSS per pipeline stage — from the
        committed benchmark JSONs. Regressions across releases are
        detectable by diffing successive `benchmarks/results/scale_*.json`.
        """
    )
    return


@app.cell
def _(REPO_ROOT, json, mo, plt):
    _bench_dir = REPO_ROOT / "benchmarks" / "results"
    _stages = {}
    for _name in ("scale_ingest", "scale_align", "scale_cluster"):
        _f = _bench_dir / f"{_name}.json"
        if _f.exists():
            _stages[_name] = json.loads(_f.read_text())
    if not _stages:
        _fig = None
        mo.md("> No benchmark JSONs in `benchmarks/results/`.")
    else:
        _fig, _axes = plt.subplots(1, 2, figsize=(9, 3.2))
        _colours = {
            "scale_ingest": "#0072B2",
            "scale_align": "#D55E00",
            "scale_cluster": "#009E73",
        }
        for _stage, _data in _stages.items():
            _xs = [m["n_sequences"] for m in _data["measurements"]]
            _walls = [m["wall_seconds"] for m in _data["measurements"]]
            _rss = [m["peak_mb"] for m in _data["measurements"]]
            _axes[0].loglog(_xs, _walls, "o-", color=_colours[_stage], label=_stage)
            _axes[1].semilogx(_xs, _rss, "o-", color=_colours[_stage], label=_stage)
        _axes[0].set_xlabel("n sequences")
        _axes[0].set_ylabel("wall seconds")
        _axes[0].set_title("Wall time vs n")
        _axes[0].grid(alpha=0.3, which="both")
        _axes[1].set_xlabel("n sequences")
        _axes[1].set_ylabel("peak RSS (MB)")
        _axes[1].set_title("Peak RSS vs n")
        _axes[1].grid(alpha=0.3, which="both")
        _axes[0].legend(fontsize=8)
        _fig.tight_layout()
    _fig
    return


# ------------------------------------------------------------------
# §I Reproducibility footer
# ------------------------------------------------------------------


@app.cell
def _(mo):
    mo.md("## §I Reproducibility footer")
    return


@app.cell
def _(REPO_ROOT, mo, platform, subprocess, sys):
    try:
        _sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        _sha = "unknown"
    try:
        _dirty = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=REPO_ROOT, text=True
        ).strip()
        _dirty_note = f"⚠️ working tree dirty ({len(_dirty.splitlines())} files)" if _dirty else "✅ clean working tree"
    except Exception:
        _dirty_note = "unknown"
    try:
        _tag = subprocess.check_output(
            ["git", "describe", "--tags", "--always"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        _tag = "unknown"
    try:
        import polars as _pl_mod
        import matplotlib as _mpl_mod

        _pkg_versions = f"polars={_pl_mod.__version__}, matplotlib={_mpl_mod.__version__}"
    except Exception:
        _pkg_versions = "unknown"

    mo.md(
        f"""
        ### Environment

        * **git SHA**: `{_sha}` — tag / describe: `{_tag}` — {_dirty_note}
        * **Python**: {sys.version.split()[0]} · **Platform**: {platform.platform()}
        * **Key packages**: {_pkg_versions}
        * **Regenerate this notebook**: `make eda`
        * **Related notebooks**: `notebooks/dtg_deep_dive.py` (Fig 4 supplement)
        """
    )
    return


@app.cell
def _(Path, REPO_ROOT, mo, pl):
    # ---------- Data-artifact provenance table ----------
    import datetime as _dt

    _artifacts = [
        ("raw ingest (LANL+GenBank merged)", "data/interim/east_africa_pol.parquet"),
        ("curated (aligned + QC + dedup)", "data/interim/curated/aligned.parquet"),
        ("Sierra subtype calls", "data/interim/sierra_subtypes.parquet"),
        ("Sierra DR calls", "data/interim/dr_results.parquet"),
        ("cluster assignments", "data/gold/cluster_assignments.parquet"),
        ("outbreak features (Paper 2 DAG input)", "data/gold/outbreak_features.parquet"),
        ("tree manifest", "data/gold/trees/manifest.parquet"),
        ("tree partition report", "data/gold/trees/partition_report.parquet"),
    ]
    _rows = []
    for _name, _relpath in _artifacts:
        _p = REPO_ROOT / _relpath
        if _p.exists():
            _stat = _p.stat()
            _mtime = _dt.datetime.fromtimestamp(_stat.st_mtime).isoformat(
                timespec="seconds"
            )
            # Pull schema_version from the parquet if the column exists.
            try:
                _sv_series = pl.read_parquet(_p, columns=["schema_version"])
                _sv = _sv_series["schema_version"].unique().to_list()
                _sv_str = ", ".join(str(s) for s in _sv if s is not None) or "—"
            except Exception:
                _sv_str = "—"
            _rows.append(
                {
                    "artifact": _name,
                    "path": _relpath,
                    "size_kb": round(_stat.st_size / 1024, 1),
                    "mtime": _mtime,
                    "schema_version": _sv_str,
                }
            )
        else:
            _rows.append(
                {
                    "artifact": _name,
                    "path": _relpath,
                    "size_kb": 0.0,
                    "mtime": "MISSING",
                    "schema_version": "—",
                }
            )
    _prov_table = pl.from_dicts(_rows)
    mo.vstack([
        mo.md(
            "### Data-artifact provenance\n\n"
            "Every downstream claim in this notebook (and in the manuscripts) "
            "traces back to one of these files. `schema_version` lets you "
            "verify the parser you're using matches the file you're reading — "
            "a version bump means the semantics changed. `mtime` is your "
            "quickest check for whether a stage has been re-run since the "
            "last upstream stage."
        ),
        _prov_table,
    ])
    return


@app.cell
def _(REPO_ROOT, mo, pl):
    # ---------- Tool versions from the tree manifest ----------
    _mp = REPO_ROOT / "data" / "gold" / "trees" / "manifest.parquet"
    if not _mp.exists():
        _tools_out = mo.md("> No tree manifest — tool-version block skipped.")
    else:
        _m = pl.read_parquet(_mp).select(
            ["subtype", "kind", "tool_version", "random_seed", "built_at"]
        )
        _tools_out = mo.vstack([
            mo.md(
                "### Tool versions (from tree manifest)\n\n"
                "Every tree row records the exact IQ-TREE / TreeTime version "
                "and RNG seed used. A reviewer wanting to reproduce a tree "
                "run needs only this row + the frozen curated Parquet."
            ),
            _m,
        ])
    _tools_out
    return


@app.cell
def _(REPO_ROOT, mo):
    # ---------- DVC pointer visibility ----------
    _artifacts = [
        ("raw ingest", "data/interim/east_africa_pol.parquet"),
        ("curated", "data/interim/curated/aligned.parquet"),
        ("Sierra DR", "data/interim/dr_results.parquet"),
        ("clusters", "data/gold/cluster_assignments.parquet"),
        ("outbreak features", "data/gold/outbreak_features.parquet"),
        ("trees dir", "data/interim/trees"),
    ]
    _tracked = []
    _untracked = []
    for _name, _relpath in _artifacts:
        _pointer = REPO_ROOT / f"{_relpath}.dvc"
        _target = REPO_ROOT / _relpath
        if _pointer.exists():
            _tracked.append((_name, _relpath, _pointer))
        elif _target.exists():
            _untracked.append((_name, _relpath))

    _lines = ["### DVC pointer status"]
    if _tracked:
        _lines.append("")
        _lines.append("**Tracked** (cite these hashes in the manuscript):")
        _lines.append("")
        for _name, _relpath, _pointer in _tracked:
            _txt = _pointer.read_text().strip()
            _lines.append(f"* `{_name}` → `{_relpath}.dvc`")
            _lines.append(f"    ```yaml")
            _lines.append(f"    {_txt}")
            _lines.append(f"    ```")
    if _untracked:
        _lines.append("")
        _lines.append(
            "> ⚠️ **Not yet DVC-tracked** — vintage cannot be cited by hash. "
            "Run **`make dvc-track`** (or `bash scripts/setup_dvc.sh`) to "
            "initialize DVC and register the following:"
        )
        _lines.append("")
        for _name, _relpath in _untracked:
            _lines.append(f"* `{_name}` — `{_relpath}`")
        _lines.append("")
        _lines.append(
            "After the script runs, commit the printed pointer files + "
            "`.dvc/config` + `.dvcignore`."
        )
    if not _tracked and not _untracked:
        _lines.append("")
        _lines.append(
            "> No DVC-trackable artifacts found on disk. Run the pipeline "
            "(`make pipeline`) to produce them, then `make dvc-track`."
        )
    _dvc_out = mo.md("\n".join(_lines))
    _dvc_out
    return


if __name__ == "__main__":
    app.run()
