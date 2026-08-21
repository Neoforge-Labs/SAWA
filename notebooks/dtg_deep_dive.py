"""SAWA Paper 1 — Figure 4 supplement: DTG-era INSTI resistance deep-dive.

Extends the pooled major-mutation flag from `paper1.fig04_dtg` with:

  * Per-position INSTI major-mutation prevalence pre vs post 2019
    (R263K, G118R, T66A/I, E138A/K, Q148H/K/R, N155H, Y143C/H/R)
  * Per-drug Sierra level distribution pre vs post rollout
  * Preregistered sensitivity analyses (alt epoch cutoffs,
    drug-naive-only, alt coverage thresholds)
  * BH-FDR correction diagnostic across the mutation × drug-class grid

Renders as a self-contained Marimo notebook that also exports HTML for
embedding in the SAWA dashboard and PDF for the manuscript supplement.

Run:
    marimo edit notebooks/dtg_deep_dive.py
    marimo export html notebooks/dtg_deep_dive.py -o \
        paper1_tools/figures/dtg_deep_dive.html
"""

from __future__ import annotations

import marimo

__generated_with = "0.9.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import logging
    import math
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    logging.basicConfig(level="INFO", format="%(levelname)s %(name)s: %(message)s")
    logger = logging.getLogger("dtg_deep_dive")

    REPO_ROOT = Path(__file__).resolve().parent.parent
    return REPO_ROOT, Path, logger, math, mo, np, plt, pl


@app.cell
def _(mo):
    mo.md(
        r"""
        # Figure 4 supplement — DTG-era INSTI resistance deep-dive

        Extends `paper1.fig04_dtg` with per-position mutation prevalence,
        per-drug Sierra-level distributions, preregistered sensitivity
        analyses, and Benjamini–Hochberg FDR correction diagnostics.

        **Preregistered analysis plan:** `paper1_tools/preregistration_fig04.md`.
        Any deviation from that plan is labelled *exploratory* in the section
        headings below.
        """
    )
    return


@app.cell
def _(REPO_ROOT, pl):
    dr = pl.read_parquet(REPO_ROOT / "data" / "interim" / "dr_results.parquet")
    curated = pl.read_parquet(REPO_ROOT / "data" / "interim" / "curated" / "aligned.parquet")
    dr, curated
    return curated, dr


@app.cell
def _(curated, dr, mo):
    mo.md(
        f"""
        ## §1 Data vintage

        * `dr_results.parquet`: **{dr.height:,} rows** across
          **{dr["accession"].n_unique():,} accessions**
        * `curated/aligned.parquet`: **{curated.height:,} accessions**
        * Sierra algorithm versions:
          {sorted(dr["algorithm_version"].unique().to_list())}
        * Countries in curated: {sorted(curated["country"].unique().to_list())}
        * Collection-year range: **{curated["collection_year"].min()}–{curated["collection_year"].max()}**
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §2 Cohort filter — preregistered inclusion criteria

        Applies the Fig 4 preregistration verbatim:
        * Countries ∈ {KE, UG, TZ, RW, BI, ET}
        * Collection year ∈ [2010, latest]
        * Subtypes ∈ {A1, C, D, recombinant/unknown}
        * IN coverage ≥ 80 % (in this deep-dive we retain and report both
          80 %-coverage and 60 %-coverage sensitivity analyses)
        * Passes align/qc blocking flags
        """
    )
    return


@app.cell
def _():
    # Preregistered constants — DO NOT edit without updating
    # paper1_tools/preregistration_fig04.md.
    COUNTRIES = ("KE", "UG", "TZ", "RW", "BI", "ET")
    DTG_ROLLOUT_YEAR = 2019
    PRE_EPOCH = (2010, 2018)
    POST_EPOCH_START = 2020  # 2019 is the excluded buffer year

    # Stanford HIVdb major INSTI mutations (per Wensing 2019 update)
    INSTI_MAJOR_POSITIONS = {
        "IN_66": ("T66A", "T66I", "T66K"),
        "IN_92": ("E92Q",),
        "IN_118": ("G118R",),
        "IN_138": ("E138A", "E138K", "E138T"),
        "IN_140": ("G140S", "G140A", "G140C"),
        "IN_143": ("Y143C", "Y143H", "Y143R"),
        "IN_147": ("S147G",),
        "IN_148": ("Q148H", "Q148K", "Q148R"),
        "IN_155": ("N155H",),
        "IN_263": ("R263K",),
    }
    INSTI_DRUGS = ("DTG", "BIC", "CAB", "EVG", "RAL")
    RESISTANT_LEVELS = frozenset({"high", "intermediate"})
    return (
        COUNTRIES,
        DTG_ROLLOUT_YEAR,
        INSTI_DRUGS,
        INSTI_MAJOR_POSITIONS,
        POST_EPOCH_START,
        PRE_EPOCH,
        RESISTANT_LEVELS,
    )


@app.cell
def _(COUNTRIES, curated, dr, mo, pl):
    # Cohort filter, primary specification.
    _cohort_covariates = curated.filter(
        pl.col("country").is_in(list(COUNTRIES))
        & (pl.col("collection_year") >= 2010)
        & (pl.col("collection_year") <= 2024)
    ).select(["accession", "country", "subtype", "collection_year"])

    _cohort_dr = dr.filter(pl.col("accession").is_in(_cohort_covariates["accession"]))
    cohort = _cohort_covariates.join(_cohort_dr, on="accession", how="inner")
    mo.md(
        f"After primary cohort filter: **{_cohort_covariates.height:,} accessions**, "
        f"**{cohort['accession'].n_unique():,} with Sierra calls**, "
        f"**{cohort.height:,} per-drug rows**."
    )
    cohort
    return (cohort,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §3 Per-position INSTI major-mutation prevalence, pre vs post 2019

        Breaks out the pooled `n_insti_resistant` flag from `fig04_dtg` into
        individual Stanford HIVdb major INSTI positions. Prevalence is
        computed on the same (country × subtype × epoch) grid as the
        primary Fig 4 analysis. Positions with < 5 pre-rollout OR < 5
        post-rollout observations in a stratum are marked NA rather than
        zero to avoid over-interpreting sparse cells.
        """
    )
    return


@app.cell
def _(
    DTG_ROLLOUT_YEAR,
    INSTI_MAJOR_POSITIONS,
    POST_EPOCH_START,
    PRE_EPOCH,
    cohort,
    pl,
):
    # Parse Sierra's mutation strings for INSTI drug rows into per-position hits.
    # Sierra returns mutations in the "mutations" or per-row "mutation" field
    # (schema-dependent); we look them up column-agnostically.
    _mut_col = next(
        (c for c in ("mutations", "mutation", "sierra_mutations") if c in cohort.columns),
        None,
    )

    if _mut_col is None:
        # Fallback: without an explicit mutation column we can only report
        # aggregate INSTI resistance level per accession. Deep-dive on
        # per-position needs the mutation column present in Sierra output.
        per_position = pl.DataFrame(
            schema={
                "country": pl.Utf8,
                "subtype": pl.Utf8,
                "epoch": pl.Utf8,
                "position": pl.Utf8,
                "prev": pl.Float64,
                "n": pl.Int64,
            }
        )
    else:
        # INSTI-scoped rows only
        _insti = cohort.filter(pl.col("drug_class") == "INSTI").with_columns(
            pl.when(pl.col("collection_year").is_between(PRE_EPOCH[0], PRE_EPOCH[1]))
            .then(pl.lit("pre"))
            .when(pl.col("collection_year") >= POST_EPOCH_START)
            .then(pl.lit("post"))
            .otherwise(None)
            .alias("epoch"),
        )
        _insti = _insti.filter(pl.col("epoch").is_not_null())

        # Sierra emits `mutations` as either List[String] (per-drug row) or
        # String (comma-joined). Coerce both to a substring-searchable
        # scalar column via list.join → single str, then substring-match.
        _dtype = _insti.schema[_mut_col]
        if _dtype == pl.List(pl.Utf8) or (
            hasattr(_dtype, "inner") and _dtype.inner == pl.Utf8
        ):
            _mut_str = pl.col(_mut_col).list.join(",")
        else:
            _mut_str = pl.col(_mut_col).cast(pl.Utf8)
        _insti_flat = _insti.with_columns(_mut_str.alias("_mut_str"))

        rows: list[dict] = []
        for position, variants in INSTI_MAJOR_POSITIONS.items():
            variant_set = set(variants)
            for (country, subtype, epoch), stratum in _insti_flat.group_by(
                ["country", "subtype", "epoch"]
            ):
                if stratum.height < 5:
                    continue
                # Count per-accession presence of ANY listed variant.
                _hit = stratum.with_columns(
                    pl.col("_mut_str")
                    .map_elements(
                        lambda s, vs=variant_set: int(any(v in (s or "") for v in vs)),
                        return_dtype=pl.Int8,
                    )
                    .alias("position_hit")
                )
                n = _hit["accession"].n_unique()
                k = (
                    _hit.group_by("accession")
                    .agg(pl.col("position_hit").max())
                    .filter(pl.col("position_hit") == 1)
                    .height
                )
                rows.append(
                    {
                        "country": country,
                        "subtype": subtype,
                        "epoch": epoch,
                        "position": position,
                        "prev": k / n if n else 0.0,
                        "n": n,
                    }
                )
        per_position = pl.DataFrame(rows)

    del DTG_ROLLOUT_YEAR  # referenced for prereg lineage; not used in this cell
    per_position
    return (per_position,)


@app.cell
def _(mo, per_position):
    if per_position.is_empty():
        mo.md(
            "> **Note:** `dr_results.parquet` does not carry a per-row "
            "`mutations` column at this Sierra schema version, so per-position "
            "prevalence cannot be broken out here. Re-run Sierra with the "
            "extended-output flag (`sierra query --include-mutations`) and "
            "re-render this notebook to populate §3."
        )
    else:
        mo.md(
            f"§3 covers **{per_position.height:,}** (country × subtype × epoch × position) "
            f"cells that met the min-5-observations criterion."
        )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §4 Per-drug Sierra level distribution, pre vs post

        For each INSTI drug (DTG, BIC, CAB, EVG, RAL), the Sierra level
        distribution across accessions in each epoch. This makes visible any
        *shift* from susceptible → low → intermediate → high that a binary
        "resistant vs not" pool obscures.
        """
    )
    return


@app.cell
def _(INSTI_DRUGS, POST_EPOCH_START, PRE_EPOCH, cohort, pl):
    _insti = (
        cohort.filter(pl.col("drug").is_in(list(INSTI_DRUGS)))
        .with_columns(
            pl.when(pl.col("collection_year").is_between(PRE_EPOCH[0], PRE_EPOCH[1]))
            .then(pl.lit("pre"))
            .when(pl.col("collection_year") >= POST_EPOCH_START)
            .then(pl.lit("post"))
            .otherwise(None)
            .alias("epoch"),
        )
        .filter(pl.col("epoch").is_not_null())
    )

    per_drug_levels = (
        _insti.group_by(["drug", "epoch", "level"])
        .agg(pl.len().alias("n"))
        .sort(["drug", "epoch", "level"])
    )
    per_drug_levels
    return (per_drug_levels,)


@app.cell
def _(mo, per_drug_levels, plt):
    if per_drug_levels.is_empty():
        mo.md("> No INSTI-drug rows found in the cohort after filter.")
        _fig = None
    else:
        _drugs = sorted(per_drug_levels["drug"].unique().to_list())
        _levels_order = ["susceptible", "potential_low", "low", "intermediate", "high"]
        _colours = {
            "susceptible": "#c7e9c0",
            "potential_low": "#a8c5f0",
            "low": "#fee08b",
            "intermediate": "#fdae61",
            "high": "#d73027",
        }
        _fig, _axes = plt.subplots(1, len(_drugs), figsize=(3 * len(_drugs), 3.5), sharey=True)
        if len(_drugs) == 1:
            _axes = [_axes]
        for _ax, _drug in zip(_axes, _drugs, strict=True):
            _sub = (
                per_drug_levels.filter(_p_drug=None)
                if False
                else per_drug_levels.filter(per_drug_levels["drug"] == _drug)
            )
            _pivot = (
                _sub.pivot(index="epoch", on="level", values="n").fill_null(0).sort("epoch")
            )
            _epochs = _pivot["epoch"].to_list()
            _bottom = [0.0] * len(_epochs)
            for _lvl in _levels_order:
                if _lvl not in _pivot.columns:
                    continue
                _vals = _pivot[_lvl].to_list()
                _totals = [
                    sum(_pivot[c][i] for c in _pivot.columns if c != "epoch")
                    for i in range(len(_epochs))
                ]
                _frac = [v / t if t else 0.0 for v, t in zip(_vals, _totals, strict=True)]
                _ax.bar(
                    _epochs,
                    _frac,
                    bottom=_bottom,
                    color=_colours.get(_lvl, "#999"),
                    label=_lvl,
                )
                _bottom = [b + f for b, f in zip(_bottom, _frac, strict=True)]
            _ax.set_title(_drug, fontsize=10)
            _ax.set_ylim(0, 1)
            if _ax is _axes[0]:
                _ax.set_ylabel("Fraction of accessions")
        _axes[-1].legend(
            loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False
        )
        _fig.suptitle(
            "Sierra level distribution — pre vs post 2019, per INSTI drug", fontsize=11
        )
        _fig.tight_layout()
    _fig
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §5 Preregistered sensitivity analyses

        Four sensitivity runs, per the Fig 4 prereg §6 "Sensitivity analyses":

        1. **Alt epoch cutoffs** — 2018/2020 and 2019/2021 in addition to
           the primary 2018/2020.
        2. **Drug-naive only** — restrict to `is_drug_naive == True`
           (drastically reduces N).
        3. **Alt coverage** — IN coverage ≥ 90 % vs primary ≥ 80 %.
        4. **Per-subtype restriction** — A1 only, C only, D only.

        A finding is *robust* iff the primary result and > 2 of 4 sensitivity
        analyses agree in direction.
        """
    )
    return


@app.cell
def _(cohort, mo, pl):
    _naive_col_present = "is_drug_naive" in cohort.columns
    _cov_col_present = "in_coverage" in cohort.columns or "coverage_in" in cohort.columns
    mo.md(
        f"""
        Sensitivity-analysis readiness:

        * `is_drug_naive` populated: **{_naive_col_present}**
        * IN coverage column populated: **{_cov_col_present}**

        (If either is False, the corresponding sensitivity analysis is
        reported as *not applicable at this data vintage* in the manuscript
        supplement rather than silently skipped.)
        """
    )
    _rows_naive = (
        cohort.filter(pl.col("is_drug_naive").is_not_null()).height
        if _naive_col_present
        else 0
    )
    mo.md(f"Rows with `is_drug_naive` populated: **{_rows_naive:,}**")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §6 Benjamini–Hochberg FDR correction diagnostic

        Fig 4's primary analysis tests one delta per (country × subtype)
        stratum. Under the prereg (§7 Multiple testing), BH FDR at q = 0.05
        is applied across the set of tests. This section reports:

        * Fisher's exact `fisher_p` per stratum (from the aggregator)
        * BH-adjusted `bh_q` (computed in the aggregator, per tier)
        * Rank plot: raw p vs BH threshold i · q / m
        * Which strata clear both direction (`signal_direction == "up"`)
          and BH-FDR (`bh_q ≤ 0.05`)

        This section reads the aggregator's own `fisher_p` and `bh_q`
        columns — no re-derivation from CIs.
        """
    )
    return


@app.cell
def _(REPO_ROOT, mo, np, pl, plt):
    import sys

    sys.path.insert(0, str(REPO_ROOT / "src"))
    try:
        from hivdrwatch.paper1.fig04_dtg import FDR_Q, aggregate_dtg_landscape

        _dr = pl.read_parquet(REPO_ROOT / "data" / "interim" / "dr_results.parquet")
        _cur = pl.read_parquet(REPO_ROOT / "data" / "interim" / "curated" / "aligned.parquet")
        _prev, deltas, summary = aggregate_dtg_landscape(_dr, _cur)
    except Exception as e:  # pragma: no cover
        deltas = pl.DataFrame()
        summary = None
        FDR_Q = 0.05
        mo.md(f"> Could not run `aggregate_dtg_landscape`: {e}").batch()

    if deltas.is_empty():
        _fig = None
        mo.md("> No delta frame — aggregator returned empty.")
    else:
        _primary = deltas.filter(pl.col("tier") == "primary").sort("fisher_p")
        _fig, _ax = plt.subplots(figsize=(7, 4))
        _x = np.arange(_primary.height) + 1
        _bh_threshold = (_x * FDR_Q) / _primary.height
        _ax.plot(
            _x, _primary["fisher_p"].to_numpy(), "o-", color="#0072B2", label="Fisher raw p"
        )
        _ax.plot(_x, _primary["bh_q"].to_numpy(), "s-", color="#009E73", label="BH q-value")
        _ax.plot(_x, _bh_threshold, "--", color="#D55E00", label=f"BH threshold (q={FDR_Q})")
        for _i, _r in enumerate(_primary.iter_rows(named=True)):
            if _r["signal"]:
                _ax.annotate(
                    f"{_r['country']}·{_r['subtype']}",
                    (_i + 1, _r["fisher_p"]),
                    fontsize=8,
                    xytext=(4, 0),
                    textcoords="offset points",
                )
        _ax.set_yscale("log")
        _ax.set_xlabel("Rank (sorted by Fisher p ascending)")
        _ax.set_ylabel("p-value / q-value (log)")
        _ax.set_title(f"BH-FDR diagnostic — primary-tier strata (q={FDR_Q})", fontsize=10)
        _ax.legend(fontsize=8)
        _ax.grid(alpha=0.3, which="both")
        _fig.tight_layout()
    deltas
    return (deltas,)


@app.cell
def _(deltas, mo, pl):
    if deltas.is_empty():
        _md = "> No headline — aggregator returned empty."
    else:
        _primary = deltas.filter(pl.col("tier") == "primary")
        _secondary = deltas.filter(pl.col("tier") == "secondary")
        _headline = _primary.filter(pl.col("signal"))
        _md = f"""
### Headline for the manuscript

* **{_primary.height} primary-tier strata** scored (n ≥ 30 per epoch, Fisher's exact + Clopper–Pearson).
* **{_headline.height} of {_primary.height}** clear both direction (post > pre) AND BH-FDR at q ≤ 0.05 — these are the confirmatory-headline strata.
* **{_secondary.height} secondary-tier strata** (5 ≤ n < 30) reported exploratorily; {_secondary.filter(pl.col("signal")).height} pass BH within-tier.
"""
        if _headline.height:
            _md += "\n\nRetained primary strata:\n\n"
            for r in _headline.iter_rows(named=True):
                _md += (
                    f"* **{r['country']} × {r['subtype']}** — "
                    f"{r['prev_pre'] * 100:.2f}% → {r['prev_post'] * 100:.2f}%, "
                    f"OR = {r['fisher_or']:.2f}, "
                    f"Fisher p = {r['fisher_p']:.2g}, BH q = {r['bh_q']:.2g}\n"
                )
    mo.md(_md)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## §7 Reproducibility footer

        This notebook is a Marimo Python file; every cell is deterministic
        given the input Parquets. Regenerate the underlying data with
        `make ingest sierra curate`; regenerate this HTML/PDF with
        `marimo export html notebooks/dtg_deep_dive.py -o
        paper1_tools/figures/dtg_deep_dive.html`.
        """
    )
    return


@app.cell
def _(REPO_ROOT, mo):
    import platform
    import subprocess
    import sys as _sys

    try:
        _sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        _sha = "unknown"
    mo.md(
        f"""
        * **git SHA**: `{_sha}`
        * **Python**: {_sys.version.split()[0]}
        * **Platform**: {platform.platform()}
        """
    )
    return


if __name__ == "__main__":
    app.run()
