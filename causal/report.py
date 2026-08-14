"""Render a human-readable HTML report of causal effect estimates.

Kept intentionally simple — one HTML file with a table per data scope.
The point is to have something MLflow-loggable that a reviewer can open
without setting up a notebook.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .schema import CausalEffect


def effects_to_frame(effects: list[CausalEffect]) -> pl.DataFrame:
    """Flatten a list of CausalEffect into a wide parquet-friendly frame."""
    if not effects:
        return pl.DataFrame()
    rows = []
    for e in effects:
        row = {
            "treatment": e.treatment,
            "outcome": e.outcome,
            "adjustment_set": ",".join(e.adjustment_set),
            "data_scope": e.data_scope,
            "n": e.n_observations,
            "method": e.method.value,
            "estimate": e.estimate,
            "ci_low": e.ci_low,
            "ci_high": e.ci_high,
            "p_value": e.p_value,
            "identified": e.identified,
            "criterion": e.identification_criterion,
            "dag_version": e.dag_version,
            "estimated_at": e.estimated_at.isoformat(),
        }
        # Refutation summary: pass/fail/inconclusive counts.
        from collections import Counter

        counter = Counter(r.verdict.value for r in e.refutations)
        row["ref_pass"] = counter.get("passes", 0)
        row["ref_fail"] = counter.get("fails", 0)
        row["ref_inconclusive"] = counter.get("inconclusive", 0)
        rows.append(row)
    return pl.from_dicts(rows)


def write_html_report(effects: list[CausalEffect], out_path: Path) -> Path:
    """Emit a single self-contained HTML page. No JS, no external CSS."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame = effects_to_frame(effects)

    scopes = frame["data_scope"].unique().to_list() if frame.height else []

    parts: list[str] = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>HIV-DR-Watch causal effects</title>",
        "<style>",
        "body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;color:#222}",
        "table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:0.9rem}",
        "th,td{padding:0.35rem 0.6rem;border-bottom:1px solid #ddd;text-align:left}",
        "th{background:#f6f6f6;font-weight:600}",
        ".pass{color:#166534}.fail{color:#991b1b}.inc{color:#6b7280}",
        ".unident{opacity:0.55}",
        "code{font-size:0.85em;background:#f6f6f6;padding:0.1em 0.3em;border-radius:3px}",
        "h2{margin-top:2.5rem;border-bottom:2px solid #eee;padding-bottom:0.3rem}",
        "</style></head><body>",
        "<h1>Causal effect estimates</h1>",
        "<p>Every row is a backdoor-adjusted causal effect estimate against the "
        "committed domain DAG. Refutation columns show how many placebo / "
        "random-common-cause / subset challenges each estimate passed.</p>",
    ]

    for scope in scopes:
        parts.append(f"<h2>{scope}</h2>")
        sub = frame.filter(pl.col("data_scope") == scope).sort(
            "identified", "estimate", descending=[True, True]
        )
        parts.append("<table><thead><tr>")
        headers = [
            "treatment",
            "outcome",
            "n",
            "estimate",
            "adjust",
            "criterion",
            "ref: pass / fail / inc",
        ]
        for h in headers:
            parts.append(f"<th>{h}</th>")
        parts.append("</tr></thead><tbody>")
        for row in sub.iter_rows(named=True):
            klass = "" if row["identified"] else "unident"
            ref = (
                f"<span class='pass'>{row['ref_pass']}</span>"
                f" / <span class='fail'>{row['ref_fail']}</span>"
                f" / <span class='inc'>{row['ref_inconclusive']}</span>"
            )
            parts.append(
                f"<tr class='{klass}'>"
                f"<td><code>{row['treatment']}</code></td>"
                f"<td><code>{row['outcome']}</code></td>"
                f"<td>{row['n']}</td>"
                f"<td>{row['estimate']:.4f}</td>"
                f"<td>{row['adjustment_set']}</td>"
                f"<td>{row['criterion'] or '-'}</td>"
                f"<td>{ref}</td>"
                f"</tr>"
            )
        parts.append("</tbody></table>")

    parts.append("</body></html>")
    out_path.write_text("".join(parts), encoding="utf-8")
    return out_path
