"""Row-order regression tests for the Paper 2 causal path.

Three things downstream of a Polars `unique()` reach a committed artifact:
the row order of each causal frame (which `_bootstrap_ci` resamples *by
position*, so a moving order moves the CI under a fixed seed), the row
order of `causal_effects.parquet`, and the section order of
`causal_report.html`.

`causal/flows/causal.py`'s drug loop is the fourth such site, but that
module has been unimportable since the Phase-0 package move (`from
..causal import` resolves to `causal.causal`), so it has no test here —
see docs/evaluations/unique-row-ordering.md §6.

See docs/evaluations/unique-row-ordering.md.
"""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path

import polars as pl

import causal.estimate as estimate_mod
from causal.data_prep import build_dr_causal_frames
from causal.estimate import estimate_dr_effects_per_subtype
from causal.report import write_html_report
from causal.schema import CausalEffect, EstimationMethod

N_ACCESSIONS = 60
N_RUNS = 8
SUBTYPES = ("A1", "C", "D", "B", "G", "A2", "F1", "CRF01_AE")
DRUGS = ("3TC", "AZT", "EFV", "NVP", "DTG", "RAL", "TDF", "ABC")


def _toy() -> tuple[pl.DataFrame, pl.DataFrame]:
    dr_rows: list[dict] = []
    cov_rows: list[dict] = []
    for i in range(N_ACCESSIONS):
        acc = f"A{i:05d}"
        dr_rows.append(
            {
                "accession": acc,
                "gene": "RT",
                "mutations": ["M184V"] if i % 2 else ["K103N"],
                "drug_class": "NRTI",
                "drug": "3TC",
                "level": "high" if i % 3 == 0 else "susceptible",
            }
        )
        cov_rows.append(
            {
                "accession": acc,
                "subtype": SUBTYPES[i % len(SUBTYPES)],
                "country": "KE",
                "collection_year": 2019,
            }
        )
        # A duplicate curated row, which is what the dedup in data_prep is for.
        cov_rows.append(cov_rows[-1])
    return pl.from_dicts(dr_rows), pl.from_dicts(cov_rows)


class TestCausalFrameOrder:
    def test_frame_row_order_is_stable_across_builds(self) -> None:
        """`_bootstrap_ci` indexes by position, so this decides the CIs."""
        dr, curated = _toy()
        orders = {
            tuple(build_dr_causal_frames(dr, curated)[0]["NRTI"]["accession"].to_list())
            for _ in range(N_RUNS)
        }
        assert len(orders) == 1, (
            f"build_dr_causal_frames produced {len(orders)} distinct row orders "
            f"in {N_RUNS} builds; a seeded bootstrap over this frame draws a "
            "different resample each time"
        )

    def test_frame_still_holds_one_row_per_accession(self) -> None:
        dr, curated = _toy()
        frame = build_dr_causal_frames(dr, curated)[0]["NRTI"]
        assert frame.height == N_ACCESSIONS
        assert frame["accession"].n_unique() == N_ACCESSIONS


class TestPerSubtypeEffectOrder:
    """`estimate_dr_effects_per_subtype` appends to the effects Parquet."""

    def _frame(self) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "subtype": [SUBTYPES[i % len(SUBTYPES)] for i in range(800)],
                "RT_184_V": [i % 2 for i in range(800)],
                "label_binary": [i % 3 == 0 for i in range(800)],
            }
        )

    def test_subtypes_are_visited_in_sorted_order(self, monkeypatch) -> None:
        seen: list[str] = []

        def _fake(df, mutation, *, outcome, data_scope="", **kw):
            seen.append(data_scope)
            return CausalEffect(
                treatment=mutation,
                outcome=outcome,
                data_scope=data_scope,
                n_observations=df.height,
                method=EstimationMethod.LOGISTIC_REGRESSION,
                estimate=0.5,
                identified=True,
                dag_version="0.1.0",
                estimated_at=date(2026, 8, 12),
            )

        monkeypatch.setattr(estimate_mod, "estimate_dr_effect", _fake)
        frame = self._frame()
        orders = set()
        for _ in range(N_RUNS):
            seen.clear()
            estimate_dr_effects_per_subtype(frame, "RT_184_V", min_samples=1)
            orders.add(tuple(seen))

        assert len(orders) == 1, (
            f"per-subtype estimation visited subtypes in {len(orders)} distinct "
            f"orders in {N_RUNS} runs; the effects Parquet's row order follows it"
        )
        visited = next(iter(orders))
        assert list(visited) == sorted(visited)


class TestHtmlReportSectionOrder:
    def _effects(self) -> list[CausalEffect]:
        return [
            CausalEffect(
                treatment="RT_184_V",
                outcome="label_binary",
                data_scope=f"DR/drug={d}/subtype={s}",
                n_observations=500,
                method=EstimationMethod.LOGISTIC_REGRESSION,
                estimate=0.5,
                identified=True,
                dag_version="0.1.0",
                estimated_at=date(2026, 8, 12),
            )
            for d in DRUGS
            for s in SUBTYPES
        ]

    def test_section_order_follows_the_effects_list(self, tmp_path: Path) -> None:
        effects = self._effects()
        expected = [e.data_scope for e in effects]

        orders = set()
        for run in range(N_RUNS):
            out = tmp_path / f"report_{run}.html"
            write_html_report(effects, out)
            orders.add(tuple(re.findall(r"<h2>(.*?)</h2>", out.read_text())))

        assert len(orders) == 1, (
            f"the HTML report emitted {len(orders)} distinct section orders in "
            f"{N_RUNS} runs on one effects list"
        )
        assert list(next(iter(orders))) == expected
