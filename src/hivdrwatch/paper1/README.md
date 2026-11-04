# `paper1` — Paper 1 figure aggregators + style

Pure-Python building blocks that Paper 1 figure notebooks import. Keeping aggregation here (not in the notebook) means the semantics are unit-testable without a matplotlib import.

## Contents

| File | Purpose |
|---|---|
| `style.py` | Okabe-Ito colorblind palette, canonical drug-class / subtype orderings, `apply_paper_style()` (Nature-journal rcParams with `pdf.fonttype=42` for selectable-text PDFs), `resistance_cmap()`. |
| `fig02_aggregate.py` | Fig 2 — resistance landscape (country × year × drug class × subtype). Aggregation locked with tests; min-cell-N masking. |
| `fig03_parity.py` | Fig 3 — cluster parity aggregator: loads the extended `cluster_validate` JSON with both our clusterer and HIV-TRACE results, produces a comparison frame and the provenance note (parity ran, HIV-TRACE not requested, or requested and skipped) shown on the parity panel. |
| `fig04_dtg.py` | Fig 4 — DTG case study aggregator: INSTI major-mutation prevalence pre/post-2019, stratified by country and subtype. |
| `fig04_its.py` | Fig 4 secondary — per-country segmented regression (interrupted time series) of yearly INSTI prevalence, `prevalence ~ t + post + t:post`, centred on the 2019 rollout. OLS hand-rolled on numpy/scipy (no statsmodels); `_ols` is cross-checked against `scipy.stats.linregress` in the tests. Shares the outcome definition with `fig04_dtg.build_insti_indicator`. |
| `fig05_scalability.py` | Fig 5 — reads `benchmarks/results/scale_*.json`, produces long-format frame + power-law fit per stage. |

## Figures

The renderers live in `notebooks/fig0N_*.py` (per CLAUDE.md: "no figure lives outside `notebooks/`"). Each notebook imports its aggregator + `style.py`, renders to `paper1_tools/figures/fig0N_*.{pdf,png}`, and prints a summary of what it built.

Regenerate:

```bash
make fig2
make fig3
make fig4
make fig5
```

Or all of Paper 1:

```bash
make paper1
```

## Tests

```bash
pytest tests/paper1 -v
```
