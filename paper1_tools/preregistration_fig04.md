# Preregistration — Paper 1 Figure 4 (DTG-era INSTI resistance case study)

**Registered:** *pending OSF submission — locked here in git until then*
**Analysis lead:** Kelyn Njeri (NeoForge Labs)
**Repository:** https://github.com/Neoforge-Labs/SAWA
**Frozen analysis code:** `notebooks/fig04_dtg_case_study.py` at git SHA `<TAG_ON_SUBMISSION>`
**Frozen data vintage:** DVC pointer `data/interim/dr_results.parquet.dvc` at same SHA

---

## 1. Background

In 2019 the WHO recommended dolutegravir (DTG) as the preferred first-line antiretroviral component for adults living with HIV in low- and middle-income countries. East African national programs (Kenya, Uganda, Tanzania, Rwanda, Burundi, Ethiopia) rolled out DTG-based regimens over 2019–2021. Pre-rollout resistance surveillance in the region focused on NRTI/NNRTI classes; the INSTI class was rarely genotyped because integrase inhibitors were not standard of care.

With DTG now the dominant first-line drug, INSTI major mutations — particularly R263K, G118R, T66A/I, E138A/K, Q148H/K/R, N155H — are the ones worth watching. Prior work has documented rare DTG resistance events under monotherapy or with pre-existing INSTI failure; whether population-level emergence is happening in East African subtypes A1, C, D under the current programmatic rollout is an open question.

## 2. Research question

**Is there evidence of emerging INSTI major-mutation prevalence in HIV-1 pol sequences from East Africa (KE, UG, TZ, RW, BI, ET) after the 2019 WHO DTG rollout, compared with the 2010–2018 baseline?**

## 3. Hypotheses

- **H1 (primary):** In post-rollout years (2020 onward), INSTI major-mutation prevalence in East African subtype A1, C, or D sequences exceeds the pre-rollout (2010–2018) baseline by ≥ 0.5 percentage points, with a 95 % confidence interval excluding zero.
- **H0:** Post-rollout prevalence is within the 95 % CI of the pre-rollout baseline.

## 4. Cohort definition (frozen before running)

**Inclusion criteria:**
- Sampling country ∈ {KE, UG, TZ, RW, BI, ET}
- Collection year ∈ [2010, latest available] (must include ≥ 3 pre-rollout and ≥ 2 post-rollout years to be scored)
- Subtype ∈ {A1, C, D, recombinant}
- Sequence covers HIV-1 integrase (IN region, HXB2 4230–5096) at ≥ 80 % coverage
- Passes `align/qc.py` blocking flags

**Exclusion criteria:**
- Duplicate patient (dedup by `patient_hash` where present)
- Sequence submitted under a clinical-trial identifier (would bias for treated patients)
- Countries with fewer than 30 sequences per pre/post period stratum

## 5. Exposure and outcome

**Exposure:** epoch ∈ {pre-rollout: 2010–2018, post-rollout: 2020–latest}. Rollout year 2019 is excluded from both epochs as an implementation-year buffer.

**Outcome:** binary indicator per sequence — presence of at least one Stanford HIVdb "major" INSTI mutation at Sierra level ∈ {intermediate, high} against any INSTI drug (DTG, BIC, CAB, EVG, RAL).

## 6. Statistical analysis (locked before running)

**Primary:** two-proportion test (Wilson score interval) per (country × subtype) stratum. Effect = post-rollout prevalence − pre-rollout prevalence, with 95 % CI. A stratum "shows a signal" iff CI excludes zero AND post > pre.

**Secondary:** interrupted time-series (ITS) with 2019 as the intervention point. Model: `prevalence ~ year + post + year:post` fit as segmented regression per country stratum. The `post` and interaction coefficients are the effect estimates of interest.

**Sensitivity analyses:**
1. Alternate epoch cutoffs (2018/2020, 2019/2021)
2. Restrict to drug-naive sequences (`is_drug_naive == True`) — likely reduces N drastically
3. Restrict to sequences with `IN` coverage ≥ 90 % (vs 80 %)
4. Restrict to subtype A1 only, C only, D only (per-subtype consistency)

If the primary result and > 2 of 4 sensitivity analyses agree in direction, we call it robust. If only the primary agrees, we report the primary with an explicit caveat that sensitivity analyses were inconclusive.

## 7. Multiple testing

The primary test runs across 6 countries × 4 subtypes = 24 strata. We apply Benjamini-Hochberg FDR correction at q = 0.05. Only BH-corrected significant strata contribute to the paper's headline claim.

## 8. Analysis code

- `src/hivdrwatch/paper1/fig04_dtg.py` — aggregator (unit-tested)
- `notebooks/fig04_dtg_case_study.py` — figure + statistical test
- `tests/paper1/test_fig04_dtg.py` — semantic tests

Any deviation from this plan must be documented as an exploratory finding, clearly labelled in the manuscript, and not counted toward the primary hypothesis test.

## 9. Data provenance

- Sequences: LANL Table exports + GenBank Entrez pulls, frozen at DVC pointer above
- Resistance calls: Stanford HIVdb Sierra, algorithm version pinned in `dr_results.parquet` per record
- Any manual curation is out-of-scope; the pipeline is fully deterministic given the frozen data vintage

## 10. Timeline

1. Preregister on OSF (this document as the submission) — before running the final analysis
2. Freeze `data/` vintage via DVC tag
3. Run `make fig4` on the frozen vintage — produce the figure + statistical test output
4. Amend `paper1_tools/manuscript.qmd` with results
5. Post preprint on bioRxiv

## 11. Authors

- **Kelyn Njeri** (NeoForge Labs, Nairobi) — lead analyst, corresponding author
- Additional co-authors added on manuscript submission per the standing MODULES.md authorship contract
