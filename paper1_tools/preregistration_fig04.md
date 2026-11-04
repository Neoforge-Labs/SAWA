# Preregistration — Paper 1 Figure 4 (DTG-era INSTI resistance case study)

**Registered:** *pending OSF submission — locked here in git until then*
**Version:** v3 (2026-11-04) — see §12. Amended **before** OSF filing and
before the analysis has been run against the frozen vintage.
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

- **H1 (primary):** In post-rollout years (2020 onward), INSTI major-mutation prevalence in East African subtype A1, C, or D sequences exceeds the pre-rollout (2010–2018) baseline. A stratum supports H1 only when **all three** of the following hold (see §6 for the operational definition):
  1. post-rollout prevalence > pre-rollout prevalence,
  2. the increase is at least **0.5 percentage points** (a clinical-relevance floor, so a statistically significant but trivial shift is not reported as a public-health signal),
  3. Fisher's exact *p*, after Benjamini–Hochberg correction within the primary tier, is ≤ 0.05.
- **H0:** Post-rollout prevalence is not distinguishable from the pre-rollout baseline under the criterion above.

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
- Strata (country × subtype) with fewer than **5** sequences in either epoch are dropped entirely and are not analysed or reported.

**Tiered inclusion (amended v2, see §12):**

| Tier | Threshold | Role |
|---|---|---|
| **Primary** | n ≥ 30 in *both* epochs | Confirmatory. Carries the preregistered FDR budget and the paper's headline claim. |
| **Secondary** | n ≥ 5 in both epochs, but not both ≥ 30 | Exploratory. Reported with its own separate BH correction, labelled exploratory in the manuscript, and never counted toward H1. |

The v1 text excluded sub-30 strata outright. Reporting them as an explicitly exploratory tier is strictly more informative than discarding them, provided they cannot contaminate the confirmatory result — which the separate FDR budget guarantees.

## 5. Exposure and outcome

**Exposure:** epoch ∈ {pre-rollout: 2010–2018, post-rollout: 2020–latest}. Rollout year 2019 is excluded from both epochs as an implementation-year buffer.

**Outcome:** binary indicator per sequence — presence of at least one Stanford HIVdb "major" INSTI mutation at Sierra level ∈ {intermediate, high} against any INSTI drug (DTG, BIC, CAB, EVG, RAL).

## 6. Statistical analysis (locked before running)

**Primary test — Fisher's exact** (amended v2, see §12) on the 2 × 2
pre-vs-post × resistant-vs-susceptible contingency table, per
(country × subtype) stratum. Fisher's exact makes no normal-approximation
assumption and stays valid down to expected counts of 1 per cell, which
the smaller East African strata (Rwanda, Burundi, Ethiopia) require. The
odds ratio uses a Haldane–Anscombe 0.5 correction when any cell is zero
so it stays finite.

**Interval estimates.** Per-epoch prevalence is reported with
**Clopper–Pearson exact binomial** 95 % CIs — conservative but valid at
any *n*, including the extreme p̂ where the Wilson score interval named
in v1 becomes anti-conservative. The pre→post difference is reported
with a **Newcombe hybrid-score** 95 % CI. Both are reported for every
stratum; neither is the significance criterion.

**Signal criterion.** A stratum "shows a signal" iff **all three** hold:

1. post-rollout prevalence > pre-rollout prevalence;
2. the difference is ≥ **0.005** (0.5 percentage points), per H1;
3. the BH-adjusted Fisher *q* within that stratum's tier is ≤ 0.05.

Criteria (1) and (3) are the statistical test; (2) is the
clinical-relevance floor from H1. The Newcombe CI is reported and
plotted but is deliberately *not* part of the criterion — using both a
CI-excludes-zero rule and an FDR rule would be two tests of the same
hypothesis at two different error rates.


**Secondary — interrupted time-series (ITS)** with 2019 as the intervention
point. Per-country segmented regression of yearly INSTI-resistance
prevalence:

    prevalence ~ t + post + t:post

Operationalised as follows (specified in v3; v1 and v2 left these details open):

* The unit of observation is one (country, year) with `prevalence = n_resistant / n_sequences`, using the same outcome definition as §5 and the same cohort as the primary analysis.
* `t = year − 2019`, so `t = 0` is the interruption and the `post` coefficient reads as the level change *at* the interruption rather than an extrapolation to year zero.
* The rollout year 2019 is excluded from the series, exactly as §5 excludes it from both epochs. `post = 1` for years after it.
* The `post` coefficient (level change) and the `t:post` interaction (slope change) are the effect estimates of interest, each with a 95 % CI and a t-test on the residual degrees of freedom.
* A country with only one segment, or with fewer than 5 usable years, cannot identify the four parameters and is reported as not fittable rather than fit.
* Estimation is ordinary least squares. This is the standard segmented-regression formulation, but it is OLS on a proportion, so the homoscedastic standard errors are approximate — the residual variance of a proportion depends on both prevalence and the per-year denominator. Per-year `n_sequences` is reported alongside every fit. This is why ITS is secondary and the Fisher test in the primary tier carries the confirmatory claim.

Implemented and unit-tested in `src/hivdrwatch/paper1/fig04_its.py`.

**Sensitivity analyses:**
1. Alternate epoch cutoffs (2018/2020, 2019/2021)
2. Restrict to drug-naive sequences (`is_drug_naive == True`) — likely reduces N drastically
3. Restrict to sequences with `IN` coverage ≥ 90 % (vs 80 %)
4. Restrict to subtype A1 only, C only, D only (per-subtype consistency)

If the primary result and > 2 of 4 sensitivity analyses agree in direction, we call it robust. If only the primary agrees, we report the primary with an explicit caveat that sensitivity analyses were inconclusive.

## 7. Multiple testing

The primary test runs across at most 6 countries × 4 subtypes = 24
strata; strata that fail the §4 thresholds are never tested, so the
realised family is smaller.

Benjamini–Hochberg FDR correction at q = 0.05 is applied **separately
within each tier** (amended v2, see §12). The primary tier carries the
preregistered FDR budget and is the only tier that contributes to the
headline claim. Secondary-tier q-values are computed over the secondary
family alone and are reported for context only — pooling the two
families would spend the confirmatory budget on exploratory strata.

## 8. Analysis code

- `src/hivdrwatch/paper1/fig04_dtg.py` — primary aggregator, epoch comparison (unit-tested)
- `src/hivdrwatch/paper1/fig04_its.py` — secondary aggregator, segmented regression (unit-tested)
- `notebooks/fig04_dtg_case_study.py` — figure + statistical test
- `tests/paper1/test_fig04_dtg.py` — semantic tests for the primary
- `tests/paper1/test_fig04_its.py` — semantic tests for the ITS, including a cross-check of the hand-rolled OLS against `scipy.stats.linregress`

Any deviation from this plan must be documented as an exploratory finding, clearly labelled in the manuscript, and not counted toward the primary hypothesis test.

Deviations found **before** OSF filing are corrected by amending this
document, with the change and its rationale recorded in §12. Deviations
found **after** filing are reported as deviations in the manuscript and
are never resolved by editing this document.

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

## 12. Amendment log

### v2 — 2026-11-04 (pre-filing, pre-analysis)

A repository audit found that `src/hivdrwatch/paper1/fig04_dtg.py`
implemented a different statistical plan from the one this document
locked in v1. Because the document has not been filed on OSF and the
analysis has not been run against the frozen vintage, the divergence is
resolved here, in advance, rather than reported as a deviation later.

| # | v1 said | v2 says | Why |
|---|---|---|---|
| 1 | Primary test: two-proportion test with Wilson score interval | Primary test: Fisher's exact; Clopper–Pearson per-epoch CIs; Newcombe CI for the difference | Fisher's exact needs no normal approximation and stays valid at the cell counts the smaller strata actually have. Clopper–Pearson is valid at extreme p̂ where Wilson is anti-conservative. |
| 2 | Signal iff "CI excludes zero AND post > pre" | Signal iff direction up **and** ≥ 0.5 pp **and** BH-adjusted Fisher *q* ≤ 0.05 | v1 mixed a CI rule with the §7 FDR rule, which tests the same hypothesis twice at two error rates. v2 makes the FDR rule the single significance criterion and states H1's 0.5 pp floor as an explicit conjunct — it had no operational definition in the code before this amendment. |
| 3 | Excluded strata with < 30 sequences per epoch outright | Two tiers: primary n ≥ 30 (confirmatory), secondary n ≥ 5 (exploratory, separate FDR family); < 5 still dropped | Reporting small strata as explicitly exploratory is more informative than discarding them, and a separate FDR family keeps them out of the confirmatory result. |
| 4 | BH correction over one family of 24 | BH applied within each tier separately | Follows from #3. |

**Direction of the change.** Every amendment above is neutral or
*more* conservative than v1: Fisher's exact and Clopper–Pearson are
both conservative relative to their normal-approximation counterparts,
the 0.5 pp floor was in v1's H1 but had no operational definition and
now gates the signal, and the tier split can only remove strata from
the confirmatory family, never add them. No amendment makes a positive
finding easier to obtain.

**Outstanding at v2.** The §6 secondary ITS analysis is still
unimplemented. It must be built and unit-tested before the analysis
runs, or withdrawn in a dated v3 amendment before filing.

### v3 — 2026-11-04 (pre-filing, pre-analysis)

Closes the item left outstanding at v2: the §6 secondary ITS analysis is
now implemented (`src/hivdrwatch/paper1/fig04_its.py`) and unit-tested
(`tests/paper1/test_fig04_its.py`). The commitment is kept, not withdrawn.

v1 and v2 specified the ITS only as a model formula. v3 writes down the
operational detail that formula left open — the unit of observation, the
centring of `t` on 2019, the exclusion of 2019 from the series, which
coefficients are the effect estimates, the not-fittable rule, and the
estimator. None of this changes the hypothesis or the primary analysis;
it removes the latitude that would otherwise let the implementation
choose after seeing the data.

Two points recorded deliberately:

* **2019 is excluded from the ITS series**, matching §5's treatment of it as an implementation buffer, so the secondary analysis runs on the same cohort as the primary. v1 said "2019 as the intervention point" without saying which segment the year itself belonged to.
* **The estimator is OLS on a proportion**, with the heteroscedasticity caveat stated in §6 rather than left for a reviewer to raise. Per-year denominators are reported with every fit. ITS remains secondary; the confirmatory claim rests on the primary tier.

No new dependency was added — the regression is least squares on numpy
with a scipy t reference distribution, and the tests check it against
`scipy.stats.linregress` on the single-predictor case it reduces to.
