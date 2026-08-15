# Preregistration — CausalPhylo v2 (real-data application against SAWA v0.1.0)

**Registered:** *pending OSF submission — locked here in git until then*
**Analysis lead:** Kelyn Njeri (NeoForge Labs)
**Repository:** https://github.com/Neoforge-Labs/SAWA
**Frozen analysis code:** `causal/` at git SHA `<TAG_ON_SUBMISSION>` (`paper2-v0.1.0`)
**Frozen data vintage:** DVC pointers under `data/interim/*.dvc` at the same SHA
**Preregistered manuscript:** `paper2_causal/manuscript.qmd` at the same SHA

---

## 1. Background

Public HIV-1 drug-resistance surveillance derives per-drug Sierra
levels from *pol* sequences without observed individual treatment
history. Naive associations between mutation features and Sierra
levels conflate causal drug pressure with subtype-specific
polymorphism and non-random sampling. CausalPhylo v1
(`paper2_causal/manuscript.qmd`) specifies a DAG-based
identification framework and validates it on synthetic ground truth.
This preregistration locks the real-data application against SAWA
v0.1.0 in advance of running the analysis.

## 2. Research question

**For each mutation × drug-class combination in the SAWA v0.1.0 East
African cohort, what is the identified causal effect on the Sierra
resistance level, adjusting for the observed backdoor set
`{subtype, country, collection_year}` and (where populated)
`is_drug_naive`? Which of the identified effects survive placebo,
random-common-cause, and subset refuters?**

## 3. Hypotheses

* **H1 (positive control):** Effects for well-established
  resistance mutations in the Stanford HIVdb major-mutation list
  (e.g. `RT_184_V`, `RT_103_N`, `IN_148_H`, `IN_155_H`) are
  identified with 95 % bootstrap CIs excluding OR = 1 and pass
  ≥ 2 of 3 refuters.
* **H2 (subtype-confounded decoy):** Mutations known to be
  subtype-associated polymorphisms (Stanford HIVdb non-resistance
  polymorphism list) that show significant *unadjusted* marginal
  associations with Sierra levels *lose significance* after
  adjustment or fail placebo refutation.
* **H3 (novel-effects reporting):** A set of currently-unclassified
  mutation × drug-class edges may be identified with confidence.
  These are reported as candidate effects to be adjudicated by
  clinical / bench collaborators, not as causal claims in the
  headline.

## 4. Cohort definition (frozen before running)

**Inclusion criteria:**

* Sampling country ∈ `{KE, UG, TZ, RW, BI, ET}`
* Collection year ∈ `[1994, 2024]` (as of SAWA v0.1.0 freeze)
* Subtype ∈ `{A1, B, C, D, CRF02_AG}` — other/recombinant reported
  as a supplement
* Sequence covers HIV-1 *pol* at ≥ 60 % coverage
* Passes `src/hivdrwatch/align/qc.py` blocking flags

**Exclusion criteria:**

* Duplicate patient (dedup by `patient_hash` where present)
* Sequences submitted under a clinical-trial identifier
* Countries × subtypes with fewer than 100 accessions after QC
  (reported as supplement but not used for primary estimates)

## 5. Adjustment set (locked)

* Backdoor adjustment set for each `mutation_i → sierra_level`
  edge: `{subtype, country, collection_year}` as derived by DoWhy
  from `causal/dag.py::DRDag`.
* Additional observed proxy: `is_drug_naive`, added when
  non-null. Per-accession completeness is population-dependent and
  reported.

## 6. Estimation (locked)

* Outcome model: logistic regression conditioning on the
  adjustment set (DoWhy `backdoor.generalized_linear_model` with
  `glm_family=Binomial()`).
* Effect summary: adjusted odds ratio with 95 % bootstrap
  percentile CI (200 accession-level resamples).
* Aggregation: per drug class (NRTI, NNRTI, PI, INSTI) in the
  primary specification. Per-drug aggregation is a preregistered
  secondary specification.

## 7. Refutation (locked)

Three refuters per estimate:

1. **Placebo treatment** — mutation column replaced by a permuted
   binary vector; refuter PASSES if
   `|new_effect_estimate - 0| < 0.1` (on the log-odds scale).
2. **Random common cause** — random confounder added; refuter
   PASSES if `|new / original - 1| < 0.20`.
3. **Subset validity** — re-estimate on random 80 % subset;
   refuter PASSES if `|new / original - 1| < 0.20`.

An effect is *retained in the headline* iff it passes ≥ 2 of 3
refuters AND its 95 % CI excludes OR = 1. Effects that fail ≥ 2 of 3
refuters are *flagged as candidate spurious associations* in the
manuscript, not silently dropped.

## 8. Feature pruning

`causal/pruning.py::derive_prune_config` derives a conservative
prune config from the effect frame. A feature is dropped only if
BOTH (i) it is identified in the DAG, AND (ii) it either (a) has a
CI that spans 0 with `|effect| < 0.5 × (ci_high − ci_low)`
(`SMALL_EFFECT_THRESHOLD = 0.5`) OR (b) failed at least 2 of 3
refuters (`MIN_FAILED_REFUTATIONS = 2`). Effects identified in any
per-subtype scope survive as separate rows and are not erased by
pooling, so subtype-heterogeneous effects remain visible in the
final report.

## 9. Multiple testing

Benjamini–Hochberg FDR at `q = 0.05` is applied across the set of
mutations that survive both the CI criterion and refutation. Only
BH-corrected significant effects contribute to the manuscript's
headline. The correction is applied at the **mutation** level after
per-mutation aggregation across drugs within a class.

## 10. Sensitivity analyses (planned, secondary)

1. Alternate CI width (99 % percentile CI).
2. Restriction to `is_drug_naive == True` sequences (likely reduces
   `N` substantially).
3. Alternative outcome model (random-forest propensity + IPW).
4. Per-country stratification (report where per-country effects
   differ from pooled by > 20 % on the log-odds scale).
5. Alternate DAG variants:
   * Adding `viral_load` (when available) as an observed proxy for
     drug pressure.
   * Removing `country` from the adjustment set (national programme
     as source of confounding vs. as mediator).

If the primary result and > 2 of 4 sensitivity analyses agree in
direction, the finding is called *robust*. If only the primary
agrees, it is reported with an explicit caveat.

## 11. Analysis code

* `causal/data_prep.py` — adapter (Sprint 1, unit-tested)
* `causal/dag.py::DRDag` — DAG definition
* `causal/dag.py::to_dowhy_gml(graph)` — GML serialiser for DoWhy
* `causal/estimate.py` — identify → estimate → refute driver;
  `_run_refutations` dispatches the three refuters and
  `_score_refutation` scores each PASSES / FAILS / INCONCLUSIVE
* `causal/flows/causal.py` — Prefect flow orchestrating end-to-end
* `causal/pruning.py::derive_prune_config` — conservative post-hoc
  pruning (drop only if identified AND (small-effect-CI-includes-zero
  OR ≥ 2/3 refuters failed))
* `benchmarks/paper2_synth_recovery.py` — synthetic-ground-truth
  benchmark (CI-verified)
* `notebooks/paper2_*.py` — Marimo notebooks for each figure

Any deviation from this plan must be documented as an exploratory
finding, clearly labelled in the manuscript, and not counted toward
the primary hypothesis test.

## 12. Data provenance

* Sequences: LANL Table exports + GenBank Entrez pulls, frozen at
  DVC pointers under `data/interim/*.dvc` at the tag SHA.
* Resistance calls: Stanford HIVdb Sierra, algorithm version pinned
  per row in `dr_results.parquet`.
* Any manual curation is out of scope; the pipeline is fully
  deterministic given the frozen data vintage.

## 13. Timeline

1. Preregister on OSF (this document as submission) — before running
   the final analysis.
2. Tag `paper2-v0.1.0`. Zenodo DOI issued for the frozen state.
3. Run the primary analysis exactly as specified above; commit
   result JSON and figures under `paper2_causal/figures/`.
4. Run sensitivity analyses; commit under `paper2_causal/figures/sensitivity/`.
5. Compile CausalPhylo v2 manuscript; submit to PLOS Computational
   Biology (or Epidemics as fallback).

## 14. Deviations from this plan (to be filled during v2 write-up)

<!--
Any deviation observed during the v2 analysis is documented here
with its date, rationale, and effect on the primary hypothesis test.
An empty section at v2 submission is the desired outcome.
-->
