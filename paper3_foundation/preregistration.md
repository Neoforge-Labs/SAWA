# Preregistration — Paper 3 v2 (foundation-model comparison against SAWA v0.1.0)

**Registered:** *pending OSF submission — locked here in git until then*
**Analysis lead:** Kelyn Njeri (NeoForge Labs)
**Repository:** https://github.com/Neoforge-Labs/SAWA
**Frozen analysis code:** `foundation_models/` at git SHA `<TAG_ON_SUBMISSION>` (`paper3-v0.1.0`)
**Frozen data vintage:** DVC pointers under `data/interim/*.dvc` at the same SHA
**Preregistered manuscript:** `paper3_foundation/manuscript.qmd` at the same SHA

---

## 1. Background

Foundation models trained on public protein or nucleotide corpora
are dominated by HIV-1 subtype B while the global disease burden is
dominated by non-B subtypes. Paper 3 v1
(`paper3_foundation/manuscript.qmd`) specifies a subtype-stratified
benchmark and reports non-FM baseline performance on synthetic
ground truth. This preregistration locks the v2 foundation-model
comparison against SAWA v0.1.0 in advance of running the fine-tune.

## 2. Research question

**Under a subtype-holdout split (train on subtype B, test on non-B),
what is the AUROC / PR-AUC gap between the training-distribution
subtype and the held-out subtypes for fine-tuned ESM-2 and Evo on
HIV-1 *pol* drug-resistance prediction, and does subtype-stratified
fine-tuning reduce that gap?**

## 3. Hypotheses

* **H1 (headline).** Under the subtype-holdout split, both ESM-2
  and Evo (vanilla fine-tune) exhibit a per-drug median AUROC gap
  between subtype B (train) and non-B (test) whose 95 % bootstrap CI
  excludes zero.
* **H2 (mitigation).** Subtype-stratified fine-tuning (curriculum:
  B, then A1, then C, then D) reduces the H1 gap by at least 50 %
  of its magnitude under vanilla fine-tune.
* **H3 (parity floor).** On the stratified split, at least one FM
  variant achieves per-drug median AUROC greater than or equal to
  the k-mer logistic regression baseline. If both FMs fail this
  test, the FM added value is null on this benchmark and the paper
  is a null result.

## 4. Cohort definition (frozen before running)

**Inclusion criteria:**

* Sampling country ∈ `{KE, UG, TZ, RW, BI, ET}`
* Collection year ∈ `[1994, 2024]` (SAWA v0.1.0 freeze)
* Subtype ∈ `{A1, B, C, D, CRF02_AG}` with at least 500 accessions
  per subtype after QC. Subtypes below the 500 threshold are
  reported in a supplement, not the primary comparison.
* Sequence covers HIV-1 *pol* at ≥ 60 % coverage
* Passes `src/hivdrwatch/align/qc.py` blocking flags

**Exclusion criteria:**

* Duplicate patient (dedup by `patient_hash` where present)
* Sequences submitted under a clinical-trial identifier
* Ambiguous-base fraction > 5 %

## 5. Splits (locked)

* **Stratified split** (`strategy="stratified"`): 70 / 15 / 15
  train / val / test stratified on `(subtype, drug, label_binary)`,
  seed = 42.
* **Subtype-holdout split** (`strategy="subtype_holdout"`,
  `holdout_train_subtypes=("B",)`): 80 / 20 train / val on B; test
  on the union of A1, C, D, CRF02_AG (subject to the 500-accession
  threshold), seed = 42.
* Splits are produced by `foundation_models.data.build_fm_dataset`
  and are deterministic given the seed. Split composition
  (`SplitCounts`) is committed to `paper3_foundation/figures/split_counts.json`.

## 6. Models (locked)

* **Baseline 1**: `KmerLogisticRegression(k=3, C=1.0, seed=42)`.
* **Baseline 2**: `OneHotRandomForest(n_estimators=200, seed=42)`.
* **ESM-2**: `esm2_t33_650M_UR50D` checkpoint. Two variants:
  frozen backbone with linear classification head; fully fine-tuned
  backbone + head. Both with cross-entropy loss, AdamW,
  learning-rate 1e-5, batch size 8, 3 epochs, seed 42.
* **Evo**: `evo-1.5-8k-base` checkpoint. Two variants: frozen
  backbone + linear head; full fine-tune. Same optimiser and seed.
* **Subtype-stratified fine-tune** (Paper 3's mitigation candidate):
  curriculum ordering B → A1 → C → D with equal per-subtype epoch
  budget; applied to both FMs. Total epoch budget matched to the
  vanilla fine-tune for a fair comparison.

All fine-tunes run under determinism flags
(`torch.use_deterministic_algorithms(True)`, cuDNN deterministic
mode).

## 7. Metrics (locked)

Per (drug × subtype × model) cell:

* AUROC, PR-AUC, Brier score, ECE at 10 bins
  (`foundation_models.evaluate.compute_metrics`).
* Cells with `n_positive < 5` are omitted.
* Bootstrap CIs: 500 accession-level percentile bootstrap resamples.

**Primary reporting metric:** per-drug median AUROC across subtypes,
with the B-vs-non-B gap as the derived quantity of interest.

## 8. Multiple testing

Benjamini–Hochberg FDR at `q = 0.05` across the drug × subtype cell
grid. Reported at the per-drug aggregate level; per-cell effects
reported without correction as supplementary.

## 9. Sensitivity analyses (planned, secondary)

1. Alternative aggregation: per-drug mean AUROC instead of median.
2. Alternate coverage threshold (≥ 80 % vs ≥ 60 %).
3. Restriction to `is_drug_naive == True` sequences.
4. Alternate FM checkpoints (`esm2_t12_35M`, `evo-1.0-8k`).
5. LoRA / QLoRA fine-tune with 1/10 parameter budget.

If the H1 result and > 2 of 5 sensitivity analyses agree in
direction, the finding is called *robust*.

## 10. Analysis code

* `foundation_models/data.py` — dataset adapter (Sprint 3)
* `foundation_models/kmers.py` — k-mer featurisation (Sprint 3)
* `foundation_models/baselines.py` — non-FM baselines (Sprint 3)
* `foundation_models/evaluate.py` — metrics (Sprint 3)
* `foundation_models/models/esm2.py`, `models/evo.py` — FM adapters
  (v2)
* `foundation_models/train.py` — fine-tune loops (v2)
* `foundation_models/failure_analysis.py` — mutation-class breakdown
  (v2)
* `benchmarks/paper3_baselines_synth.py` — synthetic-ground-truth
  baseline benchmark (Sprint 3, CI-verified)

Any deviation from this plan is documented as an exploratory finding
and excluded from the headline.

## 11. Data provenance

* Sequences: LANL Table exports + GenBank Entrez pulls, DVC pointers
  at the tag SHA.
* Resistance calls: Stanford HIVdb Sierra, algorithm version pinned
  per row in `dr_results.parquet`.
* FM checkpoints: fetched from the public model hubs at the tag SHA
  and cached under `data/models/` (DVC-tracked).

## 12. Compute plan

* Baselines: laptop-class, ≤ 30 min end-to-end.
* FM fine-tunes: single A100-class GPU, expected ≤ 24 h per variant.
  Full compute budget for the paper is ≤ 200 GPU-hours.
* Compute credit source is disclosed in Acknowledgements at v2
  submission time. No commercial dependency.

## 13. Timeline

1. Preregister on OSF (this document) — before running the final
   analysis.
2. Tag `paper3-v0.1.0`. Zenodo DOI issued.
3. Run baseline analysis on SAWA v0.1.0 real data; commit results
   under `paper3_foundation/figures/`.
4. Run FM fine-tunes; commit checkpoints (DVC) and metrics.
5. Compile Paper 3 v2; submit to Bioinformatics or the NeurIPS
   Datasets & Benchmarks track.

## 14. Deviations from this plan (to be filled during v2 write-up)

<!--
Any deviation observed during the v2 analysis is documented here
with its date, rationale, and effect on the headline. An empty
section at v2 submission is the desired outcome.
-->
