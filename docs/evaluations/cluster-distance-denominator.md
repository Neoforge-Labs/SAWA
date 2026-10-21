# The cluster distance denominator, and the FASTA it was handed

**Follows from** `docs/evaluations/align-coverage.md` §8 items 1 and 2,
which recorded both defects and deliberately left them: one is a
methodological change to a Paper 1 number. This closes both and
measures what they are worth on the real cohort.

**Needs Kelyn before any figure moves.** The distance metric is a
Paper 1 claim. Nothing in `data/gold/` has been regenerated here.

---

## 1. The two defects

**`cluster/distance.py::compute_pairwise_distances_python`** computed
`(row_i != row_j).sum() / seq_len` — the mismatch count over the full
3,012 alignment columns. A curated row is 3,012 columns wide whether
or not any sequence was placed on them, so dividing by the width
treats a column where **neither** sequence has a base as evidence that
the two agree there. Nobody sequenced it. Measured at the 0.015
threshold, on sequences built to the shape of real curated rows:

| pair | old: reported | linked at 0.015? | fixed | reference tn93 (§7) |
|---|---|---|---|---|
| two all-gap rows (3,012 dashes each) | 0.000000 | **yes** | no pair | no pair |
| 300 shared nt, 10% divergent there | 0.009960 | **yes** | 0.100000 | 0.127706 |
| 99 shared nt, identical there | 0.000000 | **yes** | no pair | no pair |
| two 500 nt rows covering different columns | 0.332005 | no | no pair | no pair |

Its docstring said it was the validation-harness path. The pipeline
disagreed: `pipeline/cluster.py` defaults `--distance-engine` to
`python`, and has since 452bfda. §6 and §7 settle which was
intended.

**`pipeline/cluster.py::materialize_aligned_fasta`** applied no filter
at all — it wrote every row whose `seq` was truthy, and a string of
3,012 dashes is truthy. All 21,358 all-gap rows went into the FASTA
handed to the distance engine, while `phylo/build.py` had been gating
its tree input on `passes_qc` and `MIN_ALIGNED_COVERAGE = 0.5` since it
was written. The tree and the network were built from different
cohorts.

## 2. Method

Exact over the whole shipped corpus, no sampling. `aligned.parquet`
holds 42,106 rows / 36,510 distinct accessions (the duplication is
`align-coverage.md` §7 item 3) behind **15,274 distinct `aligned_pol`
strings**; every accession has exactly one sequence and one QC verdict,
verified. Identical strings are distance 0 and always co-cluster, so
pair relations were computed over the 15,274 representatives with BLAS
— for each character `c`, `P_c @ P_c.T` counts the columns where both
sequences carry `c` — and then expanded back to accessions.

One exception, and it is the whole point of the fix: under the overlap
denominator a representative whose own unambiguous coverage is below
`min_overlap` cannot pair with **anything**, including an identical
copy of itself, so its accessions expand to singletons instead of
collapsing. Scoring that wrongly is what resurrects the all-gap block.

Two eligibility counts appear below and they differ by 7 records.
`phylo.build.filter_by_coverage` measures coverage as the non-`-`
fraction, so an `N` or an `R` counts as covered; the distance engine's
overlap counts only A/C/G/T. The pipeline's own gate (phylo's
definition, since that is the helper being reused) admits **3,397**
records; the probe, which used the stricter ACGT definition
throughout, admits 3,390. The 7 extra records are singletons under
either denominator and move nothing: end-to-end the engine produces
1,469 clusters where the probe predicted 1,462, and the two agree
exactly on the 7 extra singletons, on the largest component and on the
whole top-ten. That cross-check is what licenses the probe's numbers
for the columns too large to run end-to-end.

Scripts are in this change's scratch, not committed — the measurement
is a one-off against a vintage that is about to be replaced (§9).

## 3. The fixes

`compute_pairwise_distances_python` now divides by the count of
columns where **both** sequences carry an unambiguous base (A/C/G/T),
and refuses to emit a pair whose overlap is below `min_overlap`
(default **100**, matching tn93's own `-l` default: refuse the pair,
do not clamp its distance). Gap characters and IUPAC ambiguity codes
leave the denominator — an uncalled site is not evidence of agreement
either. This is one place where it departs from tn93's `-a resolve`,
which resolves an ambiguous site against its partner rather than
dropping it; dropping is the conservative direction.

`materialize_aligned_fasta` now takes the eligible record set from
`eligible_accessions`, which applies `passes_qc` and
`phylo.build.filter_by_coverage` at `MIN_ALIGNED_COVERAGE` — phylo's
own helper, not a second copy. It applies to the pre-existing-FASTA
branch too: a FASTA the curate stage wrote carries no QC or coverage
information, so the curated Parquet decides which of its records may be
used. It also writes one record per accession, since a FASTA with
repeated identifiers makes the distance CSV ambiguous about which
record a pair refers to.

Both gates are parameters (`--min-overlap`, `--min-coverage`) with
explicit defaults.

## 4. What it is worth on the real cohort

Shipped vintage, 0.015 threshold, 36,510 accessions. "width" is the
old denominator, "overlap" the new one at `min_overlap=100`; "filter"
is `passes_qc` + 50% coverage, which leaves **3,390 of 36,510 (9.3%)**
eligible on this vintage.

| | A: as the code stands | B: filter only | C: denominator only | D: both (new default) |
|---|---|---|---|---|
| sequences in | 36,510 | 3,390 | 36,510 | 3,390 |
| distinct-sequence pairs below threshold | 175,917 | 15,926 | 110,578 | 21,775 |
| clusters | 8,103 | 1,596 | 27,060 | **1,462** |
| singletons | 7,200 | 1,393 | 25,947 | 1,269 |
| sequences in a cluster of ≥2 | 29,310 (80.3%) | 1,997 | 10,563 | 2,121 |
| clusters of ≥10 | 116 | 62 | 47 | 63 |
| clusters of ≥100 | 16 | 1 | 4 | 2 |
| **largest component** | **19,164** | 150 | 6,307 | **195** |

The pair row counts pairs of *distinct sequences*, which is what the
probe computes; the engine emits one row per pair of accessions, so on
the real FASTA for column D it writes 35,852 rows for the same 21,775
distinct-sequence pairs — 4,029 of them at distance exactly 0, the
identical records behind `align-coverage.md` §7 item 3. Run end to end
on the filtered FASTA the engine gives **3,397 sequences → 1,469
clusters, 1,276 singletons, largest component 195**, 63 clusters of
≥10 and 2 of ≥100.

Largest ten components, A: 19164, 1212, 1192, 476, 449, 438, 379, 256,
211, 210. D: 195, 105, 62, 59, 50, 47, 36, 34, 34, 32.

**The 19,164-member component in column A is exactly the set of
accessions with no alignment at all** (`align-coverage.md` §2: "no
alignment at all — 19,164"). 52.5% of the corpus in one cluster,
mutually linked at distance 0.000000 because none of them had been
sequenced. Four of the next five components are the same failure at
smaller coverage.

**The correction is not one-directional.** Compare B with D: on the
same 3,390 sequences, the corrected denominator *adds* 37% more pairs
(15,926 → 21,775) and grows the largest component from 150 to 195. The
width denominator was wrong in both directions. It manufactured links
out of shared absence, and it also *suppressed* real ones: for two
sequences covering different parts of pol, every column covered by one
and not the other counts as a mismatch, so a pair that matches closely
where it actually overlaps was pushed above threshold by the regions
where only one of them exists. Any revision of figure 3 should expect
clusters to appear as well as disappear.

## 5. `min_overlap` is a safety net, not a knob

Sweeping it over the same vintage (clusters / largest component):

| `min_overlap` | no coverage filter | with coverage filter |
|---|---|---|
| 100 | 27,060 / 6,307 | 1,462 / 195 |
| 250 | 27,362 / 5,248 | 1,462 / 195 |
| 500 | 28,563 / 2,832 | 1,462 / 195 |
| 750 | 29,359 / 2,143 | 1,463 / 195 |
| 1,000 | 32,271 / 471 | 1,463 / 195 |
| 1,500 | 34,506 / 197 | 1,463 / 195 |

Two things follow. First, without the coverage filter a 100-site floor
is weak: a 6,307-member giant component survives it, assembled from
chains of sequences that barely overlap but are near-identical where
they do. Second, **with** the coverage filter the floor is inert —
every surviving record carries at least 1,506 unambiguous sites, so
every surviving pair overlaps far more than 100 and the whole sweep
returns the same answer. The coverage gate, not `min_overlap`, is what
controls this. `min_overlap` earns its place for the unfiltered path
and for a lower `--min-coverage`, not as something to tune.

## 6. The published artifact is a zero-cluster run

`data/gold/cluster_assignments.parquet` holds 36,510 rows in **36,510
clusters — every sequence a singleton**, and
`data/interim/cluster/distances.csv` is 17 bytes, a bare header. Both
are dated 2026-08-17, four days before 452bfda added the Python engine.
So the published cluster artifact is the *docker tn93* run, the one
that emitted zero pairs. None of column A above was ever shipped; it is
what the current code would produce on this vintage if it were run.

## 7. Why tn93 emitted zero pairs, and what it says once it works

`make images` had never been run here, so this is the first time the
pinned reference has been executed at all. Two things came out of it.

**The image is mis-tagged.** `docker/hivtrace/Dockerfile` clones
`--branch v1.0.13` and tags the image `hivdrwatch/tn93:1.0.13`, but
the binary it builds reports **v1.0.12**.

**Its single-file FASTA parser reads only the first record.** Given a
three-record, gap-free, LF-terminated FASTA of 600 nt sequences it
reports `"Sequences": 1`, and with `-0` (report every sequence against
itself) it emits exactly one row. Same from stdin, same wrapped at 60
columns, same with CRLF line endings, same writing to stdout with an
explicit `-f csv`, same with two 10 nt sequences, same on the real
3,397-record filtered FASTA. One sequence means no pairs, which is the
whole of "tn93 silently emits zero pairs on both synthetic and real
cohorts": not flags, not overlap, not divergence, and nothing the
FASTA writer did.

The defect sits in the build's shared FASTA reader, not in tn93's
argument handling: `seqcoverage`, a sibling tool in the same image
built from the same `tn93_shared.cc`, reports
`{"unique": 1, "total": 1}` on the same three-record file. The
distance code itself is fine — `-s SECOND_FASTA` (compare file A
against file B) reads one record from each and computes correct
distances.

This retires the hypothesis that tn93's silence was tn93 correctly
refusing low-overlap pairs. It was not. It also means **the reference
implementation has never validated our distance engine**, which is the
part that matters for a tools paper: for as long as the pipeline has
existed there has been no independent check on the numbers figure 3
rests on. Fixing the image belongs with the regeneration (§9).

Through the two-file path, which does work, the reference agrees with
the fix on every fixture in §1 — the all-gap pair, the disjoint pair
and the 99-site pair all emit nothing; the 300-nt 10%-divergent pair
is reported at 0.127706, far above threshold, where the old denominator
said 0.009960 and linked it. Two differences are worth recording:

* **`-l` is strictly greater than, not at least.** tn93's help says
  "overlap over at least OVERLAP nucleotides", but at the default the
  binary refuses a 100-site overlap and accepts 101. This
  implementation uses `overlap >= min_overlap`, the documented
  reading. With the coverage gate in force every surviving record
  carries ≥1,506 sites, so the one-site boundary cannot bind.
* **We apply the threshold to the uncorrected p-distance.** TN93 is a
  monotone transform of it, so memberships are identical *at a
  matched cut* — but the cut is not matched. Measured: a full-length
  pair at p = 0.014940 is reported by tn93 as 0.0151563, so it falls
  inside our 0.015 and outside tn93's. On the filtered cohort **759 of
  35,852 pairs sit in p ∈ [0.0147, 0.015]**, the band where the two
  conventions can disagree; 14 sit in the last 0.0001. Paper 1 claims
  the HIV-TRACE convention, so this is a decision for Kelyn: either
  apply the correction, or state in the methods that the threshold is
  applied to the uncorrected p-distance. It is not changed here —
  changing it changes the number again.

## 8. Figure 3 and the HIV-TRACE parity numbers

**Both sides of the parity comparison shared the same denominator, so
the parity numbers are internally comparable.** `cluster_validate`
runs `compute_pairwise_distances_python` **once** and hands the one CSV
to both our connected-components clusterer and `hivnetworkcsv`.
`cluster/hivtrace.py` computes no distances of its own — it parses a
network JSON. There is no asymmetry to correct.

**And the parity numbers do not move.** `cluster/simulate.py` emits
uniformly random ACGT with no gaps, no ambiguity codes and equal
length, so the pairwise overlap *is* the full alignment for every pair
and the two denominators coincide exactly. Verified: on the Fig 3
default simulation (20 clusters, sizes 5–15, 50 singletons, seed 42)
the old and new engines emit **1,066 pairs each, byte-identical**. P/R/F1,
ARI and pair-agreement are unchanged.

That is the problem with it. The harness could not have caught this
defect in any configuration, because the regime it validates is the one
regime where the defect vanishes. Figure 3 certifies the clusterer on
gap-free full-coverage input and says nothing about the input the
pipeline actually gets, where the median curated row had **0 nt** of
aligned coverage. Giving the simulator ragged coverage — draw a
start/end per sequence from the real coverage distribution — is what
would make figure 3 a test of the thing being claimed. It is not done
here.

## 9. What is left, and the order figure 3 needs

Every number in §4 and §5 is measured against the **asm10** vintage of
`aligned.parquet`, the one `align-coverage.md` showed to be wrong. They
are the like-for-like before/after of the metric change and nothing
else. They are **not** a preview of the figure. The 9.3% eligibility in
particular is an artefact of that vintage: at asm20 coverage ≥800 nt
rises from 42.5% to 94.1% and the median from 0 to 1,081 nt, so the
50%-coverage gate will admit far more — how much more cannot be read
off `align-coverage.md`'s published table and has to be measured on the
rebuilt vintage. **Do not mix vintages.** Regenerating figure 3 against
the asm10 parquet with the new metric would compare two different
things.

Order:

1. **Kelyn approves the metric change** (this document, §3). It changes
   a Paper 1 number. Nothing below should start first.
2. **Re-curate** `data/interim/curated/aligned.parquet` with the
   asm20 aligner and the `NO_ALIGNMENT` QC gate, and update DVC —
   `align-coverage.md` §8 item 4. Everything downstream derives from
   it.
3. **Re-measure eligibility** on the new vintage: how many of the
   36,510 clear `passes_qc` + 50% coverage. If that number is still
   small, `--min-coverage` is the next decision, and §5 says
   `min_overlap` will not substitute for it.
4. **Decide the `ingest/` bugs** — `align-coverage.md` §7. The CRF→C
   subtype collapse and the fabricated-country duplication both feed
   `outbreak_features` and any country- or subtype-stratified cluster
   count. Cluster *membership* is unaffected; what the clusters are
   *labelled* is not.
5. **Re-run `pipeline/cluster`** on the new vintage, then
   `outbreak_rollup`, then regenerate figure 3's cohort panel.
6. **Then** the open questions: fix the tn93 image so the reference
   can read a multi-record FASTA, re-tag it at the version it actually
   builds, and compare its pair set with ours on the filtered cohort
   (§7); decide the p-distance-vs-TN93 cut (§7); and give the
   simulator ragged coverage so the harness can see this class of
   defect (§8).

## 10. Tests

`tests/cluster/test_distance.py` — the all-gap pair, the 300-nt
10%-divergent pair and the 99-site refusal, each with the number the
old denominator reported; plus the ambiguity exclusion, the
`min_overlap` parameter, and two guards pinning the gap-free synthetic
case that must not move.

`tests/cluster/test_materialize_fasta.py` — the QC gate, the coverage
gate at phylo's threshold, all-gap rows carrying `passes_qc=True`,
duplicate-accession collapse, and the pre-existing-FASTA branch being
filtered rather than passed through.

`tests/cluster/test_features.py`'s 6 pre-existing failures were not
mine; they are fixed on `claude/sawa-audit-rename-3fc1d2` as b291210
and that commit is cherry-picked here rather than re-fixed, so
`pytest tests/cluster/` is green.
