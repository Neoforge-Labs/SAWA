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

## 11. The reference implementation, fixed — and what it says about ours

§9 item 6 deferred three things. This closes the first two: the image is
fixed and correctly tagged, and its pair set has been compared with ours
on the filtered cohort. **No unexplained disagreement survived.**

Nothing in `data/gold/` was regenerated and the metric is unchanged.
§3 still needs Kelyn.

### 11.1 Neither defect was what §7 took it for

**The mis-tag was upstream's, not ours.** `docker/hivtrace/Dockerfile`
cloned `--branch v1.0.13` and got exactly that commit — upstream simply
never bumped `VERSION.txt` at that tag, so it still reads `v1.0.12`. The
build was faithful; the binary's self-report was the tag's own slip.
Upstream appears to have noticed: a second tag, `v.1.0.13.actual`,
points at the same commit `a40aa62`.

**The reader bug is a known upstream defect, and it is architectural,
not version-specific.** It is veg/tn93 issue #46, fixed in **v1.0.15**
by commit `5c14f35`, which changed one declaration:

```c
-char validFlags [256];
+signed char validFlags [256];
```

`initAlphabets` fills that table with `-1` for every byte that is not a
sequence character, and `readFASTA`'s state machine asks
`validFlags[currentC] >= 0` to decide whether to append a byte to the
current sequence. Plain `char` is **unsigned on aarch64** and signed on
x86-64. Unsigned, `-1` reads back as 255, so the test is true for every
byte — including `>` (62) — and the reader, once past the first header,
never recognises another record boundary. Confirmed directly in the
image's own base: in `debian:12-slim` on linux/arm64, `((char)-1) < 0`
is false and `f[62] = -1` reads back as 255.

The record count §7 measured is the arithmetic signature of exactly
this. A three-record file of 600 nt sequences named `seq0..seq2` was
read as "1 sequences of length **1815**", and 1815 is what you get by
parsing the first header correctly and then appending every subsequent
byte, newlines and later headers included:

```
600 + 1 + len(">seq1") + 1 + 600 + 1 + len(">seq2") + 1 + 600 + 1 = 1815
```

That is why it reproduced on CRLF, on stdin, wrapped at 60 columns, on
10 nt sequences and on the real cohort: none of those are the variable.
The variable was the build host. It also explains the two survivors —
`-s SECOND_FASTA` and the §1 fixtures worked because each file held one
record, where "concatenate everything after the first header" and
"parse correctly" coincide.

So "tn93 v1.0.13 silently emits zero pairs" was never a property of
tn93, of our FASTA writer, or of the cohort. Those claims have been
corrected in `cluster/distance.py`, `pipeline/cluster.py`,
`pipeline/cluster_validate.py` and `cluster/simulate.py`.

### 11.2 The fix

Pinned to **v1.0.17** (`ef82c1b`), the last 1.0.x patch release:
it carries the v1.0.15 reader fix, its `VERSION.txt` matches its tag, and
unlike v1.1.0 it adds no new subquadratic distance path to a binary whose
only job here is to be the reference. The image is
`hivdrwatch/tn93:1.0.17` and `tn93 -v` reports `v1.0.17`; the tag,
`DEFAULT_IMAGE` and the binary now agree by construction. The tag was
updated in `Makefile`, `.github/workflows/benchmarks.yml`,
`cluster/distance.py` and `cluster/README.md`.

`python3` joins the image's build-only dependencies: from v1.0.15
tn93's `CMakeLists.txt` calls `find_package(PythonInterp 3 REQUIRED)`
unconditionally to wire up its ctest suite, with no option to skip it.
It is purged with the rest of the toolchain and does not reach the final
image.

**The regression check.** `RUN tn93 --help` passed throughout — it never
reads a sequence, which is why this went unnoticed. It was also written
as `… && ldd … | grep -q "not found" && exit 1 || true`, and that
trailing `|| true` swallowed any failure earlier in the chain, so the
layer could not have failed loudly as its comment claimed. Three
assertions replace it, under `set -eux` and with no `|| true`:

1. `ldd` resolves every library (the original libgomp regression).
2. `tn93 -v` equals the tag that was cloned — the exact drift that
   shipped a `v1.0.12` binary as `1.0.13`.
3. A two-record FASTA, 300 shared sites differing at 4, must report
   `"Sequences" : 2` and emit exactly one pair row. Under the aarch64
   bug this file parsed as one concatenated sequence and produced an
   empty CSV.

Two further guards need no Docker daemon and so run in the ordinary unit
pass (`tests/cluster/test_distance.py::TestPinnedReferenceImage`):
`DEFAULT_IMAGE` must match the Dockerfile's `ARG TN93_VERSION`, and the
pin must not regress below v1.0.15.

### 11.3 Engine versus reference, 3,397 records at threshold 0.015

`materialize_aligned_fasta` over the asm10 `aligned.parquet` yields the
3,397-record filtered FASTA; both engines read that one file. Pairs are
compared unordered. tn93 was also run at `-t 0.05` so that a pair it
scores *above* threshold can be told apart from one it never scored.

A **third** convention difference governs this comparison, beyond §7's
two. `run_tn93` passes `-a resolve`, tn93's default, which resolves an
IUPAC code to whichever base minimises the distance; our engine drops
ambiguous sites from numerator and denominator alike. tn93's `-a skip`
is the setting that matches ours. Both are reported.

**Matched ambiguity handling (`-a skip`) — the like-for-like run:**

| | pairs |
|---|---|
| `compute_pairwise_distances_python` @ 0.015 | 35,852 |
| reference tn93 `-a skip` @ 0.015 | 35,655 |
| intersection | **35,655** |
| ours only | 197 |
| tn93 only | **0** |
| **unexplained** | **0** |

tn93's pair set is a strict **subset** of ours, and all 197 extra pairs
sit in the one band §7 predicted: p ∈ [0.014784, 0.014983] against a
TN93-corrected [0.0150004, 0.0152116] — inside our cut, outside tn93's,
every one of them. §7 measured 759 pairs in p ∈ [0.0147, 0.015] as the
band where the conventions *can* disagree; 197 of them actually do.

The `-l` boundary binds on **nothing**, as §7 expected: the coverage
gate leaves every surviving record with ≥1,506 unambiguous sites, so no
pair comes near the 100/101 edge.

On the 35,655 shared pairs the two distances differ by
−3.4 × 10⁻⁷ to +2.7 × 10⁻⁴, median +2.3 × 10⁻⁵ — the TN93 correction,
which is monotone in p and adds more the more divergent the pair. The
small negative tail is three orders of magnitude below the positive one,
near tn93's printed precision, and it changed no link decision: it
cannot, since tn93 emitted no pair we did not.

**Production ambiguity handling (`-a resolve`):** 35,682 pairs,
intersection 35,668, 184 ours-only, **14 tn93-only**. The 184 are again
all in the correction band. The 14 are all pairs where our p-distance
is just above 0.015 and tn93's is just below — and all 14 carry between
83 and 117 sites where one sequence has an IUPAC code and the other a
called base. Those are the sites `resolve` resolves favourably and ours
discards, so they are the ambiguity convention, not a defect: the same
pairs vanish under `-a skip`, where tn93-only is 0.

### 11.4 What this does and does not establish

It establishes that our denominator, our overlap rule and our mismatch
count agree with the CDC reference exactly, on the real cohort, once
ambiguity handling is matched — the independent check the pipeline has
never had. Every disagreement is attributable to a stated convention,
and no disagreement is attributable to our arithmetic.

It does **not** settle the conventions themselves. Two decisions remain
for Kelyn, and neither is changed here:

1. **The cut** (§7). We threshold the uncorrected p-distance; Paper 1
   claims the HIV-TRACE convention. Measured cost at 0.015: 197 pairs,
   0.55% of ours.
2. **Ambiguity handling.** `run_tn93` resolves, our engine discards.
   Either match the reference, or state the rule in the methods. Only
   `-a skip` makes the two engines comparable, so an A/B run that wants
   to measure *our* metric should pass it.

One operational note: tn93's CSV row order is not stable between runs
on the same input (OpenMP scheduling); the pair set is. Compare pair
sets, never file bytes.

### 11.5 Still open from §9

Item 6's third part is untouched: `cluster/simulate.py` still emits
gap-free equal-length sequences, so the Fig 3 harness still cannot see
this class of defect (§8). Giving the simulator ragged coverage remains
the change that would make figure 3 a test of the input the pipeline
actually gets.

## 12. The cut: what changes if the correction is applied first

§11.4 left this open. Here is the measured cost, at the level that
matters. **Nothing is changed in the code — §3 and this both still need
Kelyn.**

Paper 1 claims the HIV-TRACE convention, and HIV-TRACE thresholds the
**TN93-corrected** distance. We threshold the uncorrected p-distance.
So this is not a free choice between two defensible conventions: one of
them is the one the manuscript says we use.

### 12.1 Three options, measured

Same 3,397-record filtered FASTA, ambiguity handling matched throughout
(`-a skip`), connected components via `cluster/assembly.py` at 0.015.

| | pairs | clusters ≥2 | accessions clustered | largest |
|---|---|---|---|---|
| **A.** p-distance cut at 0.015 *(status quo)* | 35,852 | 193 | 2,121 (62.44%) | 195 |
| **B.** TN93 cut at 0.015 *(HIV-TRACE)* | 35,655 | 192 | 2,118 (62.35%) | 195 |
| **C.** recalibrated p cut at 0.014796 | 35,654 | 192 | 2,118 (62.35%) | 195 |

**The pair-level difference mostly does not survive to membership.**
A→B drops 197 pairs, 0.55% — but 195 of those 197 sit *inside* a cluster
that survives on its other edges and change nothing at all. Only **2 are
bridges**. Net effect of adopting the reference cut:

* 193 → 192 clusters of size ≥ 2
* 2,121 → 2,118 accessions in a cluster — **3 of 3,397, 0.09%**
* largest cluster unchanged at 195
* ARI(A, B) = **0.999336**

The three accessions, and this is the honest reason to care:

| accession | under A | under B |
|---|---|---|
| `OP039478.1` | in a 59-member cluster | singleton |
| `MN789809.1` | paired with `MN650438.1` | singleton |
| `MN650438.1` | paired with `MN789809.1` | singleton |

Aggregate counts barely move, but `OP039478.1` is a single sequence
joining or not joining a 58-member transmission cluster. For a
surveillance tool that is the per-record call it exists to make, not a
rounding error. "0.09% of the cohort" and "we report this person as part
of that outbreak" are the same number read two ways.

### 12.2 A recalibrated p-threshold gets the right answer for the wrong reason

If TN93 were strictly monotone in p, some p-threshold would reproduce
the TN93 pair set exactly and we could keep the cheap in-process engine
untouched. **It is not, and none does.**

Sweeping every candidate, the best single p-threshold is **p ≤ 0.014796**,
which cuts the pair-level disagreement from 197 to **17** (8 we-link-they-
don't, 9 they-link-we-don't). It cannot do better, because the p-ranges
of the two classes *overlap*:

* p of TN93-**accepted** pairs: 0.000000 .. 0.014815
* p of TN93-**rejected** pairs: 0.014784 .. 0.014983

**73 pairs sit in [0.014784, 0.014815], a band where p alone cannot
decide.** TN93 weights purine (A↔G) and pyrimidine (C↔T) transitions
separately and against the pair's own base composition, so two pairs
with the same mismatch fraction get different corrected distances. In
p ∈ (0.0145, 0.015] the correction ranges +1.62 × 10⁻⁴ to +2.55 × 10⁻⁴ —
a 1.6× spread at essentially fixed p. That spread is the irreducible
error of *any* p-threshold, and it is why §11's "monotone transform"
framing is too generous. `compute_pairwise_distances_python`'s docstring
has been corrected accordingly.

And yet on this cohort option C reproduces option B's partition
**exactly** — ARI(B, C) = 1.000000, zero accessions differ — because all
17 of its mismatched pairs are the membership-neutral within-cluster
kind. That is luck, not a property:

* 0.014796 is a constant fitted to *this* cohort at *this* threshold. §9
  item 2 re-curates `aligned.parquet` at asm20, which changes the
  cohort, the coverage distribution and the base composition. The
  constant would have to be refitted, and nothing guarantees the next 17
  mismatches land inside clusters.
* It is unpublishable as stated. "We threshold the p-distance at
  0.014796, fitted to approximate TN93 at 0.015" invites exactly the
  reviewer question a tools paper should not invite.

### 12.3 What it would take to do it properly

**Option D: compute TN93 in-process.** Per pair the formula needs three
counts — purine transitions, pyrimidine transitions, transversions — and
the four base frequencies over that pair's overlap, in place of the
single mismatch count the engine computes now. Each is the same shape of
vectorized boolean reduction over the same candidate block, so this is a
contained change to one function: roughly 5–7× the inner-loop column
work, no change to the structure, plus a domain guard for pairs where
the logarithm's argument goes non-positive (tn93 returns −1 there).

The reason this is now worth doing, and was not before: **we can verify
it.** §11 established a working reference that agrees with us pair for
pair under matched ambiguity handling. An in-process TN93 can be checked
against `hivdrwatch/tn93:1.0.17` exactly, on this cohort, as a test.
Before the image was fixed there was nothing to check it against.

### 12.4 Recommendation

**Option D, then re-measure.** B is the convention the manuscript
claims, D is how to get B without shelling out to Docker on every run,
and the cost of getting it wrong is not the 0.09% aggregate — it is
per-record cluster membership in a surveillance tool.

Ordering against §9: D is cheap and independent of the re-curation, and
the engine change should land **before** step 5 re-runs
`pipeline/cluster` on the new vintage, so figure 3 is generated once
from the metric the paper claims rather than twice. If D is judged too
much for this paper, then B — call tn93 with `-a skip` — is correct and
available today at the cost of a Docker dependency in the pipeline path.

**C is not recommended at any price.** It happens to be exactly right
here and is not right in principle, which is the worst combination to
put in a methods section.

What is *not* on the table is leaving A in place while the manuscript
claims the HIV-TRACE convention. Either the metric changes or the
methods section does.
