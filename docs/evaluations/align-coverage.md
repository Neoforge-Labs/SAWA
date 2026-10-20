# Why the curated alignment covered less pol than Sierra interpreted

**Question (from `ligase-hivdb.md` §9 item 1):** why does
`data/interim/curated/aligned.parquet` place a large fraction of sequences
in a span of pol narrower than Stanford's Sierra interprets from the same
raw sequence?

**Answer: the minimap2 preset was set below the divergence of the data.**
`align/align_mappy.py` used `preset="asm10"`, which tolerates ~10%
divergence. HXB2 is subtype B; SAWA's cohort is A1/C/D and the CRFs, and
the pol it recovers sits at **13.05% divergence** from HXB2. minimap2
discarded every block it could not chain inside the asm10 budget. The
result was not only truncation: **50.7% of the shipped curated rows are
3,012 dashes and nothing else**, and every one of them carries
`passes_qc = True`.

The leading hypothesis going in — that `_best_hit` kept one mappy
alignment and dropped secondary blocks — **is refuted as the cause**. It
is a real defect, but a latent one: at asm10 mappy returned a single hit
for 99.4% of the corpus, so there was usually nothing for `_best_hit` to
discard. It becomes load-bearing once the preset is raised, which is why
it is fixed here too.

Scope: this is a SAWA finding, not a Ligase one. Nothing here is about
resistance *scoring*. `pipeline/dr_ingest.py` sends Sierra the raw
sequence, so the resistance calls SAWA publishes today are unaffected.
`phylo/` and `cluster/` read `aligned_pol`, and those are Paper 1
figures 3 and 4.

---

## 1. Method

- Raw sequences recovered from the local GenBank flat-file cache
  (`data/raw/genbank`, 1,450 batches) — **all 36,510** unique accessions
  in the curated table, no sampling.
- Every recovered sequence verified against the curated row by
  `sequence_sha256`, so what the aligner was handed is known rather than
  assumed. **100% match**, and the harness reproduces the shipped asm10
  span exactly for **36,510 / 36,510** rows. The pipeline is being
  measured, not approximated.
- Sierra side: the 37,709 cached `HIVDB_10.2` payloads under
  `data/raw/hivdb`. Mutation positions converted to pol-relative nt via
  `POL_REGIONS` (PR 169, RT 466, IN 2146).
- Script: `scripts/align_coverage_probe.py`. It keeps a local copy of the
  pre-fix single-block projection, so its `shipped` column is the old
  pipeline exactly and the comparison stays like-for-like after the fix
  landed in `align_mappy.py`. Its `asm10+merge` column isolates the
  preset change from the block-merge change.

      python scripts/align_coverage_probe.py scan-cache --out /tmp/raw.parquet
      python scripts/align_coverage_probe.py compare --raw /tmp/raw.parquet \
          --out /tmp/coverage.parquet
      python scripts/align_coverage_probe.py report --compared /tmp/coverage.parquet

A note on row counts: the curated table has 42,106 rows but only 36,510
unique accessions. That duplication is a separate ingest bug — §6.

## 2. The measurement

Coverage = ungapped nt placed on the 3,012-column reference.

| | shipped | asm10+merge | asm20+merge |
|---|---|---|---|
| no alignment at all | **19,164 (52.5%)** | 19,164 (52.5%) | 1,144 (3.1%) |
| coverage ≥ 800 nt | 15,506 (42.5%) | 15,515 (42.5%) | 34,372 (94.1%) |
| coverage ≥ 2,900 nt | 688 (1.9%) | 694 (1.9%) | 2,622 (7.2%) |
| mean coverage | 612 nt | 617 nt | 1,302 nt |
| median coverage | **0 nt** | 0 nt | 1,081 nt |

More than one alignment block: 224 sequences (0.6%) at asm10, 593 (1.6%)
at asm20.

Against the Sierra cache, per sequence. The middle column is the preset
change held back, so the two defects can be told apart: raising the
preset is what moves these numbers, and merging blocks barely registers.

| | shipped | asm10+merge | asm20+merge |
|---|---|---|---|
| ≥1 Sierra mutation outside the aligned span | 33,015 (90.4%) | 32,925 (90.2%) | 13,003 (35.6%) |
| …outside by more than 150 nt | 22,382 (61.3%) | 22,159 (60.7%) | 801 (2.2%) |
| a **scoring** (non-"Other") mutation outside | 8,280 (22.7%) | 8,244 (22.6%) | 291 (0.8%) |

These are higher than `ligase-hivdb.md` §6's 87.1% / 23% / 4.5% because
that sample was drawn from rows with ≥800 aligned nt, which excludes the
all-gap half of the table by construction.

## 3. Why asm10 was the wrong tolerance

Over a 1,500-sequence sample, the pol that asm20 recovers and asm10 did
not is **13.05% mismatched against HXB2** (141,075 nt newly placed,
18,411 mismatching). That is ordinary HIV-1 inter-subtype divergence —
not spurious alignment of non-pol sequence, which would sit near 75%.

Three independent lines agree that the missing sequence was real:

1. **Preset sweep.** On the `KY580541.1` cluster, asm10 gives pol
   1470..3012 (1,543 nt); asm20, `map-ont`, `map-pb` and `splice` all
   give 1..3012. `asm5` finds nothing at all.
2. **The block exists and was thrown away by a score floor, not by
   `_best_hit`.** Re-running asm10 with `min_dp_score=10` surfaces a
   *second* hit at ref 229..1232 that default asm10 never reports. The
   sequence was there; the preset's DP-score threshold refused it.
3. **A second backend disagrees.** Biopython's `PairwiseAligner` places
   **3,009 / 3,012 nt** of those same sequences where mappy placed 1,543.

Point 3 also falsifies a claim in `align/README.md`: "All three implement
`--keeplength` semantics so outputs share HXB2 coordinates and are
drop-in concatenable." They share a coordinate system, but they disagreed
about which coordinates a sequence occupies by up to 1,466 nt. The README
is corrected in this change.

## 4. The fix, and what it costs

Three changes in `align/`:

1. `DEFAULT_PRESET` → `asm20`.
2. `_best_hit` → `_usable_hits`, keeping **every** forward-strand block.
   `_build_keeplength` takes the list and lays blocks down best-first,
   writing only into columns still gapped, so the best hit wins any
   overlap and weaker blocks fill around it. 600 sequences (1.64%)
   genuinely split into more than one block under asm20, so this is not
   hypothetical.
3. Reverse-strand hits are now dropped rather than projected. The
   existing CIGAR walk indexes the forward query, but a minus-strand
   CIGAR walks the reverse complement, so projecting one wrote garbage.

Re-running the fixed align + QC path over all 36,510 sequences:

| | before | after |
|---|---|---|
| all-gap rows | 19,164 | 1,144 |
| all-gap rows with `passes_qc=True` | 19,164 | **0** |
| coverage ≥ 800 nt | 15,506 | 34,372 |
| mean / median coverage | 612 / 0 nt | 1,302 / 1,081 nt |
| sequences gaining coverage | — | 33,901 |
| sequences losing coverage | — | **1** |

**The one loss is a gain.** `OQ686534.1` is a GenBank record flagged
`UNVERIFIED:` whose pol matches HXB2 partly on each strand. The old code
took the minus-strand hit as "best" and projected it forward, placing
1,413 bases at **76.9% mismatch against HXB2** — random-sequence level —
into the curated table as `passes_qc=True`. It now contributes 815 nt of
correct plus-strand sequence instead.

## 5. QC: `SHORT_COVERAGE` was the wrong gate

`SHORT_COVERAGE` should stay non-blocking. The reasoning in `qc.py` is
sound: most GenBank HIV deposits are RT-only, and requiring PR+RT+IN
would discard most real data.

But it was carrying a case it was never meant to cover. A row with
*nothing* on the reference is not a short sequence, it is the absence of
one, and its 3,012 dashes are indistinguishable from data to any consumer
that only checks the string is non-empty. So this change adds a
**blocking** `NO_ALIGNMENT` flag, raised when not one PR, RT or IN site
is covered. After the preset fix it blocks 1,144 rows rather than 19,164.

## 6. What this changes for figures 3 and 4

**Figure 4 (DTG case study) — not affected.** It reads resistance calls,
and `dr_ingest` calls Sierra on the raw sequence.

**Figure 3 (cluster validation) — affected, and the alignment is only
half of it.** `phylo/build.py` already filters on
`MIN_ALIGNED_COVERAGE = 0.5`, so all-gap rows never reached tree
inference — but half the corpus was being dropped there silently, and
`partition_report`'s "predominantly single-gene" exclusions were partly
an artefact of the preset rather than of the data.

`pipeline/cluster.py` has no such filter. `materialize_aligned_fasta`
writes every row where `seq` is truthy, and a string of 3,012 dashes is
truthy, so all 21,358 all-gap rows were written into the FASTA handed to
the distance engine. What happens next is a property of the distance
engine, which divides by alignment width rather than by pairwise overlap
(`cluster/distance.py::compute_pairwise_distances_python`, the **default**
for the real pipeline via `--distance-engine python`):

| pair | reported distance | ≤ 0.015 threshold? |
|---|---|---|
| all-gap vs all-gap | 0.000000 | **yes — linked** |
| all-gap vs 600 nt real | 0.199 | no |
| two sequences sharing 300 nt, 10% divergent there | 0.00996 | **yes — linked** |

So the all-gap rows formed one mutually-linked block, and any two
sequences with a short overlap are linked on the strength of how little
they overlap. The `NO_ALIGNMENT` gate removes the first problem at the
source, and raising coverage shrinks the second by making overlaps
larger — but the denominator is a real methodological defect that the
alignment fix does not cure. **Figure 3's numbers should be regenerated
after that denominator is settled, not just after this change.** That is
a Paper 1 claim, so it is Kelyn's call, not a side effect of this work.

## 7. Provenance bugs found on the way — not fixed here

These are `ingest/`, not `align/`, and each is a one-line pure-function
bug with a visible corpus-wide effect. They are *why* the 11-sequence
cluster in `ligase-hivdb.md` §6 exists at all: those sequences should
never have been in an East Africa cohort.

All 14 accessions in that cluster are **near-full-genome HIV-1 from
Thailand, CRF01_AE**, 8.6–8.8 kb of gag→nef, from the RV217 study.

1. **`_normalize_subtype` maps every CRF to subtype C.** It tests
   `"SUBTYPE: C" in "SUBTYPE: CRF01_AE"`, which is true. Verified:
   `CRF01_AE → C`, `CRF02_AG → C`, `CRF35_AD → C`. The comment says
   "Order matters: check CRFs before bare letters", but the CRF check sits
   *after* the loop. Subtype-stratified anything — figure 2, the
   per-subtype trees, Paper 3's `subtype_holdout` split — is reading
   7,109 rows labelled C of which an unknown share are CRFs.
2. **`_normalize_country` fabricates geography.** A country outside
   `EAST_AFRICA_COUNTRIES` returns `None`, and `_parse` then falls back
   to `requested_country`. A Thai sequence is stored as Kenyan, Ugandan
   or Tanzanian depending on which query returned it.
3. **The same record is stored once per country query, and dedup cannot
   see it.** 3,156 accessions appear 2–3× (42,106 rows, 36,510 unique
   accessions, 35,616 unique sha256). `dedup.py` groups on
   `(country, cohort, token)`, and bug 2 gives each copy a different
   country, so each gets its own `patient_group_id`. Any
   country-stratified count double- or triple-counts these.

Also worth noting: `ingest/genbank.py::_parse` stores
`sequence=str(record.seq)` with `gene_region="pol"` hardcoded, so a
whole-genome deposit enters the pipeline as 8.7 kb labelled "pol". The
aligner is the only thing that ever restricts it to pol.

`aligned_length` deserves a mention of its own: it is `len(aligned_seq)`,
which is 3,012 for every row including the all-gap ones. It looks like a
coverage column and is not one. That is the single biggest reason this
went unnoticed for so long.

## 8. What is left

1. **The distance denominator** (`cluster/distance.py`) — divide by
   pairwise non-gap overlap, and require a minimum overlap before
   emitting a pair. Gates figure 3. Needs Kelyn: it changes a Paper 1
   number.
2. **`pipeline/cluster.py` has no coverage filter.** Even with
   `NO_ALIGNMENT` blocking, `materialize_aligned_fasta` ignores
   `passes_qc` entirely. It should respect it.
3. **The three `ingest/` bugs in §6**, and a decision on whether
   non-East-African records should be dropped rather than relabelled.
4. **Regenerate the curated vintage.** Every number above is measured
   against the shipped `aligned.parquet`; it has not been rebuilt, and
   DVC still points at the asm10 output. `phylo/` and `cluster/` outputs
   derive from it.
5. **`benchmarks/scale_align.py`** should be re-run for the record, but
   not because of speed: measured over 3,000 real sequences × 3 reps, the
   medians are 1,984 seq/s at asm10 and 1,868 seq/s at asm20, with
   overlapping ranges — no material cost, and both well above the
   ~700 seq/s the module docstring claims. Figure 5's shape does not
   change.
