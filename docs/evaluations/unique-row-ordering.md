# Polars `unique()` row ordering, and which of the 22 call sites it reaches

**Scope.** Row ordering only. Nothing here changes which rows survive a
dedup, and nothing here touches the subtype or country normalisation that
`ingest-provenance.md` reserves for Kelyn. The 22 `.unique(` call sites
outside `src/hivdrwatch/dashboard/` were each read and classified; five
sites were changed.

**One-line summary.** `df.unique()` returns rows in a different order on
every call unless told otherwise, and at four places that order escaped:
into the ingest Parquet whose md5 DVC cites, into the tree manifest, into
the causal frame that a seeded bootstrap resamples *by position*, and into
Paper 3's `subtype_holdout` split — where **seed 42 selected a different
set of training sequences on all 8 of 8 runs**. A fixed seed does not
reproduce a split whose input order moves.

> **This is a row-ordering defect, not a data-content defect.** The
> retained *row* was measured separately and never varied: for both
> `keep="last"` and `keep="first"`, over 12 runs, the set of surviving
> (key, payload) pairs was identical every time (§2). `keep=` picks
> correctly with respect to the input order whether or not
> `maintain_order` is set. Nothing was selecting the wrong row, no stored
> value is wrong because of this, and no figure's *numbers* move. Please
> do not go hunting for a corruption bug that does not exist.

---

## 1. Method

- Everything below is a count of distinct outputs over repeated calls on
  one unchanged input, measured on **polars 1.43.2 / CPython 3.13.12 /
  macOS arm64, 8 cores**. Nothing is inferred from the documentation.
- The reordering comes from Polars' parallel hash path, so it is
  **sensitive to frame size and to the host**, and small frames are not a
  test. The 6-row/3-key case below gave 4, 6, 5, 6 and 5 distinct orders
  over five repeats of 12 runs here — unstable, but far from the 12/12 a
  larger frame gives — while the session that first reported this issue
  saw a 3-row case pass cleanly with the fix removed. Neither result is
  something to build a regression test on, so every fixture in the suite
  below uses at least 11 distinct keys and all of them give 8/8.
- Each fix was mutation-tested: the flag or sort was removed, the test was
  confirmed to fail, and the fix was restored. Counts for the removed
  state are quoted throughout as the "before" column.

Reproduce:

```bash
python scripts/unique_ordering_probe.py polars --runs 12
python scripts/unique_ordering_probe.py sites  --runs 8
pytest tests/ingest/test_row_ordering.py tests/phylo/test_row_ordering.py \
       tests/causal/test_row_ordering.py \
       tests/foundation_models/test_row_ordering.py -q
```

---

## 2. What the primitive actually does

12 runs each, on `n` accessions repeated so every key has a duplicate.

| call | distinct row orders |
|---|---|
| 12 rows / 6 keys, `df.unique()` | **12** |
| 12 rows / 6 keys, `+ maintain_order=True` | 1 |
| 120 rows / 40 keys, `unique(subset, keep="last")` | **12** |
| 120 rows / 40 keys, `unique(subset, keep="first")` | **12** |
| 120 rows / 40 keys, `+ maintain_order=True` | 1 |
| 6 rows / 3 keys, `unique(subset, keep="last")` | 4–6 across repeats — *small frames are not a test* |
| `Series.unique()` | **12** |
| `Series.unique(maintain_order=True)` | 1 |

And the negative result, measured by giving each duplicate a
distinguishable payload and comparing the surviving (key, payload) pairs
as a set:

| call | distinct retained-value sets in 12 runs |
|---|---|
| `unique(subset, keep="last")` | **1** |
| `unique(subset, keep="first")` | **1** |

### Two things found while measuring, that the audit depended on

**`group_by()` has the same defect.** `df.group_by(col).agg(...)` produced
12 distinct orders in 12 runs; `maintain_order=True` gives 1. This is not
in scope and nothing here changes a `group_by`, but it is why two of the
fixes below are an explicit `.sort()` rather than a flag: a `group_by`
upstream of a fixed `unique()` would reintroduce the problem.

**A join does not scramble a deterministic input — but it propagates a
scrambled one, and which side it follows depends on `how`.** This matters
because several sites feed a `unique()` into a join:

| join | distinct orders in 12 runs |
|---|---|
| ordered L × ordered R, `how="inner"` | 1 |
| ordered L × ordered R, `how="left"` | 1 |
| ordered L × **scrambled** R, `how="inner"` | **12** |
| ordered L × **scrambled** R, `how="left"` | 1 |
| **scrambled** L × ordered R, `how="inner"` | 1 |

So an `inner` join's output order follows its **right** frame and a `left`
join's follows its **left**. That is an observed default, not a promise —
Polars documents the default as "might differ across Polars versions or
even between different runs" — which is the argument for sorting at the
point of output rather than pinning each operation along the way.

### Why it matters at the two places it escapes into something committed

| | distinct results in 12 runs |
|---|---|
| md5 of a Parquet written from `unique(...)`, zstd-6 | **12** |
| same, with `maintain_order=True` | 1 |
| `seed=42` train set drawn by position from the deduped frame | **12** |
| same, with `maintain_order=True` | 1 |

---

## 3. Per-site findings

22 sites outside `src/hivdrwatch/dashboard/`. One
(`ingest/merge.py:76`) already passed `maintain_order=True`. Of the
remaining 21: **15 benign, 1 numpy, 5 changed.**

### Changed

| site | what escapes | fix |
|---|---|---|
| `ingest/lanl.py:165` | the row order of `east_africa_pol.parquet` on the LANL-only path | `maintain_order=True` |
| `phylo/build.py:143` | the row order of `tree_manifest.parquet` | `sorted(...)` |
| `causal/data_prep.py:126,230,249,259` | the row order of each causal frame, which `_bootstrap_ci` resamples by position | one `.sort("accession")` downstream |
| `causal/estimate.py:141` | the row order of `causal_effects.parquet` | `sorted(...)` |
| `causal/report.py:55` | the section order of `causal_report.html` | `maintain_order=True` |
| `foundation_models/data.py:155,167` | the *membership* of Paper 3's holdout split | one `.sort(["accession","drug"])` downstream |
| `causal/flows/causal.py:75` | the row order of `causal_effects.parquet` | `sorted(...)` — **untested, see §6** |

Measured with `scripts/unique_ordering_probe.py sites --runs 8`:

| site | before | after |
|---|---|---|
| `ingest/lanl.py` — `load_directory` row order | **8** distinct orders | 1 |
| `ingest/merge.py` — merged Parquet md5 | 1 *(flag already present)* | 1 |
| `phylo/build.py` — partition key order | **8** distinct orders | 1 |
| `causal/data_prep.py` — frame row order | **8** distinct orders | 1 |
| `foundation_models/data.py` — `seed=42` train set | **8** distinct memberships | 1 |

#### `ingest/lanl.py:165` — the DVC hash

`load_directory` merges each country's LANL export and dedups on
accession; `merge_sources` then either promotes that frame directly
(`_promote_lanl`, filter + `with_columns` only, so order-preserving) or
left-joins GenBank against it. On the **LANL-only path** — which is the
path the shipped `east_africa_pol.parquet` was built by, per
`ingest-provenance.md` §7 — the dedup's order reaches the file unchanged,
and the file's md5 therefore changed on every re-run of identical input.

On the GenBank path the frame is the *right* operand of a `how="left"`
join, so by §2 its order never reached the output anyway; what pins that
path is the pre-existing `maintain_order=True` in
`_collapse_duplicate_accessions`, the last operation before the write.
Removing *that* flag makes the merged Parquet produce 8 distinct md5s in 8
runs, which `tests/ingest/test_row_ordering.py` now guards. No
`maintain_order` was added to the join itself: it is not measurably
load-bearing today, and adding it would be the noise this audit was asked
to avoid. It is worth revisiting only if Polars changes the default.

#### `foundation_models/data.py:155,167` — the strongest finding

`_subtype_holdout_split` draws `rng.permutation(train_pool.height)` from a
fixed seed and gathers rows **by position**. Positions only mean the same
thing across runs if the pool arrives in the same order, and it did not:
`labeled` and `covs` are each deduped, then inner-joined — and by §2 an
inner join follows its right operand, which is the deduped `covs`. With
the fix removed, `seed=42` produced **8 distinct train sets in 8 runs** on
a 60-accession cohort. The seed was reported in the logs and the results
JSON throughout.

Three things about this are worth stating precisely:

- **The `stratified` strategy was never affected.** `_stratified_split`
  already sorts its input (`df.sort(sort_cols)`) with a comment saying
  exactly why — "Polars' `group_by` iteration order is not guaranteed
  stable". With the fix removed, the stratified tests still pass. Only
  `subtype_holdout` — Paper 3's thesis split — was broken.
- **The fix is a sort, not a flag.** `maintain_order=True` on the two
  `unique()` calls would be discarded by the join between them and the
  split. `joined` is unique on `(accession, drug)`, so one
  `.sort(["accession","drug"])` is a total order and is immune to the
  dedups, the join and any future `group_by` alike. It also gives
  `_subtype_holdout_split` the same footing `_stratified_split` already
  had.
- **No stored label is wrong.** The `keep="last"` dedup picked the correct
  Sierra row throughout (§2). What moved is which correctly-labelled
  sequences landed in train versus val.

#### `causal/data_prep.py:126,230,249,259` — a seeded bootstrap over a moving frame

All four feed, through two inner joins, the per-scope frames that
`estimate.py` hands to DoWhy. An OLS point estimate is order-invariant, so
the *estimates* are safe — but `_bootstrap_ci` does
`pandas_df.iloc[rng.integers(0, n, size=n)]`, which is positional. A frame
arriving in a different order gives the same seed a different resample,
and therefore different percentile CIs and p-values in
`causal_effects.parquet`. One `.sort("accession")` on the joined frame
fixes all four sites; none of the four needed a flag of its own.

#### `causal/estimate.py:141`, `flows/causal.py:75`, `report.py:55`

Each iterates a `unique()` list and appends in that order to something
committed: the per-subtype loop and the per-drug loop both extend the
effects list that becomes `causal_effects.parquet`, and `report.py`'s
scopes become the `<h2>` order of `causal_report.html`.

The two loops are **sorted** rather than order-preserved: subtypes and
drugs have no meaningful frame order, and sorting is additionally stable
across *vintages*, which `maintain_order=True` would not be. `report.py`
is the exception and keeps `maintain_order=True` — now that the estimator
emits effects in a deterministic order, first-appearance order groups each
drug's scopes with its own per-subtype strata, where sorting would
interleave the drugs.

#### `phylo/build.py:143`

`build_subtype_trees` iterates `partition_by_subtype`'s dict and extends
one manifest list, so the dict's insertion order is the row order of
`tree_manifest.parquet`. Sorted, not merely pinned — which also matches
`build_partition_report` two hundred lines below, whose docstring already
says it sorts "so downstream diffs are stable". The per-subtype tree
*files* were never affected; each is written to its own directory keyed by
subtype.

### Benign — confirmed by reading, and three corrections to the triage

| site | why order cannot escape |
|---|---|
| `paper1/fig02_aggregate.py:187` | `tuple(sorted(...))` |
| `phylo/build.py:378` | `sorted(...)` |
| `foundation_models/data.py:307` | `sorted(...)`, inside an error message |
| `foundation_models/data.py:323` | `sorted(...)`, inside a log line |
| `foundation_models/compare_baselines.py:191` | `sorted(...)` |
| `causal/data_prep.py:132` | `sorted(...)` |
| `cluster/assembly.py:99` | feeds `group_by().len().sort("cluster_size")` — aggregated, then sorted |
| `foundation_models/evaluate.py:79` | `np.unique`, not Polars; returns sorted, and only its length is used |

Three sites were expected to be problems and are not:

- **`phylo/auspice.py:141` is benign.** The result is immediately wrapped
  in `set(...)` and used only for membership tests. Order cannot escape a
  set.
- **`cluster/features.py:244` is benign, twice over.** Its frame is the
  *left* operand of a chain of `how="left"` joins, so by §2 it would set
  the output order — except that `build_outbreak_features` ends with
  `.select(sorted(...)).sort(["country","year"])`, and `(country, year)`
  is unique in that frame, so the sort is total. `outbreak_features.parquet`
  was already byte-stable.
- **`ingest/merge.py:76` needed nothing** — it already had the flag, and
  the audit's contribution there is a test that keeps it.

---

## 4. Regression tests

One test module per affected package, each mutation-tested. Fixtures use
≥ 11 distinct keys and 8 repeats, because small frames pass by luck (§1).

| file | guards | failure message when the fix is removed |
|---|---|---|
| `tests/ingest/test_row_ordering.py` | `lanl.py:165`, `merge.py:76` | "load_directory returned 8 distinct row orders in 8 calls"; "merge_sources produced 8 distinct Parquet digests in 8 runs on identical input" |
| `tests/phylo/test_row_ordering.py` | `build.py:143` | "partition_by_subtype returned 8 distinct key orders in 8 calls" |
| `tests/causal/test_row_ordering.py` | `data_prep.py`, `estimate.py:141`, `report.py:55` | "build_dr_causal_frames produced 8 distinct row orders in 8 builds"; "per-subtype estimation visited subtypes in 8 distinct orders"; "the HTML report emitted 8 distinct section orders" |
| `tests/foundation_models/test_row_ordering.py` | `data.py` split | "seed=42 produced 8 distinct train sets in 8 runs" |

Each module also asserts the *content* is unchanged — one row per
accession, every eligible subtype still partitioned, train and test
subtypes still disjoint — so a future "fix" that stabilises the order by
dropping rows fails too.

---

## 5. Benchmarks

Two benchmarks exercise a changed site.

**`benchmarks/scale_ingest.py`** (→ `LANLIngester.load_directory`). Seven
runs at n=5,000 before and after:

| | median | min | max |
|---|---|---|---|
| before | 49,078 seq/s | 25,504 | 52,574 |
| after | 47,828 seq/s | 27,534 | 56,686 |

Median delta **−2.5%**, against a **55%** spread within the unchanged
baseline alone. The benchmark does ~0.12 s of measured work at this size
and cannot resolve a difference of that magnitude. No regression is
claimed in either direction and `benchmarks/results/scale_ingest.json` is
left untouched.

**`benchmarks/paper3_baselines_synth.py`** (→ `build_fm_dataset`)
**does not run, before or after.** On a pristine checkout it raises
`ValueError: no rows survived the (accession, drug) join`. That is
pre-existing and unrelated to this change — verified by running it on the
unmodified tree — but it means its committed
`results/paper3_baselines_synth.json` cannot currently be regenerated, and
no before/after comparison of its numbers is possible. Note also that
there is no stable "before" number to regress *from*: the split membership
it measured varied run to run (§3), so whatever is in that JSON is one
draw out of many. See §6.

---

## 6. What is left

1. **`causal/flows/causal.py` has not imported since the Phase-0 package
   move.** `from ..causal import (...)` inside `causal.flows.causal`
   resolves to `causal.causal`, which does not exist; `from ..tracking
   import tracked_run` has the same problem and the module it wants is
   `hivdrwatch.tracking`. Both were written for the old
   `src/hivdrwatch/causal/` location that MODULES.md records moving. The
   one-line ordering fix at line 75 is correct but **cannot be tested
   today**, and nothing else in the suite imports this module, which is
   why it has gone unnoticed. MODULES.md marks `causal/flows/causal.py`
   as ✓ shipped. Left alone deliberately: the right import for `tracking`
   is a cross-package decision, not a row-ordering one.
2. **`benchmarks/paper3_baselines_synth.py` is broken** (§5), so Paper 3's
   committed baseline numbers cannot be regenerated.
3. **Kelyn's call: do not regenerate the ingest vintage for this.** The
   recommendation is to leave it. The data content is unchanged — this
   fix moves rows, it does not change them — so regenerating would churn
   a DVC hash that may already be cited, in exchange for nothing a reader
   can see. The right moment to pick up the stable ordering is the re-run
   that `ingest-provenance.md` §8 already requires for the subtype and
   country fixes, which has to happen anyway.
4. **`group_by()` is the same defect and has not been audited** (§2). It
   is used far more widely than `unique()`. The two sorts added here make
   the causal and FM paths immune to it, but nothing else is. Worth its
   own pass.
5. **The join-order defaults relied on in §2 are observed, not
   promised.** Polars documents them as liable to change between
   releases. The two explicit sorts are safe; the three
   `maintain_order=True` flags and the two `sorted()` calls are safe; what
   would need rechecking on a Polars upgrade is the "benign" verdict for
   `cluster/features.py:244`, which rests on a total sort and so is also
   safe, and for `ingest/merge.py`'s join, which rests on left-join-
   follows-left and is **not**. A Polars bump should re-run
   `scripts/unique_ordering_probe.py polars`.
