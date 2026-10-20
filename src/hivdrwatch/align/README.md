# `align` — pairwise-to-HXB2, trimming, QC, dedup

Aligns every input sequence to the HXB2 K03455 *pol* reference, trims to PR + RT + IN, runs quality control, and deduplicates by patient hash before handing off to phylogenetics and clustering.

## Backends

| Aligner | Throughput | When to use |
|---|---|---|
| **`mappy`** *(default)* | ~700 seq/s single-thread, in-process | Everything unless you have a specific reason not to. |
| `biopython` | ~5 seq/s | Pure-Python fallback for environments where `mappy` won't install. |
| `mafft` (Docker) | ~50 seq/s, chunked | When you need joint-guide-tree MSA across newly added sequences. Rare for reference-based work. |

Dispatch is via `align_sequences(..., aligner="mappy" | "biopython" | "mafft")`. All three implement `--keeplength` semantics, so every output is 3,012 columns in HXB2 coordinates and is drop-in concatenable.

Sharing a coordinate system is not the same as agreeing on it. The backends are not interchangeable to the base, and they were measurably far apart until recently: at `preset="asm10"` mappy placed non-B pol across 1,543 columns where Biopython placed 3,009. The default preset is now `asm20`, which closes that gap. If you change the preset, re-check it against Biopython on a divergent sequence first — `tests/align/test_align_mappy.py` does exactly that, and `docs/evaluations/align-coverage.md` has the measurements.

## Trim and QC

After alignment we trim to the standard pol subregions using HXB2 pol-relative coordinates from `reference.py`:

- PR — 169..465 (99 codons)
- RT — 466..1785 (440 codons, polymerase domain)
- IN — 2146..3012 (288 codons)

QC (`qc.py`) applies flags per sequence. **Blocking** flags (sequence dropped): `HIGH_AMBIGUITY`, `APOBEC_HYPERMUT`, `STOP_CODON`, `LOW_COMPLEXITY`, `NO_ALIGNMENT`. **Non-blocking** flags (retained with warning): `SHORT_COVERAGE`.

`SHORT_COVERAGE` stays non-blocking on purpose — most GenBank HIV deposits are RT-only, and demanding PR+RT+IN would throw away most of the corpus. `NO_ALIGNMENT` blocks: it fires when not one PR, RT or IN site is covered, which is not a short sequence but the absence of one. Watch for it in a run's flag counts; a rising count means the aligner is refusing inputs, not that the deposits got shorter.

**`aligned_length` is the padded width, not the coverage.** It is 3,012 for every row, including an all-gap one. Use the non-gap count (`phylo.build.filter_by_coverage`) when you mean coverage.

Dedup (`dedup.py`) collapses to one record per patient hash where the hash is present; otherwise per accession.

## Run it

```bash
python -m hivdrwatch.pipeline.curate \
    --input data/interim/east_africa_pol.parquet \
    --lanl-dir data/raw/lanl \
    --work-dir data/interim/curate \
    --out data/interim/curated/aligned.parquet
```

Output: silver Parquet with `(accession, aligned_sequence, subtype, country, collection_year, qc_flags, passes_qc, patient_hash, schema_version)` and a FASTA of QC-passing aligned sequences alongside.

## Tests

```bash
pytest tests/align -v
```
