# `align` — pairwise-to-HXB2, trimming, QC, dedup

Aligns every input sequence to the HXB2 K03455 *pol* reference, trims to PR + RT + IN, runs quality control, and deduplicates by patient hash before handing off to phylogenetics and clustering.

## Backends

| Aligner | Throughput | When to use |
|---|---|---|
| **`mappy`** *(default)* | ~700 seq/s single-thread, in-process | Everything unless you have a specific reason not to. |
| `biopython` | ~5 seq/s | Pure-Python fallback for environments where `mappy` won't install. |
| `mafft` (Docker) | ~50 seq/s, chunked | When you need joint-guide-tree MSA across newly added sequences. Rare for reference-based work. |

Dispatch is via `align_sequences(..., aligner="mappy" | "biopython" | "mafft")`. All three implement `--keeplength` semantics so outputs share HXB2 coordinates and are drop-in concatenable.

## Trim and QC

After alignment we trim to the standard pol subregions using HXB2 pol-relative coordinates from `reference.py`:

- PR — 169..465 (99 codons)
- RT — 466..1785 (440 codons, polymerase domain)
- IN — 2146..3012 (288 codons)

QC (`qc.py`) applies flags per sequence. **Blocking** flags (sequence dropped): `HIGH_AMBIGUITY`, `APOBEC_HYPERMUT`, `STOP_CODON`, `LOW_COMPLEXITY`. **Non-blocking** flags (retained with warning): `SHORT_COVERAGE`.

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
