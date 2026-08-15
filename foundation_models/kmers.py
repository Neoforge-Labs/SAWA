"""K-mer featurizer for classical baselines.

Two entry points:
  * `kmer_counts(seq, k)` — one sequence → counter dict, useful for
    exploratory work and Naive Bayes-style baselines.
  * `kmer_matrix(seqs, k)` — many sequences → scipy sparse CSR matrix
    plus the vocabulary array. This is what `baselines.py` feeds to
    scikit-learn.

## Design

Nucleotide-only alphabet (A, C, G, T) by default with ambiguity codes
optional. Any character outside the alphabet is silently dropped (no
`N`-containing k-mers survive). Case-folded to uppercase.

We deliberately don't hash-featurize: with N=25 000 sequences × k=6
the exact vocabulary is ~4⁶ = 4 096 columns — a sparse matrix fits
in memory trivially. Hashing would only pay off past k=8.
"""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Sequence

import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_ALPHABET = "ACGT"


def kmer_counts(
    seq: str,
    k: int,
    *,
    alphabet: str = DEFAULT_ALPHABET,
) -> Counter[str]:
    """Return the multiset of length-`k` substrings that lie entirely
    within `alphabet`. Case-folded.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    allowed = set(alphabet.upper())
    seq = seq.upper()
    counts: Counter[str] = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if all(c in allowed for c in kmer):
            counts[kmer] += 1
    return counts


def kmer_matrix(
    seqs: Sequence[str],
    k: int,
    *,
    alphabet: str = DEFAULT_ALPHABET,
    normalize: bool = True,
    vocabulary: Sequence[str] | None = None,
):
    """Vectorize a sequence set into a (n_seqs × |vocab|) sparse matrix.

    Returns `(matrix, vocabulary)` where `vocabulary` is a numpy array
    of k-mer strings in column order. If `vocabulary` is provided
    explicitly, only those columns are emitted (use this for the
    val/test transform after fitting on train).

    `normalize=True` divides each row by the sequence's total kmer
    count so long/short sequences are comparable. Off by default only
    for raw-count downstream uses.
    """
    try:
        from scipy import sparse
    except ImportError as e:  # pragma: no cover - scipy is a hard dep
        raise ImportError("scipy is required for kmer_matrix (install `scipy`)") from e

    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")

    per_seq_counts: list[Counter[str]] = [kmer_counts(s, k, alphabet=alphabet) for s in seqs]

    if vocabulary is None:
        vocab_set: set[str] = set()
        for c in per_seq_counts:
            vocab_set.update(c)
        vocab = np.array(sorted(vocab_set), dtype=object)
    else:
        vocab = np.array(list(vocabulary), dtype=object)

    kmer_to_col: dict[str, int] = {km: i for i, km in enumerate(vocab)}

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for row_ix, c in enumerate(per_seq_counts):
        row_total = sum(c.values()) if normalize else 1
        if row_total == 0:
            continue
        for kmer, n in c.items():
            col_ix = kmer_to_col.get(kmer)
            if col_ix is None:
                continue  # kmer not in explicit vocab; skip
            rows.append(row_ix)
            cols.append(col_ix)
            data.append(n / row_total if normalize else float(n))

    matrix = sparse.csr_matrix(
        (data, (rows, cols)),
        shape=(len(seqs), len(vocab)),
        dtype=np.float32,
    )
    logger.info(
        "kmer_matrix: %d seqs × %d k-mers (k=%d, alphabet=%s, normalize=%s, density=%.4f)",
        matrix.shape[0],
        matrix.shape[1],
        k,
        alphabet,
        normalize,
        matrix.nnz / max(1, matrix.shape[0] * matrix.shape[1]),
    )
    return matrix, vocab


__all__ = ["DEFAULT_ALPHABET", "kmer_counts", "kmer_matrix"]
