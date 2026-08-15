"""Unit tests for foundation_models.kmers."""

from __future__ import annotations

import numpy as np
import pytest

from foundation_models.kmers import DEFAULT_ALPHABET, kmer_counts, kmer_matrix


class TestKmerCounts:
    def test_basic(self) -> None:
        counts = kmer_counts("ACGTACGT", k=3)
        assert counts["ACG"] == 2
        assert counts["CGT"] == 2
        assert counts["GTA"] == 1
        assert counts["TAC"] == 1

    def test_case_folded(self) -> None:
        assert kmer_counts("acgt", k=2) == kmer_counts("ACGT", k=2)

    def test_ambiguity_dropped(self) -> None:
        # N is not in default alphabet — kmers containing N are excluded.
        counts = kmer_counts("ACNGT", k=3)
        assert "ACN" not in counts
        assert "CNG" not in counts
        assert "NGT" not in counts

    def test_empty_sequence(self) -> None:
        assert kmer_counts("", k=3) == {}

    def test_k_larger_than_seq(self) -> None:
        assert kmer_counts("AC", k=5) == {}

    def test_bad_k_raises(self) -> None:
        with pytest.raises(ValueError, match="k must be"):
            kmer_counts("ACGT", k=0)


class TestKmerMatrix:
    def test_shape_and_vocab(self) -> None:
        seqs = ["ACGT", "ACGT", "TTTT"]
        matrix, vocab = kmer_matrix(seqs, k=2)
        assert matrix.shape[0] == 3
        assert matrix.shape[1] == len(vocab)
        assert set(vocab) == {"AC", "CG", "GT", "TT"}

    def test_normalize_defaults_to_true(self) -> None:
        matrix, _ = kmer_matrix(["ACGT"], k=2)
        # 3 kmers total (AC, CG, GT); each 1/3 after normalize
        assert np.isclose(matrix.sum(), 1.0)

    def test_normalize_off(self) -> None:
        matrix, _ = kmer_matrix(["ACGT"], k=2, normalize=False)
        assert matrix.sum() == 3

    def test_explicit_vocab(self) -> None:
        # Fit-time vocab from ACGT sequences
        train_matrix, vocab = kmer_matrix(["ACGT", "ACGA"], k=2)
        # Transform new sequence using the same vocab
        test_matrix, _ = kmer_matrix(["TTTT"], k=2, vocabulary=vocab)
        assert test_matrix.shape[1] == train_matrix.shape[1]
        # TTTT's TT kmer isn't in the ACGT-derived vocab → all zeros
        assert test_matrix.sum() == 0

    def test_alphabet_uses_default(self) -> None:
        assert DEFAULT_ALPHABET == "ACGT"

    def test_empty_seq_list(self) -> None:
        matrix, vocab = kmer_matrix([], k=3)
        assert matrix.shape[0] == 0
        assert len(vocab) == 0
