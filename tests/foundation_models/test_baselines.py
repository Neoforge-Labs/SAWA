"""Unit tests for foundation_models.baselines."""

from __future__ import annotations

import numpy as np
import pytest

from foundation_models.baselines import KmerLogisticRegression, OneHotRandomForest


def _binary_toy(
    n_per_class: int = 20, seq_len: int = 100, seed: int = 42
) -> tuple[list[str], list[int]]:
    """Toy dataset where class-1 sequences carry a signature motif."""
    rng = np.random.default_rng(seed)
    alphabet = np.array(list("ACGT"))
    seqs: list[str] = []
    labels: list[int] = []
    for _ in range(n_per_class):
        # Class 0: purely random
        seqs.append("".join(alphabet[rng.integers(0, 4, size=seq_len)]))
        labels.append(0)
    for _ in range(n_per_class):
        # Class 1: same random background but with 'AAAAAAAA' inserted
        # at a fixed position — a signal both baselines should catch.
        s = list(alphabet[rng.integers(0, 4, size=seq_len)])
        s[20:28] = list("AAAAAAAA")
        seqs.append("".join(s))
        labels.append(1)
    return seqs, labels


class TestKmerLogisticRegression:
    def test_fit_and_predict_recovers_signal(self) -> None:
        seqs, labels = _binary_toy()
        model = KmerLogisticRegression(k=3, C=1.0, seed=42).fit(seqs, labels)
        probs = model.predict_proba(seqs)
        assert probs.shape == (len(seqs),)
        # Class-1 sequences should get higher probability on avg
        avg_class1 = float(probs[len(seqs) // 2 :].mean())
        avg_class0 = float(probs[: len(seqs) // 2].mean())
        assert avg_class1 > avg_class0

    def test_predict_binary(self) -> None:
        seqs, labels = _binary_toy()
        model = KmerLogisticRegression(k=3, seed=42).fit(seqs, labels)
        preds = model.predict(seqs)
        # Should get most training-set predictions right
        acc = float((preds == np.asarray(labels)).mean())
        assert acc > 0.85, f"training accuracy suspiciously low: {acc}"

    def test_predict_before_fit_raises(self) -> None:
        with pytest.raises(RuntimeError, match="not fitted"):
            KmerLogisticRegression().predict_proba(["ACGT"])

    def test_vocab_locked_at_fit_time(self) -> None:
        train, y = _binary_toy(n_per_class=10)
        model = KmerLogisticRegression(k=3, seed=0).fit(train, y)
        # New sequences with novel k-mers still transform to fit-time width
        novel = ["N" * 100]  # entirely N → no valid k-mers
        probs = model.predict_proba(novel)
        assert probs.shape == (1,)


class TestOneHotRandomForest:
    def test_fit_and_predict_recovers_signal(self) -> None:
        seqs, labels = _binary_toy(n_per_class=15)
        model = OneHotRandomForest(n_estimators=50, seed=42).fit(seqs, labels)
        probs = model.predict_proba(seqs)
        avg_class1 = float(probs[len(seqs) // 2 :].mean())
        avg_class0 = float(probs[: len(seqs) // 2].mean())
        assert avg_class1 > avg_class0

    def test_mixed_length_raises(self) -> None:
        with pytest.raises(ValueError, match="equal-length"):
            OneHotRandomForest().fit(["ACGT", "ACG"], [0, 1])

    def test_predict_before_fit_raises(self) -> None:
        with pytest.raises(RuntimeError, match="not fitted"):
            OneHotRandomForest().predict_proba(["ACGT"])

    def test_predict_length_mismatch_raises(self) -> None:
        seqs = ["ACGT" * 25, "TGCA" * 25]
        model = OneHotRandomForest(n_estimators=5, seed=0).fit(seqs, [0, 1])
        with pytest.raises(ValueError, match="!= fit seq length"):
            model.predict_proba(["ACGT" * 20])  # wrong length
