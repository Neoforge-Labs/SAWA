"""Non-FM baselines for the Paper 3 benchmark.

Two shipped baselines:

  * `KmerLogisticRegression` — k-mer counts → L2-regularized logistic
    regression. Fast (~seconds for 25 k sequences), interpretable,
    strong on tasks where composition dominates position.
  * `OneHotRandomForest` — per-position one-hot encoding of the
    aligned sequence → random forest. Captures positional signal
    that k-mers wash out.

Both wrap scikit-learn (already a project dep). Both fit on train,
predict on val + test. Model + featurizer state is kept in the class
so `predict_proba(new_seqs)` works after fit.

## Fit / predict contract

    model = KmerLogisticRegression(k=3, C=1.0, seed=42)
    model.fit(train_seqs, train_labels)
    val_probs = model.predict_proba(val_seqs)
    test_probs = model.predict_proba(test_seqs)

Each `predict_proba` returns a 1-D numpy array of P(class=1). The
`evaluate.compute_metrics(y_true, y_prob)` helper then produces the
per-drug × per-subtype metric row.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .kmers import DEFAULT_ALPHABET, kmer_matrix

logger = logging.getLogger(__name__)


@dataclass
class KmerLogisticRegression:
    """L2-regularized logistic regression on k-mer counts."""

    k: int = 3
    C: float = 1.0
    max_iter: int = 500
    alphabet: str = DEFAULT_ALPHABET
    seed: int = 42
    # Fit-time state
    _model: object | None = field(default=None, init=False, repr=False)
    _vocab: np.ndarray | None = field(default=None, init=False, repr=False)

    def fit(self, seqs: Sequence[str], y: Sequence[int]) -> KmerLogisticRegression:
        from sklearn.linear_model import LogisticRegression

        x_train, vocab = kmer_matrix(seqs, self.k, alphabet=self.alphabet)
        self._vocab = vocab
        self._model = LogisticRegression(
            C=self.C,
            max_iter=self.max_iter,
            solver="liblinear",  # good for sparse binary problems
            random_state=self.seed,
            n_jobs=1,
        )
        self._model.fit(x_train, np.asarray(y))
        logger.info(
            "KmerLogisticRegression fitted: %d train seqs × %d k-mers (k=%d, C=%.2f)",
            x_train.shape[0],
            x_train.shape[1],
            self.k,
            self.C,
        )
        return self

    def predict_proba(self, seqs: Sequence[str]) -> np.ndarray:
        if self._model is None or self._vocab is None:
            raise RuntimeError("model not fitted; call .fit() first")
        x, _ = kmer_matrix(
            seqs,
            self.k,
            alphabet=self.alphabet,
            vocabulary=self._vocab,
        )
        return self._model.predict_proba(x)[:, 1]

    def predict(self, seqs: Sequence[str], *, threshold: float = 0.5) -> np.ndarray:
        return (self.predict_proba(seqs) >= threshold).astype(np.int8)


@dataclass
class OneHotRandomForest:
    """Per-position one-hot encoded aligned sequence → random forest.

    All sequences MUST be the same length (aligned). Ambiguous / gap
    positions become an "other" channel; this keeps the shape stable.
    """

    n_estimators: int = 200
    max_depth: int | None = None
    alphabet: str = DEFAULT_ALPHABET
    n_jobs: int = -1
    seed: int = 42
    _model: object | None = field(default=None, init=False, repr=False)
    _seq_len: int | None = field(default=None, init=False, repr=False)

    def fit(self, seqs: Sequence[str], y: Sequence[int]) -> OneHotRandomForest:
        from sklearn.ensemble import RandomForestClassifier

        x_train, seq_len = self._encode(seqs, expected_len=None)
        self._seq_len = seq_len
        self._model = RandomForestClassifier(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            n_jobs=self.n_jobs,
            random_state=self.seed,
        )
        self._model.fit(x_train, np.asarray(y))
        logger.info(
            "OneHotRandomForest fitted: %d seqs × %d features (seq_len=%d, n_estimators=%d)",
            x_train.shape[0],
            x_train.shape[1],
            seq_len,
            self.n_estimators,
        )
        return self

    def predict_proba(self, seqs: Sequence[str]) -> np.ndarray:
        if self._model is None or self._seq_len is None:
            raise RuntimeError("model not fitted; call .fit() first")
        x, _ = self._encode(seqs, expected_len=self._seq_len)
        return self._model.predict_proba(x)[:, 1]

    def predict(self, seqs: Sequence[str], *, threshold: float = 0.5) -> np.ndarray:
        return (self.predict_proba(seqs) >= threshold).astype(np.int8)

    # ---- internal ----

    def _encode(
        self, seqs: Sequence[str], *, expected_len: int | None
    ) -> tuple[np.ndarray, int]:
        """Return (n_seqs × (seq_len * n_channels)) dense uint8 matrix."""
        upper_seqs = [s.upper() for s in seqs]
        lengths = {len(s) for s in upper_seqs}
        if len(lengths) != 1:
            raise ValueError(
                f"OneHotRandomForest requires aligned (equal-length) sequences; "
                f"got lengths {sorted(lengths)}"
            )
        seq_len = lengths.pop()
        if expected_len is not None and seq_len != expected_len:
            raise ValueError(f"predict seq length {seq_len} != fit seq length {expected_len}")

        channels = [*self.alphabet.upper(), "OTHER"]
        n_channels = len(channels)
        idx_of = {ch: i for i, ch in enumerate(channels)}
        other_ix = idx_of["OTHER"]

        x = np.zeros((len(upper_seqs), seq_len * n_channels), dtype=np.uint8)
        for row, s in enumerate(upper_seqs):
            for pos, ch in enumerate(s):
                col = pos * n_channels + idx_of.get(ch, other_ix)
                x[row, col] = 1
        return x, seq_len


__all__ = ["KmerLogisticRegression", "OneHotRandomForest"]
