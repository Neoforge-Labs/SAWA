"""Unit tests for the cluster validation harness.

Docker-invoked tn93 isn't tested here — that's an integration test
against the built image. This file tests the pure helpers:
  * simulate_clusters: shape + determinism + write_fasta
  * score_clusters: pair-based P/R/F1/ARI arithmetic on known cases
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hivdrwatch.cluster import (
    ClusterValidationResult,
    simulate_clusters,
    score_clusters,
)
from hivdrwatch.cluster.simulate import _mutate, _random_sequence
import numpy as np


# ---------- simulate ----------


class TestSimulateClusters:
    def test_shape(self) -> None:
        sim = simulate_clusters(
            n_clusters=5, cluster_size=(3, 3),
            n_singletons=4, seq_length=100, seed=1,
        )
        assert sim.n_clusters == 5
        assert len(sim.singletons) == 4
        # 5 clusters * 3 members + 4 singletons = 19
        assert sim.n_sequences == 5 * 3 + 4
        # Every sequence is exactly seq_length nt
        for seq in sim.accession_to_seq.values():
            assert len(seq) == 100
            assert set(seq) <= set("ACGT")

    def test_deterministic(self) -> None:
        a = simulate_clusters(
            n_clusters=3, cluster_size=(4, 4),
            n_singletons=2, seq_length=50, seed=7,
        )
        b = simulate_clusters(
            n_clusters=3, cluster_size=(4, 4),
            n_singletons=2, seq_length=50, seed=7,
        )
        assert a.accession_to_seq == b.accession_to_seq
        assert a.ground_truth == b.ground_truth
        assert a.singletons == b.singletons

    def test_different_seeds_differ(self) -> None:
        a = simulate_clusters(n_clusters=2, cluster_size=(3, 3), n_singletons=1, seed=1)
        b = simulate_clusters(n_clusters=2, cluster_size=(3, 3), n_singletons=1, seed=2)
        assert a.accession_to_seq != b.accession_to_seq

    def test_truth_membership_shape(self) -> None:
        sim = simulate_clusters(
            n_clusters=2, cluster_size=(3, 3), n_singletons=2, seed=1
        )
        truth = sim.truth_membership()
        assert len(truth) == 2 * 3 + 2
        # Members of the same simulated cluster share a truth label
        for cluster_id, accessions in sim.ground_truth.items():
            labels = {truth[a] for a in accessions}
            assert len(labels) == 1
            assert next(iter(labels)) == cluster_id
        # Every singleton has a unique label
        singleton_labels = {truth[a] for a in sim.singletons}
        assert len(singleton_labels) == len(sim.singletons)

    def test_write_fasta_roundtrip(self, tmp_path: Path) -> None:
        sim = simulate_clusters(
            n_clusters=2, cluster_size=(3, 3), n_singletons=2, seed=1
        )
        fasta = sim.write_fasta(tmp_path / "sim.fasta")
        text = fasta.read_text()
        # Every accession appears as a header line
        for accession in sim.accession_to_seq:
            assert f">{accession}\n" in text

    def test_cluster_size_int_and_tuple(self) -> None:
        # Fixed int
        a = simulate_clusters(n_clusters=3, cluster_size=5, n_singletons=0, seed=1)
        for accs in a.ground_truth.values():
            assert len(accs) == 5

        # Range tuple
        b = simulate_clusters(n_clusters=5, cluster_size=(2, 4), n_singletons=0, seed=1)
        for accs in b.ground_truth.values():
            assert 2 <= len(accs) <= 4

    def test_bad_cluster_size_raises(self) -> None:
        with pytest.raises(ValueError):
            simulate_clusters(n_clusters=1, cluster_size=1)  # too small
        with pytest.raises(ValueError):
            simulate_clusters(n_clusters=1, cluster_size=(5, 3))  # inverted range

    def test_mutation_produces_close_but_different_sequences(self) -> None:
        rng = np.random.default_rng(0)
        founder = _random_sequence(rng, 3000)
        member = _mutate(rng, founder, rate=0.005)
        # ~0.5% mutation → sequences share most positions
        diffs = sum(1 for a, b in zip(founder, member) if a != b)
        assert 5 < diffs < 40


# ---------- validate ----------


class TestScoreClusters:
    def test_perfect_recovery(self) -> None:
        truth = {"a": "C1", "b": "C1", "c": "C2", "d": "C2"}
        predicted = {"a": "X", "b": "X", "c": "Y", "d": "Y"}
        r = score_clusters(predicted, truth)
        # Same-pair set is {(a,b), (c,d)}; both correct → 100% P/R/F1
        assert r.precision == 1.0
        assert r.recall == 1.0
        assert r.f1 == 1.0
        assert r.tp == 2 and r.fp == 0 and r.fn == 0 and r.tn == 4

    def test_all_singletons(self) -> None:
        # Predicted no clusters at all — recall must be 0 since truth has clusters
        truth = {"a": "C1", "b": "C1", "c": "C2", "d": "C2"}
        predicted = {"a": "A", "b": "B", "c": "C", "d": "D"}
        r = score_clusters(predicted, truth)
        assert r.tp == 0
        assert r.recall == 0.0
        # Precision is undefined (0/0); we set it to 0 by convention
        assert r.precision == 0.0

    def test_over_merging(self) -> None:
        # Predicted lumps everything into one cluster
        truth = {"a": "C1", "b": "C1", "c": "C2", "d": "C2"}
        predicted = {"a": "X", "b": "X", "c": "X", "d": "X"}
        r = score_clusters(predicted, truth)
        # 6 pairs; 2 are same-truth → tp=2, fp=4
        assert r.tp == 2
        assert r.fp == 4
        assert r.precision == pytest.approx(2 / 6)
        assert r.recall == 1.0

    def test_partial_recovery(self) -> None:
        # a & b caught; c misassigned to a new cluster
        truth = {"a": "C1", "b": "C1", "c": "C1", "d": "C2"}
        predicted = {"a": "X", "b": "X", "c": "Y", "d": "Y"}
        r = score_clusters(predicted, truth)
        # Pairs: (a,b) same_pred=T same_true=T → tp
        #        (a,c) same_pred=F same_true=T → fn
        #        (a,d) same_pred=F same_true=F → tn
        #        (b,c) same_pred=F same_true=T → fn
        #        (b,d) same_pred=F same_true=F → tn
        #        (c,d) same_pred=T same_true=F → fp
        assert r.tp == 1
        assert r.fp == 1
        assert r.fn == 2
        assert r.tn == 2
        assert r.precision == 0.5
        assert r.recall == pytest.approx(1 / 3)

    def test_returns_ari(self) -> None:
        truth = {f"s{i}": "C1" if i < 5 else "C2" for i in range(10)}
        predicted = dict(truth)  # perfect
        r = score_clusters(predicted, truth)
        assert r.adjusted_rand_index == pytest.approx(1.0)

    def test_end_to_end_simulator_then_scorer(self) -> None:
        # Score the simulator's truth against itself — perfect result
        sim = simulate_clusters(
            n_clusters=5, cluster_size=(4, 4), n_singletons=5, seed=42
        )
        truth = sim.truth_membership()
        r = score_clusters(truth, truth)
        assert r.f1 == 1.0
        assert r.adjusted_rand_index == pytest.approx(1.0)

    def test_no_overlap_raises(self) -> None:
        with pytest.raises(ValueError, match="no overlap"):
            score_clusters({"a": "X"}, {"b": "Y"})


class TestClusterValidationResult:
    def test_as_dict(self) -> None:
        r = ClusterValidationResult(
            n_sequences=10, n_pairs=45,
            tp=10, fp=2, fn=3, tn=30,
            precision=10/12, recall=10/13, f1=0.75,
            adjusted_rand_index=0.8,
            n_predicted_clusters=3, n_truth_clusters=4,
        )
        d = r.as_dict()
        assert d["precision"] == 10/12
        assert d["adjusted_rand_index"] == 0.8
