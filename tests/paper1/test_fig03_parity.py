"""Unit tests for the Fig 3 cluster-parity aggregator."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hivdrwatch.paper1.fig03_parity import (
    CLUSTERER_ORDER,
    FIG3_METRICS,
    confusion_frame,
    load_parity_payload,
    metrics_frame,
    summarize,
)


def _make_payload(
    with_hivtrace: bool = True,
    ours_f1: float = 0.95,
    hivtrace_f1: float = 0.93,
    agreement: float = 0.99,
) -> dict:
    def _mk_metrics(f1: float) -> dict:
        return {
            "n_sequences": 200,
            "n_pairs": 19900,
            "tp": int(f1 * 100),
            "fp": 5,
            "fn": 3,
            "tn": 19800,
            "precision": 0.98,
            "recall": 0.97,
            "f1": f1,
            "adjusted_rand_index": 0.94,
            "n_predicted_clusters": 20,
            "n_truth_clusters": 20,
        }

    results = [
        {
            "clusterer": "tn93+connected_components",
            "threshold": 0.015,
            "n_predicted_clusters": 20,
            "metrics": _mk_metrics(ours_f1),
        }
    ]
    parity = None
    if with_hivtrace:
        results.append(
            {
                "clusterer": "hivtrace-hivnetworkcsv",
                "threshold": 0.015,
                "n_predicted_clusters": 20,
                "metrics": _mk_metrics(hivtrace_f1),
            }
        )
        parity = {
            "n_pairs": 19900,
            "n_agree": int(19900 * agreement),
            "fraction_agree": agreement,
        }

    return {
        "clusterer": "tn93+connected_components",
        "threshold": 0.015,
        "metrics": _mk_metrics(ours_f1),
        "results": results,
        "parity": parity,
        "sim": {
            "n_clusters": 20,
            "n_sequences": 200,
            "n_singletons": 50,
            "seq_length": 3012,
            "seed": 42,
        },
        "run_date": "2026-08-13",
    }


class TestLoadPayload:
    def test_new_format_passthrough(self, tmp_path: Path) -> None:
        p = tmp_path / "cv.json"
        p.write_text(json.dumps(_make_payload()))
        payload = load_parity_payload(p)
        assert "results" in payload
        assert len(payload["results"]) == 2

    def test_legacy_format_backfilled(self, tmp_path: Path) -> None:
        legacy = {
            "clusterer": "tn93+connected_components",
            "threshold": 0.015,
            "metrics": {
                "precision": 0.9,
                "recall": 0.9,
                "f1": 0.9,
                "adjusted_rand_index": 0.8,
                "tp": 10,
                "fp": 1,
                "fn": 1,
                "tn": 100,
                "n_predicted_clusters": 5,
            },
            "sim": {
                "n_clusters": 5,
                "n_sequences": 50,
                "n_singletons": 10,
                "seq_length": 3012,
                "seed": 42,
            },
        }
        p = tmp_path / "cv.json"
        p.write_text(json.dumps(legacy))
        payload = load_parity_payload(p)
        assert "results" in payload
        assert len(payload["results"]) == 1
        assert payload["parity"] is None
        assert payload["results"][0]["clusterer"] == "tn93+connected_components"


class TestMetricsFrame:
    def test_one_row_per_clusterer_metric(self) -> None:
        payload = _make_payload()
        frame = metrics_frame(payload)
        assert frame.height == 2 * len(FIG3_METRICS)
        assert set(frame["clusterer"].unique().to_list()) == {
            "tn93+connected_components",
            "hivtrace-hivnetworkcsv",
        }
        assert set(frame["metric"].unique().to_list()) == set(FIG3_METRICS)

    def test_ordering_columns(self) -> None:
        payload = _make_payload()
        frame = metrics_frame(payload)
        assert frame.columns == ["clusterer", "metric", "value"]

    def test_empty_results_raises(self) -> None:
        with pytest.raises(ValueError, match="no scored metrics"):
            metrics_frame({"results": []})


class TestConfusionFrame:
    def test_shape(self) -> None:
        frame = confusion_frame(_make_payload())
        assert set(frame.columns) == {"clusterer", "tp", "fp", "fn", "tn"}
        assert frame.height == 2


class TestSummarize:
    def test_with_hivtrace(self) -> None:
        payload = _make_payload(
            with_hivtrace=True, ours_f1=0.95, hivtrace_f1=0.93, agreement=0.99
        )
        s = summarize(payload)
        assert s.n_clusterers == 2
        assert s.ours_f1 == pytest.approx(0.95)
        assert s.hivtrace_f1 == pytest.approx(0.93)
        assert s.pair_agreement == pytest.approx(0.99)
        assert s.n_sequences == 200
        assert s.n_truth_clusters == 20

    def test_without_hivtrace(self) -> None:
        payload = _make_payload(with_hivtrace=False, ours_f1=0.95)
        s = summarize(payload)
        assert s.n_clusterers == 1
        assert s.hivtrace_f1 is None
        assert s.pair_agreement is None


class TestConstants:
    def test_clusterer_order_locked(self) -> None:
        assert CLUSTERER_ORDER == (
            "tn93+connected_components",
            "hivtrace-hivnetworkcsv",
        )

    def test_fig3_metrics_locked(self) -> None:
        assert FIG3_METRICS == ("precision", "recall", "f1", "adjusted_rand_index")
