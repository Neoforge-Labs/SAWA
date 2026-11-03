"""Unit tests for the Fig 3 cluster-parity aggregator."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hivdrwatch.paper1.fig03_parity import (
    CLUSTERER_ORDER,
    FIG3_METRICS,
    UNRECORDED_SKIP_REASON,
    ParitySummary,
    confusion_frame,
    load_parity_payload,
    metrics_frame,
    provenance_note,
    provenance_state,
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


def _deliberate_tn93_only() -> dict:
    """What cluster_validate writes when --with-hivtrace was not passed."""
    return {
        **_make_payload(with_hivtrace=False),
        "with_hivtrace": False,
        "hivtrace_skipped": None,
    }


def _requested_and_skipped(reason: str = "image_unavailable") -> dict:
    """What cluster_validate writes when --with-hivtrace was passed but the run degraded."""
    return {
        **_make_payload(with_hivtrace=False),
        "with_hivtrace": True,
        "hivtrace_skipped": reason,
    }


def _parity_ran() -> dict:
    return {
        **_make_payload(with_hivtrace=True, agreement=0.99),
        "with_hivtrace": True,
        "hivtrace_skipped": None,
    }


class TestSummarizeProvenance:
    def test_reads_skipped_run(self) -> None:
        s = summarize(_requested_and_skipped("image_unavailable"))
        assert s.with_hivtrace is True
        assert s.hivtrace_skipped == "image_unavailable"

    def test_reads_deliberate_run(self) -> None:
        s = summarize(_deliberate_tn93_only())
        assert s.with_hivtrace is False
        assert s.hivtrace_skipped is None

    def test_reads_successful_run(self) -> None:
        s = summarize(_parity_ran())
        assert s.with_hivtrace is True
        assert s.hivtrace_skipped is None

    def test_results_json_without_provenance_keys_loads_as_not_requested(
        self, tmp_path: Path
    ) -> None:
        old = _make_payload(with_hivtrace=False)
        assert "with_hivtrace" not in old
        assert "hivtrace_skipped" not in old
        p = tmp_path / "old.json"
        p.write_text(json.dumps(old))
        s = summarize(load_parity_payload(p))
        assert s.with_hivtrace is False
        assert s.hivtrace_skipped is None
        assert provenance_state(s) == "not_requested"

    def test_existing_construction_sites_keep_working(self) -> None:
        s = ParitySummary(
            n_clusterers=1,
            n_sequences=200,
            n_truth_clusters=20,
            pair_agreement=None,
            ours_f1=0.95,
            hivtrace_f1=None,
        )
        assert s.with_hivtrace is False
        assert s.hivtrace_skipped is None

    def test_loader_keeps_provenance_when_normalising_single_clusterer_payload(
        self, tmp_path: Path
    ) -> None:
        legacy_shape = _requested_and_skipped("image_unavailable")
        del legacy_shape["results"]
        del legacy_shape["parity"]
        p = tmp_path / "cv.json"
        p.write_text(json.dumps(legacy_shape))
        payload = load_parity_payload(p)
        assert len(payload["results"]) == 1
        s = summarize(payload)
        assert s.with_hivtrace is True
        assert s.hivtrace_skipped == "image_unavailable"


class TestProvenanceNote:
    def test_skipped_run_reads_differently_from_deliberate_tn93_only(self) -> None:
        deliberate = _deliberate_tn93_only()
        skipped = _requested_and_skipped()
        # Same scored data on both sides: only the provenance may tell them apart.
        assert metrics_frame(deliberate).equals(metrics_frame(skipped))

        deliberate_note = provenance_note(summarize(deliberate))
        skipped_note = provenance_note(summarize(skipped))

        assert deliberate_note != skipped_note
        assert provenance_state(summarize(deliberate)) == "not_requested"
        assert provenance_state(summarize(skipped)) == "requested_skipped"

    def test_skipped_note_says_requested_skipped_and_names_reason(self) -> None:
        note = provenance_note(summarize(_requested_and_skipped("image_unavailable")))
        assert "requested but skipped" in note
        assert "image_unavailable" in note
        assert "does not support a parity claim" in note

    def test_skipped_note_names_whatever_reason_was_recorded(self) -> None:
        note = provenance_note(summarize(_requested_and_skipped("container_timeout")))
        assert "container_timeout" in note
        assert "image_unavailable" not in note

    def test_deliberate_note_does_not_suggest_a_failure(self) -> None:
        note = provenance_note(summarize(_deliberate_tn93_only()))
        assert "not requested" in note
        for word in ("skipped", "DEGRADED", "requested but"):
            assert word not in note

    def test_parity_ran_note_reports_agreement_and_no_degradation(self) -> None:
        note = provenance_note(summarize(_parity_ran()))
        assert "0.990" in note
        for word in ("skipped", "DEGRADED", "not requested"):
            assert word not in note

    def test_three_states_yield_three_distinct_notes(self) -> None:
        notes = {
            provenance_note(summarize(p))
            for p in (_deliberate_tn93_only(), _requested_and_skipped(), _parity_ran())
        }
        assert len(notes) == 3

    def test_hivtrace_result_without_flag_is_parity_ran_not_not_requested(self) -> None:
        # A parity run written before the flag existed has the bars but no flag.
        payload = _make_payload(with_hivtrace=True)
        assert "with_hivtrace" not in payload
        s = summarize(payload)
        assert provenance_state(s) == "parity_ran"
        assert "not requested" not in provenance_note(s)

    def test_requested_with_no_result_and_no_reason_is_not_read_as_deliberate(self) -> None:
        payload = {**_make_payload(with_hivtrace=False), "with_hivtrace": True}
        s = summarize(payload)
        assert provenance_state(s) == "requested_skipped"
        note = provenance_note(s)
        assert UNRECORDED_SKIP_REASON in note
        assert "not requested" not in note

    def test_skip_reason_alone_marks_the_run_degraded(self) -> None:
        payload = {
            **_make_payload(with_hivtrace=False),
            "hivtrace_skipped": "image_unavailable",
        }
        assert provenance_state(summarize(payload)) == "requested_skipped"


class TestConstants:
    def test_clusterer_order_locked(self) -> None:
        assert CLUSTERER_ORDER == (
            "tn93+connected_components",
            "hivtrace-hivnetworkcsv",
        )

    def test_fig3_metrics_locked(self) -> None:
        assert FIG3_METRICS == ("precision", "recall", "f1", "adjusted_rand_index")
