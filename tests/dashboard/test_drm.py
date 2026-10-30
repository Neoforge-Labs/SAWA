"""Tests for the surveillance drug-resistance mutation (SDRM) matcher.

Sierra reports every variant against consensus, so its most frequent mutation
strings are background polymorphisms (PR R41K, RT V35T, ...). These tests pin
that the matcher keeps only the Stanford surveillance list, and that the
mixture, insertion and malformed-input conventions are applied consistently.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hivdrwatch.dashboard.export.drm import (
    SDRM_PANELS,
    MissingPanelError,
    SdrmPanel,
    load_sdrm_panels,
    match_sdrm,
    parse_mutation,
)

REFERENCE_PATH = (
    Path(__file__).resolve().parents[2] / "data" / "reference" / "drm_positions.json"
)


def _panel(name: str) -> SdrmPanel:
    panels, _ = load_sdrm_panels(REFERENCE_PATH)
    return next(p for p in panels if p.name == name)


def test_loads_only_pr_and_rt_panels():
    panels, version = load_sdrm_panels(REFERENCE_PATH)
    assert tuple(p.name for p in panels) == ("PR", "RT_NRTI", "RT_NNRTI")
    assert all(p.gene in {"PR", "RT"} for p in panels)
    assert version == "SDRM-2019"


def test_integrase_is_not_a_surveillance_panel():
    # Bennett et al. 2009 has no integrase component, and the reference file's
    # IN block predates dolutegravir-era mutations.
    assert "IN" not in SDRM_PANELS
    panels, _ = load_sdrm_panels(REFERENCE_PATH)
    assert "IN" not in {p.name for p in panels}


def test_panel_positions_match_reference_file():
    raw = json.loads(REFERENCE_PATH.read_text())["genes"]
    for panel in load_sdrm_panels(REFERENCE_PATH)[0]:
        assert set(panel.positions) == {int(k) for k in raw[panel.name]["positions"]}
    assert 220 not in _panel("RT_NNRTI").positions


def test_missing_panel_raises_a_specific_error(tmp_path: Path):
    path = tmp_path / "drm.json"
    path.write_text(
        json.dumps(
            {
                "algorithm_version_pinned": "SDRM-2019",
                "genes": {"PR": {"positions": {"46": ["I"]}}},
            }
        )
    )
    with pytest.raises(MissingPanelError, match="RT_NRTI"):
        load_sdrm_panels(path)


@pytest.mark.parametrize("raw", ["R41K", "M36I", "V35T", "K122E", "D177E"])
def test_known_polymorphisms_are_not_sdrms(raw):
    panels, _ = load_sdrm_panels(REFERENCE_PATH)
    assert all(match_sdrm(raw, p) is None for p in panels)


@pytest.mark.parametrize(
    ("raw", "panel_name", "expected"),
    [("M184V", "RT_NRTI", "M184V"), ("K103N", "RT_NNRTI", "K103N"), ("M46I", "PR", "M46I")],
)
def test_known_drms_match_their_panel(raw, panel_name, expected):
    panel = _panel(panel_name)
    assert match_sdrm(raw, panel) == expected


def test_mixture_counts_as_carrying_any_matching_state():
    # Surveillance convention (spec §7.2): any state in the mixture counts.
    assert match_sdrm("M184MV", _panel("RT_NRTI")) == "M184V"
    assert match_sdrm("T215TFILMS", _panel("RT_NRTI")) is not None


def test_mixture_matching_several_states_resolves_deterministically():
    # F, I and S are all SDRM states at 215; the label must not depend on the
    # order the states were written in or on set iteration order.
    panel = _panel("RT_NRTI")
    labels = {match_sdrm(raw, panel) for raw in ("T215TFIS", "T215TSIF", "T215SIFT")}
    assert labels == {"T215F"}


def test_mixture_with_no_matching_state_does_not_match():
    # RT_NRTI 184 states are {V, I}; T is neither.
    assert match_sdrm("M184MT", _panel("RT_NRTI")) is None


def test_both_sdrm_states_at_a_position_match():
    # 184 lists V and I: M184I is an SDRM in its own right, not a near-miss.
    assert match_sdrm("M184I", _panel("RT_NRTI")) == "M184I"


def test_stop_codons_and_deletions_are_excluded():
    # Both sit at real SDRM positions (RT_NRTI 184 and 69), so these assert the
    # exclusion rule rather than a position miss.
    assert match_sdrm("M184*", _panel("RT_NRTI")) is None
    assert match_sdrm("T69del", _panel("RT_NRTI")) is None


def test_rt69_insertion_matches_the_ins_state():
    assert match_sdrm("T69S_SS", _panel("RT_NRTI")) == "T69ins"


def test_sdrm_position_with_non_sdrm_state_does_not_match():
    assert match_sdrm("K103R", _panel("RT_NNRTI")) is None


def test_same_mutation_does_not_match_in_the_wrong_panel():
    assert match_sdrm("M184V", _panel("PR")) is None
    assert match_sdrm("K103N", _panel("RT_NRTI")) is None


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("M184V", (184, frozenset({"V"}))),
        ("T215TFILMS", (215, frozenset("TFILMS"))),
        ("T69S_SS", (69, frozenset({"ins"}))),
    ],
)
def test_parse_mutation(raw, expected):
    assert parse_mutation(raw) == expected


@pytest.mark.parametrize("raw", ["X", "", "T215", "...", "69"])
def test_malformed_strings_are_rejected_without_raising(raw):
    assert parse_mutation(raw) is None
    assert match_sdrm(raw, _panel("RT_NRTI")) is None
