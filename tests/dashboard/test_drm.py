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
    # F, I and S are all 215 states; sorted order picks F.
    assert match_sdrm("T215TFILMS", _panel("RT_NRTI")) == "T215F"


@pytest.mark.parametrize(
    ("raws", "panel_name"),
    [
        (("T215TFIS", "T215TSIF", "T215SIFT"), "RT_NRTI"),
        (("T215TFILMS", "T215YFTSCD", "T215DSCYE"), "RT_NRTI"),
        (("A82AFTSLMC", "A82CMLSTF", "A82TFA"), "PR"),
        (("K65KRNE", "K65ENRK", "K65NK"), "RT_NRTI"),
    ],
)
def test_mixture_label_is_the_smallest_matching_state(raws, panel_name):
    # The label must not depend on the order the states were written in or on
    # set iteration order, so it is pinned to the smallest matching state.
    panel = _panel(panel_name)
    for raw in raws:
        parsed = _parse_with_reference(raw)
        reference, position, observed = parsed
        expected = f"{reference}{position}{min(panel.positions[position] & observed)}"
        assert match_sdrm(raw, panel) == expected


@pytest.mark.parametrize(
    ("raw", "panel_name", "attributed", "other_matching_state"),
    [
        ("M46IL", "PR", "M46I", "L"),
        ("M46LI", "PR", "M46I", "L"),
        ("V82AF", "PR", "V82A", "F"),
        ("V82FA", "PR", "V82A", "F"),
    ],
)
def test_a_mixture_is_attributed_to_the_alphabetically_first_matching_state(
    raw: str, panel_name: str, attributed: str, other_matching_state: str
):
    # The page tells readers that a multi-state mixture goes to the alphabetically first
    # matching state, which inflates the earlier label at the later one's expense. That
    # sentence is only true while this holds, so a matcher change (counting every state, or
    # the last one) must fail here and send someone to the copy.
    panel = _panel(panel_name)
    position = int(raw[1:3])
    assert other_matching_state in panel.positions[position]
    assert match_sdrm(raw, panel) == attributed
    assert attributed[-1] < other_matching_state


def _parse_with_reference(raw: str) -> tuple[str, int, frozenset[str]]:
    parsed = parse_mutation(raw)
    assert parsed is not None
    position, states = parsed
    return raw[0], position, states


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


def test_substitution_beside_an_insertion_is_still_counted():
    assert match_sdrm("M46I_E", _panel("PR")) == "M46I"
    assert match_sdrm("M230L_C", _panel("RT_NNRTI")) == "M230L"


def test_insertion_beside_a_non_sdrm_residue_is_not_an_sdrm():
    # 65 lists E/N/R, so K65K_K is neither a substitution nor an insertion SDRM.
    assert match_sdrm("K65K_K", _panel("RT_NRTI")) is None


def test_ins_letters_in_a_mixture_are_ordinary_residues():
    panel = _panel("RT_NRTI")
    # Matched set is {F, I, S, Y}; the label is the smallest.
    assert match_sdrm("T215TFINSY", panel) == "T215F"
    # 215 lists I and S, both carried by the I-N-S run.
    assert match_sdrm("T215TINS", panel) == "T215I"


@pytest.mark.parametrize(
    "raw",
    ["A299L_LACAL*PH", "A538ADEVG_X", "D549N_TWMVLQASTS*SR*SR", "T69_S", "T69SS_S"],
)
def test_messy_underscore_forms_stay_rejected(raw):
    assert parse_mutation(raw) is None


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
        # The residue before the underscore is a substitution at the position and
        # the part after it is an insertion; both are carried.
        ("M46I_E", (46, frozenset({"I", "ins"}))),
        ("T69S_SS", (69, frozenset({"S", "ins"}))),
        ("K65K_K", (65, frozenset({"K", "ins"}))),
        # Sierra writes a mixture as the retained consensus residue (if any) followed by
        # the variant residues in alphabetical order, so the I-N-S run in these strings is
        # three adjacent letters, not an insertion marker.
        ("A554INS", (554, frozenset({"I", "N", "S"}))),
        ("I47INST", (47, frozenset({"I", "N", "S", "T"}))),
        ("T215TFINSY", (215, frozenset({"T", "F", "I", "N", "S", "Y"}))),
    ],
)
def test_parse_mutation(raw, expected):
    assert parse_mutation(raw) == expected


@pytest.mark.parametrize("raw", ["X", "", "T215", "...", "69"])
def test_malformed_strings_are_rejected_without_raising(raw):
    assert parse_mutation(raw) is None
    assert match_sdrm(raw, _panel("RT_NRTI")) is None


GOLDEN_PANELS_PATH = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "fixtures"
    / "dashboard"
    / "sdrm_panels_golden.json"
)


def test_surveillance_panels_are_exactly_the_positions_and_states_pinned_in_the_golden():
    # The page attributes these panels to Bennett et al. 2009 (updated 2019), while the data
    # file's own notes call them candidate model features, and nothing here checks them
    # against a published list. Until that is verified, any edit to a position or state
    # changes what the page publishes under that citation, so it must be a deliberate edit
    # to the golden in the same commit and read as one in review. Deliberately not wired to
    # `--update-golden`: regenerating this file wholesale would defeat the point.
    golden = json.loads(GOLDEN_PANELS_PATH.read_text())
    panels, version = load_sdrm_panels(REFERENCE_PATH)
    assert version == golden["algorithm_version_pinned"]
    assert [p.name for p in panels] == list(golden["panels"])
    for panel in panels:
        pinned = {int(pos): states for pos, states in golden["panels"][panel.name].items()}
        loaded = {pos: sorted(states) for pos, states in panel.positions.items()}
        assert loaded == pinned, (
            f"{panel.name} differs from the pinned golden; if the change is intended, edit "
            f"{GOLDEN_PANELS_PATH.name} in the same commit"
        )
