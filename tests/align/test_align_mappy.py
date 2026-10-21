"""Coverage regression tests for the mappy backend.

These pin the property the curated layer actually depends on: an input
sequence that *contains* HXB2 pol must be placed across the whole of it,
not across whichever block minimap2 could chain at the configured
divergence tolerance.

The fixture is a frozen GenBank pol region (CRF01_AE, ~13% diverged from
HXB2 subtype B) — the divergence regime the whole corpus sits in, since
HXB2 is subtype B and SAWA's cohort is not.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pytest

from hivdrwatch.align.align_biopython import align_sequences_biopython
from hivdrwatch.align.align_mappy import (
    DEFAULT_PRESET,
    _build_keeplength,
    align_sequences_mappy,
)
from hivdrwatch.align.qc import passes_qc, qc_flags
from hivdrwatch.align.reference import HXB2_POL_LENGTH
from hivdrwatch.align.schema import QCFlag


@dataclass
class _FakeHit:
    """Stands in for a `mappy.Alignment` — only the fields we project from."""

    r_st: int
    q_st: int
    mlen: int
    cigar: list[tuple[int, int]] = field(default_factory=list)
    strand: int = 1


FIXTURE = Path(__file__).parent.parent / "fixtures" / "align" / "genbank_KY580541.1_pol.fasta"
HXB2_POL = Path(__file__).parents[2] / "data" / "reference" / "hxb2_pol.fasta"


def _read_fasta(path: Path) -> tuple[str, str]:
    header, buf = "", []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].split()[0]
        else:
            buf.append(line.strip())
    return header, "".join(buf).upper()


def _covered(aligned: str) -> int:
    return sum(1 for c in aligned if c != "-")


def _span(aligned: str) -> tuple[int, int]:
    idx = [i for i, c in enumerate(aligned) if c != "-"]
    return (idx[0] + 1, idx[-1] + 1) if idx else (0, 0)


@pytest.fixture(scope="module")
def divergent_pol() -> tuple[str, str]:
    return _read_fasta(FIXTURE)


@pytest.fixture(scope="module")
def reference() -> Path:
    if not HXB2_POL.exists():
        pytest.skip(f"HXB2 pol reference not present at {HXB2_POL}")
    return HXB2_POL


class TestDivergentSubtypeCoverage:
    """A non-B pol that spans the reference must be aligned across the reference."""

    def test_covers_almost_all_of_the_reference(
        self, divergent_pol: tuple[str, str], reference: Path
    ) -> None:
        acc, seq = divergent_pol
        # The fixture is the pol region itself, so there is nothing in it
        # that could legitimately fall outside HXB2 pol.
        assert len(seq) == HXB2_POL_LENGTH

        aligned = align_sequences_mappy([(acc, seq)], reference)[acc]
        assert len(aligned) == HXB2_POL_LENGTH
        covered = _covered(aligned)
        assert covered >= int(0.95 * HXB2_POL_LENGTH), (
            f"{acc} covers only {covered}/{HXB2_POL_LENGTH} nt "
            f"(span {_span(aligned)}) — the aligner dropped sequence it was given"
        )

    def test_protease_is_not_silently_lost(
        self, divergent_pol: tuple[str, str], reference: Path
    ) -> None:
        """PR is pol 169..465. Losing it costs every PI call and half of Fig 2."""
        acc, seq = divergent_pol
        aligned = align_sequences_mappy([(acc, seq)], reference)[acc]
        pr = aligned[168:465]
        assert _covered(pr) >= int(0.95 * len(pr)), (
            f"{acc} protease coverage {_covered(pr)}/{len(pr)}"
        )


class TestBackendAgreement:
    """`align/README.md` promises all backends share HXB2 coordinates.

    Sharing a coordinate system is worth nothing if the backends disagree
    about which coordinates the same sequence occupies.
    """

    def test_mappy_matches_biopython_span(
        self, divergent_pol: tuple[str, str], reference: Path
    ) -> None:
        acc, seq = divergent_pol
        recs = [(acc, seq)]
        m = align_sequences_mappy(recs, reference)[acc]
        b = align_sequences_biopython(recs, reference, n_workers=1)[acc]
        assert abs(_covered(m) - _covered(b)) <= int(0.05 * HXB2_POL_LENGTH), (
            f"mappy covers {_covered(m)} nt (span {_span(m)}), "
            f"biopython covers {_covered(b)} nt (span {_span(b)})"
        )


class TestMultiBlockProjection:
    """A query minimap2 splits into several blocks must keep all of them.

    1.64% of the corpus splits under the default preset. Keeping only
    the largest block silently discards the rest of the sequence.
    """

    def test_second_block_fills_around_the_first(self) -> None:
        # Two blocks, disjoint on both query and reference.
        query = "A" * 10 + "C" * 10
        big = _FakeHit(r_st=0, q_st=0, cigar=[(10, 0)], mlen=10)
        small = _FakeHit(r_st=20, q_st=10, cigar=[(10, 0)], mlen=10)
        out = _build_keeplength(query, [big, small], 40)
        assert out == "A" * 10 + "-" * 10 + "C" * 10 + "-" * 10

    def test_best_hit_wins_overlapping_columns(self) -> None:
        query = "A" * 10 + "C" * 10
        best = _FakeHit(r_st=0, q_st=0, cigar=[(10, 0)], mlen=10)
        weaker = _FakeHit(r_st=5, q_st=10, cigar=[(10, 0)], mlen=4)
        out = _build_keeplength(query, [best, weaker], 20)
        # columns 0..9 stay with the best hit; 10..14 come from the weaker one
        assert out == "A" * 10 + "C" * 5 + "-" * 5

    def test_single_block_is_unchanged(self) -> None:
        query = "ACGTACGTAC"
        hit = _FakeHit(r_st=2, q_st=0, cigar=[(10, 0)], mlen=10)
        assert _build_keeplength(query, [hit], 14) == "--ACGTACGTAC--"

    def test_insertions_are_dropped_and_deletions_become_gaps(self) -> None:
        # M3, I2 (dropped), M3, D2 (gap), M2
        query = "AAAGGCCCTT"
        hit = _FakeHit(r_st=0, q_st=0, cigar=[(3, 0), (2, 1), (3, 0), (2, 2), (2, 0)], mlen=8)
        assert _build_keeplength(query, [hit], 12) == "AAACCC--TT--"


class TestNoAlignmentBlocks:
    """An all-gap row must not reach downstream as a QC pass."""

    def test_all_gap_regions_flag_and_block(self) -> None:
        flags = qc_flags(
            ambiguity=0.0,
            apobec_ratio=None,
            region_coverage={"PR": 0.0, "RT": 0.0, "IN": 0.0},
            region_stops={"PR": False, "RT": False, "IN": False},
        )
        assert QCFlag.NO_ALIGNMENT in flags
        assert not passes_qc(flags)

    def test_partial_coverage_still_passes(self) -> None:
        """RT-only genotyping sequences are the common case — keep them."""
        flags = qc_flags(
            ambiguity=0.0,
            apobec_ratio=None,
            region_coverage={"PR": 0.0, "RT": 0.98, "IN": 0.0},
            region_stops={"PR": False, "RT": False, "IN": False},
        )
        assert QCFlag.NO_ALIGNMENT not in flags
        assert QCFlag.SHORT_COVERAGE in flags
        assert passes_qc(flags)


class TestPresetTolerance:
    """The default preset must tolerate inter-subtype divergence.

    HXB2 is subtype B; the cohort is A1/C/D/CRF. Measured divergence of
    recovered pol against HXB2 across the corpus is ~13%, so a preset
    advertising ~10% is below the data.
    """

    def test_default_preset_tolerates_at_least_twenty_percent(self) -> None:
        assert DEFAULT_PRESET == "asm20"
