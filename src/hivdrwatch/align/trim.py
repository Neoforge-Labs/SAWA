"""Slice HXB2-aligned sequences into PR / RT / IN gene regions.

Alignment coordinates are HXB2 pol-relative (1-based inclusive). After
MAFFT --keeplength, every aligned sequence has exactly HXB2_POL_LENGTH
columns, so we can slice by pol-relative range without per-sequence
coordinate mapping.

If a sequence doesn't cover a region (all gaps in that range), the slice
comes back as all-gap and downstream QC flags SHORT_COVERAGE.
"""

from __future__ import annotations

from .reference import HXB2_POL_LENGTH, POL_REGIONS


def slice_region(aligned_sequence: str, region: str) -> str:
    """Return the gap-preserving slice for the named region ('PR', 'RT', 'IN')."""
    if region not in POL_REGIONS:
        raise ValueError(f"unknown region {region!r}; expected one of {list(POL_REGIONS)}")
    if len(aligned_sequence) != HXB2_POL_LENGTH:
        raise ValueError(
            f"aligned sequence length {len(aligned_sequence)} != expected "
            f"{HXB2_POL_LENGTH}; was MAFFT run with --keeplength?"
        )
    start, stop = POL_REGIONS[region]
    return aligned_sequence[start - 1 : stop]


def slice_all_regions(aligned_sequence: str) -> dict[str, str]:
    """Slice into every named region at once."""
    return {region: slice_region(aligned_sequence, region) for region in POL_REGIONS}


def coverage_fraction(region_seq: str) -> float:
    """Fraction of non-gap positions in a region slice.

    Used to flag SHORT_COVERAGE when < 0.8 for the trimmed region.
    """
    if not region_seq:
        return 0.0
    non_gap = sum(1 for c in region_seq if c != "-")
    return non_gap / len(region_seq)
