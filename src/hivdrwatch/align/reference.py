"""HXB2 pol reference sequence — the coordinate system everything else lives in.

HXB2 (GenBank K03455.1) is the field-standard HIV-1 reference. Every DRM
position, alignment coordinate, and Sierra call is expressed relative to it.
We pin the pol region as K03455 positions 2085..5096 (3012 bp, 1-based
inclusive), which spans PR + RT + IN — everything Sierra scores.

We don't ship the HXB2 sequence as a raw file in the repo. Instead we
extract it from any LANL alignment on hand (LANL always ships HXB2 as the
alignment reference — its header contains "HXB2" and its accession is
K03455). This keeps our reference downstream of the *actual* alignment
context the sequences come from.

Standard HIV region boundaries in HXB2 K03455 nucleotide coordinates:

    Region     K03455 nt      pol-relative nt   codons (aa)
    ------     ------------   ---------------   -----------
    PR          2253..2549       169..465           99
    RT          2550..3869       466..1785         440   (polymerase domain)
    RT (p66)    2550..4229       466..2145         560   (full p66 subunit)
    IN          4230..5096      2146..3012         288

We use `RT` = 2550..3869 by default because that covers all Stanford-scored
NRTI/NNRTI positions; extend to 4229 only if you need the RNase-H domain.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


# HXB2 K03455.1 nucleotide coordinates (1-based, inclusive).
HXB2_POL_START = 2085
HXB2_POL_STOP = 5096
HXB2_POL_LENGTH = HXB2_POL_STOP - HXB2_POL_START + 1  # 3012

# All ranges below in **pol-relative** 1-based inclusive coordinates
# (subtract HXB2_POL_START - 1 = 2084 from an HXB2 coordinate to convert).
POL_REGIONS: dict[str, tuple[int, int]] = {
    "PR": (169, 465),      # HXB2 2253..2549, 99 codons
    "RT": (466, 1785),     # HXB2 2550..3869, 440 codons (polymerase domain)
    "IN": (2146, 3012),    # HXB2 4230..5096, 288 codons
}


@dataclass(frozen=True)
class HXB2Reference:
    """Extracted HXB2 pol reference sequence."""

    sequence: str            # ungapped nucleotide sequence
    length: int              # sequence length
    source_path: Path        # where we extracted it from

    def write_fasta(self, out_path: Path, *, header: str = "HXB2_pol_K03455_2085_5096") -> Path:
        """Write as a single-record FASTA — used as MAFFT --add reference."""
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(f">{header}\n{self.sequence}\n", encoding="utf-8")
        return out_path


def extract_hxb2_from_lanl_fasta(lanl_fasta: Path) -> HXB2Reference:
    """Extract the HXB2 record from a LANL alignment FASTA.

    LANL always ships HXB2 as the first record; its header contains the
    substring 'HXB2'. We keep the first matching record and strip alignment
    gaps to get the ungapped reference.
    """
    header: str | None = None
    buffer: list[str] = []
    capturing = False
    with lanl_fasta.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if capturing:
                    break  # done — we've collected the HXB2 record
                if "HXB2" in line.upper():
                    header = line[1:].strip()
                    capturing = True
                continue
            if capturing:
                buffer.append(line.strip())

    if not capturing or header is None:
        raise ValueError(f"no HXB2 record found in {lanl_fasta}")

    ungapped = "".join(buffer).replace("-", "").upper()
    logger.info(
        "extracted HXB2 pol reference (%d nt) from %s [header=%s]",
        len(ungapped),
        lanl_fasta.name,
        header,
    )
    return HXB2Reference(sequence=ungapped, length=len(ungapped), source_path=lanl_fasta)


def ensure_hxb2_reference(
    reference_path: Path, lanl_raw_dir: Path
) -> HXB2Reference:
    """Return the HXB2 pol reference, extracting it if not yet cached.

    Idempotent: if `reference_path` already exists, we read it back.
    Otherwise we walk `lanl_raw_dir` for any hiv-db.fasta, extract HXB2,
    and cache it as a FASTA at `reference_path`.
    """
    if reference_path.exists():
        seq = _read_single_fasta(reference_path)
        return HXB2Reference(
            sequence=seq, length=len(seq), source_path=reference_path
        )

    fasta_candidates = sorted(lanl_raw_dir.rglob("hiv-db.fasta"))
    if not fasta_candidates:
        raise FileNotFoundError(
            f"no LANL FASTA under {lanl_raw_dir}; can't extract HXB2 reference. "
            f"Run ingestion first or place HXB2 pol at {reference_path}."
        )
    ref = extract_hxb2_from_lanl_fasta(fasta_candidates[0])

    if ref.length != HXB2_POL_LENGTH:
        logger.warning(
            "extracted HXB2 length %d != expected %d — LANL alignment may span "
            "a different region than K03455 2085..5096",
            ref.length,
            HXB2_POL_LENGTH,
        )
    ref.write_fasta(reference_path)
    return ref


def _read_single_fasta(path: Path) -> str:
    """Read a single-record FASTA and return the ungapped uppercase sequence."""
    buffer: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            buffer.append(line.strip())
    return "".join(buffer).replace("-", "").upper()
