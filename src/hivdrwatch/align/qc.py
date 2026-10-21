"""Quality control for aligned HIV pol sequences.

Three checks that matter for downstream modeling:

1. **Ambiguity fraction** — count of non-ACGT-gap positions. High ambiguity
   sequences (>5%) usually mean bulk sequencing of a mixed population, or
   base-calling problems. Both bias mutation calling.

2. **APOBEC3G hypermutation signature** — a heavily G→A-mutated sequence
   isn't representing the transmitted virus; it's a defective form the
   host's innate immunity edited. Detect via GG→AG bias vs GA→AA
   (Rose & Korber 2000; Hypermut algorithm). Not a full Hypermut re-impl
   here — we compute the ratio and flag anything above a conservative
   threshold. Rule out with a proper Hypermut run before publishing.

3. **Premature stop codons in PR/RT/IN** — a real sequence has open reading
   frames across all three. Stops mid-gene often mean a hypermutated or
   frameshifted read; drop from resistance calling.

Everything here is a pure function on strings — trivial to unit-test.
"""

from __future__ import annotations

from .schema import QCFlag

# Thresholds — deliberately conservative. Downstream can override via a
# QCConfig object if we ever want per-project tuning.
MAX_AMBIGUITY_FRACTION = 0.05
MIN_REGION_COVERAGE = 0.80
APOBEC_RATIO_THRESHOLD = 2.0  # GG→AG / GA→AA — hypermut algorithm hint

_UNAMBIG = frozenset("ACGT")


def ambiguity_fraction(sequence: str) -> float:
    """Fraction of non-ACGT positions among non-gap characters.

    Ignores gap columns entirely (gaps are alignment artifacts, not
    sequence quality). Returns 0.0 for all-gap sequences.
    """
    non_gap = [c for c in sequence.upper() if c != "-"]
    if not non_gap:
        return 0.0
    ambig = sum(1 for c in non_gap if c not in _UNAMBIG)
    return ambig / len(non_gap)


def apobec_signature_ratio(sequence: str, reference: str) -> float | None:
    """Compare G→A substitutions in GG vs GA dinucleotide contexts.

    APOBEC3G preferentially deaminates G in a GG context (5'-GG → GA
    on the reverse strand appears as GG → AG on the query). Real
    variation is roughly context-neutral, so a query:reference ratio
    of (GG→AG changes) / (GA→AA changes) noticeably above 1 flags
    hypermutation.

    Returns None when the reference has no eligible context sites
    (denominator would be 0) — that's ambiguous rather than clean.
    """
    q = sequence.upper()
    r = reference.upper()
    if len(q) != len(r):
        raise ValueError(f"length mismatch: query={len(q)} ref={len(r)}")

    gg_to_ag = 0
    gg_sites = 0
    ga_to_aa = 0
    ga_sites = 0

    # Walk positions where reference has G with a downstream context nt.
    for i in range(len(r) - 1):
        if r[i] != "G" or r[i + 1] not in "ACGT":
            continue
        # Skip if query is gapped or ambiguous here.
        if q[i] not in "AG" or q[i + 1] not in "ACGT":
            continue
        context = r[i + 1]
        if context == "G":
            gg_sites += 1
            if q[i] == "A" and q[i + 1] == "G":
                gg_to_ag += 1
        elif context == "A":
            ga_sites += 1
            if q[i] == "A" and q[i + 1] == "A":
                ga_to_aa += 1

    if gg_sites == 0 or ga_sites == 0:
        return None
    gg_rate = gg_to_ag / gg_sites
    ga_rate = ga_to_aa / ga_sites
    if ga_rate == 0:
        # No control substitutions but GG→AG present → strong signal;
        # return a large sentinel so the caller flags hypermut.
        return float("inf") if gg_to_ag > 0 else 0.0
    return gg_rate / ga_rate


def translate_in_frame(nucleotide_seq: str) -> str:
    """Translate a gap-preserving nucleotide slice to amino acids.

    Gaps are handled per-codon: if any position in the codon is '-'
    or ambiguous, the amino acid is 'X'. Stop codons render as '*'.
    This is intentionally lenient because we only care whether stops
    show up in the *sequence-derived* codons, not what X represents.
    """
    seq = nucleotide_seq.upper()
    aa_out: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        aa_out.append(_CODON_TABLE.get(codon, "X"))
    return "".join(aa_out)


def has_premature_stop(nucleotide_seq: str) -> bool:
    """True if the translation contains a '*' anywhere but the final codon."""
    protein = translate_in_frame(nucleotide_seq)
    if not protein:
        return False
    # Ignore a stop at the very end — some regions include the natural
    # terminator (though PR/RT/IN slices typically don't).
    return "*" in protein[:-1]


# Standard genetic code (extended: gaps → X). Only the standard codons
# are listed here; anything else falls back to 'X' in translate_in_frame.
_CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def qc_flags(
    *,
    ambiguity: float,
    apobec_ratio: float | None,
    region_coverage: dict[str, float],
    region_stops: dict[str, bool],
) -> list[QCFlag]:
    """Aggregate all QC signals into a list of flags."""
    flags: list[QCFlag] = []
    if region_coverage and not any(c > 0.0 for c in region_coverage.values()):
        flags.append(QCFlag.NO_ALIGNMENT)
    if ambiguity > MAX_AMBIGUITY_FRACTION:
        flags.append(QCFlag.HIGH_AMBIGUITY)
    if apobec_ratio is not None and apobec_ratio > APOBEC_RATIO_THRESHOLD:
        flags.append(QCFlag.APOBEC_HYPERMUT)
    for _region, coverage in region_coverage.items():
        if coverage < MIN_REGION_COVERAGE:
            flags.append(QCFlag.SHORT_COVERAGE)
            break
    for _region, has_stop in region_stops.items():
        if has_stop:
            flags.append(QCFlag.STOP_CODON)
            break
    return flags


# Flags that mean the sequence itself is unreliable — these block
# admission to feature building. SHORT_COVERAGE is intentionally NOT
# here: most GenBank HIV sequences cover only RT (for genotyping) or
# only part of pol; requiring full coverage across PR + RT + IN would
# throw away >95% of real data. Downstream feature builders already
# emit per-position `covered` flags, so a sequence that only covers
# RT still contributes to NRTI/NNRTI features with zero on PR/IN.
#
# NO_ALIGNMENT *does* block. A row with no PR, RT or IN site covered is
# not a short sequence, it is the absence of one, and the 3012 dashes it
# carries are indistinguishable from data to anything that only checks
# the string is non-empty — which is how 21,358 all-gap rows reached
# tn93 as `passes_qc=True`. See `docs/evaluations/align-coverage.md`.
BLOCKING_QC_FLAGS = frozenset(
    {
        QCFlag.HIGH_AMBIGUITY,
        QCFlag.APOBEC_HYPERMUT,
        QCFlag.STOP_CODON,
        QCFlag.LOW_COMPLEXITY,
        QCFlag.NO_ALIGNMENT,
    }
)


def passes_qc(flags: list[QCFlag]) -> bool:
    """A sequence passes QC iff no *blocking* flags are set.

    Non-blocking flags (SHORT_COVERAGE) are recorded in `qc_flags` for
    downstream consumers but don't gate ingestion. See BLOCKING_QC_FLAGS.
    """
    return not any(f in BLOCKING_QC_FLAGS for f in flags)
