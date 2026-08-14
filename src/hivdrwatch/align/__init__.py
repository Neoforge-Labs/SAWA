"""Silver-layer curation: alignment, gene-region trim, QC, patient dedup."""

from .align import align_sequences, read_fasta, run_mafft, write_fasta
from .dedup import assign_patient_groups, extract_patient_token
from .qc import (
    ambiguity_fraction,
    apobec_signature_ratio,
    has_premature_stop,
    passes_qc,
    qc_flags,
    translate_in_frame,
)
from .reference import (
    HXB2_POL_LENGTH,
    HXB2_POL_START,
    HXB2_POL_STOP,
    POL_REGIONS,
    HXB2Reference,
    ensure_hxb2_reference,
    extract_hxb2_from_lanl_fasta,
)
from .schema import CURATE_SCHEMA_VERSION, CuratedSequence, QCFlag
from .trim import coverage_fraction, slice_all_regions, slice_region

__all__ = [
    "CURATE_SCHEMA_VERSION",
    "HXB2_POL_LENGTH",
    "HXB2_POL_START",
    "HXB2_POL_STOP",
    "POL_REGIONS",
    "CuratedSequence",
    "HXB2Reference",
    "QCFlag",
    "align_sequences",
    "ambiguity_fraction",
    "apobec_signature_ratio",
    "assign_patient_groups",
    "coverage_fraction",
    "ensure_hxb2_reference",
    "extract_hxb2_from_lanl_fasta",
    "extract_patient_token",
    "has_premature_stop",
    "passes_qc",
    "qc_flags",
    "read_fasta",
    "run_mafft",
    "slice_all_regions",
    "slice_region",
    "translate_in_frame",
    "write_fasta",
]
