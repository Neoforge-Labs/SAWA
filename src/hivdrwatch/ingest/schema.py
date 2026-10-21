"""Canonical schema for HIV sequence records.

Every ingestion source (GenBank, LANL, KEMRI) must produce records that
validate against SequenceRecord. This is the contract that keeps downstream
pipelines source-agnostic.

Schema versioning: bump SCHEMA_VERSION on any breaking change. The version
is written into every Parquet output so future analyses can detect drift.

Migration notes
---------------
**0.1.0 → 0.2.0 — `subtype` and `country` semantics.**

Three ingest provenance defects changed what these two fields mean, so
0.1.0 and 0.2.0 records must not be pooled.

1. `subtype` — `genbank._normalize_subtype` matched the declared subtype
   by substring, so `"SUBTYPE: C"` matched inside `"SUBTYPE: CRF01_AE"`
   and every CRF was stored as subtype C. Matching is now token-based,
   and the free-text forms the depositors actually use are parsed
   (`genotype:`/`clade:` keys, mosaic shorthands, COMET confidences).
2. `subtype` — `A` and `F` are new members. The enum had A1/A2 and F1/F2
   but no parent clade, so the ~1,500 records whose note says only
   `subtype: A` were stored as `unknown`. Rolling A1/A2 up to A is a
   downstream decision; ingest no longer has to choose.
3. `country` — it was the country of the *query* that returned the
   record, not the country the record states. The source feature is now
   read (including NCBI's `/geo_loc_name=`, which replaced `/country=`
   and is the only form present in SAWA's corpus) and there is no
   fallback to the query. Two new optional fields record where the value
   came from: `country_raw` (verbatim, keeping any subdivision) and
   `country_source`.

A 0.1.0 vintage cannot be migrated in place: the correct subtype and
country are only recoverable from the source records, so re-run ingest
against `data/raw/genbank`.
"""

from __future__ import annotations

import re
from datetime import date
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, Field, StringConstraints, field_validator

SCHEMA_VERSION = "0.2.0"


class Subtype(StrEnum):
    """HIV-1 subtypes and common East African CRFs.

    Not exhaustive — extend as needed. Unknown falls through to `UNKNOWN`
    rather than raising, so ingestion never blocks on annotation gaps.
    """

    A = "A"  # parent clade, used when the depositor did not say A1 or A2
    A1 = "A1"
    A2 = "A2"
    B = "B"
    C = "C"
    D = "D"
    F = "F"  # parent clade, as above
    F1 = "F1"
    F2 = "F2"
    G = "G"
    H = "H"
    J = "J"
    K = "K"
    CRF01_AE = "CRF01_AE"
    CRF02_AG = "CRF02_AG"
    CRF10_CD = "CRF10_CD"
    CRF16_A2D = "CRF16_A2D"
    CRF21_A2D = "CRF21_A2D"
    CRF35_AD = "CRF35_AD"
    RECOMBINANT = "recombinant"  # catch-all for other CRFs/URFs
    UNKNOWN = "unknown"


class Source(StrEnum):
    GENBANK = "genbank"
    LANL = "lanl"
    KEMRI = "kemri"


# Subtype letters that have appeared in a mosaic shorthand ("AD", "A1C",
# "AE"). Wider than the enum: E and I are not standalone subtypes but the
# shorthands that name them are in the literature and in GenBank notes.
_MOSAIC_LETTERS = "ABCDEFGHIJKL"

_SUBTYPE_BY_VALUE: dict[str, Subtype] = {st.value.upper(): st for st in Subtype}

# Tokens that sit inside a subtype declaration without naming a clade.
_DESIGNATION_NOISE = frozenset({"HIV", "GROUP", "AND", "OR", "SUBTYPE", "GENOTYPE", "CLADE"})

# Words that assert a mosaic genome without naming an enumerated form.
_MOSAIC_WORDS = frozenset(
    {"RECOMBINANT", "RECOMBINANTS", "MOSAIC", "INTERSUBTYPE", "URF", "URFS"}
)

# "CRF 10_CD", "CRF_10_CD" and "CRF10_CD" are the same designation.
_CRF_SPACING = re.compile(r"\b(CRF|URF)[_ ]*(\d+)[_ ]*", re.IGNORECASE)
_DESIGNATION_TOKEN = re.compile(r"[A-Za-z0-9_]+")
# "HIV-1C" and "1C" both mean subtype C.
_LEADING_HIV1 = re.compile(r"^(?:HIV[-_]?)?1(?=[A-K])")
_MOSAIC_SHORTHAND = re.compile(rf"^(?:[{_MOSAIC_LETTERS}][0-9]?){{2,}}$")


def subtype_from_designation(text: str | None) -> Subtype:
    """Classify a subtype *designation* — the bit after `subtype:`.

    Shared by both ingesters so a GenBank note and a LANL Subtype column
    cannot disagree about the same string. Scans left to right and stops
    at the first token that is not part of a designation, which is what
    keeps trailing prose ("subtype C from Kenya") out of the answer.

    Two or more named parents, or any mosaic marker, resolve to
    `RECOMBINANT`; a named CRF resolves to itself because the enum
    carries more information than the catch-all does.
    """
    if not text:
        return Subtype.UNKNOWN

    normalized = _CRF_SPACING.sub(lambda m: f"{m.group(1).upper()}{m.group(2)}_", text)
    named: list[Subtype] = []
    mosaic = False

    for raw_token in _DESIGNATION_TOKEN.findall(normalized):
        token = _LEADING_HIV1.sub("", raw_token.upper())
        if not token or token.isdigit() or token in _DESIGNATION_NOISE:
            continue
        if token in _MOSAIC_WORDS:
            mosaic = True
            continue
        member = _SUBTYPE_BY_VALUE.get(token)
        if member is Subtype.UNKNOWN:
            break  # the depositor explicitly said they do not know
        if member is not None:
            named.append(member)
            continue
        if token.startswith(("CRF", "URF")):
            # A CRF with no enum member of its own, or a truncated one.
            mosaic = True
            continue
        # Shorthand mosaics are written in upper case; the same letters in
        # lower case are prose ("subtype had not been determined").
        if raw_token.isupper() and _MOSAIC_SHORTHAND.match(token):
            mosaic = True
            continue
        break

    # A parent clade alongside one of its own sub-clades ("A (A1)") is the
    # depositor being more specific, not a list of mosaic parents.
    distinct = {
        st
        for st in named
        if not any(
            other.value != st.value
            and other.value.startswith(st.value)
            and other.value[len(st.value) :].isdigit()
            for other in named
        )
    }

    if len(distinct) == 1:
        only = next(iter(distinct))
        if only.value.startswith("CRF") or not mosaic:
            return only
        return Subtype.RECOMBINANT
    if mosaic or len(distinct) > 1:
        return Subtype.RECOMBINANT
    return Subtype.UNKNOWN


# GenBank's `/country=` (now `/geo_loc_name=`) and LANL's Country column,
# lower-cased, → ISO 3166-1 alpha-2. Covers every name present in SAWA's
# corpus plus the neighbours a sub-Saharan pull is likely to reach; extend
# as new geographies arrive. An unmapped name yields no country rather
# than a guess.
COUNTRY_NAME_TO_CODE: dict[str, str] = {
    # East Africa — the priority cohort
    "kenya": "KE",
    "uganda": "UG",
    "tanzania": "TZ",
    "united republic of tanzania": "TZ",
    "rwanda": "RW",
    "burundi": "BI",
    "ethiopia": "ET",
    "south sudan": "SS",
    "somalia": "SO",
    "eritrea": "ER",
    "djibouti": "DJ",
    "sudan": "SD",
    # Rest of sub-Saharan Africa
    "south africa": "ZA",
    "zambia": "ZM",
    "zimbabwe": "ZW",
    "malawi": "MW",
    "mozambique": "MZ",
    "botswana": "BW",
    "namibia": "NA",
    "angola": "AO",
    "nigeria": "NG",
    "cameroon": "CM",
    "cote d'ivoire": "CI",
    "côte d'ivoire": "CI",
    "ivory coast": "CI",
    "ghana": "GH",
    "senegal": "SN",
    "gabon": "GA",
    "republic of the congo": "CG",
    "congo": "CG",
    "democratic republic of the congo": "CD",
    "madagascar": "MG",
    # Elsewhere — present in the corpus because the country search
    # over-matches, and kept so the records say where they are from.
    "thailand": "TH",
    "india": "IN",
    "israel": "IL",
    "brazil": "BR",
    "peru": "PE",
    "argentina": "AR",
    "cuba": "CU",
    "haiti": "HT",
    "dominican republic": "DO",
    "spain": "ES",
    "italy": "IT",
    "finland": "FI",
    "sweden": "SE",
    "united kingdom": "GB",
    "uk": "GB",
    "united states": "US",
    "usa": "US",
}


def country_name_to_code(raw: str | None) -> str | None:
    """Country name → ISO 3166-1 alpha-2, or None when unmapped."""
    if not raw:
        return None
    return COUNTRY_NAME_TO_CODE.get(raw.strip().lower())


# ISO 3166-1 alpha-2. Not validated against the full list here — the ingester
# normalizes country names to codes and warns on unknowns.
CountryCode = Annotated[str, StringConstraints(min_length=2, max_length=2, to_upper=True)]

# GenBank accession pattern: 1-2 letters + 5-8 digits, optional .version
AccessionStr = Annotated[
    str, StringConstraints(pattern=r"^[A-Z]{1,2}\d{5,8}(\.\d+)?$", to_upper=True)
]


class SequenceRecord(BaseModel):
    """One HIV sequence with the metadata phylodynamics actually needs."""

    # --- identity ---
    accession: AccessionStr
    source: Source
    schema_version: str = SCHEMA_VERSION

    # --- sequence ---
    sequence: str = Field(..., min_length=200)  # anything shorter is noise
    length: int
    gene_region: str = "pol"  # pol only for HIV-DR-Watch v1

    # --- provenance ---
    country: CountryCode | None = None
    # Verbatim country as the source states it, subdivision included
    # ("Kenya: Nairobi"), so an unmapped name is still on the record.
    country_raw: str | None = None
    # Where `country` came from: "genbank_source_feature", "lanl_curated",
    # or None when no source asserted one. Never the search query.
    country_source: str | None = None
    collection_date: date | None = None
    collection_year: int | None = Field(None, ge=1980, le=2100)
    subtype: Subtype = Subtype.UNKNOWN
    subtype_source: str | None = None  # "genbank", "lanl_curated", "comet", "scueal"

    # --- clinical context (rarely available in public data) ---
    is_drug_naive: bool | None = None
    risk_group: str | None = None
    patient_id_hashed: str | None = None  # for deduping same-patient sequences

    # --- publication trail ---
    pubmed_id: str | None = None
    study_accession: str | None = None  # e.g. PANGEA cohort ID

    # --- ingestion trail ---
    ingested_at: date
    raw_record_hash: str  # SHA-256 of the source record for change detection

    @field_validator("sequence")
    @classmethod
    def _uppercase_and_check_alphabet(cls, v: str) -> str:
        v = v.upper().replace("U", "T")
        # IUPAC nucleotide codes including ambiguities and gap
        allowed = set("ACGTRYSWKMBDHVN-")
        bad = set(v) - allowed
        if bad:
            raise ValueError(f"sequence contains non-IUPAC characters: {sorted(bad)}")
        return v

    @field_validator("length")
    @classmethod
    def _length_matches_sequence(cls, v: int, info) -> int:
        seq = info.data.get("sequence")
        if seq is not None and len(seq) != v:
            raise ValueError(f"length {v} does not match sequence length {len(seq)}")
        return v
