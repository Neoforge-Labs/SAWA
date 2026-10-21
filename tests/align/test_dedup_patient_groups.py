"""Patient grouping for sequences that carry no name.

`dedup.py`'s docstring promises "a per-accession group_id" for sequences
it cannot tie to a patient. GenBank rows have no `sequence_name` column at
all, so the shared "no_name" token put the whole corpus into one group per
country: the shipped `curated/aligned.parquet` has 6 distinct
`patient_group_id` values across 42,106 rows, and
`is_earliest_per_patient` is true for 6 of them.
"""

from __future__ import annotations

from hivdrwatch.align.dedup import assign_patient_groups


class TestPatientGroupsForUnnamedSequences:
    def test_unnamed_rows_are_independent_patients(self) -> None:
        rows = [
            {"accession": "MK000001", "country": "KE", "collection_year": 2019},
            {"accession": "MK000002", "country": "KE", "collection_year": 2020},
            {"accession": "MK000003", "country": "KE", "collection_year": 2021},
        ]
        out = assign_patient_groups(rows)
        assert len({r["patient_group_id"] for r in out}) == 3
        assert all(r["is_earliest_per_patient"] for r in out)

    def test_named_longitudinal_sequences_still_share_a_group(self) -> None:
        rows = [
            {
                "accession": "AF000001",
                "sequence_name": "KER2003_1998",
                "country": "KE",
                "collection_year": 1998,
            },
            {
                "accession": "AF000002",
                "sequence_name": "KER2003_2001",
                "country": "KE",
                "collection_year": 2001,
            },
        ]
        out = assign_patient_groups(rows)
        assert out[0]["patient_group_id"] == out[1]["patient_group_id"]
        assert out[0]["is_earliest_per_patient"] is True
        assert out[1]["is_earliest_per_patient"] is False
