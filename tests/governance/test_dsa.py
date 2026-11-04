"""Tests for the KEMRI Data Sharing Agreement gate.

The gate is the only thing standing between a future KEMRI adapter and
PHI leaving the partner lab's control, so it is tested at both levels:
the bare `require_dsa()` call and the `SequenceRecord` choke point that
every ingestion path must cross.
"""

from __future__ import annotations

from datetime import date

import pytest
from pydantic import ValidationError

from hivdrwatch.governance import (
    DSA_ENV_VAR,
    DSANotSignedError,
    dsa_is_signed,
    require_dsa,
)
from hivdrwatch.ingest.schema import SequenceRecord, Source


@pytest.fixture(autouse=True)
def _clear_gate(monkeypatch: pytest.MonkeyPatch) -> None:
    """Every test starts with the gate closed, whatever the host env says."""
    monkeypatch.delenv(DSA_ENV_VAR, raising=False)


def _record(**overrides: object) -> dict:
    return {
        "accession": "AB123456",
        "source": Source.LANL,
        "sequence": "ACGT" * 100,
        "length": 400,
        "ingested_at": date(2026, 9, 22),
        "raw_record_hash": "0" * 64,
        **overrides,
    }


class TestDsaIsSigned:
    def test_false_when_unset(self) -> None:
        assert dsa_is_signed() is False

    @pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes", "on"])
    def test_true_for_affirmative_values(
        self, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        monkeypatch.setenv(DSA_ENV_VAR, value)
        assert dsa_is_signed() is True

    @pytest.mark.parametrize("value", ["", "0", "false", "no", "off", "maybe"])
    def test_false_for_everything_else(
        self, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        monkeypatch.setenv(DSA_ENV_VAR, value)
        assert dsa_is_signed() is False


class TestRequireDsa:
    def test_raises_when_unsigned(self) -> None:
        with pytest.raises(DSANotSignedError):
            require_dsa("kemri ingest")

    def test_error_names_the_context(self) -> None:
        with pytest.raises(DSANotSignedError, match="kemri ingest"):
            require_dsa("kemri ingest")

    def test_error_names_the_env_var_so_the_fix_is_obvious(self) -> None:
        with pytest.raises(DSANotSignedError, match=DSA_ENV_VAR):
            require_dsa("kemri ingest")

    def test_passes_when_signed(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv(DSA_ENV_VAR, "1")
        require_dsa("kemri ingest")

    def test_logs_a_warning_when_the_gate_is_open(
        self, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ) -> None:
        monkeypatch.setenv(DSA_ENV_VAR, "1")
        with caplog.at_level("WARNING"):
            require_dsa("kemri ingest")
        assert "kemri ingest" in caplog.text


class TestSequenceRecordGate:
    def test_kemri_source_rejected_when_unsigned(self) -> None:
        with pytest.raises(DSANotSignedError, match=DSA_ENV_VAR):
            SequenceRecord(**_record(source=Source.KEMRI))

    def test_kemri_source_rejected_when_given_as_string(self) -> None:
        with pytest.raises(DSANotSignedError):
            SequenceRecord(**_record(source="kemri"))

    def test_block_is_not_a_validation_error(self) -> None:
        """A governance block must not look like a malformed record.

        Ingestion loops catch `ValidationError` (and CLAUDE.md allows one
        broad `except` there) to skip bad records and keep going. If the
        DSA block arrived as a ValidationError it would be swallowed
        per-record and logged as a parse failure instead of halting the
        run. Pydantic only converts ValueError, so a RuntimeError
        subclass propagates — that is load-bearing, not incidental.
        """
        with pytest.raises(DSANotSignedError) as excinfo:
            SequenceRecord(**_record(source=Source.KEMRI))
        assert not isinstance(excinfo.value, ValidationError)

    def test_kemri_source_accepted_when_signed(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv(DSA_ENV_VAR, "1")
        rec = SequenceRecord(**_record(source=Source.KEMRI))
        assert rec.source is Source.KEMRI

    @pytest.mark.parametrize("source", [Source.LANL, Source.GENBANK])
    def test_public_sources_unaffected(self, source: Source) -> None:
        rec = SequenceRecord(**_record(source=source))
        assert rec.source is source
