"""Tests for the GenBank incremental-ingest watermark file."""

from __future__ import annotations

from datetime import date
from pathlib import Path

from hivdrwatch.ingest.watermark import (
    WATERMARK_FILENAME,
    bump_watermark,
    get_since,
    read_watermarks,
    write_watermarks,
)


class TestReadWrite:
    def test_roundtrip(self, tmp_path: Path) -> None:
        watermarks = {"KE": date(2026, 8, 11), "UG": date(2026, 8, 10)}
        write_watermarks(tmp_path, watermarks)
        assert (tmp_path / WATERMARK_FILENAME).exists()
        assert read_watermarks(tmp_path) == watermarks

    def test_missing_file_returns_empty(self, tmp_path: Path) -> None:
        assert read_watermarks(tmp_path) == {}

    def test_corrupt_file_returns_empty(self, tmp_path: Path) -> None:
        (tmp_path / WATERMARK_FILENAME).write_text("{not json")
        assert read_watermarks(tmp_path) == {}

    def test_bad_date_skipped(self, tmp_path: Path) -> None:
        (tmp_path / WATERMARK_FILENAME).write_text(
            '{"KE": "not-a-date", "UG": "2026-08-11"}'
        )
        result = read_watermarks(tmp_path)
        assert result == {"UG": date(2026, 8, 11)}

    def test_country_codes_uppercased(self, tmp_path: Path) -> None:
        (tmp_path / WATERMARK_FILENAME).write_text('{"ke": "2026-08-11"}')
        assert read_watermarks(tmp_path) == {"KE": date(2026, 8, 11)}


class TestGetSince:
    def test_returns_watermark_when_present(self) -> None:
        wm = {"KE": date(2026, 8, 1)}
        assert get_since(wm, "KE") == date(2026, 8, 1)
        # Case-insensitive
        assert get_since(wm, "ke") == date(2026, 8, 1)

    def test_returns_fallback_when_missing(self) -> None:
        floor = date(1980, 1, 1)
        assert get_since({}, "KE", fallback=floor) == floor

    def test_returns_none_without_fallback(self) -> None:
        assert get_since({}, "KE") is None


class TestBumpWatermark:
    def test_bumps_named_country_only(self) -> None:
        before = {"KE": date(2026, 8, 1), "UG": date(2026, 8, 1)}
        after = bump_watermark(before, "KE", date(2026, 8, 11))
        assert after == {"KE": date(2026, 8, 11), "UG": date(2026, 8, 1)}
        # Doesn't mutate input
        assert before == {"KE": date(2026, 8, 1), "UG": date(2026, 8, 1)}

    def test_case_insensitive(self) -> None:
        after = bump_watermark({}, "ke", date(2026, 8, 11))
        assert after == {"KE": date(2026, 8, 11)}
