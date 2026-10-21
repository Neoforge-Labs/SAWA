"""Pairwise-distance tests: the denominator must be the pairwise overlap.

Every case here is a pair that the alignment-width denominator reports
as closely related and the overlap denominator does not. They are the
regression guards for the defect recorded in
`docs/evaluations/cluster-distance-denominator.md`: dividing the
mismatch count by the full 3,012 alignment columns makes gap-vs-gap
count as agreement, so a pair is linked in proportion to how little of
it was ever sequenced.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from hivdrwatch.cluster import distance
from hivdrwatch.cluster.distance import (
    DEFAULT_MIN_OVERLAP,
    DEFAULT_THRESHOLD,
    compute_pairwise_distances_python,
    read_distances,
)

ALIGNMENT_WIDTH = 3012  # HXB2 pol, the width every curated row carries


def _write_fasta(tmp_path: Path, records: dict[str, str]) -> Path:
    fasta = tmp_path / "aligned.fasta"
    with fasta.open("w", encoding="utf-8") as fh:
        for accession, seq in records.items():
            fh.write(f">{accession}\n{seq}\n")
    return fasta


def _pairs(tmp_path: Path, records: dict[str, str], **kwargs) -> dict[tuple[str, str], float]:
    """Run the engine over `records` and return {(a, b): distance}."""
    fasta = _write_fasta(tmp_path, records)
    out = compute_pairwise_distances_python(fasta, tmp_path / "distances.csv", **kwargs)
    df = read_distances(out)
    return {
        (row["accession_a"], row["accession_b"]): row["distance"]
        for row in df.iter_rows(named=True)
    }


def _padded(prefix: str, width: int = ALIGNMENT_WIDTH) -> str:
    """A sequence occupying the first `len(prefix)` columns, gapped after."""
    return prefix + "-" * (width - len(prefix))


def _window(seed: str, n: int) -> str:
    """`n` nucleotides of a repeating, non-degenerate pattern."""
    return (seed * ((n // len(seed)) + 1))[:n]


class TestAllGapPairs:
    """A row of 3,012 dashes is the absence of a sequence, not a sequence."""

    def test_all_gap_pair_is_not_emitted(self, tmp_path: Path) -> None:
        # Under the alignment-width denominator both rows are identical
        # at all 3,012 columns, so the pair reports distance 0.000000
        # and links into the same cluster.
        records = {
            "GAP_A": "-" * ALIGNMENT_WIDTH,
            "GAP_B": "-" * ALIGNMENT_WIDTH,
        }
        assert _pairs(tmp_path, records) == {}

    def test_all_gap_row_links_to_nothing(self, tmp_path: Path) -> None:
        real = _window("ACGT", 2000)
        records = {
            "GAP_A": "-" * ALIGNMENT_WIDTH,
            "REAL_A": _padded(real),
            "REAL_B": _padded(real),
        }
        pairs = _pairs(tmp_path, records)
        # The two real sequences are identical over 2,000 shared columns.
        assert pairs == {("REAL_A", "REAL_B"): pytest.approx(0.0)}


class TestShortOverlap:
    """Divergence is measured where both sequences exist, or not at all."""

    def test_300nt_overlap_at_10pct_divergence_is_not_linked(self, tmp_path: Path) -> None:
        # 300 shared columns, 30 of them mismatching: true divergence
        # 10%, nearly 7x the 1.5% threshold. The alignment-width
        # denominator reports 30/3012 = 0.00996 and links the pair.
        base = list(_window("ACGT", 300))
        diverged = list(base)
        for i in range(0, 300, 10):
            diverged[i] = "G" if base[i] != "G" else "A"
        records = {
            "SHORT_A": _padded("".join(base)),
            "SHORT_B": _padded("".join(diverged)),
        }
        assert _pairs(tmp_path, records) == {}

    def test_overlap_divergence_is_reported_when_above_threshold(self, tmp_path: Path) -> None:
        base = list(_window("ACGT", 300))
        diverged = list(base)
        for i in range(0, 300, 10):
            diverged[i] = "G" if base[i] != "G" else "A"
        records = {
            "SHORT_A": _padded("".join(base)),
            "SHORT_B": _padded("".join(diverged)),
        }
        pairs = _pairs(tmp_path, records, threshold=0.2)
        assert pairs == {("SHORT_A", "SHORT_B"): pytest.approx(0.1, abs=1e-6)}


class TestMinimumOverlap:
    """Below `min_overlap` shared sites, no distance is emitted at all."""

    def test_pair_below_default_min_overlap_is_refused(self, tmp_path: Path) -> None:
        # 99 shared columns, identical in all of them. Both denominators
        # call this distance 0.0; tn93's own `-l` default refuses the
        # pair rather than reporting a distance estimated from 99 sites.
        shared = _window("ACGT", DEFAULT_MIN_OVERLAP - 1)
        records = {
            "THIN_A": _padded(shared),
            "THIN_B": _padded(shared),
        }
        assert _pairs(tmp_path, records) == {}

    def test_pair_at_min_overlap_is_emitted(self, tmp_path: Path) -> None:
        shared = _window("ACGT", DEFAULT_MIN_OVERLAP)
        records = {
            "THIN_A": _padded(shared),
            "THIN_B": _padded(shared),
        }
        pairs = _pairs(tmp_path, records)
        assert pairs == {("THIN_A", "THIN_B"): pytest.approx(0.0)}

    def test_min_overlap_is_a_parameter(self, tmp_path: Path) -> None:
        shared = _window("ACGT", 50)
        records = {
            "THIN_A": _padded(shared),
            "THIN_B": _padded(shared),
        }
        assert _pairs(tmp_path, records, min_overlap=50) != {}
        assert _pairs(tmp_path, records, min_overlap=51) == {}

    def test_min_overlap_must_be_positive(self, tmp_path: Path) -> None:
        records = {"A": _padded(_window("ACGT", 500)), "B": _padded(_window("ACGT", 500))}
        fasta = _write_fasta(tmp_path, records)
        with pytest.raises(ValueError, match="min_overlap"):
            compute_pairwise_distances_python(fasta, tmp_path / "distances.csv", min_overlap=0)


class TestAmbiguityHandling:
    """Ambiguous sites are not evidence of agreement either."""

    def test_ambiguous_sites_leave_the_denominator(self, tmp_path: Path) -> None:
        # 200 unambiguous shared columns with 6 mismatches (3%), plus
        # 400 columns where one side is N. Counting the N columns as
        # agreement would report 6/600 = 0.01 and link the pair.
        base = list(_window("ACGT", 200))
        diverged = list(base)
        for i in range(0, 180, 30):  # exactly 6 mismatching sites
            diverged[i] = "G" if base[i] != "G" else "A"
        filler = _window("ACGT", 400)
        records = {
            "AMB_A": _padded("".join(base) + filler),
            "AMB_B": _padded("".join(diverged) + "N" * 400),
        }
        assert _pairs(tmp_path, records) == {}
        pairs = _pairs(tmp_path, records, threshold=0.2)
        assert pairs == {("AMB_A", "AMB_B"): pytest.approx(6 / 200, abs=1e-6)}


class TestHarnessRegression:
    """The gap-free synthetic case the validation harness runs must not move."""

    def test_full_length_identical_pair(self, tmp_path: Path) -> None:
        seq = _window("ACGTTGCA", ALIGNMENT_WIDTH)
        pairs = _pairs(tmp_path, {"SIM_A": seq, "SIM_B": seq})
        assert pairs == {("SIM_A", "SIM_B"): pytest.approx(0.0)}

    def test_full_length_pair_at_one_percent(self, tmp_path: Path) -> None:
        base = list(_window("ACGTTGCA", ALIGNMENT_WIDTH))
        diverged = list(base)
        n_diff = 30  # 30/3012 ≈ 0.00996, below the 1.5% threshold
        for i in range(n_diff):
            pos = i * 100
            diverged[pos] = "G" if base[pos] != "G" else "A"
        pairs = _pairs(tmp_path, {"SIM_A": "".join(base), "SIM_B": "".join(diverged)})
        # Gap-free and unambiguous: overlap is the full alignment, so
        # the two denominators agree exactly.
        assert pairs == {("SIM_A", "SIM_B"): pytest.approx(n_diff / ALIGNMENT_WIDTH, abs=1e-6)}

    def test_threshold_default_is_unchanged(self) -> None:
        assert DEFAULT_THRESHOLD == 0.015


class TestEngineGuards:
    """Input the engine must reject rather than guess at."""

    def test_missing_fasta_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="aligned FASTA missing"):
            compute_pairwise_distances_python(
                tmp_path / "absent.fasta", tmp_path / "distances.csv"
            )

    def test_single_sequence_emits_only_a_header(self, tmp_path: Path) -> None:
        fasta = _write_fasta(tmp_path, {"ONLY": _padded(_window("ACGT", 2000))})
        out = compute_pairwise_distances_python(fasta, tmp_path / "distances.csv")
        assert out.read_text() == "ID1,ID2,Distance\n"

    def test_mixed_lengths_raise(self, tmp_path: Path) -> None:
        fasta = _write_fasta(
            tmp_path, {"A": _window("ACGT", 3012), "B": _window("ACGT", 2000)}
        )
        with pytest.raises(ValueError, match="mixed lengths"):
            compute_pairwise_distances_python(fasta, tmp_path / "distances.csv")

    def test_lowercase_sequence_is_read_as_sequence(self, tmp_path: Path) -> None:
        seq = _window("acgt", ALIGNMENT_WIDTH)
        pairs = _pairs(tmp_path, {"LOWER_A": seq, "LOWER_B": seq.upper()})
        assert pairs == {("LOWER_A", "LOWER_B"): pytest.approx(0.0)}


class TestReferenceInvocation:
    """`run_tn93` is the reference path; pin the command line it builds.

    Never runs Docker — `subprocess.run` is replaced. What matters here
    is that no `-l` is passed, so tn93 applies its own 100-site
    minimum-overlap default, which is the semantics
    `DEFAULT_MIN_OVERLAP` matches.
    """

    @staticmethod
    def _fake_run(returncode: int = 0, write: str = "ID1,ID2,Distance\nA,B,0.001\n"):
        calls: list[list[str]] = []

        def runner(cmd: list[str], **kwargs: object) -> SimpleNamespace:
            calls.append(cmd)
            out_flag = cmd.index("-o")
            name = Path(cmd[out_flag + 1]).name
            mount = cmd[cmd.index("-v") + 1].split(":")[0]
            (Path(mount) / name).write_text(write)
            return SimpleNamespace(returncode=returncode, stdout=b"", stderr=b"")

        return runner, calls

    def test_builds_the_documented_command(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fasta = _write_fasta(tmp_path, {"A": _window("ACGT", 500)})
        runner, calls = self._fake_run()
        monkeypatch.setattr(distance.subprocess, "run", runner)

        invocation = distance.run_tn93(fasta, tmp_path / "distances.csv", threshold=0.015)

        assert invocation.args == (
            "-t",
            "0.015",
            "-a",
            "resolve",
            "-o",
            "/work/distances.csv",
            "/work/aligned.fasta",
        )
        assert "-l" not in invocation.args  # tn93's own default applies
        assert calls[0][:3] == ["docker", "run", "--rm"]

    def test_nonzero_exit_raises(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        fasta = _write_fasta(tmp_path, {"A": _window("ACGT", 500)})
        runner, _ = self._fake_run(returncode=2)
        monkeypatch.setattr(distance.subprocess, "run", runner)
        with pytest.raises(RuntimeError, match="tn93 failed"):
            distance.run_tn93(fasta, tmp_path / "distances.csv")

    def test_missing_fasta_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="aligned FASTA missing"):
            distance.run_tn93(tmp_path / "absent.fasta", tmp_path / "distances.csv")

    def test_output_must_share_the_mounted_directory(self, tmp_path: Path) -> None:
        fasta = _write_fasta(tmp_path, {"A": _window("ACGT", 500)})
        with pytest.raises(ValueError, match="share a parent"):
            distance.run_tn93(fasta, tmp_path / "elsewhere" / "distances.csv")


class TestPinnedReferenceImage:
    """The reference image's tag, its source pin and its reader must agree.

    Both of the defects recorded in
    `docs/evaluations/cluster-distance-denominator.md` §11 were
    invisible to the test suite. The image was tagged `1.0.13` while the
    binary inside reported `v1.0.12`, and the binary could not read a
    multi-record FASTA at all. The Dockerfile now asserts both at build
    time; these two guards are the parts that need no Docker daemon, so
    they run in the ordinary unit pass.
    """

    # First upstream release carrying veg/tn93 commit 5c14f35, which
    # changed `char validFlags[256]` to `signed char`. Before it, the
    # -1 "not a sequence character" sentinel read back as 255 wherever
    # plain `char` is unsigned (every aarch64 host), so `>` tested as a
    # valid base and the FASTA reader never saw a record boundary.
    FIRST_FIXED_VERSION = (1, 0, 15)

    @staticmethod
    def _dockerfile_version() -> str:
        dockerfile = Path(__file__).resolve().parents[2] / "docker" / "hivtrace" / "Dockerfile"
        for line in dockerfile.read_text().splitlines():
            if line.startswith("ARG TN93_VERSION="):
                return line.split("=", 1)[1].strip()
        raise AssertionError(f"no ARG TN93_VERSION in {dockerfile}")

    def test_default_image_tag_matches_the_dockerfile_pin(self) -> None:
        pinned = f"hivdrwatch/tn93:{self._dockerfile_version()}"
        assert pinned == distance.DEFAULT_IMAGE

    def test_pin_is_at_or_past_the_aarch64_fasta_reader_fix(self) -> None:
        pinned = tuple(int(p) for p in self._dockerfile_version().split("."))
        assert pinned >= self.FIRST_FIXED_VERSION, (
            f"tn93 {self._dockerfile_version()} predates the signed-char fix; "
            "the reference binary would read only the first FASTA record "
            "on any aarch64 build host and silently emit zero pairs"
        )
