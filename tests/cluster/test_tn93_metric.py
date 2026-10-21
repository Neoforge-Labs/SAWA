"""The opt-in Tamura-Nei 1993 metric, and the default it must not disturb.

`compute_pairwise_distances_python` scores a pair either way over the
same overlap: `p-distance` (the default, the uncorrected mismatch
fraction) or `tn93` (the correction HIV-TRACE thresholds). The cut
Paper 1 claims is still Kelyn's decision, so half of what these tests
guard is that asking for nothing changed anything
(`docs/evaluations/cluster-distance-denominator.md` §12, §13).

The correction itself is checked three ways: against a scalar
transcription of the published formula, against two values small enough
to work out by hand, and — behind the `reference_image` marker, which
the default run deselects — against the pinned tn93 image itself.
"""

from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

import pytest

from hivdrwatch.cluster import distance
from hivdrwatch.cluster.distance import (
    DEFAULT_DISTANCE_METRIC,
    DEFAULT_IMAGE,
    DISTANCE_METRICS,
    compute_pairwise_distances_python,
)

# The corrected path writes six significant figures, which is tn93's
# own `printf("%g")`, so a value read back from the CSV sits within
# 5e-6 of the computed one in relative terms.
PRINTED_PRECISION = 1e-5

PURINES = {"A": "G", "G": "A"}
PYRIMIDINES = {"C": "T", "T": "C"}


def _write_fasta(tmp_path: Path, records: dict[str, str], name: str = "aligned.fasta") -> Path:
    fasta = tmp_path / name
    with fasta.open("w", encoding="utf-8") as fh:
        for accession, seq in records.items():
            fh.write(f">{accession}\n{seq}\n")
    return fasta


def _rows(csv_path: Path) -> dict[frozenset[str], float]:
    with csv_path.open(encoding="utf-8") as fh:
        reader = csv.reader(fh)
        assert [h.lower() for h in next(reader)] == ["id1", "id2", "distance"]
        return {frozenset((a, b)): float(d) for a, b, d in reader}


def _distances(
    tmp_path: Path, records: dict[str, str], **kwargs: object
) -> dict[frozenset[str], float]:
    fasta = _write_fasta(tmp_path, records)
    out = compute_pairwise_distances_python(fasta, tmp_path / "distances.csv", **kwargs)  # type: ignore[arg-type]
    return _rows(out)


def _one(tmp_path: Path, a: str, b: str, **kwargs: object) -> float | None:
    """The distance reported for a single pair, or None if it got no row."""
    rows = _distances(tmp_path, {"A": a, "B": b}, **kwargs)
    return rows.get(frozenset(("A", "B")))


def _tn93_scalar(a: str, b: str) -> float:
    """Tamura-Nei 1993 for one pair, transcribed from the published formula.

    Deliberately an independent implementation — scalar, dict-counted,
    no numpy — so that agreement with the vectorized engine is evidence
    about the engine rather than about shared code. Ambiguous and gap
    columns are skipped, which is tn93's `-a skip`.
    """
    cells: dict[tuple[str, str], int] = {}
    for x, y in zip(a.upper(), b.upper(), strict=True):
        if x in "ACGT" and y in "ACGT":
            cells[(x, y)] = cells.get((x, y), 0) + 1

    sites = sum(cells.values())
    if sites == 0:
        return math.inf

    counts = {base: 0 for base in "ACGT"}
    for (x, y), count in cells.items():
        counts[x] += count
        counts[y] += count
    freq = {base: counts[base] / (2 * sites) for base in "ACGT"}

    purine = sum(cells.get(pair, 0) for pair in (("A", "G"), ("G", "A"))) / sites
    pyrimidine = sum(cells.get(pair, 0) for pair in (("C", "T"), ("T", "C"))) / sites
    identical = sum(cells.get((base, base), 0) for base in "ACGT") / sites
    transversion = 1.0 - identical - purine - pyrimidine

    if min(counts.values()) == 0:  # a base the pair never shows: Kimura 2-parameter
        w_transition = 1.0 - 2.0 * (purine + pyrimidine) - transversion
        w_transversion = 1.0 - 2.0 * transversion
        if w_transition <= 0.0 or w_transversion <= 0.0:
            return math.inf
        return -0.5 * math.log(w_transition) - 0.25 * math.log(w_transversion)

    f_r = freq["A"] + freq["G"]
    f_y = freq["C"] + freq["T"]
    k1 = 2.0 * freq["A"] * freq["G"] / f_r
    k2 = 2.0 * freq["C"] * freq["T"] / f_y
    k3 = 2.0 * (
        f_r * f_y - freq["A"] * freq["G"] * f_y / f_r - freq["C"] * freq["T"] * f_r / f_y
    )
    w_purine = 1.0 - purine / k1 - 0.5 * transversion / f_r
    w_pyrimidine = 1.0 - pyrimidine / k2 - 0.5 * transversion / f_y
    w_transversion = 1.0 - 0.5 * transversion / (f_y * f_r)
    if w_purine <= 0.0 or w_pyrimidine <= 0.0 or w_transversion <= 0.0:
        return math.inf
    return (
        -k1 * math.log(w_purine) - k2 * math.log(w_pyrimidine) - k3 * math.log(w_transversion)
    )


def _mutate(seq: str, positions: dict[int, str]) -> str:
    chars = list(seq)
    for position, base in positions.items():
        chars[position] = base
    return "".join(chars)


def _composition(pattern: str, width: int) -> str:
    return (pattern * ((width // len(pattern)) + 1))[:width]


class TestDefaultMetricUnchanged:
    """Asking for nothing must get exactly what it got before.

    The cut is a Paper 1 claim awaiting Kelyn (§12.4), so the
    correction is reachable only by naming it.
    """

    def test_default_is_the_uncorrected_p_distance(self) -> None:
        assert DEFAULT_DISTANCE_METRIC == "p-distance"
        assert DISTANCE_METRICS == ("p-distance", "tn93")

    def test_default_matches_the_named_default(self, tmp_path: Path) -> None:
        records = {
            "A": _composition("ACGTGCAT", 900),
            "B": _mutate(_composition("ACGTGCAT", 900), {7: "G", 99: "A", 451: "C"}),
        }
        implicit = compute_pairwise_distances_python(
            _write_fasta(tmp_path, records), tmp_path / "implicit.csv"
        )
        explicit = compute_pairwise_distances_python(
            _write_fasta(tmp_path, records), tmp_path / "explicit.csv", metric="p-distance"
        )
        assert implicit.read_bytes() == explicit.read_bytes()

    def test_default_reports_the_mismatch_fraction_not_the_correction(
        self, tmp_path: Path
    ) -> None:
        """Nine transitions in 900 sites is p = 0.01 and TN93 > 0.01."""
        reference = _composition("ACGTGCAT", 900)
        mutated = _mutate(
            reference,
            {
                position: PURINES.get(reference[position], "A")
                for position in range(0, 900, 100)
            },
        )
        p_distance = _one(tmp_path, reference, mutated)
        corrected = _one(tmp_path, reference, mutated, metric="tn93")
        assert p_distance == pytest.approx(0.01, abs=5e-7)
        assert corrected is not None
        assert corrected > p_distance > 0.0

    def test_unknown_metric_is_refused(self, tmp_path: Path) -> None:
        fasta = _write_fasta(tmp_path, {"A": _composition("ACGT", 500)})
        with pytest.raises(ValueError, match="unknown distance metric"):
            compute_pairwise_distances_python(
                fasta,
                tmp_path / "distances.csv",
                metric="tamura-nei",  # type: ignore[arg-type]
            )


class TestCorrectionOnHandComputableFixtures:
    """Two cases whose arithmetic fits in a comment, plus the formula itself."""

    def test_identical_sequences_are_zero(self, tmp_path: Path) -> None:
        seq = _composition("ACGTGCAT", 900)
        assert _one(tmp_path, seq, seq, metric="tn93") == pytest.approx(0.0, abs=1e-12)

    def test_purine_only_alignment_reduces_to_kimura(self, tmp_path: Path) -> None:
        """A pair with no pyrimidine at all: TN93's own K2P fallback.

        1,000 columns of A and G, two of them an A↔G transition and no
        transversions, so P = 0.002, Q = 0 and
        d = −½ ln(1 − 2(0.002)) = −½ ln(0.996) = 0.0020040...
        """
        reference = _composition("AG", 1000)
        mutated = _mutate(reference, {10: "G", 11: "A"})
        expected = -0.5 * math.log(0.996)
        assert expected == pytest.approx(0.00200401, abs=1e-8)
        assert _one(tmp_path, reference, mutated, metric="tn93") == pytest.approx(
            expected, rel=PRINTED_PRECISION
        )

    def test_transversion_only_alignment_is_the_published_k2p(self, tmp_path: Path) -> None:
        """A/C only, four transversions in 1,000: P = 0, Q = 0.004.

        d = −½ ln(1 − 0.004) − ¼ ln(1 − 0.008) = 0.0040120...
        """
        reference = _composition("AC", 1000)
        mutated = _mutate(reference, {0: "C", 1: "A", 2: "C", 3: "A"})
        expected = -0.5 * math.log(0.996) - 0.25 * math.log(0.992)
        assert expected == pytest.approx(0.00401205, abs=1e-8)
        assert _one(tmp_path, reference, mutated, metric="tn93") == pytest.approx(
            expected, rel=PRINTED_PRECISION
        )

    @pytest.mark.parametrize(
        ("pattern", "edits"),
        [
            ("ACGTGCAT", {3: "C", 17: "A", 400: "G"}),
            ("AAACGGGT", {5: "A", 6: "T", 700: "C", 701: "G"}),
            ("ACGT", {12: "T", 13: "G", 14: "C", 15: "A", 888: "A"}),
            ("ACCGGTTA", {1: "N", 2: "R", 300: "-", 301: "A", 500: "T"}),
        ],
    )
    def test_agrees_with_the_scalar_formula(
        self, tmp_path: Path, pattern: str, edits: dict[int, str]
    ) -> None:
        reference = _composition(pattern, 900)
        mutated = _mutate(reference, edits)
        reported = _one(tmp_path, reference, mutated, metric="tn93")
        assert reported is not None
        assert reported == pytest.approx(
            _tn93_scalar(reference, mutated), rel=PRINTED_PRECISION
        )

    def test_the_correction_is_not_a_function_of_p_alone(self, tmp_path: Path) -> None:
        """§12.2: two pairs, same mismatch fraction, different corrections.

        This is why no recalibrated p-threshold reproduces the TN93 pair
        set, and so why the correction had to be computed rather than
        approximated.
        """
        reference = _composition("ACGTGCAT", 800)
        transitions = _mutate(
            reference,
            {
                position: PURINES.get(reference[position], "A")
                for position in range(0, 800, 100)
            },
        )
        transversions = _mutate(
            reference,
            {
                position: ("C" if reference[position] in PURINES else "A")
                for position in range(0, 800, 100)
            },
        )
        by_transition = _one(tmp_path, reference, transitions, metric="tn93")
        by_transversion = _one(tmp_path, reference, transversions, metric="tn93")
        p_transition = _one(tmp_path, reference, transitions)
        p_transversion = _one(tmp_path, reference, transversions)

        assert p_transition == p_transversion
        assert by_transition is not None and by_transversion is not None
        assert by_transition != pytest.approx(by_transversion, rel=1e-4)


class TestNonPositiveLogarithmGuard:
    """A pair whose logarithm leaves the domain gets no row, not a number.

    tn93 reports its `TN93_MAX_DIST` sentinel (1000.0) for these, which
    is above any surveillance threshold; refusing the row is the same
    observable outcome and the same treatment a `min_overlap` refusal
    gets.
    """

    def test_every_site_a_transversion_is_refused(self, tmp_path: Path) -> None:
        """A-vs-C throughout: K2P's 1 − 2Q = −1, outside the domain."""
        assert _one(tmp_path, "A" * 900, "C" * 900, metric="tn93", threshold=10.0) is None

    def test_every_site_a_transition_is_refused(self, tmp_path: Path) -> None:
        """All four bases present, so the TN93 branch runs: 1 − P/k1 = −1."""
        reference = _composition("ACGT", 900)
        swapped = "".join(PURINES.get(base) or PYRIMIDINES[base] for base in reference)
        assert _one(tmp_path, reference, swapped, metric="tn93", threshold=10.0) is None

    def test_a_refused_pair_does_not_suppress_its_neighbours(self, tmp_path: Path) -> None:
        """One out-of-domain pair in a trio leaves the other two reported."""
        reference = _composition("ACGTGCAT", 900)
        swapped = "".join(PURINES.get(base) or PYRIMIDINES[base] for base in reference)
        near = _mutate(reference, {44: PURINES.get(reference[44], "A")})
        rows = _distances(
            tmp_path,
            {"REF": reference, "SWAP": swapped, "NEAR": near},
            metric="tn93",
            threshold=10.0,
        )
        assert frozenset(("REF", "SWAP")) not in rows
        assert frozenset(("REF", "NEAR")) in rows
        assert frozenset(("SWAP", "NEAR")) not in rows


@pytest.mark.reference_image
@pytest.mark.skipif(shutil.which("docker") is None, reason="docker not on PATH")
class TestAgreementWithTheReferenceImage:
    """The in-process correction against tn93 itself, pair for pair.

    Opt-in: the default run deselects `reference_image`, because this
    needs a Docker daemon and the built `hivdrwatch/tn93` image. Run it
    with `pytest tests/cluster -m reference_image`.

    `-a skip` is not optional here. tn93's default `-a resolve`
    resolves an IUPAC code against its partner where this engine drops
    the column, which is a third convention difference on top of the
    metric (§11.3) — it would make a correct implementation look wrong.
    """

    WIDTH = 900
    RECORDS = 24

    @staticmethod
    def _synthetic_alignment() -> dict[str, str]:
        """A seeded alignment with ragged coverage and a few IUPAC codes."""
        import numpy as np

        rng = np.random.default_rng(93_2026)
        bases = np.array(list("ACGT"))
        founder = rng.choice(
            bases, size=TestAgreementWithTheReferenceImage.WIDTH, p=[0.36, 0.17, 0.24, 0.23]
        )
        records: dict[str, str] = {}
        for index in range(TestAgreementWithTheReferenceImage.RECORDS):
            seq = founder.copy()
            sites = rng.random(seq.size) < (0.002 + 0.02 * rng.random())
            for position in np.flatnonzero(sites):
                base = str(seq[position])
                transition = PURINES.get(base) or PYRIMIDINES[base]
                seq[position] = transition if rng.random() < 0.7 else rng.choice(bases)
            seq[: rng.integers(0, 50)] = "-"
            tail = int(rng.integers(0, 50))
            if tail:
                seq[seq.size - tail :] = "-"
            for position in rng.choice(seq.size, size=int(rng.integers(0, 8)), replace=False):
                seq[position] = rng.choice(list("NRYW"))
            records[f"_synthetic_{index:03d}"] = "".join(seq)
        return records

    def test_pair_set_and_values_match(self, tmp_path: Path) -> None:
        records = self._synthetic_alignment()
        fasta = _write_fasta(tmp_path, records)
        ours = compute_pairwise_distances_python(
            fasta, tmp_path / "ours.csv", threshold=1.0, metric="tn93"
        )
        reference = distance.run_tn93(
            fasta,
            tmp_path / "reference.csv",
            threshold=1.0,
            image=DEFAULT_IMAGE,
            ambiguity_handling="skip",
        )

        mine = _rows(ours)
        theirs = _rows(reference.output_csv)
        assert mine.keys() == theirs.keys()
        assert mine
        # tn93 prints `%g`, six significant figures, and so does the
        # corrected path — so agreement is exact at the printed digit.
        for pair, value in theirs.items():
            assert mine[pair] == pytest.approx(value, rel=PRINTED_PRECISION, abs=1e-12), (
                sorted(pair)
            )
