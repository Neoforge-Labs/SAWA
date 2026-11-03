"""Repo-wide pytest plumbing: the `--update-golden` switch for golden-file tests.

The option is registered here, not beside the tests that use it, because pytest
only loads a conftest before parsing options when the run targets its directory,
so `pytest tests/ --update-golden` would otherwise be rejected.
"""

from __future__ import annotations

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Rewrite golden files from the current output instead of comparing against them.",
    )


@pytest.fixture
def golden_update(request: pytest.FixtureRequest) -> bool:
    return bool(request.config.getoption("--update-golden"))
