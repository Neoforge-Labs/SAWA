"""Pytest plumbing for the dashboard tests: the `--update-golden` switch.

Pytest only loads this conftest before parsing options when the run targets
this directory, so regenerate with `pytest tests/dashboard --update-golden`.
"""

from __future__ import annotations

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Rewrite the dashboard golden files from the current output instead of comparing.",
    )


@pytest.fixture
def golden_update(request: pytest.FixtureRequest) -> bool:
    return bool(request.config.getoption("--update-golden"))
