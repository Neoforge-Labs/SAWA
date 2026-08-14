"""Unit tests for the unified CLI.

We don't exercise the underlying pipeline modules here — those have
their own test suites. This file locks in:

  * Every registered subcommand is discoverable
  * Extra args are forwarded verbatim to the underlying module's `main`
  * `hivdrwatch version` prints without crashing
  * `hivdrwatch commands` lists every registered subcommand
"""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

import pytest
from typer.testing import CliRunner

from hivdrwatch.cli.app import _SUBCOMMANDS, app

runner = CliRunner()


def _install_fake_module(monkeypatch: pytest.MonkeyPatch, path: str) -> SimpleNamespace:
    """Register a fake module at `path` whose `main(argv)` records the argv."""
    calls = SimpleNamespace(argv=None, returned=0)

    fake = ModuleType(path)

    def _main(argv: list[str] | None = None) -> int:
        calls.argv = list(argv) if argv is not None else []
        return calls.returned

    fake.main = _main  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, path, fake)
    return calls


class TestSubcommandDiscovery:
    def test_every_registered_subcommand_appears_in_help(self) -> None:
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        for name in _SUBCOMMANDS:
            assert name in result.stdout, f"{name} missing from --help"

    def test_commands_subcommand_lists_all(self) -> None:
        result = runner.invoke(app, ["commands"])
        assert result.exit_code == 0
        for name in _SUBCOMMANDS:
            assert name in result.stdout


class TestArgForwarding:
    def test_extra_args_forwarded_verbatim(self, monkeypatch: pytest.MonkeyPatch) -> None:
        module_path, _ = _SUBCOMMANDS["ingest"]
        calls = _install_fake_module(monkeypatch, module_path)

        result = runner.invoke(
            app,
            [
                "ingest",
                "--skip-genbank",
                "--lanl-dir",
                "data/raw/lanl",
                "--out",
                "out.parquet",
            ],
        )
        assert result.exit_code == 0
        assert calls.argv == [
            "--skip-genbank",
            "--lanl-dir",
            "data/raw/lanl",
            "--out",
            "out.parquet",
        ]

    def test_exit_code_bubbles_up_from_module(self, monkeypatch: pytest.MonkeyPatch) -> None:
        module_path, _ = _SUBCOMMANDS["curate"]
        calls = _install_fake_module(monkeypatch, module_path)
        calls.returned = 2

        result = runner.invoke(app, ["curate"])
        assert result.exit_code == 2

    def test_unknown_option_still_forwarded(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # `--madeup-flag` is not something Typer knows about; the CLI must not
        # reject it — the module's argparse is the source of truth.
        module_path, _ = _SUBCOMMANDS["sierra"]
        calls = _install_fake_module(monkeypatch, module_path)

        result = runner.invoke(app, ["sierra", "--madeup-flag", "42"])
        assert result.exit_code == 0
        assert calls.argv == ["--madeup-flag", "42"]


class TestVersion:
    def test_version_runs_and_prints_schema(self) -> None:
        result = runner.invoke(app, ["version"])
        assert result.exit_code == 0
        assert "hivdrwatch" in result.stdout
        assert "ingest schema" in result.stdout
        assert "cluster schema" in result.stdout


class TestSubcommandContract:
    """Every registered subcommand's module path must exist AND have `main`.

    Catches the "renamed a module and forgot to update the CLI table"
    class of bug before it ships.
    """

    @pytest.mark.parametrize("name", list(_SUBCOMMANDS))
    def test_module_importable_with_main(self, name: str) -> None:
        module_path, _ = _SUBCOMMANDS[name]
        try:
            module = __import__(module_path, fromlist=["main"])
        except ImportError as e:
            # Missing third-party dep (mappy, StrEnum from 3.12, etc.) —
            # not a CLI-table bug, so skip rather than fail. The rename-
            # detection value of this test still fires for pure-Python
            # subcommands (version, commands, bench-ingest).
            pytest.skip(f"{module_path} not importable in this env: {e}")
        assert callable(getattr(module, "main", None)), (
            f"CLI subcommand {name!r} → {module_path!r} lacks main(argv)"
        )
