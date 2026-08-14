"""Entry point for `python -m hivdrwatch.cli`.

Also invoked by the `hivdrwatch` script defined in pyproject.toml
under [project.scripts]; both hit `app()`.
"""

from __future__ import annotations

from .app import app


def main() -> None:  # pragma: no cover - trivial dispatch
    app()


if __name__ == "__main__":
    main()
