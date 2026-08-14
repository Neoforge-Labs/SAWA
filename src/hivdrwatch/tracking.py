"""Thin MLflow adapter.

The whole codebase talks to MLflow through this one module. That gives us
three things at low cost:

1. **Env-driven URI.** `MLFLOW_TRACKING_URI` (usually `http://localhost:5000`
   in dev, set by the compose stack via `.env`) is the only knob. Nothing
   hardcoded.
2. **No-op when unset.** If `MLFLOW_TRACKING_URI` is empty, every call
   here becomes a no-op — tests, CI, and offline runs don't need a live
   MLflow server. This is the biggest reason tracking gets shelved
   elsewhere; we're avoiding that trap.
3. **One place to add auth, tags, or a second backend later.**

Usage inside a Prefect flow (or anywhere):

    from hivdrwatch.tracking import tracked_run

    with tracked_run("train-dr", tags={"drug": "3TC"}) as run:
        run.log_params({"C": 0.5, "l1_ratio": 0.5})
        model, metrics = ...
        run.log_metrics(metrics)
        run.log_artifact(model_path)
        run.log_dict(feature_importance, "feature_importance.json")

If MLflow is unreachable the context manager still enters and every call
is silently swallowed. Nothing crashes.
"""

from __future__ import annotations

import contextlib
import logging
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def _mlflow_uri() -> str | None:
    """Env override wins; no default so tests are hermetic."""
    uri = os.environ.get("MLFLOW_TRACKING_URI")
    return uri.strip() if uri and uri.strip() else None


def is_enabled() -> bool:
    """True iff MLFLOW_TRACKING_URI is set. Cheap probe for callers that
    want to skip expensive artifact prep when tracking is off."""
    return _mlflow_uri() is not None


# ------------------------------------------------------------------
# Public API — the two shapes callers use
# ------------------------------------------------------------------


class TrackedRun:
    """Thin wrapper over an MLflow active run.

    All methods swallow exceptions and log at WARNING. Missing tracking
    server never breaks a training run — dropped metrics are visible in
    logs but don't propagate.
    """

    def __init__(self, mlflow_module: Any | None) -> None:
        self._mlflow = mlflow_module

    def log_params(self, params: dict[str, Any]) -> None:
        if self._mlflow is None:
            return
        try:
            self._mlflow.log_params(_stringify(params))
        except Exception as e:
            logger.warning("mlflow log_params failed: %s", e)

    def log_metrics(self, metrics: dict[str, float], *, step: int | None = None) -> None:
        if self._mlflow is None:
            return
        try:
            flat = {k: float(v) for k, v in metrics.items() if _is_numeric(v)}
            if step is not None:
                self._mlflow.log_metrics(flat, step=step)
            else:
                self._mlflow.log_metrics(flat)
        except Exception as e:
            logger.warning("mlflow log_metrics failed: %s", e)

    def log_artifact(self, path: Path, *, artifact_path: str | None = None) -> None:
        if self._mlflow is None:
            return
        try:
            self._mlflow.log_artifact(str(path), artifact_path=artifact_path)
        except Exception as e:
            logger.warning("mlflow log_artifact failed for %s: %s", path, e)

    def log_dict(self, payload: dict[str, Any], filename: str) -> None:
        if self._mlflow is None:
            return
        try:
            self._mlflow.log_dict(payload, filename)
        except Exception as e:
            logger.warning("mlflow log_dict failed for %s: %s", filename, e)

    def set_tags(self, tags: dict[str, str]) -> None:
        if self._mlflow is None:
            return
        try:
            self._mlflow.set_tags({k: str(v) for k, v in tags.items()})
        except Exception as e:
            logger.warning("mlflow set_tags failed: %s", e)


@contextlib.contextmanager
def tracked_run(
    experiment: str,
    *,
    run_name: str | None = None,
    tags: dict[str, str] | None = None,
) -> Iterator[TrackedRun]:
    """Context manager that starts (and finishes) an MLflow run.

    When `MLFLOW_TRACKING_URI` isn't set, yields a no-op `TrackedRun`
    so caller code stays uniform.
    """
    uri = _mlflow_uri()
    if uri is None:
        logger.info("MLflow tracking disabled (MLFLOW_TRACKING_URI unset); no-op")
        yield TrackedRun(None)
        return

    try:
        import mlflow  # type: ignore[import-not-found]
    except ImportError:
        logger.warning("mlflow not installed; tracking disabled for this run")
        yield TrackedRun(None)
        return

    try:
        mlflow.set_tracking_uri(uri)
        mlflow.set_experiment(experiment)
        with mlflow.start_run(run_name=run_name) as _run:
            wrapper = TrackedRun(mlflow)
            if tags:
                wrapper.set_tags(tags)
            wrapper.set_tags(_default_tags())
            yield wrapper
    except Exception as e:
        logger.warning(
            "mlflow run setup failed (uri=%s experiment=%s): %s — running untracked",
            uri,
            experiment,
            e,
        )
        yield TrackedRun(None)


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------


def _stringify(d: dict[str, Any]) -> dict[str, str]:
    """MLflow params must be strings ≤ 6000 chars."""
    return {k: str(v)[:6000] for k, v in d.items()}


def _is_numeric(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _default_tags() -> dict[str, str]:
    """Tags applied to every run — git SHA, host, user."""
    tags: dict[str, str] = {}
    for env_key, tag_key in (
        ("USER", "user"),
        ("HOSTNAME", "host"),
        ("GITHUB_SHA", "git_sha"),
    ):
        v = os.environ.get(env_key)
        if v:
            tags[tag_key] = v
    # Fallback: local git SHA if we're inside a repo.
    if "git_sha" not in tags:
        sha = _git_sha()
        if sha:
            tags["git_sha"] = sha
    return tags


def _git_sha() -> str | None:
    """Read the current git SHA if we're in a checkout. No git dep required."""
    try:
        # Traverse up from cwd looking for .git/HEAD
        cur = Path.cwd().resolve()
        for _ in range(6):
            head = cur / ".git" / "HEAD"
            if head.exists():
                ref = head.read_text(encoding="utf-8").strip()
                if ref.startswith("ref: "):
                    ref_path = cur / ".git" / ref[5:]
                    if ref_path.exists():
                        return ref_path.read_text(encoding="utf-8").strip()[:12]
                    return None  # ref file missing (packed refs, etc.)
                return ref[:12]  # detached HEAD — raw SHA
            if cur.parent == cur:
                break
            cur = cur.parent
    except OSError:
        return None
    return None
