"""Live streaming of subprocess stdout/stderr to the module logger.

Bio tools like IQ-TREE 2 and TreeTime run for minutes-to-hours; buffering
their output until the process exits gives operators zero visibility
into progress. This helper runs a command with `Popen`, spawns two
drain threads (one per pipe), and forwards each line to the logger as
soon as the tool flushes it.

Design choices:
  * **Per-line, unbuffered**: `bufsize=1` + `text=True` gets line
    buffering; combined with pipe-side flushes from the child, lines
    appear near-real-time. IQ-TREE and TreeTime both flush regularly.
  * **Return the captured text** so callers can grep the tail for
    error hints, same as `subprocess.run(capture_output=True)`.
  * **Cap what we hold in memory** at `MAX_CAPTURED_CHARS` per stream —
    a multi-hour IQ-TREE run can spit tens of MB into stderr; keep
    only the last chunk so the return value stays bounded.
"""

from __future__ import annotations

import logging
import subprocess
import threading
from collections import deque

# Keep at most this many characters of captured stdout/stderr around
# for the return value. Older lines are discarded but were already
# streamed to the logger, so nothing is truly lost.
MAX_CAPTURED_CHARS = 200_000


def stream_subprocess(
    cmd: list[str],
    *,
    prefix: str,
    logger: logging.Logger,
    stderr_prefix: str | None = None,
    timeout_s: float | None = None,
) -> tuple[int, str, str]:
    """Run `cmd`, stream stdout+stderr to `logger` line-by-line.

    Returns `(returncode, stdout_tail, stderr_tail)`. Tails are capped
    at `MAX_CAPTURED_CHARS` per stream — full output already went to
    the logger.

    `prefix` labels stdout lines in the log (e.g. `iqtree/A1`);
    `stderr_prefix` defaults to `prefix + "!"` so stderr is visually
    distinguishable in the Prefect UI.
    """
    if stderr_prefix is None:
        stderr_prefix = f"{prefix}!"

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
        text=True,
    )
    if proc.stdout is None or proc.stderr is None:
        raise RuntimeError("Popen returned pipes as None")

    stdout_buf: deque[str] = deque()
    stderr_buf: deque[str] = deque()
    stdout_chars = [0]  # mutable counter for the closure
    stderr_chars = [0]

    def _drain(stream, buf: deque[str], counter: list[int], line_prefix: str) -> None:
        for line in stream:
            line = line.rstrip("\n")
            if not line:
                continue
            logger.info("[%s] %s", line_prefix, line)
            buf.append(line)
            counter[0] += len(line) + 1
            # Trim from the front if the buffer got too fat.
            while counter[0] > MAX_CAPTURED_CHARS and buf:
                dropped = buf.popleft()
                counter[0] -= len(dropped) + 1
        stream.close()

    t_out = threading.Thread(
        target=_drain,
        args=(proc.stdout, stdout_buf, stdout_chars, prefix),
        daemon=True,
    )
    t_err = threading.Thread(
        target=_drain,
        args=(proc.stderr, stderr_buf, stderr_chars, stderr_prefix),
        daemon=True,
    )
    t_out.start()
    t_err.start()

    try:
        returncode = proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        proc.kill()
        t_out.join(timeout=2)
        t_err.join(timeout=2)
        raise

    t_out.join()
    t_err.join()

    return returncode, "\n".join(stdout_buf), "\n".join(stderr_buf)
