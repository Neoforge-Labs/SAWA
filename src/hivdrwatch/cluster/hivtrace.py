"""HIV-TRACE clustering wrapper — parity baseline for our clusterer.

Wraps `hivnetworkcsv` (from `hivclustering`) via the pinned Docker
image so Paper 1 Fig 3 can score our distance clusterer AND HIV-TRACE
against the same simulated ground truth on the same TN93 distances.

The parity claim: given identical TN93 distances and identical
threshold, both clusterers should produce identical connected
components. Any divergence is either (a) a bug in ours worth catching,
or (b) a HIV-TRACE-specific filter (contam filter, fringe filter) worth
reporting.

## Output shape

`hivnetworkcsv -f plain -J out.json` writes a JSON with `Nodes` (list
of {id: accession, cluster: int, ...}) and `Edges` (list of
{sequences: [i, j], length, ...}). We parse `Nodes` directly into
`{accession: cluster_id}` — singletons (no edges below threshold) get
`SINGLETON_<accession>` IDs matching our ground-truth convention in
`SimulationResult.truth_membership()`.
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_IMAGE = "hivdrwatch/hivtrace:0.1"
DEFAULT_THRESHOLD = 0.015          # HIV-TRACE convention
DEFAULT_FORMAT = "plain"           # tn93's CSV output shape


@dataclass(frozen=True)
class HIVTraceInvocation:
    """Provenance record for one HIV-TRACE run."""

    image: str
    args: tuple[str, ...]
    distances_csv: Path
    network_json: Path
    threshold: float


class HIVTraceNotAvailable(RuntimeError):
    """Raised when Docker or the HIV-TRACE image isn't available."""


def is_available(image: str = DEFAULT_IMAGE) -> bool:
    """Return True iff Docker is on PATH and the image exists locally."""
    if shutil.which("docker") is None:
        return False
    try:
        result = subprocess.run(
            ["docker", "image", "inspect", image],
            capture_output=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


def run_hivtrace(
    distances_csv: Path,
    network_json: Path,
    all_accessions: list[str],
    *,
    threshold: float = DEFAULT_THRESHOLD,
    image: str = DEFAULT_IMAGE,
    fmt: str = DEFAULT_FORMAT,
) -> tuple[dict[str, str], HIVTraceInvocation]:
    """Cluster with HIV-TRACE; return ({accession: cluster_id}, invocation).

    `all_accessions` is required because `hivnetworkcsv` only lists
    accessions that appear in at least one below-threshold edge —
    everything else (singletons) has to be added back explicitly with
    a `SINGLETON_<accession>` cluster ID, matching our
    `SimulationResult.truth_membership()` convention.
    """
    distances_csv = distances_csv.resolve()
    network_json = network_json.resolve()

    if not distances_csv.exists():
        raise FileNotFoundError(f"distances CSV missing: {distances_csv}")
    if distances_csv.parent != network_json.parent:
        raise ValueError(
            "distances CSV and network JSON must share a parent for the docker mount"
        )
    if not is_available(image):
        raise HIVTraceNotAvailable(
            f"HIV-TRACE image '{image}' not available. Run `make images` first."
        )

    work_dir = distances_csv.parent
    # `-J` in hivnetworkcsv is a BOOLEAN flag ("emit JSON") — the actual
    # path can't be passed positionally. Some versions accept `-O <path>`;
    # others only write to stdout. Capture stdout and write it ourselves
    # for maximum portability across hivclustering versions.
    #
    # `--no-degree-fit` skips the power-law degree-distribution fitting
    # step that hivclustering's `describe_network` calls unconditionally.
    # That step needs the `hppy` (HyPhy Python bindings) optional dep,
    # which requires building HyPhy. We don't need the degree fit for
    # cluster-membership parity — only for HIV-TRACE's own network
    # visualization — so skipping it removes a heavy dep with no cost.
    hn_args = (
        "-i", f"/work/{distances_csv.name}",
        "-t", str(threshold),
        "-f", fmt,
        "--no-degree-fit",
        "-J",
    )
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{work_dir}:/work",
        image,
        *hn_args,
    ]
    logger.info("running HIV-TRACE: %s", " ".join(docker_cmd))

    result = subprocess.run(docker_cmd, capture_output=True, check=False)
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")[-2000:]
        raise RuntimeError(f"hivnetworkcsv failed (exit {result.returncode}):\n{stderr}")

    network_json.write_bytes(result.stdout)
    if network_json.stat().st_size == 0:
        raise RuntimeError(
            "hivnetworkcsv exited 0 but produced no JSON on stdout — "
            "check tn93 distances CSV shape."
        )

    membership = parse_network_json(network_json, all_accessions)
    invocation = HIVTraceInvocation(
        image=image,
        args=hn_args,
        distances_csv=distances_csv,
        network_json=network_json,
        threshold=threshold,
    )
    logger.info(
        "HIV-TRACE parsed %d accessions into %d clusters",
        len(membership),
        len({v for v in membership.values() if not v.startswith("SINGLETON_")}),
    )
    return membership, invocation


def parse_network_json(
    network_json: Path,
    all_accessions: list[str],
) -> dict[str, str]:
    """Parse hivnetworkcsv's JSON output into `{accession: cluster_id}`.

    Handles two shapes hivclustering ships depending on `-J` vs `-j`
    and version:

      1. **Compact `-J`** (hivclustering 1.9.x): parallel arrays with
         indirected cluster labels::

             {"Nodes": {
                 "id":      ["S000102", "S000103", ...],
                 "cluster": {
                     "keys":   {"0": 1, "1": 2, ...},   # compact idx → real label
                     "values": [0, 0, 0, 1, 1, ...],    # per-node compact idx
                 }
             }}

      2. **Legacy list-of-dicts** (older builds and `-j` pretty mode)::

             {"Nodes": [{"id": "S000102", "cluster": 1}, ...]}

    Accessions not present in Nodes are singletons; every input
    accession lands in the returned mapping (singletons get their own
    `SINGLETON_<accession>` id), matching the ground-truth convention
    in `SimulationResult.truth_membership()`.
    """
    payload = json.loads(network_json.read_text())
    nodes = payload.get("Nodes", [])

    membership: dict[str, str] = {}

    if isinstance(nodes, dict) and "id" in nodes and "cluster" in nodes:
        # Compact shape.
        ids = nodes["id"]
        cluster_field = nodes["cluster"]
        if isinstance(cluster_field, dict) and "values" in cluster_field:
            keys_map = cluster_field.get("keys") or {}
            values = cluster_field["values"]
            for acc, compact_idx in zip(ids, values, strict=False):
                # `keys_map` is {str(compact_idx): real_cluster_id}. Fall
                # back to compact_idx itself if the mapping is absent.
                real_label = keys_map.get(str(compact_idx), compact_idx)
                membership[str(acc)] = f"HT-{int(real_label):04d}"
        elif isinstance(cluster_field, list):
            for acc, label in zip(ids, cluster_field, strict=False):
                membership[str(acc)] = f"HT-{int(label):04d}"
    elif isinstance(nodes, list):
        # Legacy list-of-dicts shape.
        for node in nodes:
            if not isinstance(node, dict):
                continue
            acc = node.get("id")
            cluster = node.get("cluster")
            if acc is None or cluster is None:
                continue
            membership[str(acc)] = f"HT-{int(cluster):04d}"

    # Add singletons — every input accession must land somewhere.
    for acc in all_accessions:
        if acc not in membership:
            membership[acc] = f"SINGLETON_{acc}"

    return membership


__all__ = [
    "DEFAULT_IMAGE",
    "DEFAULT_THRESHOLD",
    "HIVTraceInvocation",
    "HIVTraceNotAvailable",
    "is_available",
    "parse_network_json",
    "run_hivtrace",
]
