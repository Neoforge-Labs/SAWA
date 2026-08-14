"""Phylogenetic clustering: TN93 distances → connected components → dynamics.

Also includes the ground-truth validation harness for Paper 1 Figure 3:
`simulate.simulate_clusters` produces a synthetic dataset with known
membership; `validate.score_clusters` computes pair-based P/R/F1/ARI
of a predicted labeling against that truth.
"""

from .assembly import (
    ClusterAssignments,
    assign_clusters,
    build_cluster_graph,
    cluster_size_distribution,
)
from .distance import (
    DEFAULT_THRESHOLD,
    TN93Invocation,
    compute_pairwise_distances_python,
    read_distances,
    run_tn93,
)
from .dynamics import (
    cluster_geographic_spread,
    cluster_growth_rate,
    cluster_subtype_mix,
    cluster_time_series,
)
from .hivtrace import (
    HIVTraceInvocation,
    HIVTraceNotAvailable,
    is_available as hivtrace_is_available,
    run_hivtrace,
)
from .schema import CLUSTER_SCHEMA_VERSION, ClusterMembership, PairwiseDistance
from .simulate import (
    DEFAULT_MUTATION_RATE,
    DEFAULT_SEQ_LENGTH,
    SimulationResult,
    simulate_clusters,
)
from .validate import ClusterValidationResult, score_clusters

__all__ = [
    "CLUSTER_SCHEMA_VERSION",
    "ClusterAssignments",
    "ClusterMembership",
    "ClusterValidationResult",
    "DEFAULT_MUTATION_RATE",
    "DEFAULT_SEQ_LENGTH",
    "DEFAULT_THRESHOLD",
    "HIVTraceInvocation",
    "HIVTraceNotAvailable",
    "PairwiseDistance",
    "SimulationResult",
    "TN93Invocation",
    "assign_clusters",
    "build_cluster_graph",
    "compute_pairwise_distances_python",
    "cluster_geographic_spread",
    "cluster_growth_rate",
    "cluster_size_distribution",
    "cluster_subtype_mix",
    "cluster_time_series",
    "hivtrace_is_available",
    "read_distances",
    "run_hivtrace",
    "run_tn93",
    "score_clusters",
    "simulate_clusters",
]
