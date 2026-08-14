"""Utilities for pinning the HIVdb algorithm version across a run.

Every DRResult carries `algorithm_version`. When training a model, you want
to filter to a single version — mixing versions in the label column means
mixing labeling functions, which invalidates the training loss.

Usage:
    algo = pin_algorithm(client)   # queries Sierra once, caches result
    ...
    if algo.version != old_version:
        log.warning("HIVdb algorithm bumped %s → %s", old_version, algo.version)
"""

from __future__ import annotations

import logging

from .sierra import AlgorithmVersion, SierraClient

logger = logging.getLogger(__name__)


def pin_algorithm(client: SierraClient) -> AlgorithmVersion:
    """Fetch (or read from cache) the current HIVdb algorithm version.

    Idempotent: safe to call at the start of every flow run.
    """
    algo = client.current_algorithm()
    logger.info(
        "HIVdb algorithm pinned: %s (published %s)",
        algo.version,
        algo.publish_date,
    )
    return algo
