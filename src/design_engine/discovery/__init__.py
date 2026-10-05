"""Discovery strategies for the unified design engine.

Each strategy is a reaction-level candidate generator. The pipeline runs the
selected strategies (optionally in parallel), pools their proposals, and feeds
the pool to the single shared verification/ranking core.
"""

from .base import DiscoveryParams, Baseline, STRATEGIES, get_strategy

__all__ = ["DiscoveryParams", "Baseline", "STRATEGIES", "get_strategy"]
