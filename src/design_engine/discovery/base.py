"""
base.py - The discovery-strategy contract and shared knobs.

A discovery strategy turns a (discovery) model + target/biomass into a list of
reaction-level :class:`~design_engine.schema.Proposal` objects. It never
verifies, audits, ranks, or enumerates genes -- that is the shared core's job.
Strategies must be side-effect free on the model they are given (restore any
bounds they change) so the same model can be reused across strategies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Protocol

import cobra

from ..schema import Proposal


@dataclass(frozen=True)
class Baseline:
    """Full-model wild-type baseline, computed once and shared by all strategies."""
    mu_wt: float
    floor_wt: float
    ceiling_wt: float


@dataclass(frozen=True)
class DiscoveryParams:
    """Knobs forwarded to every discovery strategy (a strategy uses what it needs)."""
    min_growth_frac: float = 0.05
    min_growth_abs: Optional[float] = None
    eps_prod: float = 1e-6
    class_tol: float = 1e-6
    gamma: float = 0.001
    # OptForce
    enable_k2: bool = False
    # OptKnock
    max_knockouts: int = 1
    # Solver routing for a strategy's own LP handles.
    backend: str = "cobra"
    solver_timeout: Optional[float] = None


class DiscoveryStrategy(Protocol):
    """A named reaction-level candidate generator."""
    name: str

    def discover(self, disc_model: cobra.Model, target: str, biomass: str,
                 substrate: Optional[str], baseline: Baseline,
                 params: DiscoveryParams, notes: list[str]) -> list[Proposal]:
        ...


# Registry, populated by the strategy modules on import.
STRATEGIES: dict[str, "DiscoveryStrategy"] = {}


def register(strategy: "DiscoveryStrategy") -> "DiscoveryStrategy":
    STRATEGIES[strategy.name] = strategy
    return strategy


def get_strategy(name: str) -> "DiscoveryStrategy":
    # Import lazily so a missing optional dependency in one strategy never
    # blocks the other from loading.
    from . import optforce as _optforce  # noqa: F401
    from . import optknock as _optknock  # noqa: F401
    if name not in STRATEGIES:
        raise KeyError(name)
    return STRATEGIES[name]
