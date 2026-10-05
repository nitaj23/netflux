"""
schema.py - Result dataclasses shared by every design-engine stage.

These types were previously duplicated verbatim in ``optforce.pipeline`` and
``optknock.pipeline``. They are extracted here once so discovery strategies,
verification, ranking, and reporting all speak one vocabulary.

A :class:`Proposal` is what a discovery strategy emits (a reaction-level design
plus its discovery-model diagnostics). Everything else is produced downstream by
the shared verification/ranking core and consumed by :mod:`design_engine.report`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .gpr import GeneAnnotation


# All reaction-level actions the engine understands. OptKnock emits only
# 'knockout'; OptForce may also emit 'upregulation'/'downregulation'.
SUPPORTED_ACTIONS = ("upregulation", "downregulation", "knockout")


@dataclass(frozen=True)
class ReactionAction:
    """One normalized reaction-level intervention (the retained bound op)."""
    reaction_id: str
    action: str          # 'upregulation' | 'downregulation' | 'knockout'
    forced_flux: float

    def identifier(self) -> str:
        return f"{self.reaction_id}:{self.action}"


@dataclass
class Proposal:
    """A reaction-level design handed from a discovery strategy to verification.

    ``source`` names the strategy and cardinality (e.g. ``optforce_k1``,
    ``optknock_k2``). When two strategies independently propose the same action
    set, the pool merges them into a single Proposal whose ``source`` lists both
    (e.g. ``optforce+optknock_k1``).
    """
    actions: tuple[ReactionAction, ...]
    source: str
    discovery_floor: Optional[float]
    discovery_growth: Optional[float]


@dataclass(frozen=True)
class GeneOption:
    """A gene-deletion alternative (KO) or putative regulatory target set."""
    genes: tuple[str, ...]
    kind: str            # 'knockout' | 'putative_regulation'
    complete: bool = True


@dataclass(frozen=True)
class GeneConstructResult:
    """Full-model verification of one gene-deletion construct."""
    genes: tuple[str, ...]
    silenced_reactions: tuple[str, ...]
    verified_growth: float
    verified_floor: float
    verified_ceiling: float
    retained_mu: float
    import_loop_status: str
    buildable: bool


@dataclass(frozen=True)
class InterventionCandidate:
    """A verified, audited, (optionally ranked) intervention set."""
    source: str
    reaction_actions: tuple[ReactionAction, ...]
    cardinality: int

    # Discovery-model metrics (never used for ranking).
    discovery_floor: Optional[float]
    discovery_growth: Optional[float]

    # Full-model reaction-bound verification.
    verified_floor: float
    verified_ceiling: float
    verified_growth: float
    floor_delta: float
    ceiling_delta: float
    coupling_potential: Optional[float]
    retained_mu: float
    mandatory_production: Optional[bool]

    guaranteed_growth: float
    pinned_growth: Optional[float]
    growth_robustness_ratio: Optional[float]
    coupling_class: str

    gene_options: tuple[GeneOption, ...]
    translation_status: str
    gene_constructs: tuple[GeneConstructResult, ...]

    import_loop_status: str
    target_import_reactions: tuple[str, ...]
    closed_import_floor: Optional[float]

    annotations: tuple[GeneAnnotation, ...]
    artifact_flags: tuple[str, ...]
    passes_checks: bool
    pareto_front: Optional[int]
    rank: Optional[int]
    rejection_reason: Optional[str]

    # --- unified-engine additions ---
    # True iff every action is a knockout (so the whole design can be realized
    # as a gene-deletion construct rather than a regulatory hypothesis).
    all_knockout: bool = True
    # True iff a gene-deletion construct exists that both grows and holds its
    # production floor on the full model (i.e. at least one buildable construct).
    buildable: bool = False


@dataclass
class RankingResult:
    """Output of rank_intervention_candidates()."""
    accepted: list[InterventionCandidate]           # ranked, buildable shortlist
    rejected: list[InterventionCandidate]           # diagnostics, deterministic order
    hypotheses: list[InterventionCandidate]         # verified regulatory designs, unranked
    mu_wt: float
    floor_wt: float
    ceiling_wt: float
    substrate_rxn: Optional[str]
    methods: tuple[str, ...]
    extraction_fallback: bool
    extraction_failure_reason: Optional[str]
    pareto_objectives: tuple[str, ...]
    timings: dict[str, float]
    notes: list[str] = field(default_factory=list)
