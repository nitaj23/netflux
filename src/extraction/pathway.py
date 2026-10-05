# src/extraction/pathway.py
"""
Flux-essential bioproduction pathway extraction (FASTCC -> MFG -> GUDIE).

One call that ties the flux-consistency boundary (extraction.fastcc) to GUDIE
subgraph extraction (extraction.gudie). Given a GEM and a bioproduction target
(a product exchange/demand reaction, or the metabolite to drain), it returns the
localized subgraph H: the flux-essential production pathway for that target.

Pipeline
--------
    fastcc            -> flux-consistent universe (drop blocked reactions)
    seed + maximise   -> productive FBA flux state for the target
    build_mfg         -> mass-flow graph, restricted to its active subgraph G_F
    run_gudie         -> localized H around the target (VUDIE PPR + LUDIE FCA)

Hyper-parameters are the validated all-organism-sweep defaults (see
scripts/validate_gudie_all_organisms, test_fastcc_gudie): FASTCC eps=1e-4 and
GUDIE alpha_r=0.3, k=0.1, hops=8, lam=0.05 -- already GudieConfig's defaults, so
`extract_bioproduction_pathway(model, target)` needs no tuning to reproduce them.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import Sequence

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba
import networkx as nx

from src.extraction.fastcc import fastcc_from_cobra
from src.extraction.gudie import GudieConfig, GudieResult, run_gudie
from src.extraction.repair import RepairResult, repair_connectivity
from src.graph_ops import (
    build_mfg, fpe_similarity_to_target, graph_json, graphml_bytes)

logger = logging.getLogger(__name__)

FASTCC_EPS = 1e-4  # flux-consistency / producibility threshold (validated default)

# ── FVA design-annotation constants ──────────────────────────────────────────
COUPLING_TAU = 0.1          # c_i > tau => reaction is coupled to the target
MAX_OPERATING_POINTS = 8    # cap on |operating points| per annotation
ANNOTATION_LP_BUDGET = 4000  # ceiling on 2*|H|*points*loopless_factor LPs
LOOPLESS_FACTOR = 4         # loopless FVA multiplies the LP cost by ~this

__all__ = [
    "PathwayResult", "extract_bioproduction_pathway", "FASTCC_EPS",
    "ReactionAnnotation", "PathwayAnnotation",
    "TargetNotFound", "TargetBlocked", "TargetNotProducible", "TargetNoMassFlow",
    "InfeasibleModel", "InvalidOperatingPoint", "AnnotationBudgetExceeded",
    "COUPLING_TAU", "MAX_OPERATING_POINTS", "ANNOTATION_LP_BUDGET",
]


# ── Typed errors (subclass ValueError/RuntimeError so `except ValueError`
#    in existing callers/scripts keeps working) ────────────────────────────────
class TargetNotFound(ValueError):
    """`target` is neither a reaction nor a metabolite id in the model."""

class TargetBlocked(ValueError):
    """`target` reaction survives id-check but is dropped by FASTCC (blocked)."""

class InfeasibleModel(ValueError):
    """The consistent model is infeasible when maximising the target."""

class TargetNotProducible(ValueError):
    """Max target flux is ~0: the target cannot be produced under this medium."""

class TargetNoMassFlow(ValueError):
    """Target carries no mass flow in the productive state; nothing to extract."""

class InvalidOperatingPoint(ValueError):
    """An FVA operating-point fraction is outside (0, 1] or too many were given."""

class AnnotationBudgetExceeded(RuntimeError):
    """The FVA annotation would exceed the LP budget ceiling."""


@dataclass(slots=True)
class ReactionAnnotation:
    reaction_id: str
    gene_id: str                 # first gene / GPR string ('' if none)
    coupling: float | None       # FCA c_i in [0,1] (None if FCA unavailable)
    per_op: list[dict]           # one dict per operating point (see _annotate_fva)


@dataclass(slots=True)
class PathwayAnnotation:
    operating_points: list[float]
    v_star_target: float
    coupling_tau: float
    per_reaction: dict[str, ReactionAnnotation]


@dataclass(slots=True)
class PathwayResult:
    subgraph: nx.DiGraph        # H, the flux-essential production pathway (MFG edges)
    reactions: list[str]        # reaction ids retained in H
    target_id: str              # seed reaction (existing exchange/demand, or the added demand)
    n_consistent: int           # |flux-consistent universe| after FASTCC
    n_active: int               # |active MFG| (mass-carrying nodes feeding the flux state)
    gudie: GudieResult          # full GUDIE diagnostics (ppr, coupling, threshold, fallbacks)
    annotation: PathwayAnnotation | None = None  # FVA design-annotation (if requested)
    repair: RepairResult | None = None           # connectivity-repair record (None if repair=False)
    consistent_ids: frozenset[str] = frozenset() # F, the FASTCC flux-consistent universe
    #   under this model's medium. Carried (additively) so the strain-design engine can
    #   reuse F when chained after extraction, instead of re-running FASTCC (spec v3 §6,
    #   mode-A entry). Empty only for PathwayResults built before this field existed.
    fastcc_medium_key: str = ""  # SHA256 medium fingerprint of the model F was computed on
    #   (design.universe.medium_key). Lets the design engine PROVE the reused F matches the
    #   medium it will design under, instead of trusting it -- flux-consistency is
    #   medium-dependent, so a mismatched F is a wrong universe (spec v3 §3 Stage 1). Empty
    #   for PathwayResults built before this field existed => reuse falls back to recompute.

    def to_graphml(self) -> bytes:
        """The derived pathway H as in-memory GraphML bytes (for a downloadable
        .graphml file). Delegates to graph_ops.graphml_bytes."""
        return graphml_bytes(self.subgraph)

    def to_json(self, *, as_str: bool = False) -> dict | str:
        """The derived pathway H as node-link JSON — a dict (default) or string
        (as_str=True) for a downloadable .json file. graph_ops.graph_json."""
        return graph_json(self.subgraph, as_str=as_str)


def _seed_target(model: cobra.Model, target: str) -> str:
    """Resolve `target` to a seed reaction id. A reaction id is used as-is; a
    metabolite id gets a demand reaction (M ->), reusing an existing boundary if
    one is already present (add_boundary raises when it is)."""
    if target in model.reactions:
        return target
    if target not in model.metabolites:
        raise TargetNotFound(f"target {target!r} is neither a reaction nor a metabolite id")
    met = model.metabolites.get_by_id(target)
    try:
        return model.add_boundary(met, type="demand").id
    except ValueError:
        dm = next((r for r in model.reactions
                   if len(r.metabolites) == 1 and met in r.metabolites), None)
        if dm is None:
            raise
        return dm.id


def extract_bioproduction_pathway(
    model: cobra.Model,
    target: str,
    *,
    source_ids: Sequence[str] | None = None,
    fastcc_eps: float = FASTCC_EPS,
    config: GudieConfig | None = None,
    repair: bool = False,
    fva_operating_points: Sequence[float] | None = None,
    fva_loopless: bool = False,
) -> PathwayResult:
    """
    Extract the flux-essential bioproduction pathway for `target` from `model`.

    Parameters
    ----------
    model       cobra GEM; its bounds define the medium/condition. Not mutated
                (work is done on a copy).
    target      product exchange/demand reaction id, or a metabolite id for which
                a demand reaction (M ->) is added and maximised.
    source_ids  optional substrate/uptake reaction ids. When given, GUDIE's
                dual-anchor selects the source->target corridor and force-keeps
                the uptakes; None => single-anchor (target-seeded) extraction.
    fastcc_eps  flux-consistency threshold (default 1e-4, the validated value).
    config      GUDIE hyper-parameters; defaults are the validated sweep values.
    repair      opt-in connectivity repair (default False => byte-identical output).
                When True and source_ids are given, reconnects the source->target
                corridor if GUDIE's H dropped it, before FVA annotation so the added
                reactions are design-tagged in the same pass (see repair.py).

    Returns
    -------
    PathwayResult with H (induced MFG subgraph, original edge directions) and the
    FASTCC/MFG/GUDIE diagnostics.

    Raises
    ------
    ValueError if the target survives FASTCC but is infeasible, not producible
    (max flux ~ 0), or carries no mass flow in the productive state.
    """
    cfg = config or GudieConfig()

    # 1. FASTCC: restrict a working copy to the flux-consistent universe.
    _, consistent_ids = fastcc_from_cobra(model, fastcc_eps)
    consistent = set(consistent_ids)
    # The medium fingerprint was consumed only by the strain-design engine (F-reuse
    # validation); that engine has been removed, so the field is left empty until it
    # is rebuilt (see PathwayResult.fastcc_medium_key).
    fastcc_medium_key = ""
    reduced = model.copy()
    reduced.remove_reactions(
        [r for r in reduced.reactions if r.id not in consistent], remove_orphans=True)

    # A pre-existing target reaction that FASTCC dropped is blocked (never
    # feasible), distinct from an infeasible medium below.
    if target in model.reactions and target not in consistent:
        raise TargetBlocked(f"target {target!r} is flux-blocked (removed by FASTCC)")

    # 2. Seed the target and maximise its flux -> productive FBA state.
    target_id = _seed_target(reduced, target)
    reduced.objective = target_id
    sol = reduced.optimize()
    if sol.status != "optimal":
        raise InfeasibleModel(f"target {target_id!r} infeasible in consistent model (status={sol.status})")
    if abs(sol.objective_value) < fastcc_eps:
        raise TargetNotProducible(f"target {target_id!r} not producible (max flux {sol.objective_value:.2e} ~ 0)")

    # 3. Mass-flow graph, restricted to the active (mass-carrying) subgraph.
    G = build_mfg(reduced, sol)
    G_F = G.subgraph([n for n in G.nodes if G.degree(n) > 0]).copy()
    if target_id not in G_F:
        raise TargetNoMassFlow(f"target {target_id!r} carries no mass flow; nothing to extract")

    # 4. GUDIE localization (Phase B FPE blend consulted only when lam > 0).
    embedding_sim = (fpe_similarity_to_target(G_F, target_id, list(G_F.nodes()))
                     if cfg.lam > 0.0 else None)
    res = run_gudie(G_F, target_id, source_ids=source_ids, model=reduced,
                    config=cfg, embedding_sim=embedding_sim)

    logger.info("pathway: |consistent|=%d -> |G_F|=%d -> |H|=%d (target=%s)",
                len(consistent), G_F.number_of_nodes(), len(res.nodes), target_id)

    # 4.5 Optional connectivity repair, BEFORE annotation so the single FVA pass tags
    # the added reactions too (no double FVA; §7 of the spec). Additive only: when it
    # fires, res is augmented to H+ (same GudieResult shape, node_order/coupling/
    # threshold carried through) so all downstream code is unchanged.
    repair_info: RepairResult | None = None
    if repair and source_ids:
        repair_info = repair_connectivity(res, G_F, source_ids, target_id, config=cfg)
        if repair_info.repaired:
            h_plus = [n for n in res.node_order if n in repair_info.nodes]
            res = replace(res, subgraph=G_F.subgraph(h_plus).copy(), nodes=h_plus)

    # 5. Optional FVA design-annotation over H at the requested operating points.
    annotation = None
    if fva_operating_points is not None:
        annotation = _annotate_fva(reduced, res, target_id,
                                   float(sol.objective_value),
                                   fva_operating_points, fva_loopless)

    return PathwayResult(
        subgraph=res.subgraph,
        reactions=res.nodes,
        target_id=target_id,
        n_consistent=len(consistent),
        n_active=G_F.number_of_nodes(),
        gudie=res,
        annotation=annotation,
        repair=repair_info,
        consistent_ids=frozenset(consistent),
        fastcc_medium_key=fastcc_medium_key,
    )


def _annotate_fva(
    reduced: cobra.Model,
    res: GudieResult,
    target_id: str,
    v_star: float,
    points: Sequence[float],
    loopless: bool,
) -> PathwayAnnotation:
    """FVA design-annotation of H at each operating point (fraction of v*_target).

    For each fraction the target is pinned at frac*v* and FVA is run over H's
    reactions ONLY (never the full active graph 𝓕 — the LP-budget trap). Each
    reaction is classified backbone (0 not in [v_min,v_max], must carry flux) vs
    flexible, then composed with its FCA coupling into a strain-design tag:
        backbone+coupled   -> protected
        backbone+uncoupled -> nontarget_essential
        flexible+coupled   -> alt_route
        flexible+uncoupled -> ko_candidate
    """
    pts = [float(f) for f in points]
    if not pts:
        raise InvalidOperatingPoint("no operating points given")
    if len(pts) > MAX_OPERATING_POINTS:
        raise InvalidOperatingPoint(
            f"too many operating points ({len(pts)} > {MAX_OPERATING_POINTS})")
    for f in pts:
        if not (0.0 < f <= 1.0):
            raise InvalidOperatingPoint(f"operating point {f!r} outside (0, 1]")

    H = res.nodes
    budget = 2 * len(H) * len(pts) * (LOOPLESS_FACTOR if loopless else 1)
    if budget > ANNOTATION_LP_BUDGET:
        raise AnnotationBudgetExceeded(
            f"annotation needs ~{budget} LPs (> {ANNOTATION_LP_BUDGET}); "
            f"use fewer operating points, disable loopless, or narrow |H|={len(H)}")

    coupling_map = ({res.node_order[i]: float(res.coupling[i])
                     for i in range(len(res.node_order))}
                    if res.coupling is not None else {})

    target = reduced.reactions.get_by_id(target_id)
    per_reaction: dict[str, ReactionAnnotation] = {
        rid: ReactionAnnotation(
            reaction_id=rid,
            gene_id=(reduced.reactions.get_by_id(rid).gene_reaction_rule or ""),
            coupling=coupling_map.get(rid),
            per_op=[],
        )
        for rid in H
    }

    for frac in pts:
        level = frac * v_star
        with reduced:
            target.bounds = (level, level)
            try:
                v_op = pfba(reduced).fluxes
            except Exception:  # pFBA is illustrative only; FVA is the signal
                v_op = None
            ranges = flux_variability_analysis(
                reduced, reaction_list=list(H), processes=1,
                fraction_of_optimum=0.0,
                loopless=("cycleFreeFlux" if loopless else None))
        for rid in H:
            lo = float(ranges.at[rid, "minimum"])
            hi = float(ranges.at[rid, "maximum"])
            backbone = lo > FASTCC_EPS or hi < -FASTCC_EPS
            cls = "backbone" if backbone else "flexible"
            c = coupling_map.get(rid)
            coupled = c is not None and c > COUPLING_TAU
            if backbone:
                design = "protected" if coupled else "nontarget_essential"
            else:
                design = "alt_route" if coupled else "ko_candidate"
            per_reaction[rid].per_op.append({
                "frac": frac,
                "level": level,
                "v": (float(v_op[rid]) if v_op is not None else None),
                "v_min": lo,
                "v_max": hi,
                "width": hi - lo,
                "cls": cls,
                "design": design,
            })

    return PathwayAnnotation(
        operating_points=pts,
        v_star_target=v_star,
        coupling_tau=COUPLING_TAU,
        per_reaction=per_reaction,
    )
