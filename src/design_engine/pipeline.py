"""
pipeline.py - Unified design-engine entry point.

:func:`rank_intervention_candidates` runs one or more discovery strategies
(OptForce, OptKnock, or both), pools their reaction-level proposals, and feeds
the pool to a single shared verification/audit/ranking core. The result is one
Pareto front across whichever engines were requested.

Two guarantees distinguish this from the former per-engine pipelines:

* **Buildable by default.** With ``require_buildable=True`` (the default), every
  ranked candidate is a knockout design for which a concrete gene-deletion
  construct both grows and holds its production floor on the full model.
  Regulatory (up/down) and mixed designs are never ranked; they are returned,
  verified, in ``RankingResult.hypotheses``.
* **One ranking.** Proposals from every engine are co-ranked on identical
  full-model verified values, so a cheaper design from either engine can
  dominate a deeper one from the other.

Discovery of the selected strategies can run in parallel (``parallel=True``),
each strategy in its own process on its own model copy; results are pooled only
after all strategies finish, so ranking is unaffected by concurrency.
"""

from __future__ import annotations

import time
from dataclasses import replace
from typing import Optional

import numpy as np
import cobra

from . import errors as E
from .schema import Proposal, InterventionCandidate, RankingResult, SUPPORTED_ACTIONS
from .lpkernel import make_lp
from . import verify
from . import rank as _rank
from .discovery.base import Baseline, DiscoveryParams, get_strategy

PARETO_OBJECTIVES = _rank.PARETO_OBJECTIVES
DEFAULT_METHODS = ("optforce", "optknock")


# ---------------------------------------------------------------------------
# Proposal pooling with source-merge dedup
# ---------------------------------------------------------------------------

def _split_source(src: str) -> tuple[set[str], str]:
    """('optforce+optknock_k2') -> ({'optforce','optknock'}, 'k2')."""
    if "_k" in src:
        eng, k = src.rsplit("_k", 1)
        return set(eng.split("+")), "k" + k
    return set(src.split("+")), ""


def _merge_source(a: str, b: str) -> str:
    eng_a, ka = _split_source(a)
    eng_b, kb = _split_source(b)
    engines = "+".join(sorted(eng_a | eng_b))
    if ka == kb and ka:
        return f"{engines}_{ka}"
    # Different cardinality tags (unusual for an identical action set): keep both.
    return "|".join(sorted({a, b}))


def _pool_proposals(proposals: list[Proposal]) -> list[Proposal]:
    """Dedup identical action sets, merging their sources so provenance from
    both engines is preserved (e.g. 'optforce+optknock_k1')."""
    groups: dict[tuple, Proposal] = {}
    order: list[tuple] = []
    for p in proposals:
        key = tuple(sorted(a.identifier() for a in p.actions))
        if key not in groups:
            groups[key] = p
            order.append(key)
        else:
            cur = groups[key]
            groups[key] = Proposal(
                actions=cur.actions,
                source=_merge_source(cur.source, p.source),
                discovery_floor=max((x for x in (cur.discovery_floor, p.discovery_floor)
                                     if x is not None), default=None),
                discovery_growth=max((x for x in (cur.discovery_growth, p.discovery_growth)
                                      if x is not None), default=None),
            )
    return [groups[k] for k in order]


# ---------------------------------------------------------------------------
# Discovery dispatch (serial or parallel)
# ---------------------------------------------------------------------------

def _run_strategy(name: str, model: cobra.Model, target: str, biomass: str,
                  substrate: Optional[str], baseline: Baseline,
                  params: DiscoveryParams) -> tuple[str, list[Proposal], list[str]]:
    """Top-level (picklable) worker: run one strategy, return its proposals."""
    strat = get_strategy(name)
    notes: list[str] = []
    props = strat.discover(model, target, biomass, substrate, baseline, params, notes)
    return name, props, notes


def _discover_all(methods, disc_model, target, biomass, substrate, baseline,
                  params, parallel, notes) -> list[Proposal]:
    if parallel and len(methods) > 1:
        try:
            from concurrent.futures import ProcessPoolExecutor
            pooled: list[Proposal] = []
            with ProcessPoolExecutor(max_workers=len(methods)) as ex:
                # Each worker gets its own model copy so concurrent solves never
                # touch a shared solver.
                futures = [ex.submit(_run_strategy, m, disc_model.copy(), target,
                                     biomass, substrate, baseline, params)
                           for m in methods]
                results = [f.result() for f in futures]
            # Deterministic order: by the requested method order.
            by_name = {name: (props, ns) for name, props, ns in results}
            for m in methods:
                props, ns = by_name.get(m, ([], []))
                pooled.extend(props)
                notes.extend(ns)
            return pooled
        except Exception as exc:  # pickling / spawn / solver issue -> serial
            notes.append(f"parallel discovery failed ({exc!r}); ran serially.")

    pooled = []
    for m in methods:
        _, props, ns = _run_strategy(m, disc_model, target, biomass, substrate,
                                     baseline, params)
        pooled.extend(props)
        notes.extend(ns)
    return pooled


# ---------------------------------------------------------------------------
# Discovery-model validation (reduced-model extraction fallback)
# ---------------------------------------------------------------------------

def _validate_discovery(discovery_model, full_model, target, biomass, substrate,
                        class_tol, notes):
    for rid in (target, biomass):
        if rid not in discovery_model.reactions:
            raise E.DiscoveryModelInvalid(
                f"Discovery model is missing reaction {rid!r}.",
                "Rebuild the reduced model with the target and biomass retained.")
    with discovery_model:
        mu, status = verify.solve(discovery_model, biomass, "max")
    if status != "optimal" or not np.isfinite(mu) or mu <= class_tol:
        notes.append("extraction fallback: discovery model growth not preserved.")
        return True, "growth_not_preserved", full_model
    with discovery_model:
        cap, s2 = verify.solve(discovery_model, target, "max")
    if s2 != "optimal" or not np.isfinite(cap) or cap <= class_tol:
        notes.append("extraction fallback: discovery model producibility not preserved.")
        return True, "producibility_not_preserved", full_model
    return False, None, discovery_model


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def rank_intervention_candidates(
    full_model: cobra.Model,
    target_rxn_id: str,
    biomass_rxn_id: str,
    substrate_uptake_rxn_id: Optional[str] = None,
    *,
    methods: tuple[str, ...] = DEFAULT_METHODS,
    discovery_model: Optional[cobra.Model] = None,
    require_buildable: bool = True,
    parallel: bool = False,
    max_knockouts: int = 1,
    enable_k2: bool = False,
    min_growth_frac: float = 0.05,
    min_growth_abs: Optional[float] = None,
    gamma: float = 0.001,
    eps_prod: float = 1e-6,
    class_tol: float = 1e-6,
    backend: str = "cobra",
    top_k: Optional[int] = None,
    solver_timeout: Optional[float] = None,
    gpr_cap: int = verify.DEFAULT_GPR_CAP,
) -> RankingResult:
    """Discover (via ``methods``), verify, audit, and Pareto-rank interventions.

    ``methods`` selects the discovery strategies to run and pool (default: both
    OptForce and OptKnock). ``require_buildable`` (default True) restricts the
    ranked set to knockout designs with a buildable gene-deletion construct;
    regulatory and mixed designs are returned unranked in ``hypotheses``. Both
    supplied models are left unmodified. Ranking uses full-model verified values
    only.
    """
    methods = tuple(dict.fromkeys(methods))  # dedup, preserve order
    if not methods:
        raise E.DesignEngineError("No discovery methods selected.",
                                  "Pass methods=('optforce',), ('optknock',), or both.")
    for m in methods:
        try:
            get_strategy(m)
        except KeyError:
            raise E.DesignEngineError(
                f"Unknown discovery method {m!r}.",
                "Supported methods are 'optforce' and 'optknock'.")
    if backend != "cobra" and "highs" not in backend:
        raise E.DesignEngineError(
            f"Unsupported backend {backend!r}.",
            "Use backend='cobra' (reference) or 'highs' (accelerated HiGHS kernel).")
    if not (0.0 < gamma < 1.0):
        raise E.InvalidGamma(f"gamma must be in (0, 1); got {gamma}.")
    max_knockouts = max(1, min(3, int(max_knockouts)))

    timings: dict[str, float] = {}
    notes: list[str] = []

    for rid in (target_rxn_id, biomass_rxn_id):
        if rid not in full_model.reactions:
            raise E.MissingFullModelReaction(f"Reaction {rid!r} not found on the full model.")

    try:
        lp = make_lp(full_model, backend, timeout=solver_timeout)
    except ImportError as exc:
        raise E.DesignEngineError(str(exc), "Install highspy or use backend='cobra'.") from exc

    with verify.timeout_scope(full_model, solver_timeout):
        # --- Stage 1: baseline ---
        t0 = time.perf_counter()
        mu_wt, floor_wt, ceiling_wt = verify.wt_baseline(
            lp, target_rxn_id, biomass_rxn_id, gamma, class_tol)
        timings["baseline"] = time.perf_counter() - t0

        substrate = substrate_uptake_rxn_id or verify.detect_substrate(full_model)
        min_growth = min_growth_abs if min_growth_abs is not None else min_growth_frac * mu_wt
        baseline = Baseline(mu_wt, floor_wt, ceiling_wt)
        params = DiscoveryParams(
            min_growth_frac=min_growth_frac, min_growth_abs=min_growth_abs,
            eps_prod=eps_prod, class_tol=class_tol, gamma=gamma,
            enable_k2=enable_k2, max_knockouts=max_knockouts,
            backend=backend, solver_timeout=solver_timeout)

        # --- Stage 2: discovery (with extraction fallback) ---
        t0 = time.perf_counter()
        extraction_fallback = False
        extraction_reason: Optional[str] = None
        disc = discovery_model if discovery_model is not None else full_model
        if discovery_model is not None:
            extraction_fallback, extraction_reason, disc = _validate_discovery(
                discovery_model, full_model, target_rxn_id, biomass_rxn_id,
                substrate, class_tol, notes)
        with verify.timeout_scope(disc, solver_timeout):
            raw = _discover_all(methods, disc, target_rxn_id, biomass_rxn_id,
                                substrate, baseline, params, parallel, notes)
        timings["discovery"] = time.perf_counter() - t0

        unique = _pool_proposals(raw)

        # --- Stages 3-4: annotate, verify, audit, route ---
        accepted: list[InterventionCandidate] = []
        rejected: list[InterventionCandidate] = []
        hypotheses: list[InterventionCandidate] = []
        for p in unique:
            cand = verify.process_proposal(
                lp, target_rxn_id, biomass_rxn_id, p, mu_wt, floor_wt,
                ceiling_wt, gamma, class_tol, eps_prod, min_growth, gpr_cap, timings)
            if cand.rejection_reason is not None or cand.import_loop_status != E.IMPORT_PASSES:
                rejected.append(cand)
            elif not cand.all_knockout:
                # Regulatory or mixed design: a verified hypothesis, never ranked.
                hypotheses.append(cand)
            elif require_buildable and not cand.buildable:
                rejected.append(replace(cand, rejection_reason=E.NOT_BUILDABLE,
                                        passes_checks=False))
            else:
                accepted.append(cand)

        # --- Stage 5: Pareto ranking (buildable KO designs only) ---
        t0 = time.perf_counter()
        best_k1_floor = max((c.verified_floor for c in accepted if c.cardinality == 1),
                            default=None)
        accepted = _rank.finalize(accepted, best_k1_floor, mu_wt, class_tol, top_k)
        pareto_objs = _rank.batch_objectives(accepted)
        timings["ranking"] = time.perf_counter() - t0

    rejected.sort(key=lambda c: (c.rejection_reason or c.import_loop_status,
                                 tuple(sorted(a.identifier() for a in c.reaction_actions))))
    hypotheses.sort(key=lambda c: (-_rank._q(c.verified_floor, class_tol),
                                   tuple(sorted(a.identifier() for a in c.reaction_actions))))

    return RankingResult(
        accepted=accepted, rejected=rejected, hypotheses=hypotheses,
        mu_wt=mu_wt, floor_wt=floor_wt, ceiling_wt=ceiling_wt,
        substrate_rxn=substrate, methods=methods,
        extraction_fallback=extraction_fallback,
        extraction_failure_reason=extraction_reason,
        pareto_objectives=pareto_objs, timings=timings, notes=notes,
    )
