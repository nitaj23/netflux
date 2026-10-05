"""
verify.py - Shared full-model verification, audit, and GPR translation.

Every discovery strategy hands the pipeline reaction-level :class:`Proposal`
objects; this module turns each into a fully verified
:class:`InterventionCandidate` using the *full* model only. It contains the
machinery that used to be copied byte-for-byte between the OptForce and OptKnock
pipelines: the wild-type baseline, substrate detection, reaction-bound
verification, the import-loop audit, and the GPR knockout -> gene-construct
translation.

One behavioral change from the old pipelines lives here: the gene-construct
search is *early-exit over all minimal gene sets* rather than a fixed top-4
slice, and it records whether any construct is buildable. The pipeline uses that
flag to guarantee every ranked result is buildable.
"""

from __future__ import annotations

import time
from contextlib import contextmanager
from typing import Optional

import numpy as np
import cobra

from . import errors as E
from .gpr import annotate, DEFAULT_GPR_CAP, _dnf_terms, _minimal_hitting_sets
from .schema import (
    ReactionAction, Proposal, GeneOption, GeneConstructResult, InterventionCandidate,
)

# Hard ceiling on how many minimal gene sets the buildable search verifies for a
# single design. The search stops at the first buildable set; this only bounds
# the worst case (a design with no buildable construct). Tied to gpr_cap.
_CONSTRUCT_SEARCH_CEILING = DEFAULT_GPR_CAP


# ---------------------------------------------------------------------------
# Low-level solver helpers
# ---------------------------------------------------------------------------

@contextmanager
def timeout_scope(model: cobra.Model, seconds: Optional[float]):
    if seconds is None:
        yield
        return
    cfg = model.solver.configuration
    old = getattr(cfg, "timeout", None)
    try:
        # GLPK's time limit is an integer (ms); coerce so a float timeout from a
        # caller never trips glp_smcp_tm_lim_set. HiGHS gets the float directly
        # via make_lp(..., timeout=...), so this only affects the cobra path.
        cfg.timeout = int(round(seconds))
        yield
    finally:
        try:
            cfg.timeout = old
        except Exception:
            pass


def solve(model: cobra.Model, rxn_id: str, direction: str) -> tuple[float, str]:
    """Optimize a single reaction on a cobra model (discovery-side handoff)."""
    model.objective = model.reactions.get_by_id(rxn_id)
    model.objective_direction = direction
    value = model.slim_optimize()   # NaN when not optimal
    return float(value), model.solver.status


def _solve_or_raise(lp, rxn_id: str, direction: str, what: str) -> float:
    value, status = (lp.maximize(rxn_id) if direction == "max" else lp.minimize(rxn_id))
    if status == "time_limit":
        raise E.SolverTimeout(f"Solver timed out computing {what}.",
                              "Raise solver_timeout or reduce the model size.")
    if status != "optimal" or not np.isfinite(value):
        raise E.SolverFailure(f"Could not compute {what} (status={status}).",
                              "Check the model medium and bounds.")
    return value


def _floor_ceiling(lp, target: str, biomass: str, mu: float,
                   gamma: float) -> tuple[float, str, float, str]:
    """Min and max target flux with biomass held near ``mu`` (>= (1-gamma)*mu)."""
    lp.set_bound(biomass, lb=(1.0 - gamma) * mu)
    lo, s_lo = lp.minimize(target)
    hi, s_hi = lp.maximize(target)
    return lo, s_lo, hi, s_hi


# ---------------------------------------------------------------------------
# Stage 1: full-model wild-type baseline
# ---------------------------------------------------------------------------

def wt_baseline(lp, target: str, biomass: str, gamma: float,
                class_tol: float) -> tuple[float, float, float]:
    with lp.scope():
        mu_wt = _solve_or_raise(lp, biomass, "max", "wild-type growth")
        if mu_wt <= class_tol:
            raise E.BiomassNotGrowing(
                f"Wild-type biomass optimum is {mu_wt:g} <= {class_tol:g}.",
                "Check the medium / forced-uptake bounds.")
        target_cap = _solve_or_raise(lp, target, "max", "target producibility")
        if target_cap <= class_tol:
            raise E.TargetNotProducible(
                f"Target producibility is {target_cap:g} <= {class_tol:g}.",
                "Loosen an uptake bound or choose another target.")
        floor_wt, s_lo, ceiling_wt, s_hi = _floor_ceiling(lp, target, biomass, mu_wt, gamma)
        if s_lo != "optimal" or s_hi != "optimal":
            raise E.BaseModelInfeasible(
                "Wild-type production envelope is infeasible at the biomass floor.",
                "Check gamma and the medium.")
    if floor_wt > ceiling_wt + class_tol:
        raise E.BaseModelInfeasible(
            f"Wild-type floor {floor_wt:g} exceeds ceiling {ceiling_wt:g}.")
    return mu_wt, floor_wt, ceiling_wt


# ---------------------------------------------------------------------------
# Substrate detection (largest organic-carbon uptake)
# ---------------------------------------------------------------------------

_INORGANIC_CARBON_FORMULAS = {"CO2", "CHO3", "CO3", "HCO3"}


def detect_substrate(model: cobra.Model) -> Optional[str]:
    """Largest organic-carbon uptake route."""
    best, best_uptake = None, 0.0
    for rxn in model.boundary:
        if rxn.lower_bound >= 0:
            continue
        carbons = sum(met.elements.get("C", 0) for met in rxn.metabolites)
        if carbons <= 0:
            continue
        formulas = {(met.formula or "").replace(" ", "") for met in rxn.metabolites}
        if formulas <= _INORGANIC_CARBON_FORMULAS:
            continue
        uptake = -rxn.lower_bound
        if uptake > best_uptake:
            best, best_uptake = rxn.id, uptake
    return best


# ---------------------------------------------------------------------------
# Direction-flip detection
# ---------------------------------------------------------------------------

def _is_direction_flip(v_req: float, wt_lo: float, wt_hi: float, tol: float) -> bool:
    """True when the forced flux crosses the WT feasible direction."""
    if v_req < -tol and wt_lo >= -tol and wt_hi > tol:
        return True
    if v_req > tol and wt_hi <= tol and wt_lo < -tol:
        return True
    return False


def _wt_interval(lp, rxn_id: str, biomass: str, mu_wt: float,
                 gamma: float) -> tuple[float, float]:
    with lp.scope():
        lp.set_bound(biomass, lb=(1.0 - gamma) * mu_wt)
        lo, _ = lp.minimize(rxn_id)
        hi, _ = lp.maximize(rxn_id)
    return lo, hi


# ---------------------------------------------------------------------------
# Import-loop audit
# ---------------------------------------------------------------------------

def _target_metabolites(model: cobra.Model, target: str) -> list[cobra.Metabolite]:
    rxn = model.reactions.get_by_id(target)
    return [m for m, c in rxn.metabolites.items() if c < 0]


def _import_routes(model: cobra.Model, target_mets, target: str):
    met_ids = {m.id for m in target_mets}
    routes = []
    for rxn in model.boundary:
        if rxn.id == target:
            continue
        for met, coeff in rxn.metabolites.items():
            if met.id in met_ids and rxn.lower_bound < 0:
                routes.append((rxn.id, coeff))
                break
    return routes


def _elemental_carbon_ok(lp, target: str, tol: float) -> Optional[bool]:
    value, status, fluxes = lp.solve_primal(target, "min")
    if status != "optimal":
        return None
    c_in = c_out = 0.0
    any_formula = False
    for rxn in lp.model.boundary:
        flux = fluxes.get(rxn.id, 0.0)
        if abs(flux) < tol:
            continue
        for met, coeff in rxn.metabolites.items():
            c = met.elements.get("C", 0)
            if c:
                any_formula = True
                moved = coeff * flux * c
                if moved > 0:
                    c_in += moved
                else:
                    c_out += -moved
    if not any_formula:
        return None
    return c_out <= c_in + max(tol, 1e-6 * max(c_in, 1.0))


def import_audit(lp, target: str, biomass: str,
                 actions: tuple[ReactionAction, ...], mu_star: float,
                 verified_floor: float, gamma: float, eps_prod: float,
                 class_tol: float) -> dict:
    """Counterfactual target-import closure test."""
    target_mets = _target_metabolites(lp.model, target)
    if not target_mets:
        return {"status": E.IMPORT_AMBIGUOUS, "routes": (), "closed_floor": None,
                "observed_inward": (), "supported_fraction": None, "elemental": "unavailable"}

    with lp.scope():
        for a in actions:
            lp.apply_action(a.reaction_id, a.action, a.forced_flux)
        routes = _import_routes(lp.model, target_mets, target)
        if not routes:
            elem = _elemental_carbon_ok(lp, target, class_tol)
            return {"status": E.IMPORT_PASSES, "routes": (), "closed_floor": None,
                    "observed_inward": (), "supported_fraction": 0.0,
                    "elemental": "unavailable" if elem is None else ("ok" if elem else "carbon_exceeded")}

        lp.set_bound(biomass, lb=(1.0 - gamma) * mu_star)
        _, status, fluxes = lp.solve_primal(target, "min")
        observed_inward = []
        for rxn_id, coeff in routes:
            if status != "optimal":
                continue
            flux = float(fluxes.get(rxn_id, 0.0))
            if coeff * flux > class_tol:
                observed_inward.append((rxn_id, flux))

        for rxn_id, _ in routes:
            lb, ub = lp.bounds(rxn_id)
            new_ub = 0.0 if ub < 0 else ub
            lp.set_bound(rxn_id, lb=0.0, ub=new_ub)
        mu_closed, s_mu = lp.maximize(biomass)
        if s_mu != "optimal" or not np.isfinite(mu_closed) or mu_closed <= class_tol:
            closed_floor = 0.0
        else:
            lp.set_bound(biomass, lb=(1.0 - gamma) * mu_closed)
            closed_floor, s_cf = lp.minimize(target)
            if s_cf != "optimal" or not np.isfinite(closed_floor):
                closed_floor = 0.0
        elem = _elemental_carbon_ok(lp, target, class_tol)

    supported = (max(0.0, verified_floor - closed_floor) / verified_floor
                 if verified_floor > class_tol else None)
    used_import = bool(observed_inward)
    loop = used_import and closed_floor < eps_prod
    return {
        "status": E.IMPORT_LOOP_STATUS if loop else E.IMPORT_PASSES,
        "routes": tuple(r for r, _ in routes),
        "closed_floor": closed_floor,
        "observed_inward": tuple(observed_inward),
        "supported_fraction": supported,
        "elemental": "unavailable" if elem is None else ("ok" if elem else "carbon_exceeded"),
    }


# ---------------------------------------------------------------------------
# Reaction-bound verification
# ---------------------------------------------------------------------------

def _pinned_growth(lp, biomass: str,
                   actions: tuple[ReactionAction, ...]) -> Optional[float]:
    with lp.scope():
        for a in actions:
            lp.set_bound(a.reaction_id, lb=a.forced_flux, ub=a.forced_flux)
        val, status = lp.maximize(biomass)
    return val if status == "optimal" and np.isfinite(val) else None


def _verify_reaction_bounds(lp, target: str, biomass: str,
                            actions: tuple[ReactionAction, ...], mu_wt: float,
                            gamma: float, class_tol: float, eps_prod: float) -> dict:
    with lp.scope():
        for a in actions:
            lp.apply_action(a.reaction_id, a.action, a.forced_flux)
        bio_lb0, bio_ub0 = lp.bounds(biomass)
        mu_star, s_mu = lp.maximize(biomass)
        if s_mu == "time_limit":
            return {"reason": E.CANDIDATE_INFEASIBLE, "mu_star": None}
        if s_mu != "optimal" or not np.isfinite(mu_star):
            return {"reason": E.CANDIDATE_INFEASIBLE, "mu_star": None}
        lp.set_bound(biomass, lb=(1.0 - gamma) * mu_star)
        floor, s_lo = lp.minimize(target)
        ceiling, s_hi = lp.maximize(target)
        if s_lo != "optimal" or not np.isfinite(floor):
            return {"reason": E.CANDIDATE_INFEASIBLE, "mu_star": mu_star}
        if s_hi != "optimal" or not np.isfinite(ceiling):
            ceiling = floor
        tgt_lb, tgt_ub = lp.bounds(target)
        lp.set_bound(biomass, lb=bio_lb0, ub=bio_ub0)
        # "Target off" diagnostic: cap target flux at eps_prod. If the medium
        # forces a secretion floor (tgt_lb > eps_prod, e.g. iNF517 EX_ac_e), drop
        # the lower bound too so lb <= ub -- otherwise set_bound would raise.
        lp.set_bound(target, lb=min(tgt_lb, eps_prod), ub=eps_prod)
        mu_no_prod, s_np = lp.maximize(biomass)
        lp.set_bound(target, lb=tgt_lb, ub=tgt_ub)
    if s_np != "optimal" or not np.isfinite(mu_no_prod):
        U, mandatory = mu_star, True
    else:
        U, mandatory = max(0.0, min(mu_star, mu_star - mu_no_prod)), False
    return {"reason": None, "mu_star": mu_star, "floor": floor, "ceiling": ceiling,
            "U": U, "mandatory": mandatory}


# ---------------------------------------------------------------------------
# Gene translation & construct verification (early-exit buildable search)
# ---------------------------------------------------------------------------

def _verify_gene_construct(lp, target, biomass, gene_set, ko_rxns, actions,
                           mu_wt, gamma, class_tol, eps_prod) -> GeneConstructResult:
    genes = set(gene_set)
    silenced = []
    for rxn in lp.model.reactions:
        rule = (rxn.gene_reaction_rule or "").strip()
        if not rule:
            continue
        if genes & {g.id for g in rxn.genes} and not rxn.gpr.eval(genes):
            silenced.append(rxn.id)
    silenced_sorted = tuple(sorted(silenced))

    with lp.scope():
        for rxn_id in silenced:
            lp.knockout(rxn_id)
        for a in actions:
            if a.action != "knockout":
                lp.apply_action(a.reaction_id, a.action, a.forced_flux)
        mu_star, s_mu = lp.maximize(biomass)
        if s_mu != "optimal" or not np.isfinite(mu_star) or mu_star <= class_tol:
            return GeneConstructResult(tuple(sorted(genes)), silenced_sorted,
                                       mu_star if np.isfinite(mu_star) else 0.0,
                                       0.0, 0.0, 0.0, E.IMPORT_UNAVAILABLE, False)
        lp.set_bound(biomass, lb=(1.0 - gamma) * mu_star)
        floor, s_lo = lp.minimize(target)
        ceiling, s_hi = lp.maximize(target)
    floor = floor if s_lo == "optimal" and np.isfinite(floor) else 0.0
    ceiling = ceiling if s_hi == "optimal" and np.isfinite(ceiling) else floor
    buildable = mu_star > class_tol and floor >= eps_prod
    return GeneConstructResult(
        tuple(sorted(genes)), silenced_sorted, mu_star, floor, ceiling,
        mu_star / mu_wt, E.IMPORT_PASSES, buildable)


def translate_and_verify_genes(lp, target: str, biomass: str,
                               actions: tuple[ReactionAction, ...], mu_wt: float,
                               gamma: float, class_tol: float, eps_prod: float,
                               gpr_cap: int) -> tuple[tuple[GeneOption, ...], str,
                                                      tuple[GeneConstructResult, ...]]:
    """Translate the knockout actions to gene sets and verify them on the full
    model, searching minimal gene sets smallest-first and stopping at the first
    buildable construct (so a buildable design past the old top-4 slice is not
    lost). Regulatory actions are reported as putative targets only."""
    model = lp.model
    ko_rxns = [model.reactions.get_by_id(a.reaction_id)
               for a in actions if a.action == "knockout"]
    reg_rxns = [model.reactions.get_by_id(a.reaction_id)
                for a in actions if a.action != "knockout"]

    gene_options: list[GeneOption] = []
    for rxn in reg_rxns:
        genes = tuple(sorted(g.id for g in rxn.genes))
        if genes:
            gene_options.append(GeneOption(genes, "putative_regulation", True))

    if not ko_rxns:
        status = E.GPR_PUTATIVE if reg_rxns else E.GPR_REACTION_ONLY
        return tuple(gene_options), status, ()

    combined_terms: list = []
    empty_gpr = False
    parse_fail = False
    for rxn in ko_rxns:
        rule = (rxn.gene_reaction_rule or "").strip()
        if not rule:
            empty_gpr = True
            continue
        terms = _dnf_terms(rxn.gpr.body, gpr_cap)
        if terms is None:
            parse_fail = True
            break
        combined_terms.extend(terms)

    if parse_fail:
        return tuple(gene_options), E.GPR_FAILURE, ()
    if not combined_terms:
        return tuple(gene_options), E.GPR_NONE, ()

    gene_sets, complete = _minimal_hitting_sets(combined_terms, gpr_cap)
    verified_sets = []
    for s in gene_sets:
        if all(not rxn.gpr.eval(set(s)) for rxn in ko_rxns if (rxn.gene_reaction_rule or "").strip()):
            verified_sets.append(s)
    if not verified_sets:
        return tuple(gene_options), E.GPR_NONE if empty_gpr else E.GPR_FAILURE, ()

    for s in verified_sets:
        gene_options.append(GeneOption(tuple(sorted(s)), "knockout", complete))

    # Early-exit buildable search: verify minimal gene sets in order (they are
    # already smallest-first) and stop at the first buildable construct. Bounded
    # by the search ceiling so a design with no buildable set still terminates.
    constructs: list[GeneConstructResult] = []
    for s in verified_sets[:_CONSTRUCT_SEARCH_CEILING]:
        gc = _verify_gene_construct(
            lp, target, biomass, s, ko_rxns, actions, mu_wt, gamma, class_tol, eps_prod)
        constructs.append(gc)
        if gc.buildable:
            break

    if empty_gpr:
        status = E.GPR_NONE
    elif len([o for o in gene_options if o.kind == "knockout"]) == 1 and complete:
        status = E.GPR_EXACT_KO
    elif not complete:
        status = E.GPR_TRANSLATION_INCOMPLETE
    else:
        status = E.GPR_ALTERNATIVE_KO
    return tuple(gene_options), status, tuple(constructs)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def _coupling_class(floor_delta: float, ceiling_delta: float, class_tol: float) -> str:
    if floor_delta > class_tol:
        return "strong"
    if ceiling_delta > class_tol:
        return "weak"
    return "none"


# ---------------------------------------------------------------------------
# Per-proposal verification (produces an InterventionCandidate)
# ---------------------------------------------------------------------------

def process_proposal(lp, target, biomass, p: Proposal, mu_wt, floor_wt,
                     ceiling_wt, gamma, class_tol, eps_prod, min_growth, gpr_cap,
                     timings) -> InterventionCandidate:
    actions = p.actions
    cardinality = len(actions)
    all_knockout = all(a.action == "knockout" for a in actions)

    def _reject(reason, import_status=E.IMPORT_PASSES, **extra):
        base = dict(
            source=p.source, reaction_actions=actions, cardinality=cardinality,
            discovery_floor=p.discovery_floor, discovery_growth=p.discovery_growth,
            verified_floor=0.0, verified_ceiling=0.0, verified_growth=0.0,
            floor_delta=0.0, ceiling_delta=0.0, coupling_potential=None,
            retained_mu=0.0, mandatory_production=None, guaranteed_growth=0.0,
            pinned_growth=None, growth_robustness_ratio=None, coupling_class="none",
            gene_options=(), translation_status=E.GPR_REACTION_ONLY, gene_constructs=(),
            import_loop_status=import_status, target_import_reactions=(),
            closed_import_floor=None, annotations=(), artifact_flags=(),
            passes_checks=False, pareto_front=None, rank=None, rejection_reason=reason,
            all_knockout=all_knockout, buildable=False)
        base.update(extra)
        return InterventionCandidate(**base)

    for a in actions:
        if not lp.has_reaction(a.reaction_id):
            return _reject(E.MISSING_FULL_MODEL_REACTION)

    for a in actions:
        wt_lo, wt_hi = _wt_interval(lp, a.reaction_id, biomass, mu_wt, gamma)
        if _is_direction_flip(a.forced_flux, wt_lo, wt_hi, class_tol):
            return _reject(E.UNSUPPORTED_DIRECTION_FLIP)

    t0 = time.perf_counter()
    v = _verify_reaction_bounds(lp, target, biomass, actions, mu_wt, gamma,
                                class_tol, eps_prod)
    timings["verify"] = timings.get("verify", 0.0) + (time.perf_counter() - t0)
    if v["reason"] is not None:
        return _reject(v["reason"])
    mu_star = v["mu_star"]
    if mu_star <= max(min_growth, class_tol):
        return _reject(E.BELOW_GROWTH_FLOOR)
    verified_floor, verified_ceiling = v["floor"], v["ceiling"]
    if verified_floor < eps_prod:
        reason = (E.FULL_MODEL_BYPASS
                  if (p.discovery_floor or 0.0) > eps_prod else E.BELOW_PRODUCTION_THRESHOLD)
        return _reject(reason)

    floor_delta = verified_floor - floor_wt
    ceiling_delta = verified_ceiling - ceiling_wt
    retained_mu = mu_star / mu_wt

    pinned = _pinned_growth(lp, biomass, actions)
    robustness = pinned / mu_star if (pinned is not None and mu_star > class_tol) else None

    t0 = time.perf_counter()
    audit = import_audit(lp, target, biomass, actions, mu_star, verified_floor,
                         gamma, eps_prod, class_tol)
    timings["import_audit"] = timings.get("import_audit", 0.0) + (time.perf_counter() - t0)

    t0 = time.perf_counter()
    annotations = tuple(annotate(lp.model.reactions.get_by_id(a.reaction_id), gpr_cap)
                        for a in actions)
    gene_options, translation_status, gene_constructs = translate_and_verify_genes(
        lp, target, biomass, actions, mu_wt, gamma, class_tol, eps_prod, gpr_cap)
    timings["gpr"] = timings.get("gpr", 0.0) + (time.perf_counter() - t0)

    buildable = any(gc.buildable for gc in gene_constructs)
    coupling_class = _coupling_class(floor_delta, ceiling_delta, class_tol)

    import_status = audit["status"]
    rejection = E.IMPORT_LOOP if import_status == E.IMPORT_LOOP_STATUS else None

    flags = []
    if robustness is not None and pinned is not None and pinned < 0.5 * mu_star:
        flags.append("collapse")
    if mu_star < 0.10 * mu_wt:
        flags.append("starved")

    return InterventionCandidate(
        source=p.source, reaction_actions=actions, cardinality=cardinality,
        discovery_floor=p.discovery_floor, discovery_growth=p.discovery_growth,
        verified_floor=verified_floor, verified_ceiling=verified_ceiling,
        verified_growth=mu_star, floor_delta=floor_delta, ceiling_delta=ceiling_delta,
        coupling_potential=v["U"], retained_mu=retained_mu,
        mandatory_production=v["mandatory"], guaranteed_growth=mu_star,
        pinned_growth=pinned, growth_robustness_ratio=robustness,
        coupling_class=coupling_class, gene_options=gene_options,
        translation_status=translation_status, gene_constructs=gene_constructs,
        import_loop_status=import_status, target_import_reactions=audit["routes"],
        closed_import_floor=audit["closed_floor"], annotations=annotations,
        artifact_flags=tuple(flags), passes_checks=(rejection is None),
        pareto_front=None, rank=None, rejection_reason=rejection,
        all_knockout=all_knockout, buildable=buildable,
    )
