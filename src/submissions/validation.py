"""
validation.py - The automated curation gate.

With no user accounts and no reviewer queue, this is where curation actually
happens: a submission that passes here reaches a human already knowing the
model parses, grows, carries GPRs, and has the structure every downstream
engine assumes. Severities were calibrated by running every check over the 11
curated GEMs -- **no error-level check below rejects a model already in the
corpus**, and the warn thresholds sit clear of the worst real member
(GPR coverage 0.738 in iMM904, 14 transport reactions in iCN1361, a multi-term
objective in iAF987, unnamed compartments in iJN1026, autotrophy in iLB1027).

``validate_gem`` is pure: a parsed model in, a report out, no DB. The
DB-dependent checks (namespace overlap, duplicates, name collisions) live in
``service.validate_submission`` because they need the catalogue.

A failing check is not an exception -- it is the recorded outcome of
validation, with enough ``detail`` for a reviewer to act without opening the
model.
"""

from __future__ import annotations

import cobra

#: CO2 / bicarbonate. Excluded from the *organic* carbon test but counted as a
#: carbon source: the corpus contains autotrophs (iLB1027 grows at full rate on
#: CO2 alone), so an organic-only rule would reject a curated phototroph.
INORGANIC_CARBON = {"CO2", "CHO3", "CO3", "HCO3"}

GROWTH_EPS = 1e-6
MIN_GPR_COVERAGE = 0.5
BIOMASS_PATTERN = ("biomass", "growth")

ERROR, WARN, INFO = "error", "warn", "info"


def _check(report: dict, cid: str, name: str, severity: str, passed: bool, **detail):
    report["checks"].append({"id": cid, "name": name, "severity": severity,
                             "passed": bool(passed), "detail": detail})
    if not passed:
        if severity == ERROR:
            report["passed"] = False
            report["errors"].append(name)
        elif severity == WARN:
            report["warnings"].append(name)
    return passed


def carbon_uptakes(model: cobra.Model) -> tuple[list[str], list[str]]:
    """Open carbon-bearing uptake exchanges, split organic / inorganic."""
    organic, inorganic = [], []
    exchange_ids = {r.id for r in model.exchanges}
    for rxn in model.boundary:
        if rxn.lower_bound >= 0 or rxn.id not in exchange_ids:
            continue
        if sum(met.elements.get("C", 0) for met in rxn.metabolites) <= 0:
            continue
        formulas = {(met.formula or "").replace(" ", "") for met in rxn.metabolites}
        (inorganic if formulas <= INORGANIC_CARBON else organic).append(rxn.id)
    return organic, inorganic


def _growth_without(model: cobra.Model, rxn_ids: list[str]) -> float:
    """Growth with those uptakes closed, on a scoped copy of the bounds."""
    with model:
        for rxn_id in rxn_ids:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.bounds = (0.0, max(0.0, rxn.upper_bound))
        return float(model.slim_optimize(error_value=0.0) or 0.0)


def _objective_reactions(model: cobra.Model) -> list[str]:
    return sorted(str(r.id) for r in
                  cobra.util.solver.linear_reaction_coefficients(model))


def validate_gem(model: cobra.Model, *, conditions: list[dict] | None = None,
                 trophic_mode: str | None = None,
                 upload_format: str = "json") -> dict:
    """Run V2-V15 over a parsed model. Never raises for a bad model."""
    report: dict = {"passed": True, "checks": [], "errors": [], "warnings": []}
    conditions = conditions or []

    internal = [r for r in model.reactions if not r.boundary]

    _check(report, "V2", "metabolites_present", ERROR,
           len(model.metabolites) > 0, n=len(model.metabolites))
    _check(report, "V3", "reactions_present", ERROR,
           len(model.reactions) > 0, n=len(model.reactions))

    with_gpr = [r for r in internal if (r.gene_reaction_rule or "").strip()]
    _check(report, "V4", "genes_and_gpr_present", ERROR,
           bool(model.genes) and bool(with_gpr),
           n_genes=len(model.genes), n_gpr=len(with_gpr),
           hint=None if with_gpr else
           "no GPRs: OptForce and the KO engines cannot translate a reaction "
           "intervention into a construct")
    coverage = len(with_gpr) / len(internal) if internal else 0.0
    _check(report, "V4b", "gpr_coverage", WARN,
           coverage >= MIN_GPR_COVERAGE, coverage=round(coverage, 3))

    _check(report, "V5", "exchanges_present", ERROR,
           len(model.exchanges) > 0, n=len(model.exchanges))
    _check(report, "V5b", "boundary_gap", INFO, True,
           n_boundary=len(model.boundary), n_exchanges=len(model.exchanges))

    transport = [r for r in internal
                 if len({m.compartment for m in r.metabolites}) >= 2]
    _check(report, "V6", "transport_present", ERROR, bool(transport),
           n=len(transport), n_compartments=len(model.compartments))

    objective = _objective_reactions(model)
    biomass_like = [r.id for r in model.reactions
                    if any(p in (r.id or "").lower() or p in (r.name or "").lower()
                           for p in BIOMASS_PATTERN)]
    _check(report, "V7", "biomass_available", ERROR,
           bool(objective) and bool(biomass_like),
           objective=objective[:5], n_biomass_like=len(biomass_like))
    _check(report, "V7b", "single_reaction_objective", WARN,
           len(objective) == 1, n_objective=len(objective))

    organic, inorganic = carbon_uptakes(model)
    _check(report, "V8", "carbon_source_present", ERROR,
           bool(organic or inorganic), organic=organic[:10], inorganic=inorganic)

    baseline = float(model.slim_optimize(error_value=0.0) or 0.0)
    if organic or inorganic:
        closed_all = _growth_without(model, organic + inorganic)
        _check(report, "V8b", "carbon_source_load_bearing", WARN,
               closed_all <= GROWTH_EPS,
               growth=round(baseline, 6), growth_carbon_closed=round(closed_all, 6),
               hint=None if closed_all <= GROWTH_EPS else
               "grows with all carbon exchanges closed -- carbon enters via a "
               "sink/demand; confirm that is intentional")
        closed_organic = _growth_without(model, organic)
        detected = "autotroph" if closed_organic > GROWTH_EPS else "heterotroph"
        _check(report, "V8c", "trophic_mode", WARN,
               trophic_mode is None or trophic_mode == detected,
               detected=detected, declared=trophic_mode,
               growth_organic_closed=round(closed_organic, 6))

    feasible, infeasible = [], []
    for condition in conditions:
        with model:
            for rxn_id, bounds in (condition.get("exchange_bounds") or {}).items():
                if rxn_id in model.reactions:
                    model.reactions.get_by_id(rxn_id).bounds = tuple(bounds)
            value = float(model.slim_optimize(error_value=0.0) or 0.0)
        (feasible if value > GROWTH_EPS else infeasible).append(
            {"condition_id": condition["condition_id"], "growth": round(value, 6)})
    if conditions:
        _check(report, "V9", "growth_feasible", ERROR, bool(feasible),
               feasible=feasible, infeasible=infeasible)

    unnamed = [cid for cid, cname in model.compartments.items() if not cname]
    no_compartment = [m.id for m in model.metabolites if not m.compartment]
    _check(report, "V10", "compartment_info", WARN,
           len(model.compartments) >= 2 and not no_compartment and not unnamed,
           n_compartments=len(model.compartments), unnamed=unnamed[:10],
           metabolites_without_compartment=no_compartment[:10])

    no_formula = [m.id for m in model.metabolites if not m.formula]
    _check(report, "V11", "formula_coverage", WARN, not no_formula,
           n_missing=len(no_formula), sample=no_formula[:10])

    with_subsystem = sum(1 for r in model.reactions if (r.subsystem or "").strip())
    sub_coverage = with_subsystem / len(model.reactions) if model.reactions else 0.0
    # Error only for SBML with nothing left after the group repair: the
    # subsystems are simply gone, and reactions.subsystem is an indexed
    # catalogue column the pathway atlas colours by.
    _check(report, "V15", "subsystem_coverage",
           ERROR if (upload_format == "sbml" and with_subsystem == 0) else WARN,
           with_subsystem > 0, coverage=round(sub_coverage, 3),
           hint=None if with_subsystem else
           "no subsystems survived conversion -- re-export the SBML with groups, "
           "or submit the cobra-JSON")

    report["summary"] = {
        "n_reactions": len(model.reactions), "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes), "objective": objective[:5],
        "growth": round(baseline, 6), "n_transport": len(transport),
        "carbon_sources": {"organic": organic[:10], "inorganic": inorganic},
    }
    return report
