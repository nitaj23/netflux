"""
optforce.py - OptForce discovery strategy.

Wraps the existing OptForce candidate-discovery machinery (FVA -> MUST sets ->
FORCE) from the ``optforce`` package, which remains the single source of truth
for that algorithm. This strategy only adapts its output into the unified
:class:`~design_engine.schema.Proposal` vocabulary; verification, the import
audit, GPR translation, and ranking all happen later in the shared core.

OptForce may emit upregulation, downregulation, and knockout actions. A
knockout (flux forced to 0) is emitted explicitly by OptForce whenever zero is
attainable in the mutant, so there is no "downregulation that is really a
knockout" to reclassify here.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import cobra

from ..schema import Proposal, ReactionAction
from ..verify import solve
from .base import Baseline, DiscoveryParams, register

# The OptForce discovery algorithm (FVA -> MUST -> FORCE), vendored under this
# package so the design engine is fully self-contained. This is the algorithm
# only; the former OptForce ranking pipeline is superseded by the shared core.
from ._optforce_core.fva import fva_optforce, StrainConstraints
from ._optforce_core.must_sets import find_must_l, find_must_u, ConstrOpt
from ._optforce_core.optforce import run_optforce


class OptForceStrategy:
    name = "optforce"

    def discover(self, disc_model: cobra.Model, target: str, biomass: str,
                 substrate: Optional[str], baseline: Baseline,
                 params: DiscoveryParams, notes: list[str]) -> list[Proposal]:
        min_growth_frac = params.min_growth_frac
        eps_prod = params.eps_prod
        enable_k2 = params.enable_k2

        with disc_model:
            mu_wt, s_mu = solve(disc_model, biomass, "max")
            if s_mu != "optimal" or not np.isfinite(mu_wt) or mu_wt <= 0:
                notes.append("optforce: discovery model cannot grow; no candidates.")
                return []
            min_growth = min_growth_frac * mu_wt

            disc_model.reactions.get_by_id(biomass).lower_bound = min_growth
            ceiling, _ = solve(disc_model, target, "max")
            disc_model.reactions.get_by_id(biomass).lower_bound = 0.95 * mu_wt
            wt_max, _ = solve(disc_model, target, "max")

        if not (np.isfinite(ceiling) and np.isfinite(wt_max)):
            notes.append("optforce: target FVA non-finite; no candidates.")
            return []

        headroom = ceiling - wt_max
        epsilon_must = max(0.5 * ceiling / 155.0, 1e-4)
        if headroom <= max(0.05 * ceiling, 10.0 * epsilon_must):
            notes.append(f"optforce: native product (headroom {headroom:.4g}); "
                         "no reaction is forced to change.")
            return []
        demand = wt_max + 0.90 * headroom

        proposals: list[Proposal] = []
        with disc_model:
            constr_wt = StrainConstraints([biomass], [0.95 * mu_wt], ["b"])
            constr_mt = StrainConstraints(
                [biomass, target], [min_growth_frac * mu_wt, demand], ["l", "l"])
            fva = fva_optforce(disc_model, constr_wt, constr_mt,
                               objective_rxn=biomass, processes=1)

            constr_opt = ConstrOpt([substrate], [disc_model.reactions.get_by_id(substrate).lower_bound]) \
                if substrate else None
            constr_mt_ms = ConstrOpt([target], [demand])
            must_u = find_must_u(disc_model, fva.min_fluxes_wt, fva.max_fluxes_wt,
                                 constr_opt=constr_opt, constr_mt=constr_mt_ms,
                                 epsilon=epsilon_must)
            must_l = find_must_l(disc_model, fva.min_fluxes_wt, fva.max_fluxes_wt,
                                 constr_opt=constr_opt, constr_mt=constr_mt_ms,
                                 epsilon=epsilon_must)

            if not (must_u.must_u_set or must_l.must_l_set):
                notes.append("optforce: MUST sets empty; no reaction is forced to change.")
                return []

            cardinalities = [1] + ([2] if enable_k2 else [])
            for k in cardinalities:
                res = run_optforce(
                    disc_model, target, biomass,
                    must_u.must_u_set, must_l.must_l_set,
                    fva.min_fluxes_wt, fva.max_fluxes_wt,
                    fva.min_fluxes_mt, fva.max_fluxes_mt,
                    k=k, n_sets=50,
                    constr_opt=(ConstrOpt([substrate], [disc_model.reactions.get_by_id(substrate).lower_bound])
                                if substrate else None),
                    target_production=eps_prod,
                    min_growth=min_growth,
                )
                for s in res.opt_force_sets:
                    actions = tuple(
                        ReactionAction(r, t, float(f))
                        for r, t, f in zip(s.reactions, s.type_reg, s.fluxes)
                    )
                    proposals.append(Proposal(
                        actions=actions, source=f"optforce_k{k}",
                        discovery_floor=s.min_target_flux,
                        discovery_growth=s.max_growth_rate,
                    ))
        return proposals


register(OptForceStrategy())
