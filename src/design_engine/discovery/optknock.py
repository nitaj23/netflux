"""
optknock.py - OptKnock discovery strategy.

A greedy bilevel knockout search (k = 1..3). k=1 exhaustively screens every
gene-associated internal reaction on GLPK (cobra); k>=2 extends the promising
singles with a greedy beam on the persistent HiGHS kernel, falling back to cobra
if ``highspy`` is unavailable. Each candidate is scored by its guaranteed floor
-- the OptKnock objective -- but that discovery score is only a filter here;
the authoritative numbers come from the shared verification core downstream.

Ported from the former ``optknock.pipeline`` discovery stage; the beam itself is
unique to OptKnock (OptForce has no equivalent), so it lives here rather than
being reused from another package.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
import cobra

from ..schema import Proposal, ReactionAction
from ..lpkernel import make_lp
from .base import Baseline, DiscoveryParams, register


# Search bounds (module-level, not public knobs: the meaningful knob is
# ``max_knockouts``).
_MAX_PROPOSALS = 50        # per cardinality, handed to verification
_BEAM = 8                  # seeds carried into each pair/triple scan
_PARTNER_CAP = 150         # partner knockouts scanned per seed (by |WT flux|)
_DISCOVERY_BUDGET_S = 120  # per-target wall cap on the k>=2 beam (anti-stall)
_MARGIN_FRAC = 0.02        # a deeper knockout set must beat the shallower best by this


def _knockout_candidates(model: cobra.Model, target: str, biomass: str,
                         substrate: Optional[str]) -> list[str]:
    skip = {target, biomass} | ({substrate} if substrate else set())
    boundary = {r.id for r in model.boundary}
    return [r.id for r in model.reactions
            if r.id not in skip and r.id not in boundary
            and (r.gene_reaction_rule or "").strip()
            and not (r.lower_bound == 0.0 and r.upper_bound == 0.0)]


def _ko_floor(lp, target: str, biomass: str, ids, min_growth: float,
              gamma: float) -> tuple[float, float]:
    with lp.scope():
        for rid in ids:
            lp.set_bound(rid, lb=0.0, ub=0.0)
        mu, s_mu = lp.maximize(biomass)
        if s_mu != "optimal" or not np.isfinite(mu) or mu <= min_growth:
            return 0.0, (float(mu) if np.isfinite(mu) else 0.0)
        lp.set_bound(biomass, lb=(1.0 - gamma) * mu)
        f, s_f = lp.minimize(target)
    return (float(f) if s_f == "optimal" and np.isfinite(f) else 0.0), float(mu)


def _beam_scan(lp, target, biomass, seeds, partners, base_floor, min_growth,
               gamma, eps_prod, budget_s, t0, notes, tag):
    margin = max(_MARGIN_FRAC * abs(base_floor), 10 * eps_prod)
    out, seen = [], set()
    for seed in seeds:
        if time.perf_counter() - t0 > budget_s:
            notes.append(f"optknock: {tag} budget hit ({budget_s:g}s); partial.")
            break
        for p in partners:
            if p in seed:
                continue
            ids = seed + (p,)
            key = frozenset(ids)
            if key in seen:
                continue
            seen.add(key)
            f, mu = _ko_floor(lp, target, biomass, ids, min_growth, gamma)
            if f > base_floor + margin and f > eps_prod:
                out.append((f, mu, ids))
    out.sort(reverse=True)
    return out


class OptKnockStrategy:
    name = "optknock"

    def discover(self, disc_model: cobra.Model, target: str, biomass: str,
                 substrate: Optional[str], baseline: Baseline,
                 params: DiscoveryParams, notes: list[str]) -> list[Proposal]:
        gamma = params.gamma
        eps_prod = params.eps_prod
        max_k = max(1, min(3, int(params.max_knockouts)))
        cands = _knockout_candidates(disc_model, target, biomass, substrate)

        # ---- k=1 on GLPK (cobra reference path) ----
        lp1 = make_lp(disc_model, "cobra")
        mu_wt, s = lp1.maximize(biomass)
        if s != "optimal" or not np.isfinite(mu_wt) or mu_wt <= 0:
            notes.append("optknock: discovery model cannot grow; no candidates.")
            return []
        if params.min_growth_abs is not None:
            min_growth = params.min_growth_abs
        else:
            min_growth = params.min_growth_frac * mu_wt

        singles = []
        for rid in cands:
            f, mu = _ko_floor(lp1, target, biomass, (rid,), min_growth, gamma)
            singles.append((f, mu, rid))
        singles.sort(reverse=True)
        best1 = singles[0][0] if singles else 0.0

        props = [Proposal((ReactionAction(rid, "knockout", 0.0),), "optknock_k1", f, mu)
                 for f, mu, rid in singles[:_MAX_PROPOSALS] if f > eps_prod]

        if max_k < 2:
            if not props:
                notes.append("optknock: no single knockout couples production to growth.")
            return props

        # ---- k>=2 greedy beam on HiGHS ----
        try:
            lp = make_lp(disc_model, "highs")
        except ImportError:
            lp = lp1
            notes.append("optknock: highspy unavailable; k>=2 scan on cobra.")

        t0 = time.perf_counter()
        viable = [(rid, mu) for f, mu, rid in singles if mu > min_growth]
        _, sp, fluxes = lp.solve_primal(biomass, "max")
        wtflux = {k: abs(v) for k, v in fluxes.items()} if sp == "optimal" else {}
        by_flux = [rid for rid, _ in sorted(viable, key=lambda t: wtflux.get(t[0], 0.0), reverse=True)]
        partners = by_flux[:_PARTNER_CAP]
        seen, seed_singles = set(), []
        for rid in [r for _, _, r in singles[:_BEAM]] + by_flux[:_BEAM]:
            if rid not in seen:
                seen.add(rid); seed_singles.append(rid)

        pairs = _beam_scan(lp, target, biomass, [(s,) for s in seed_singles], partners,
                           best1, min_growth, gamma, eps_prod, _DISCOVERY_BUDGET_S, t0,
                           notes, "k2")
        props += [Proposal(tuple(ReactionAction(r, "knockout", 0.0) for r in ids),
                           "optknock_k2", f, mu) for f, mu, ids in pairs[:_MAX_PROPOSALS]]
        best2 = pairs[0][0] if pairs else best1

        if max_k >= 3:
            triples = _beam_scan(lp, target, biomass, [ids for _, _, ids in pairs[:_BEAM]],
                                 partners, best2, min_growth, gamma, eps_prod,
                                 _DISCOVERY_BUDGET_S, t0, notes, "k3")
            props += [Proposal(tuple(ReactionAction(r, "knockout", 0.0) for r in ids),
                               "optknock_k3", f, mu) for f, mu, ids in triples[:_MAX_PROPOSALS]]

        if not props:
            notes.append("optknock: no knockout set couples production to growth.")
        return props


register(OptKnockStrategy())
