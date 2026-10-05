"""
optforce.py - Main OptForce algorithm (Step 3: the FORCE problem)

Given the MUST sets (from find_must_u / find_must_l) and the FVA flux ranges
for both strains, this searches for minimal sets of k interventions
(upregulations, downregulations, knockouts) that *force* overproduction of a
target metabolite.

This follows the structure of the OptForce FORCE problem (Ranganathan et al.,
2010, PLOS Comp Bio 6(4): e1000744; optForce.m in the COBRA Toolbox):

  * Candidate interventions are drawn **from the MUST sets** -- the whole
    reason those sets are computed -- not re-derived from raw FVA.
  * Each intervention is applied as a **one-sided regulation bound**, matching
    both the biology and the MATLAB primal constraints (optForce.m:908-925):

        upregulation   -> v_j >= minFluxesM[j]   (`primal3`, raises the lower bound)
        downregulation -> v_j <= maxFluxesM[j]   (`primal4`, lowers the upper bound)
        knockout       -> v_j == 0               (`primal5`/`primal6`)

    The cell is otherwise free to distribute flux.
  * A set is accepted only if it **guarantees** target production: after
    fixing the interventions, the minimum attainable target flux still meets
    the required level. This max-min guarantee is the heart of OptForce.

Deviations from optForce.m
--------------------------
These are deliberate. Read them before comparing numbers against MATLAB.

1. **Inner problem.** MATLAB's FORCE MILP (`dual1` at optForce.m:938, outer
   objective at :951, `osense = -1` at :1105) has inner objective ``c =
   e_target`` and *nothing else* -- there is no growth maximisation anywhere
   in it. Biomass enters only if the caller pins it through ``constrOpt``.
   MATLAB therefore computes::

       max_{yu,yl,y0}  min_v { v_target : Sv=0, intervened bounds, v_ic=b_ic }

   By default this module instead evaluates the guarantee *at maximum growth*
   (maximise biomass, pin it, then minimise the target), which is the stricter
   and more biologically meaningful growth-coupled guarantee. Set
   ``growth_mode="fixed"`` to reproduce the MATLAB inner problem exactly for
   cross-validation.

2. **Search strategy.** MATLAB solves one MILP repeatedly, adding a no-good
   cut per solution found (optForce.m:1072-1080). This module enumerates all
   C(n_candidates, k) combinations, which is exhaustive rather than
   solver-order dependent, and ranks by the guaranteed production. It avoids
   the bigM = 1000 linearisation (optForce.m:788) entirely.

3. **Knockout candidates.** MATLAB leaves ``y0(j)`` unconstrained for *every*
   reaction (only ``yu``/``yl`` are restricted to MustU/MustL at
   optForce.m:866-867), so it may knock out reactions no MUST-set analysis
   implicated. This module restricts knockouts to MustL, which is closer to
   the paper's intent. Within MustL a reaction may be offered as both a
   downregulation and a knockout, mirroring MATLAB's ``yu+yl+y0 <= 1``
   (optForce.m:961) rather than pre-committing to one type.
"""

from __future__ import annotations

import itertools
import warnings
import numpy as np
from dataclasses import dataclass, field
from typing import Optional
import cobra


# Default numerical tolerance, in flux units.
DEFAULT_EPSILON = 1e-3

# Tolerance used when pinning growth at its optimum before probing the target.
DEFAULT_TOL = 1e-7


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ExcludedRxns:
    """
    Reactions to exclude from specific intervention types.

    Attributes
    ----------
    rxn_list  : list of str  Reaction IDs.
    type_reg  : list of str  One of 'U' (upregulation), 'D' (downregulation),
                             'K' (knockout), or 'A' (all three) per reaction.
                             'A' mirrors the `typeReg` option handled at
                             optForce.m:884-888.

    Example
    -------
    >>> excl = ExcludedRxns(
    ...     rxn_list=['SUCt', 'R68_b'],
    ...     type_reg=['U', 'D']
    ... )
    """
    rxn_list: list[str] = field(default_factory=list)
    type_reg: list[str] = field(default_factory=list)

    def __post_init__(self):
        assert len(self.rxn_list) == len(self.type_reg)
        for t in self.type_reg:
            assert t in ('U', 'D', 'K', 'A'), (
                f"type_reg must be U, D, K, or A. Got '{t}'."
            )


@dataclass
class OptForceSet:
    """
    A single OptForce intervention set.

    Attributes
    ----------
    reactions       : list of str   Reaction IDs in the set.
    type_reg        : list of str   'upregulation', 'downregulation', or 'knockout'.
    fluxes          : list of float Forced flux level for each reaction.
    obj_value       : float         Guaranteed (minimum) target flux.
    min_target_flux : float         Min production of target after interventions.
    max_target_flux : float         Max production of target after interventions.
    max_growth_rate : float         Max growth rate after interventions.
    """
    reactions: list[str]
    type_reg: list[str]
    fluxes: list[float]
    obj_value: float
    min_target_flux: float
    max_target_flux: float
    max_growth_rate: float


@dataclass
class OptForceResult:
    """Output of run_optforce()."""
    opt_force_sets: list[OptForceSet]          # list of intervention sets found
    n_sets_found: int
    # Cardinality actually used. With k_min set, this is the smallest k that
    # produced any accepted set (cf. `currentK` at optForce.m:380).
    k_used: int = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_candidates(
    must_u: list[str],
    must_l: list[str],
    rxn_index: dict[str, int],
    min_fluxes_mt: np.ndarray,
    max_fluxes_mt: np.ndarray,
    target_rxn: str,
    biomass_rxn: str,
    excl_u: set[str],
    excl_d: set[str],
    excl_k: set[str],
    epsilon: float,
) -> list[tuple[str, str, float]]:
    """
    Turn the MUST sets into a candidate-intervention pool.

    Each candidate is (rxn_id, intervention_type, forced_flux), where the
    forced flux comes from the mutant (overproducing) FVA range and matches
    the MATLAB primal constraints:

      * MustU reaction -> 'upregulation', forced to its mutant *minimum*
        (`primal3`, optForce.m:910: ``v_j >= minFluxesM(j)``).
      * MustL reaction -> 'downregulation', forced to its mutant *maximum*
        (`primal4`, optForce.m:916: ``v_j <= maxFluxesM(j)``).
      * MustL reaction that *can* carry zero flux in the mutant ->
        additionally offered as a 'knockout' (`primal5`/`primal6`,
        optForce.m:922,928: ``v_j == 0``).

    A MustL reaction whose mutant range collapses to zero yields only the
    knockout (the downregulation would be identical). One whose range merely
    touches zero yields both, so the search can choose -- pre-committing to
    one type would hide the better intervention.
    """
    candidates: list[tuple[str, str, float]] = []

    for rxn_id in must_u:
        if rxn_id in (target_rxn, biomass_rxn) or rxn_id in excl_u:
            continue
        j = rxn_index.get(rxn_id)
        if j is None:
            continue
        candidates.append((rxn_id, 'upregulation', float(min_fluxes_mt[j])))

    for rxn_id in must_l:
        if rxn_id in (target_rxn, biomass_rxn):
            continue
        j = rxn_index.get(rxn_id)
        if j is None:
            continue
        min_mt = float(min_fluxes_mt[j])
        max_mt = float(max_fluxes_mt[j])
        collapsed = abs(min_mt) < epsilon and abs(max_mt) < epsilon

        if not collapsed and rxn_id not in excl_d:
            candidates.append((rxn_id, 'downregulation', max_mt))
        # Zero is attainable in the mutant, so a knockout is a valid option
        # (and it is the only option when the range has collapsed to zero).
        if min_mt <= epsilon and max_mt >= -epsilon and rxn_id not in excl_k:
            candidates.append((rxn_id, 'knockout', 0.0))

    return candidates


def _apply_regulation(
    model: cobra.Model,
    rxn_id: str,
    type_reg: str,
    forced_flux: float,
) -> None:
    """
    Apply a single intervention as a one-sided regulation bound (in place).

    upregulation   -> lower_bound raised to the forced level (flux >= forced)
    downregulation -> upper_bound lowered to the forced level (flux <= forced)
    knockout       -> flux pinned to 0

    Only ever tightens. The MATLAB primal constraints (optForce.m:908-925)
    add a row on one side and leave the model's own bound on the other, so a
    forced level outside the model's bounds is not something OptForce can
    produce from self-consistent FVA input -- it means the flux arrays do not
    belong to this model, which is raised rather than silently accommodated.
    """
    rxn = model.reactions.get_by_id(rxn_id)
    if type_reg == 'upregulation':
        if forced_flux > rxn.upper_bound:
            raise ValueError(
                f"Cannot upregulate {rxn_id} to {forced_flux:g}: it exceeds the "
                f"reaction's upper bound ({rxn.upper_bound:g}). The mutant FVA "
                "arrays do not match this model."
            )
        rxn.lower_bound = max(rxn.lower_bound, forced_flux)
    elif type_reg == 'downregulation':
        if forced_flux < rxn.lower_bound:
            raise ValueError(
                f"Cannot downregulate {rxn_id} to {forced_flux:g}: it is below "
                f"the reaction's lower bound ({rxn.lower_bound:g}). The mutant "
                "FVA arrays do not match this model."
            )
        rxn.upper_bound = min(rxn.upper_bound, forced_flux)
    elif type_reg == 'knockout':
        rxn.lower_bound = 0.0
        rxn.upper_bound = 0.0
    else:
        raise ValueError(f"Unknown intervention type: {type_reg!r}")


def _guaranteed_production(
    model: cobra.Model,
    target_rxn: str,
    biomass_rxn: str,
    interventions: list[tuple[str, str, float]],
    growth_mode: str = "max_growth",
    min_growth: float = 0.0,
    tol: float = DEFAULT_TOL,
) -> tuple[float, float, float, bool]:
    """
    Apply a set of interventions and return the guaranteed production.

    Steps (the OptForce inner guarantee):
      1. Apply each intervention as a one-sided regulation bound.
      2. If growth_mode == "max_growth": maximise growth and hold it at its
         optimum. A set whose maximum growth is <= max(min_growth, tol) is
         rejected as lethal -- see the inline note. If growth_mode ==
         "fixed": skip this entirely -- growth is whatever `constr_opt`
         pinned it to, which is what optForce.m does (the tutorial pins
         biomass to 0, so zero-growth designs are expected there).
      3. Find the min/max target flux under those conditions.

    Returns (max_growth_rate, min_target, max_target, feasible). The *minimum*
    target is the value OptForce guarantees: even when the cell works against
    you, production cannot fall below it.
    """
    with model:
        for rxn_id, type_reg, forced_flux in interventions:
            _apply_regulation(model, rxn_id, type_reg, forced_flux)

        biomass = model.reactions.get_by_id(biomass_rxn)

        if growth_mode == "max_growth":
            model.objective = biomass
            model.objective_direction = 'max'
            sol = model.optimize()
            if sol.status != 'optimal':
                return 0.0, 0.0, 0.0, False
            max_growth = sol.objective_value
            # A lethal intervention set has no "maximum growth" to
            # evaluate the guarantee at. Without this check the pin below
            # becomes `biomass >= -tol`, i.e. NO growth requirement, so a
            # dead strain routes all carbon to the target and outranks
            # every viable design. Reject it instead.
            if max_growth <= max(min_growth, tol):
                return max_growth, 0.0, 0.0, False
            # Hold growth at its optimum, then probe the target range.
            biomass.lower_bound = max_growth - tol
        else:
            # MATLAB-equivalent inner problem: no growth optimisation at all.
            model.objective = biomass
            model.objective_direction = 'max'
            sol = model.optimize()
            if sol.status != 'optimal':
                return 0.0, 0.0, 0.0, False
            max_growth = sol.objective_value

        model.objective = model.reactions.get_by_id(target_rxn)
        model.objective_direction = 'min'
        s_min = model.optimize()
        if s_min.status != 'optimal':
            return max_growth, 0.0, 0.0, False
        min_target = s_min.objective_value

        model.objective_direction = 'max'
        s_max = model.optimize()
        max_target = s_max.objective_value if s_max.status == 'optimal' else min_target

    return max_growth, min_target, max_target, True


# ---------------------------------------------------------------------------
# Main OptForce function
# ---------------------------------------------------------------------------

def run_optforce(
    model: cobra.Model,
    target_rxn: str,
    biomass_rxn: str,
    must_u: list[str],
    must_l: list[str],
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    min_fluxes_mt: np.ndarray,
    max_fluxes_mt: np.ndarray,
    k: int = 1,
    n_sets: int = 1,
    constr_opt=None,
    excluded_rxns: Optional[ExcludedRxns] = None,
    target_production: Optional[float] = None,
    min_growth: float = 0.0,
    k_min: Optional[int] = None,
    growth_mode: str = "max_growth",
    epsilon: float = DEFAULT_EPSILON,
    print_level: int = 0,
) -> OptForceResult:
    """
    Run the OptForce FORCE step to find intervention sets for overproduction.

    Candidate interventions are taken from the MUST sets (`must_u`, `must_l`).
    For every combination of `k` candidates, the interventions are applied as
    one-sided regulation bounds and the set is accepted only if it *guarantees*
    target production -- i.e. the minimum target flux meets `target_production`.
    Accepted sets are ranked by that guaranteed flux, best first.

    Parameters
    ----------
    model : cobra.Model
    target_rxn : str            Reaction ID for the target metabolite.
    biomass_rxn : str           Reaction ID for the biomass/growth objective.
    must_u : list of str        MustU set (upregulation candidates).
    must_l : list of str        MustL set (down-regulation / knockout candidates).
    min_fluxes_wt, max_fluxes_wt : np.ndarray  WT FVA flux ranges.
    min_fluxes_mt, max_fluxes_mt : np.ndarray  Mutant FVA flux ranges.
    k : int                     Number of interventions per set (default 1).
    n_sets : int                Maximum number of sets to return (default 1).
    constr_opt : ConstrOpt, optional   Fixed bounds applied before the search
                                       (e.g. substrate uptake).
    excluded_rxns : ExcludedRxns, optional  Reactions barred from a given
                                            intervention type.
    target_production : float, optional
        Required minimum guaranteed target flux for a set to be accepted.
        Defaults to `epsilon`, i.e. a set must force *some* strictly positive
        production; sets that guarantee nothing are rejected. Pass a larger
        value to demand a minimum yield.
    min_growth : float
        Minimum acceptable growth rate (default 0.0). Under
        growth_mode='max_growth' a set is rejected when its maximum growth
        falls to or below max(min_growth, 1e-7), because the guarantee is
        evaluated AT maximum growth and a lethal set would otherwise be
        scored with no growth requirement at all -- ranking dead strains
        above viable ones. Raise it to demand a genuinely viable strain
        (e.g. 10%% of wild-type growth).
    k_min : int, optional
        If given, search cardinalities k_min, k_min+1, ... k and return the
        sets from the *first* cardinality that yields any solution -- the
        minimal intervention set. Mirrors the `kMin`/`loop` behaviour at
        optForce.m:277,380, which is likewise opt-in. Default None searches
        exactly `k`.
    growth_mode : {'max_growth', 'fixed'}
        'max_growth' (default) evaluates the guarantee at maximum growth.
        'fixed' skips growth maximisation and minimises the target under
        `constr_opt` alone, reproducing the MATLAB inner problem. See the
        module docstring.
    epsilon : float             Numerical tolerance (default 1e-3).
    print_level : int           0 = silent, 1 = verbose.

    Returns
    -------
    OptForceResult

    Example
    -------
    >>> result = run_optforce(
    ...     model, 'EX_succ_e', 'Biomass_Ecoli_core',
    ...     must_u=must_u_combined, must_l=must_l_combined,
    ...     min_fluxes_wt=min_wt, max_fluxes_wt=max_wt,
    ...     min_fluxes_mt=min_mt, max_fluxes_mt=max_mt,
    ...     k=2, n_sets=5,
    ... )
    >>> for s in result.opt_force_sets:
    ...     print(s.reactions, s.type_reg, s.min_target_flux)
    """
    if growth_mode not in ("max_growth", "fixed"):
        raise ValueError(
            f"growth_mode must be 'max_growth' or 'fixed', got {growth_mode!r}"
        )

    rxn_ids = [r.id for r in model.reactions]
    rxn_index = {r: i for i, r in enumerate(rxn_ids)}

    min_fluxes_mt = np.asarray(min_fluxes_mt, dtype=float)
    max_fluxes_mt = np.asarray(max_fluxes_mt, dtype=float)
    if not (np.all(np.isfinite(min_fluxes_mt)) and np.all(np.isfinite(max_fluxes_mt))):
        raise ValueError(
            "Mutant FVA ranges contain non-finite values. OptForce requires "
            "bounded flux ranges (MATLAB caps everything at +/-1000)."
        )

    # Per-type exclusion sets. 'A' excludes a reaction from all three types,
    # matching the typeReg handling at optForce.m:884-888.
    excl_u: set[str] = set()
    excl_d: set[str] = set()
    excl_k: set[str] = set()
    if excluded_rxns:
        for rxn_id, t in zip(excluded_rxns.rxn_list, excluded_rxns.type_reg):
            if t == 'A':
                excl_u.add(rxn_id)
                excl_d.add(rxn_id)
                excl_k.add(rxn_id)
            else:
                {'U': excl_u, 'D': excl_d, 'K': excl_k}[t].add(rxn_id)

    # A set must force strictly positive production. MATLAB enforces this
    # structurally: lb_bl(z) = 0 (optForce.m:872) and the outer problem
    # maximises z = v_target, so a zero-production set is never a useful
    # solution. Note the threshold is used directly -- subtracting `epsilon`
    # from it here would cancel the default and admit zero-production sets.
    if target_production is None or target_production <= 0:
        target_production = epsilon
    accept_threshold = max(target_production, epsilon)

    candidates = _build_candidates(
        must_u, must_l, rxn_index,
        min_fluxes_mt, max_fluxes_mt,
        target_rxn, biomass_rxn,
        excl_u, excl_d, excl_k,
        epsilon,
    )

    if not candidates:
        warnings.warn(
            "No candidate interventions: the MUST sets are empty or fully "
            "excluded. OptForce has nothing to search.",
            RuntimeWarning,
            stacklevel=2,
        )
        return OptForceResult(opt_force_sets=[], n_sets_found=0, k_used=0)

    k_start = k if k_min is None else k_min
    if k_start < 1:
        raise ValueError(f"k_min must be >= 1, got {k_start}")
    if k_start > k:
        raise ValueError(f"k_min ({k_start}) must not exceed k ({k})")

    # Distinct reactions available bound the achievable cardinality: a set
    # may not intervene twice on the same reaction (optForce.m:961).
    n_distinct = len({c[0] for c in candidates})
    if n_distinct < k_start:
        warnings.warn(
            f"Only {n_distinct} distinct candidate reaction(s) available but "
            f"k={k_start} interventions requested; no set can be formed.",
            RuntimeWarning,
            stacklevel=2,
        )
        return OptForceResult(opt_force_sets=[], n_sets_found=0, k_used=0)

    if print_level:
        print(f"[OptForce] {len(candidates)} candidate interventions from MUST "
              f"sets ({n_distinct} distinct reactions).")
        print(f"[OptForce] Required guaranteed target production >= "
              f"{accept_threshold:.4f}; growth_mode={growth_mode}.")

    found: list[OptForceSet] = []
    k_used = 0

    with model:
        if constr_opt:
            for rxn_id, val in zip(constr_opt.rxn_list, constr_opt.values):
                r = model.reactions.get_by_id(rxn_id)
                r.lower_bound = val
                r.upper_bound = val

        for current_k in range(k_start, k + 1):
            if current_k > n_distinct:
                break
            if print_level:
                print(f"[OptForce] searching sets of size k={current_k}...")

            for combo in itertools.combinations(candidates, current_k):
                # A reaction may carry at most one intervention type, so a
                # combination that names the same reaction twice (e.g. both
                # the downregulation and the knockout of one MustL reaction)
                # is not a valid set. Cf. `outer2` at optForce.m:961.
                if len({c[0] for c in combo}) != current_k:
                    continue

                max_growth, min_target, max_target, feasible = _guaranteed_production(
                    model, target_rxn, biomass_rxn, list(combo),
                    growth_mode=growth_mode, min_growth=min_growth,
                )

                if not feasible:
                    if print_level:
                        print(f"[OptForce] set {[c[0] for c in combo]}: "
                              "infeasible after interventions, skipped.")
                    continue
                if max_growth < min_growth - epsilon:
                    continue
                if min_target < accept_threshold:
                    continue

                found.append(OptForceSet(
                    reactions=[c[0] for c in combo],
                    type_reg=[c[1] for c in combo],
                    fluxes=[c[2] for c in combo],
                    obj_value=min_target,      # guaranteed production = objective
                    min_target_flux=min_target,
                    max_target_flux=max_target,
                    max_growth_rate=max_growth,
                ))
                if print_level:
                    print(f"[OptForce] accepted {[c[0] for c in combo]}: "
                          f"growth={max_growth:.4f}, "
                          f"guaranteed target={min_target:.4f}")

            if found:
                # Smallest cardinality that yields a solution wins, and the
                # search stops there (cf. optForce.m:380,451).
                k_used = current_k
                break

    # Rank by the strongest guaranteed production and keep the best n_sets.
    found.sort(key=lambda s: s.obj_value, reverse=True)
    opt_force_sets = found[:n_sets]

    if print_level:
        print(f"[OptForce] {len(opt_force_sets)} set(s) returned "
              f"(of {len(found)} feasible at k={k_used}).")

    return OptForceResult(
        opt_force_sets=opt_force_sets,
        n_sets_found=len(opt_force_sets),
        k_used=k_used,
    )
