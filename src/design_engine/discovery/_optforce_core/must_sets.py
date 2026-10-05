"""
must_sets.py - Find MUST sets for OptForce

Implements findMustL, findMustU, findMustLL, findMustUU, findMustUL.

Each function identifies reactions that MUST be up- or down-regulated in the
mutant strain to achieve the target overproduction.

Translated from the COBRA Toolbox MATLAB implementations
(opencobra/cobratoolbox, src/design/optForce).

Reference: Ranganathan et al. (2010) OptForce, PLOS Comp Bio 6(4): e1000744.

Relationship to the MATLAB bilevel MILPs
----------------------------------------
The MATLAB versions solve a bilevel MILP in which a binary y(j) selects one
reaction, strong duality forces the inner LP to optimality, and the outer
objective z is thresholded. Decoding those formulations, each collapses to a
single LP test per reaction (or per pair):

  findMustL.m  (outer obj line 496, z >= 0.5 line 539)
      j in MustL  <=>  max v_j  <  minFluxesW(j) - eps
  findMustU.m  (outer obj line 491, z >= 0.5 line 534)
      j in MustU  <=>  min v_j  >  maxFluxesW(j) + eps
  findMustLL.m (outer obj line 559, z >= 0.1 line 626)
      (i,j) in MustLL <=> max(v_i + v_j) < minW_i + minW_j - eps
  findMustUU.m (outer obj line 557, z >= 0.1 line 627)
      (i,j) in MustUU <=> min(v_i + v_j) > maxW_i + maxW_j + eps
  findMustUL.m (outer obj line 557, z >= 0.1 line 624)
      (i,j) in MustUL <=> min(v_i - v_j) > maxW_i - minW_j + eps

Enumerating those LPs directly is equivalent to the MATLAB MILP loop (which
re-solves, excluding previously found reactions, until infeasible) and avoids
the bigM=1000 linearisation entirely. Note that the pairwise search is
O(n^2) LPs and becomes intractable on genome-scale models.

Two candidate filters from the MATLAB source are reproduced here and must not
be dropped -- without them the MUST sets fill up with artifacts:

  * the `can` set (findMustL.m:268-270): reactions blocked in the wild type
    (minFluxesW == maxFluxesW == 0) are never eligible.
  * the constrained set (`sel_nic` in findMustL.m:419, `sel_c_nic` in
    findMustLL.m:453): reactions pinned by constrOpt/constrMT are never
    eligible, since a pinned reaction is trivially "forced" to change.
"""

from __future__ import annotations

import warnings
import numpy as np
from dataclasses import dataclass, field
from typing import Optional
import cobra


# Strict-inequality margins, in flux units. These match the thresholds
# hardcoded into the MATLAB outer objectives (z >= 0.5 for first order,
# z >= 0.1 for second order). They are NOT solver tolerances -- they express
# "the forced change must be large enough to matter". Lowering them towards
# zero admits reactions whose forced change is numerical noise.
DEFAULT_EPSILON_FIRST_ORDER = 0.5
DEFAULT_EPSILON_SECOND_ORDER = 0.1

# Tolerance for deciding whether a wild-type FVA range is identically zero
# when building the `can` set. MATLAB compares against exact zero; a small
# tolerance is used here because the ranges come back as floats.
DEFAULT_CAN_TOL = 1e-9


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ConstrOpt:
    """
    Fixed-flux constraints (both bounds equal) for MUST set problems.

    Only include reactions whose flux is pinned (lb == ub).
    Variable-bound reactions should be set directly in model.lb / model.ub.

    Reactions listed here are automatically excluded from the MUST sets,
    matching `sel_nic` in findMustL.m:419.

    Attributes
    ----------
    rxn_list : list of str
    values   : list of float
    """
    rxn_list: list[str] = field(default_factory=list)
    values: list[float] = field(default_factory=list)

    def __post_init__(self):
        assert len(self.rxn_list) == len(self.values)


@dataclass
class MustLResult:
    """Output of find_must_l()."""
    must_l_set: list[str]       # reaction IDs in MustL
    pos_must_l: list[int]       # positions in model.reactions
    # Probed extreme flux per MustL reaction (MATLAB's vmax in findMustL.m).
    flux_values: list[float] = field(default_factory=list)


@dataclass
class MustUResult:
    """Output of find_must_u()."""
    must_u_set: list[str]
    pos_must_u: list[int]
    # Probed extreme flux per MustU reaction (MATLAB's vmin in findMustU.m).
    flux_values: list[float] = field(default_factory=list)


@dataclass
class MustLLResult:
    """Output of find_must_ll()."""
    must_ll: list[tuple[str, str]]          # pairs of reaction IDs
    pos_must_ll: list[tuple[int, int]]      # pairs of positions
    must_ll_linear: list[str]               # unique reaction IDs
    pos_must_ll_linear: list[int]


@dataclass
class MustUUResult:
    """Output of find_must_uu()."""
    must_uu: list[tuple[str, str]]
    pos_must_uu: list[tuple[int, int]]
    must_uu_linear: list[str]
    pos_must_uu_linear: list[int]


@dataclass
class MustULResult:
    """Output of find_must_ul()."""
    must_ul: list[tuple[str, str]]
    pos_must_ul: list[tuple[int, int]]
    must_ul_linear: list[str]
    pos_must_ul_linear: list[int]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_rxn_index(model: cobra.Model) -> dict[str, int]:
    return {r.id: i for i, r in enumerate(model.reactions)}


def _apply_constr_opt(model: cobra.Model, constr_opt: Optional[ConstrOpt]):
    """Pin reaction bounds according to ConstrOpt (modifies model in place)."""
    if constr_opt is None:
        return
    for rxn_id, val in zip(constr_opt.rxn_list, constr_opt.values):
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = val
        rxn.upper_bound = val


def _constrained_ids(*constr_opts: Optional[ConstrOpt]) -> set[str]:
    """
    Collect the IDs of every reaction pinned by a ConstrOpt.

    MATLAB structurally forbids these from entering a MUST set (their binary
    y(j) is multiplied by `sel_nic` / `sel_c_nic`). A pinned reaction is
    trivially "forced": biomass pinned to 0 in the mutant always satisfies the
    MustL test against a nonzero wild-type minimum, so including them reports
    the caller's own constraints back as biology.
    """
    ids: set[str] = set()
    for c in constr_opts:
        if c is not None:
            ids.update(c.rxn_list)
    return ids


def _can_set(
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    tol: float = DEFAULT_CAN_TOL,
) -> np.ndarray:
    """
    Boolean mask of reactions eligible for a MUST set.

    Reproduces findMustL.m:268-270::

        can = zeros(n_rxns, 1);
        can(minFluxesW ~= 0) = 1;
        can(maxFluxesW ~= 0) = 1;

    Reactions blocked in the wild type (both FVA bounds zero) are excluded.
    On genome-scale models this removes a large fraction of reactions.
    """
    return (np.abs(min_fluxes_wt) > tol) | (np.abs(max_fluxes_wt) > tol)


def _check_finite(min_fluxes_wt: np.ndarray, max_fluxes_wt: np.ndarray) -> None:
    """
    Reject non-finite FVA ranges.

    MATLAB caps everything at H = bigM = 1000, so an unbounded range never
    reaches the comparisons. In Python an infinite bound would silently make
    every threshold comparison False (or NaN), yielding empty MUST sets with
    no error.
    """
    finite = np.isfinite(min_fluxes_wt) & np.isfinite(max_fluxes_wt)
    if not np.all(finite):
        bad = np.where(~finite)[0]
        raise ValueError(
            "Wild-type FVA ranges contain non-finite values at reaction "
            f"indices {bad.tolist()[:10]}{' ...' if len(bad) > 10 else ''}. "
            "OptForce requires bounded flux ranges; set finite bounds on the "
            "model (the MATLAB implementation caps everything at +/-1000)."
        )


def _warn_infeasible(n_skipped: int, total: int, where: str) -> None:
    if n_skipped:
        warnings.warn(
            f"{where}: {n_skipped} of {total} subproblems were not solved to "
            "optimality (infeasible or unbounded) and were skipped. An empty "
            "or small MUST set may reflect an infeasible network rather than "
            "an absence of required changes.",
            RuntimeWarning,
            stacklevel=3,
        )


def _warn_empty(name: str, n: int) -> None:
    if n == 0:
        warnings.warn(
            f"{name} is empty: no reaction was found that must change. Check "
            "that the overproduction constraints actually force the target, "
            "and that `epsilon` is not set too high.",
            RuntimeWarning,
            stacklevel=3,
        )


# ---------------------------------------------------------------------------
# Public API - first order
# ---------------------------------------------------------------------------

def find_must_l(
    model: cobra.Model,
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    constr_opt: Optional[ConstrOpt] = None,
    constr_mt: Optional[ConstrOpt] = None,
    excluded_rxns: Optional[list[str]] = None,
    epsilon: float = DEFAULT_EPSILON_FIRST_ORDER,
    can_tol: float = DEFAULT_CAN_TOL,
    print_level: int = 0,
) -> MustLResult:
    """
    Find the first-order MustL set.

    A reaction is in MustL if in the mutant strain its flux MUST decrease
    below the wild-type minimum (i.e., it must be downregulated or knocked out).

    Implemented as: for each eligible reaction j, solve

        max  v_j                        <- MAXIMISE to find the ceiling
        s.t. Sv = 0, lb <= v <= ub,
             constr_opt (e.g. glucose pinned),
             constr_mt  (e.g. succinate overproduction requirement)

    If the ceiling v_j* < minFluxesWT[j] - epsilon, then even at its best
    the mutant flux stays below the WT minimum -> reaction j in MustL.

    This is the closed form of the bilevel MILP in findMustL.m (outer
    objective line 496 with the z >= 0.5 threshold at line 539).

    Eligibility follows MATLAB: a reaction is skipped if it is blocked in the
    wild type (`can`, findMustL.m:268-270) or pinned by constr_opt/constr_mt
    (`sel_nic`, findMustL.m:419).

    Parameters
    ----------
    model : cobra.Model
    min_fluxes_wt : np.ndarray  shape (n_rxns,)
    max_fluxes_wt : np.ndarray  shape (n_rxns,)
    constr_opt : ConstrOpt, optional   Fixed constraints (e.g. glucose uptake).
    constr_mt  : ConstrOpt, optional   Overproduction constraints (e.g. succinate).
    excluded_rxns : list of str, optional  Additional reactions to skip.
    epsilon : float
        Strict-inequality margin in flux units. Default 0.5, matching the
        `z >= 0.5` constraint hardcoded at findMustL.m:539.
    can_tol : float  Tolerance for the `can` (not-blocked-in-WT) test.
    print_level : int  0 = silent, 1 = verbose.

    Returns
    -------
    MustLResult
    """
    min_fluxes_wt = np.asarray(min_fluxes_wt, dtype=float)
    max_fluxes_wt = np.asarray(max_fluxes_wt, dtype=float)
    _check_finite(min_fluxes_wt, max_fluxes_wt)

    excluded = set(excluded_rxns or []) | _constrained_ids(constr_opt, constr_mt)
    can = _can_set(min_fluxes_wt, max_fluxes_wt, can_tol)

    must_l_set: list[str] = []
    pos_must_l: list[int] = []
    flux_values: list[float] = []
    n_skipped = 0
    n_tested = 0

    with model:
        # Apply both constraint sets to define the mutant network
        _apply_constr_opt(model, constr_opt)
        _apply_constr_opt(model, constr_mt)

        for i, rxn in enumerate(model.reactions):
            if rxn.id in excluded or not can[i]:
                continue
            n_tested += 1

            with model:
                # MAXIMISE v_j to find the highest possible flux in the mutant.
                # If even this ceiling is below WT min, the reaction MUST decrease.
                model.objective = rxn
                model.objective_direction = 'max'
                sol = model.optimize()

                if sol.status != 'optimal':
                    n_skipped += 1
                    continue

                v_max_mt = sol.objective_value
                # MustL condition: mutant ceiling < WT minimum
                if v_max_mt < min_fluxes_wt[i] - epsilon:
                    must_l_set.append(rxn.id)
                    pos_must_l.append(i)
                    flux_values.append(float(v_max_mt))
                    if print_level:
                        print(
                            f"[MustL] {rxn.id}: MT_max={v_max_mt:.4f} "
                            f"< WT_min={min_fluxes_wt[i]:.4f}"
                        )

    _warn_infeasible(n_skipped, n_tested, "find_must_l")
    _warn_empty("MustL set", len(must_l_set))
    return MustLResult(
        must_l_set=must_l_set, pos_must_l=pos_must_l, flux_values=flux_values
    )


def find_must_u(
    model: cobra.Model,
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    constr_opt: Optional[ConstrOpt] = None,
    constr_mt: Optional[ConstrOpt] = None,
    excluded_rxns: Optional[list[str]] = None,
    epsilon: float = DEFAULT_EPSILON_FIRST_ORDER,
    can_tol: float = DEFAULT_CAN_TOL,
    print_level: int = 0,
) -> MustUResult:
    """
    Find the first-order MustU set.

    A reaction is in MustU if in the mutant its flux MUST exceed the
    wild-type maximum (i.e., it must be upregulated).

    For each eligible reaction j, solve

        min  v_j                        <- MINIMISE to find the floor
        s.t. Sv = 0, lb <= v <= ub,
             constr_opt (e.g. glucose pinned),
             constr_mt  (e.g. succinate overproduction requirement)

    If the floor v_j* > maxFluxesWT[j] + epsilon, then even at its lowest
    the mutant flux exceeds the WT maximum -> reaction j in MustU.

    This is the closed form of the bilevel MILP in findMustU.m (outer
    objective line 491 with the z >= 0.5 threshold at line 534).

    Parameters
    ----------
    model : cobra.Model
    min_fluxes_wt : np.ndarray
    max_fluxes_wt : np.ndarray
    constr_opt : ConstrOpt, optional   Fixed constraints (e.g. glucose uptake).
    constr_mt  : ConstrOpt, optional   Overproduction constraints.
    excluded_rxns : list of str, optional
    epsilon : float
        Strict-inequality margin in flux units. Default 0.5, matching the
        `z >= 0.5` constraint hardcoded at findMustU.m:534.
    can_tol : float
    print_level : int

    Returns
    -------
    MustUResult
    """
    min_fluxes_wt = np.asarray(min_fluxes_wt, dtype=float)
    max_fluxes_wt = np.asarray(max_fluxes_wt, dtype=float)
    _check_finite(min_fluxes_wt, max_fluxes_wt)

    excluded = set(excluded_rxns or []) | _constrained_ids(constr_opt, constr_mt)
    can = _can_set(min_fluxes_wt, max_fluxes_wt, can_tol)

    must_u_set: list[str] = []
    pos_must_u: list[int] = []
    flux_values: list[float] = []
    n_skipped = 0
    n_tested = 0

    with model:
        # Apply both constraint sets to define the mutant network
        _apply_constr_opt(model, constr_opt)
        _apply_constr_opt(model, constr_mt)

        for i, rxn in enumerate(model.reactions):
            if rxn.id in excluded or not can[i]:
                continue
            n_tested += 1

            with model:
                # MINIMISE v_j to find the lowest possible flux in the mutant.
                # If even this floor exceeds WT max, the reaction MUST increase.
                model.objective = rxn
                model.objective_direction = 'min'
                sol = model.optimize()

                if sol.status != 'optimal':
                    n_skipped += 1
                    continue

                v_min_mt = sol.objective_value
                # MustU condition: mutant floor > WT maximum
                if v_min_mt > max_fluxes_wt[i] + epsilon:
                    must_u_set.append(rxn.id)
                    pos_must_u.append(i)
                    flux_values.append(float(v_min_mt))
                    if print_level:
                        print(
                            f"[MustU] {rxn.id}: MT_min={v_min_mt:.4f} "
                            f"> WT_max={max_fluxes_wt[i]:.4f}"
                        )

    _warn_infeasible(n_skipped, n_tested, "find_must_u")
    _warn_empty("MustU set", len(must_u_set))
    return MustUResult(
        must_u_set=must_u_set, pos_must_u=pos_must_u, flux_values=flux_values
    )


# ---------------------------------------------------------------------------
# Public API - second order
# ---------------------------------------------------------------------------

def _second_order_setup(
    model: cobra.Model,
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    constr_opt: Optional[ConstrOpt],
    constr_mt: Optional[ConstrOpt],
    excluded_rxns: Optional[list[str]],
    exclude_first_order: Optional[list[str]],
    can_tol: float,
):
    """
    Shared eligibility/indexing setup for the pairwise MUST finders.

    Returns (rxn_index, eligible_reactions). Eligibility mirrors
    `sel_c_nic = can & ~sel_ic` (findMustLL.m:453) plus the caller-supplied
    exclusions, which is how the MATLAB tutorial keeps first-order reactions
    out of the second-order search.
    """
    min_fluxes_wt = np.asarray(min_fluxes_wt, dtype=float)
    max_fluxes_wt = np.asarray(max_fluxes_wt, dtype=float)
    _check_finite(min_fluxes_wt, max_fluxes_wt)

    excluded = (
        set(excluded_rxns or [])
        | set(exclude_first_order or [])
        | _constrained_ids(constr_opt, constr_mt)
    )
    can = _can_set(min_fluxes_wt, max_fluxes_wt, can_tol)

    # Build the index map once, not once per pair.
    rxn_index = _get_rxn_index(model)
    eligible = [
        r for r in model.reactions
        if r.id not in excluded and can[rxn_index[r.id]]
    ]
    return rxn_index, eligible


def _warn_pairwise_size(n: int, ordered: bool = False) -> None:
    n_pairs = n * (n - 1) if ordered else n * (n - 1) // 2
    if n_pairs > 20000:
        warnings.warn(
            f"Pairwise MUST search over {n} reactions requires {n_pairs} LP "
            "solves. This enumeration is exhaustive but does not scale to "
            "genome-scale models; narrow the search with `excluded_rxns`.",
            RuntimeWarning,
            stacklevel=3,
        )


def _add_seen(seen: list[str], seen_set: set[str], *rxn_ids: str) -> None:
    """Append reaction IDs preserving first-seen order (stable output)."""
    for rid in rxn_ids:
        if rid not in seen_set:
            seen_set.add(rid)
            seen.append(rid)


def find_must_ll(
    model: cobra.Model,
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    constr_opt: Optional[ConstrOpt] = None,
    constr_mt: Optional[ConstrOpt] = None,
    excluded_rxns: Optional[list[str]] = None,
    exclude_first_order: Optional[list[str]] = None,
    epsilon: float = DEFAULT_EPSILON_SECOND_ORDER,
    can_tol: float = DEFAULT_CAN_TOL,
    print_level: int = 0,
) -> MustLLResult:
    """
    Find the second-order MustLL set.

    A pair (i, j) is in MustLL if both reactions must be downregulated
    *jointly* (their summed flux is forced below the basal minimum), even
    though neither one is necessarily forced down on its own.

    Anchoring on the first-order logic: a reaction "must decrease" when even
    its highest attainable flux stays below the basal minimum. The joint
    analogue therefore *maximises* v_i + v_j (the best case in which both
    stay high) and tests whether that ceiling is still below the basal-min
    sum::

        max  v_i + v_j
        s.t. Sv = 0, lb <= v <= ub, (constr_opt, constr_mt)

        (i, j) in MustLL  iff  max(v_i + v_j) < minWT[i] + minWT[j] - eps

    This matches the outer objective at findMustLL.m:559 (whose `dual2`
    y-coefficients are +Icnic, making the inner LP a maximisation) with the
    z >= 0.1 threshold at line 626.

    First-order MustL reactions should be excluded via `exclude_first_order`
    so that only genuinely second-order pairs are reported -- this is what
    the MATLAB tutorial passes as `excludedRxns`.

    Parameters
    ----------
    model : cobra.Model
    min_fluxes_wt : np.ndarray
    max_fluxes_wt : np.ndarray
    constr_opt : ConstrOpt, optional
    constr_mt : ConstrOpt, optional
        Overproduction constraints. Accepted here for symmetry with the
        first-order finders; merged with constr_opt before solving.
    excluded_rxns : list of str, optional
    exclude_first_order : list of str, optional
        First-order MustL reaction IDs to remove before the pairwise search.
    epsilon : float
        Default 0.1, matching `z >= 0.1` at findMustLL.m:626.
    can_tol : float
    print_level : int

    Returns
    -------
    MustLLResult
    """
    rxn_index, rxns = _second_order_setup(
        model, min_fluxes_wt, max_fluxes_wt,
        constr_opt, constr_mt, excluded_rxns, exclude_first_order, can_tol,
    )
    min_fluxes_wt = np.asarray(min_fluxes_wt, dtype=float)
    n = len(rxns)
    _warn_pairwise_size(n)

    must_ll: list[tuple[str, str]] = []
    pos_must_ll: list[tuple[int, int]] = []
    seen: list[str] = []
    seen_set: set[str] = set()

    with model:
        _apply_constr_opt(model, constr_opt)
        _apply_constr_opt(model, constr_mt)

        for a in range(n):
            for b in range(a + 1, n):
                rxn_a, rxn_b = rxns[a], rxns[b]
                idx_a, idx_b = rxn_index[rxn_a.id], rxn_index[rxn_b.id]

                with model:
                    # Maximise the joint flux: if even the ceiling of v_a + v_b
                    # is below the basal-min sum, the pair must decrease together.
                    model.objective = {rxn_a: 1.0, rxn_b: 1.0}
                    model.objective_direction = 'max'
                    sol = model.optimize()

                    if sol.status != 'optimal':
                        continue

                    va = sol.fluxes[rxn_a.id]
                    vb = sol.fluxes[rxn_b.id]
                    basal_min_sum = min_fluxes_wt[idx_a] + min_fluxes_wt[idx_b]
                    if (va + vb) < basal_min_sum - epsilon:
                        must_ll.append((rxn_a.id, rxn_b.id))
                        pos_must_ll.append((idx_a, idx_b))
                        _add_seen(seen, seen_set, rxn_a.id, rxn_b.id)
                        if print_level:
                            print(f"[MustLL] {rxn_a.id} + {rxn_b.id}")

    return MustLLResult(
        must_ll=must_ll,
        pos_must_ll=pos_must_ll,
        must_ll_linear=seen,
        pos_must_ll_linear=[rxn_index[r] for r in seen],
    )


def find_must_uu(
    model: cobra.Model,
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    constr_opt: Optional[ConstrOpt] = None,
    constr_mt: Optional[ConstrOpt] = None,
    excluded_rxns: Optional[list[str]] = None,
    exclude_first_order: Optional[list[str]] = None,
    epsilon: float = DEFAULT_EPSILON_SECOND_ORDER,
    can_tol: float = DEFAULT_CAN_TOL,
    print_level: int = 0,
) -> MustUUResult:
    """
    Find the second-order MustUU set.

    A pair (i, j) is in MustUU if both must be upregulated *jointly* (their
    summed flux is forced above the basal maximum), even though neither one
    is necessarily forced up on its own.

    Anchoring on the first-order logic: a reaction "must increase" when even
    its lowest attainable flux exceeds the basal maximum. The joint analogue
    therefore *minimises* v_i + v_j (the worst case in which both stay low)
    and tests whether that floor is still above the basal-max sum::

        min  v_i + v_j
        s.t. Sv = 0, lb <= v <= ub, (constr_opt, constr_mt)

        (i, j) in MustUU  iff  min(v_i + v_j) > maxWT[i] + maxWT[j] + eps

    This matches the outer objective at findMustUU.m:557 (whose `dual2`
    y-coefficients are -Icnic, making the inner LP a minimisation) with the
    z >= 0.1 threshold at line 627.

    Parameters
    ----------
    model : cobra.Model
    min_fluxes_wt : np.ndarray
    max_fluxes_wt : np.ndarray
    constr_opt : ConstrOpt, optional
    constr_mt : ConstrOpt, optional
    excluded_rxns : list of str, optional
    exclude_first_order : list of str, optional
        First-order MustU reaction IDs to remove before the pairwise search.
    epsilon : float
        Default 0.1, matching `z >= 0.1` at findMustUU.m:627.
    can_tol : float
    print_level : int

    Returns
    -------
    MustUUResult
    """
    rxn_index, rxns = _second_order_setup(
        model, min_fluxes_wt, max_fluxes_wt,
        constr_opt, constr_mt, excluded_rxns, exclude_first_order, can_tol,
    )
    max_fluxes_wt = np.asarray(max_fluxes_wt, dtype=float)
    n = len(rxns)
    _warn_pairwise_size(n)

    must_uu: list[tuple[str, str]] = []
    pos_must_uu: list[tuple[int, int]] = []
    seen: list[str] = []
    seen_set: set[str] = set()

    with model:
        _apply_constr_opt(model, constr_opt)
        _apply_constr_opt(model, constr_mt)

        for a in range(n):
            for b in range(a + 1, n):
                rxn_a, rxn_b = rxns[a], rxns[b]
                idx_a, idx_b = rxn_index[rxn_a.id], rxn_index[rxn_b.id]

                with model:
                    # Minimise the joint flux: if even the floor of v_a + v_b
                    # exceeds the basal-max sum, the pair must increase together.
                    model.objective = {rxn_a: 1.0, rxn_b: 1.0}
                    model.objective_direction = 'min'
                    sol = model.optimize()

                    if sol.status != 'optimal':
                        continue

                    va = sol.fluxes[rxn_a.id]
                    vb = sol.fluxes[rxn_b.id]
                    basal_max_sum = max_fluxes_wt[idx_a] + max_fluxes_wt[idx_b]
                    if (va + vb) > basal_max_sum + epsilon:
                        must_uu.append((rxn_a.id, rxn_b.id))
                        pos_must_uu.append((idx_a, idx_b))
                        _add_seen(seen, seen_set, rxn_a.id, rxn_b.id)
                        if print_level:
                            print(f"[MustUU] {rxn_a.id} + {rxn_b.id}")

    return MustUUResult(
        must_uu=must_uu,
        pos_must_uu=pos_must_uu,
        must_uu_linear=seen,
        pos_must_uu_linear=[rxn_index[r] for r in seen],
    )


def find_must_ul(
    model: cobra.Model,
    min_fluxes_wt: np.ndarray,
    max_fluxes_wt: np.ndarray,
    constr_opt: Optional[ConstrOpt] = None,
    constr_mt: Optional[ConstrOpt] = None,
    excluded_rxns: Optional[list[str]] = None,
    exclude_first_order: Optional[list[str]] = None,
    epsilon: float = DEFAULT_EPSILON_SECOND_ORDER,
    can_tol: float = DEFAULT_CAN_TOL,
    print_level: int = 0,
) -> MustULResult:
    """
    Find the second-order MustUL set.

    A pair (i, j) is in MustUL if reaction i must be upregulated AND
    reaction j must be downregulated *jointly* -- the gap v_i - v_j is forced
    above what the basal ranges allow, even though neither reaction need move
    on its own.

    Anchoring on the first-order logic, the worst case for "i up while j down"
    is when v_i is as low and v_j as high as possible, i.e. we *minimise*
    v_i - v_j. If even that floor exceeds the basal gap maxWT[i] - minWT[j],
    the pair must move in opposite directions together::

        min  v_i - v_j
        s.t. Sv = 0, lb <= v <= ub, (constr_opt, constr_mt)

        (i, j) in MustUL  iff  min(v_i - v_j) > maxWT[i] - minWT[j] + eps

    This matches the outer objective at findMustUL.m:557 (whose `dual2`
    y-coefficients are -y1 +y2, making the inner LP min(v_i - v_j)) with the
    z >= 0.1 threshold at line 624. Unlike LL/UU the pair is ordered, so both
    (i, j) and (j, i) are tested -- mirroring `y1(j) + y2(j) <= 1` at
    findMustUL.m:677, which only forbids i == j.

    Parameters
    ----------
    model : cobra.Model
    min_fluxes_wt : np.ndarray
    max_fluxes_wt : np.ndarray
    constr_opt : ConstrOpt, optional
    constr_mt : ConstrOpt, optional
    excluded_rxns : list of str, optional
    exclude_first_order : list of str, optional
        First-order MustU/MustL reaction IDs to remove before the search.
    epsilon : float
        Default 0.1, matching `z >= 0.1` at findMustUL.m:624.
    can_tol : float
    print_level : int

    Returns
    -------
    MustULResult
    """
    rxn_index, rxns = _second_order_setup(
        model, min_fluxes_wt, max_fluxes_wt,
        constr_opt, constr_mt, excluded_rxns, exclude_first_order, can_tol,
    )
    min_fluxes_wt = np.asarray(min_fluxes_wt, dtype=float)
    max_fluxes_wt = np.asarray(max_fluxes_wt, dtype=float)
    n = len(rxns)
    _warn_pairwise_size(n, ordered=True)

    must_ul: list[tuple[str, str]] = []
    pos_must_ul: list[tuple[int, int]] = []
    seen: list[str] = []
    seen_set: set[str] = set()

    with model:
        _apply_constr_opt(model, constr_opt)
        _apply_constr_opt(model, constr_mt)

        for a in range(n):
            for b in range(n):
                if a == b:
                    continue
                rxn_a = rxns[a]   # candidate for upregulation
                rxn_b = rxns[b]   # candidate for downregulation
                idx_a, idx_b = rxn_index[rxn_a.id], rxn_index[rxn_b.id]

                with model:
                    # min v_a - v_b: worst case for "a up while b down". If even
                    # this floor exceeds the basal gap, the pair must move apart.
                    model.objective = {rxn_a: 1.0, rxn_b: -1.0}
                    model.objective_direction = 'min'
                    sol = model.optimize()

                    if sol.status != 'optimal':
                        continue

                    va = sol.fluxes[rxn_a.id]
                    vb = sol.fluxes[rxn_b.id]
                    basal_gap = max_fluxes_wt[idx_a] - min_fluxes_wt[idx_b]

                    if (va - vb) > basal_gap + epsilon:
                        must_ul.append((rxn_a.id, rxn_b.id))
                        pos_must_ul.append((idx_a, idx_b))
                        _add_seen(seen, seen_set, rxn_a.id, rxn_b.id)
                        if print_level:
                            print(f"[MustUL] U:{rxn_a.id} L:{rxn_b.id}")

    return MustULResult(
        must_ul=must_ul,
        pos_must_ul=pos_must_ul,
        must_ul_linear=seen,
        pos_must_ul_linear=[rxn_index[r] for r in seen],
    )
