"""
utils.py - Shared helper functions for OptForce

Contains input validation, model checking, printing, and tolerance
utilities used across fva.py, must_sets.py, optforce.py, and analyze.py.

There is no direct MATLAB equivalent — in the COBRA Toolbox these helpers
were scattered inline across each .m file. Here they are centralised.
"""

from __future__ import annotations

import numpy as np
import cobra
from typing import Optional


# ---------------------------------------------------------------------------
# Model validation
# ---------------------------------------------------------------------------

REQUIRED_MODEL_FIELDS = ['reactions', 'metabolites', 'objective']

def validate_model(model: cobra.Model) -> None:
    """
    Check that a COBRApy model has the minimum fields OptForce needs.

    Raises
    ------
    ValueError if any required field is missing or model is empty.
    """
    for field in REQUIRED_MODEL_FIELDS:
        if not hasattr(model, field):
            raise ValueError(f"Model is missing required field: '{field}'")
    if len(model.reactions) == 0:
        raise ValueError("Model has no reactions.")
    if len(model.metabolites) == 0:
        raise ValueError("Model has no metabolites.")


def validate_reaction_ids(model: cobra.Model, rxn_ids: list[str], label: str = "") -> None:
    """
    Verify that all reaction IDs exist in the model.

    Parameters
    ----------
    model : cobra.Model
    rxn_ids : list of str
    label : str   Label for the source of the IDs (for error messages).

    Raises
    ------
    ValueError if any ID is not found.
    """
    model_rxns = {r.id for r in model.reactions}
    missing = [r for r in rxn_ids if r not in model_rxns]
    if missing:
        src = f" in {label}" if label else ""
        raise ValueError(f"Reaction IDs not found in model{src}: {missing}")


def validate_flux_arrays(
    min_fluxes: np.ndarray,
    max_fluxes: np.ndarray,
    n_rxns: int,
    label: str = "",
) -> None:
    """
    Validate min/max flux arrays returned by FVAOptForce.

    Checks shape matches n_rxns and min <= max (within tolerance).

    Parameters
    ----------
    min_fluxes : np.ndarray
    max_fluxes : np.ndarray
    n_rxns : int
    label : str  e.g. 'wild-type' or 'mutant'
    """
    tag = f" ({label})" if label else ""
    if min_fluxes.shape != (n_rxns,):
        raise ValueError(f"min_fluxes{tag} must have shape ({n_rxns},), got {min_fluxes.shape}")
    if max_fluxes.shape != (n_rxns,):
        raise ValueError(f"max_fluxes{tag} must have shape ({n_rxns},), got {max_fluxes.shape}")
    finite = np.isfinite(min_fluxes) & np.isfinite(max_fluxes)
    if not np.all(finite):
        bad = np.where(~finite)[0]
        raise ValueError(
            f"min/max_fluxes{tag} contain non-finite values at reaction "
            f"indices: {bad.tolist()[:10]}{' ...' if len(bad) > 10 else ''}. "
            "OptForce requires bounded flux ranges; the MATLAB "
            "implementation caps everything at +/-1000."
        )
    violations = np.where(min_fluxes > max_fluxes + 1e-6)[0]
    if len(violations) > 0:
        raise ValueError(
            f"min_fluxes{tag} > max_fluxes for reaction indices: {violations.tolist()}"
        )


# ---------------------------------------------------------------------------
# Tolerance and numerical helpers
# ---------------------------------------------------------------------------

# General-purpose numerical tolerance for classifying an observed flux.
DEFAULT_EPSILON = 1e-3
DEFAULT_TOL = 1e-7

# MUST-set strict-inequality margins, in flux units, as hardcoded in the
# MATLAB outer objectives. Mirrored from must_sets for convenience.
DEFAULT_EPSILON_FIRST_ORDER = 0.5    # findMustL.m:539, findMustU.m:534
DEFAULT_EPSILON_SECOND_ORDER = 0.1   # findMustLL.m:626, findMustUU.m:627


def is_zero(value: float, tol: float = DEFAULT_TOL) -> bool:
    """Return True if value is numerically zero within tol."""
    return abs(value) < tol


def is_upregulated(
    flux: float,
    max_wt: float,
    epsilon: float = DEFAULT_EPSILON,
) -> bool:
    """Return True if flux exceeds the wild-type maximum."""
    return flux > max_wt + epsilon


def is_downregulated(
    flux: float,
    min_wt: float,
    epsilon: float = DEFAULT_EPSILON,
) -> bool:
    """Return True if flux is below the wild-type minimum."""
    return flux < min_wt - epsilon


def is_knockout(flux: float, tol: float = DEFAULT_TOL) -> bool:
    """Return True if flux is effectively zero (knockout)."""
    return abs(flux) < tol


def classify_intervention(
    flux: float,
    min_wt: float,
    max_wt: float,
    epsilon: float = DEFAULT_EPSILON,
    tol: float = DEFAULT_TOL,
) -> str:
    """
    Classify the type of intervention for a reaction given its flux
    in the mutant vs the wild-type range.

    Returns
    -------
    'upregulation'   if flux > max_wt + epsilon
    'downregulation' if flux < min_wt - epsilon (but not zero)
    'knockout'       if |flux| < tol
    'unchanged'      otherwise
    """
    if is_knockout(flux, tol):
        return 'knockout'
    elif is_upregulated(flux, max_wt, epsilon):
        return 'upregulation'
    elif is_downregulated(flux, min_wt, epsilon):
        return 'downregulation'
    else:
        return 'unchanged'


# ---------------------------------------------------------------------------
# Printing / reporting
# ---------------------------------------------------------------------------

def print_fva_summary(
    rxn_ids: list[str],
    min_wt: np.ndarray,
    max_wt: np.ndarray,
    min_mt: np.ndarray,
    max_mt: np.ndarray,
    n_show: int = 10,
) -> None:
    """
    Print a summary table of FVA results for both strains.

    Parameters
    ----------
    rxn_ids : list of str
    min_wt, max_wt : np.ndarray   Wild-type FVA bounds.
    min_mt, max_mt : np.ndarray   Mutant FVA bounds.
    n_show : int   Number of reactions to show (default 10).
    """
    header = f"{'Reaction':<20} {'WT_min':>10} {'WT_max':>10} {'MT_min':>10} {'MT_max':>10}"
    print(header)
    print("-" * len(header))
    for i, rxn in enumerate(rxn_ids[:n_show]):
        print(f"{rxn:<20} {min_wt[i]:>10.4f} {max_wt[i]:>10.4f} {min_mt[i]:>10.4f} {max_mt[i]:>10.4f}")
    if len(rxn_ids) > n_show:
        print(f"  ... and {len(rxn_ids) - n_show} more reactions.")


def print_must_sets(must_u: list[str], must_l: list[str]) -> None:
    """Print MustU and MustL sets."""
    print(f"\nMustU set ({len(must_u)} reactions): {must_u}")
    print(f"MustL set ({len(must_l)} reactions): {must_l}")


def print_optforce_results(opt_force_sets: list) -> None:
    """
    Print OptForce intervention sets in a readable table.

    Parameters
    ----------
    opt_force_sets : list of OptForceSet
    """
    if not opt_force_sets:
        print("No OptForce sets found.")
        return

    print(f"\nOptForce found {len(opt_force_sets)} intervention set(s):\n")
    for i, s in enumerate(opt_force_sets, 1):
        print(f"  Set {i}:")
        for rxn, t, f in zip(s.reactions, s.type_reg, s.fluxes):
            print(f"    {rxn:<20} {t:<16} flux = {f:>8.4f} mmol/gDW/h")
        print(f"    Min target flux : {s.min_target_flux:.4f}")
        print(f"    Max target flux : {s.max_target_flux:.4f}")
        print(f"    Max growth rate : {s.max_growth_rate:.4f} h⁻¹\n")


# ---------------------------------------------------------------------------
# Solver helpers
# ---------------------------------------------------------------------------

def get_solver_name(model: cobra.Model) -> str:
    """Return the name of the solver currently configured for the model."""
    return str(model.solver).split()[0]


def check_solver_available(preferred: str = "gurobi") -> str:
    """
    Check if a preferred solver is available; fall back to glpk.

    Parameters
    ----------
    preferred : str  e.g. 'gurobi', 'cplex', 'glpk'

    Returns
    -------
    str  Name of the available solver to use.
    """
    try:
        import optlang
        solvers = optlang.available_solvers
        if preferred in solvers and solvers[preferred]:
            return preferred
    except Exception:
        pass
    return "glpk"