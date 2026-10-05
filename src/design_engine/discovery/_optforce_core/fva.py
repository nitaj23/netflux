"""
fva.py - FVAOptForce: Flux Variability Analysis for OptForce

Calculates flux ranges for each reaction in the wild-type and mutant strain.
Translated from FVAOptForce.m (COBRA Toolbox, opencobra/cobratoolbox).
"""

from __future__ import annotations

import warnings

import numpy as np
import cobra
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Optional

from .utils import validate_reaction_ids


@dataclass
class StrainConstraints:
    """
    Constraints to apply to a strain (wild-type or mutant).

    Attributes
    ----------
    rxn_list : list of str
        Reaction IDs to constrain.
    values : list of float
        Flux values for each reaction.
    rxn_bound_type : list of str
        Type of bound per reaction:
          'b' -> set both lower and upper bound to value
          'l' -> set only lower bound
          'u' -> set only upper bound

    Example
    -------
    >>> constr = StrainConstraints(
    ...     rxn_list=['R75', 'EX_suc'],
    ...     values=[0.0, 155.55],
    ...     rxn_bound_type=['b', 'b']
    ... )
    """
    rxn_list: list[str] = field(default_factory=list)
    values: list[float] = field(default_factory=list)
    rxn_bound_type: list[str] = field(default_factory=list)

    def __post_init__(self):
        assert len(self.rxn_list) == len(self.values) == len(self.rxn_bound_type), (
            "rxn_list, values, and rxn_bound_type must all have the same length."
        )
        for bt in self.rxn_bound_type:
            assert bt in ('b', 'l', 'u'), (
                f"rxn_bound_type must be 'b', 'l', or 'u'. Got '{bt}'."
            )


@dataclass
class FVAResult:
    """
    Output of fva_optforce().

    Attributes
    ----------
    min_fluxes_wt : np.ndarray  shape (n_rxns,)
    max_fluxes_wt : np.ndarray  shape (n_rxns,)
    min_fluxes_mt : np.ndarray  shape (n_rxns,)
    max_fluxes_mt : np.ndarray  shape (n_rxns,)
    bounds_wt : list of [rxn_id, min, max]
    bounds_mt : list of [rxn_id, min, max]
    """
    min_fluxes_wt: np.ndarray
    max_fluxes_wt: np.ndarray
    min_fluxes_mt: np.ndarray
    max_fluxes_mt: np.ndarray
    bounds_wt: list
    bounds_mt: list


@contextmanager
def _apply_constraints(
    model: cobra.Model,
    constraints: Optional[StrainConstraints],
    objective_rxn: Optional[str] = None,
):
    """
    Context manager: temporarily apply StrainConstraints to the model.

    `objective_rxn` pins the FVA objective explicitly. This matters because
    both the COBRA Toolbox (fluxVariability.m:248,259, with optPercentage=0)
    and COBRApy (fraction_of_optimum=0.0) add the constraint
    `objective_expression >= 0` before varying each reaction. Leaving the
    model's current objective in place therefore silently imposes a
    constraint that depends on whatever the caller last optimised.
    """
    with model:
        if objective_rxn is not None:
            model.objective = model.reactions.get_by_id(objective_rxn)
            model.objective_direction = 'max'
        if constraints is not None:
            for rxn_id, value, bound_type in zip(
                constraints.rxn_list,
                constraints.values,
                constraints.rxn_bound_type,
            ):
                rxn = model.reactions.get_by_id(rxn_id)
                if bound_type == 'b':
                    rxn.lower_bound = value
                    rxn.upper_bound = value
                elif bound_type == 'l':
                    rxn.lower_bound = value
                elif bound_type == 'u':
                    rxn.upper_bound = value
        yield model


def fva_optforce(
    model: cobra.Model,
    constr_wt: Optional[StrainConstraints] = None,
    constr_mt: Optional[StrainConstraints] = None,
    fraction_of_optimum: float = 0.0,
    objective_rxn: Optional[str] = None,
    processes: Optional[int] = 1,
) -> FVAResult:
    """
    Run Flux Variability Analysis for OptForce on wild-type and mutant strains.

    This is a Python translation of FVAOptForce.m from the COBRA Toolbox.
    At least one of constr_wt or constr_mt should be provided.

    Parameters
    ----------
    model : cobra.Model
        A COBRApy metabolic model with fields:
        reactions, metabolites, S matrix, bounds, objective.
    constr_wt : StrainConstraints, optional
        Constraints for the wild-type strain.
    constr_mt : StrainConstraints, optional
        Constraints for the mutant (overproducing) strain.
    fraction_of_optimum : float
        Fraction of the objective optimum to require during FVA (default 0.0).
        OptForce normally encodes the growth/production requirements directly
        in `constr_wt` / `constr_mt` (e.g. biomass pinned for the wild type,
        biomass = 0 and target forced for the mutant), so the default of 0.0
        avoids double-constraining the objective. Only raise this if you are
        NOT pinning the objective through the StrainConstraints.
    objective_rxn : str, optional
        Reaction to use as the FVA objective for both strains. Strongly
        recommended: with fraction_of_optimum=0.0 COBRApy still adds
        `objective >= 0` to the problem, so leaving whatever objective the
        model happens to carry silently changes the FVA result. MATLAB's
        FVAOptForce.m:142 uses the model's `c` vector, which in the OptForce
        tutorial is the biomass reaction. Default None = use the model's
        current objective (and warn).
    processes : int, optional
        Number of parallel processes for FVA. Defaults to 1 (serial),
        matching fluxVariability in FVAOptForce.m:142, which is also serial.
        Values > 1 make COBRApy spawn a multiprocessing pool; on Windows and
        anywhere else using the 'spawn' start method the calling script MUST
        be guarded by `if __name__ == "__main__":`, or each child re-imports
        the script and spawns its own pool recursively.

    Returns
    -------
    FVAResult
        Contains min/max flux arrays and bounds tables for both strains.

    Example
    -------
    >>> import cobra.test
    >>> model = cobra.test.create_test_model('textbook')
    >>> constr_wt = StrainConstraints(
    ...     rxn_list=['EX_glc__D_e'],
    ...     values=[-10.0],
    ...     rxn_bound_type=['b']
    ... )
    >>> constr_mt = StrainConstraints(
    ...     rxn_list=['EX_succ_e'],
    ...     values=[155.55],
    ...     rxn_bound_type=['l']
    ... )
    >>> result = fva_optforce(model, constr_wt, constr_mt)
    """
    if constr_wt is None and constr_mt is None:
        raise ValueError(
            "At least one of constr_wt or constr_mt must be provided."
        )

    if objective_rxn is None:
        warnings.warn(
            "fva_optforce() called without `objective_rxn`; using the model's "
            "current objective. Because fraction_of_optimum adds an "
            "`objective >= 0` constraint, this makes the FVA result depend on "
            "whatever was last optimised. Pass objective_rxn explicitly "
            "(normally the biomass reaction).",
            RuntimeWarning,
            stacklevel=2,
        )
    else:
        validate_reaction_ids(model, [objective_rxn], "objective_rxn")

    rxn_ids = [r.id for r in model.reactions]

    def _run(constraints, label):
        with _apply_constraints(model, constraints, objective_rxn) as m:
            try:
                result = flux_variability_analysis(
                    m,
                    fraction_of_optimum=fraction_of_optimum,
                    processes=processes,
                )
            except OptimizationError as exc:
                raise OptimizationError(
                    f"FVA failed for the {label} strain: the model is "
                    f"infeasible or unbounded under its constraints ({exc})."
                ) from exc
        # Reindex to match model reaction order
        result = result.reindex(rxn_ids)
        return result["minimum"].to_numpy(), result["maximum"].to_numpy()

    min_fluxes_wt, max_fluxes_wt = _run(constr_wt, "wild-type")
    min_fluxes_mt, max_fluxes_mt = _run(constr_mt, "mutant")

    # --- Build bounds tables: [[rxn_id, min, max], ...] ---
    bounds_wt = [
        [rxn_id, mn, mx]
        for rxn_id, mn, mx in zip(rxn_ids, min_fluxes_wt, max_fluxes_wt)
    ]
    bounds_mt = [
        [rxn_id, mn, mx]
        for rxn_id, mn, mx in zip(rxn_ids, min_fluxes_mt, max_fluxes_mt)
    ]

    return FVAResult(
        min_fluxes_wt=min_fluxes_wt,
        max_fluxes_wt=max_fluxes_wt,
        min_fluxes_mt=min_fluxes_mt,
        max_fluxes_mt=max_fluxes_mt,
        bounds_wt=bounds_wt,
        bounds_mt=bounds_mt,
    )