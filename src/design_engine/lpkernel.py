"""
lpkernel.py - OptForce-local LP backends for the ranking pipeline.

The verification stage of :mod:`optforce.pipeline` solves many small LPs against
one fixed model: max biomass, then min/max target at a biomass floor, pinned
growth, the no-production diagnostic, and the import-audit re-solves. Every one
of these is a plain LP whose only mutations are single-reaction column bounds
and a single-reaction linear objective -- OptForce needs no temporary rows here,
because a biomass floor or a target demand is a column bound, not a constraint
row.

Two backends implement one small interface so the pipeline code path is
identical either way:

* :class:`CobraLP` -- the OptForce/COBRA reference path (GLPK via cobra). It is
  the correctness authority and the default.
* :class:`LPKernel` -- a persistent ``highspy.Highs`` built once from the
  model's stoichiometry (CSC-equivalent, one thread, explicit tolerances). It
  mutates column bounds and the objective in place and never rebuilds the
  matrix, so a whole candidate is verified without reconstructing the LP.

The stoichiometric model is built once and reused; bound and objective changes
are made inside :meth:`scope`, which snapshots and restores them so evaluation
order cannot leak between candidates. Neither backend mutates the caller's cobra
model: ``LPKernel`` keeps the model only for structural reads (metabolites,
genes, boundary reactions) and holds all solver state in HiGHS; ``CobraLP``
confines every change to a ``with model:`` block.
"""

from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager

import numpy as np
import cobra


# Feasibility tolerances (kept distinct from the pipeline's class_tol / gamma).
_PRIMAL_TOL = 1e-9
_DUAL_TOL = 1e-9


def _apply_action_bounds(get_bounds, set_bounds, rxn_id, action, forced_flux):
    """
    Shared one-sided-regulation semantics (mirrors optforce._apply_regulation):
    upregulation raises lb, downregulation lowers ub, knockout pins to 0. Only
    ever tightens; a forced level outside the model bounds means the FVA arrays
    do not belong to this model and is raised rather than silently accommodated.
    """
    lb, ub = get_bounds(rxn_id)
    if action == "upregulation":
        if forced_flux > ub + 1e-9:
            raise ValueError(
                f"Cannot upregulate {rxn_id} to {forced_flux:g}: exceeds upper "
                f"bound {ub:g}.")
        set_bounds(rxn_id, max(lb, forced_flux), ub)
    elif action == "downregulation":
        if forced_flux < lb - 1e-9:
            raise ValueError(
                f"Cannot downregulate {rxn_id} to {forced_flux:g}: below lower "
                f"bound {lb:g}.")
        set_bounds(rxn_id, lb, min(ub, forced_flux))
    elif action == "knockout":
        set_bounds(rxn_id, 0.0, 0.0)
    else:
        raise ValueError(f"Unknown intervention type: {action!r}")


class CobraLP:
    """cobra/GLPK reference backend. Every change lives inside a model scope."""

    def __init__(self, model: cobra.Model):
        self.model = model

    def has_reaction(self, rxn_id: str) -> bool:
        return rxn_id in self.model.reactions

    def bounds(self, rxn_id: str) -> tuple[float, float]:
        r = self.model.reactions.get_by_id(rxn_id)
        return r.lower_bound, r.upper_bound

    def set_bound(self, rxn_id, lb=None, ub=None):
        r = self.model.reactions.get_by_id(rxn_id)
        # Set atomically so cobra never sees a transient lb > ub.
        new_lb = r.lower_bound if lb is None else lb
        new_ub = r.upper_bound if ub is None else ub
        r.bounds = (new_lb, new_ub)

    def knockout(self, rxn_id):
        self.set_bound(rxn_id, 0.0, 0.0)

    def apply_action(self, rxn_id, action, forced_flux):
        _apply_action_bounds(self.bounds, self.set_bound, rxn_id, action, forced_flux)

    def _optimize(self, rxn_id, direction):
        self.model.objective = self.model.reactions.get_by_id(rxn_id)
        self.model.objective_direction = direction
        value = self.model.slim_optimize()
        return float(value), self.model.solver.status

    def maximize(self, rxn_id):
        return self._optimize(rxn_id, "max")

    def minimize(self, rxn_id):
        return self._optimize(rxn_id, "min")

    def solve_primal(self, rxn_id, direction):
        self.model.objective = self.model.reactions.get_by_id(rxn_id)
        self.model.objective_direction = direction
        try:
            sol = self.model.optimize()
        except Exception:
            return float("nan"), "error", {}
        if sol.status != "optimal":
            return float("nan"), sol.status, {}
        return float(sol.objective_value), "optimal", dict(sol.fluxes)

    @contextmanager
    def scope(self):
        with self.model:
            yield


class LPKernel:
    """Persistent HiGHS backend built once from the model's stoichiometry."""

    def __init__(self, model: cobra.Model, timeout: float | None = None):
        import highspy  # local import: highspy is an optional dependency

        self._hs = highspy
        self.model = model
        self.rxn_ids = [r.id for r in model.reactions]
        self.index = {rid: j for j, rid in enumerate(self.rxn_ids)}
        n = len(self.rxn_ids)

        self.cur_lb = np.array([r.lower_bound for r in model.reactions], dtype=float)
        self.cur_ub = np.array([r.upper_bound for r in model.reactions], dtype=float)

        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("threads", 1)
        h.setOptionValue("primal_feasibility_tolerance", _PRIMAL_TOL)
        h.setOptionValue("dual_feasibility_tolerance", _DUAL_TOL)
        if timeout is not None:
            h.setOptionValue("time_limit", float(timeout))
        h.addVars(n, self.cur_lb.copy(), self.cur_ub.copy())

        # Steady-state rows: S v = 0, built row-wise (metabolite balances).
        met_row = {met.id: i for i, met in enumerate(model.metabolites)}
        rows: dict[int, list[tuple[int, float]]] = defaultdict(list)
        for j, rxn in enumerate(model.reactions):
            for met, coeff in rxn.metabolites.items():
                rows[met_row[met.id]].append((j, float(coeff)))
        m = len(model.metabolites)
        starts, idx, val = [], [], []
        for i in range(m):
            starts.append(len(idx))
            for j, c in rows[i]:
                idx.append(j)
                val.append(c)
        zeros = np.zeros(m, dtype=float)
        h.addRows(m, zeros, zeros, len(idx),
                  np.array(starts, dtype=np.int32),
                  np.array(idx, dtype=np.int32),
                  np.array(val, dtype=float))
        self._h = h
        self._obj_col: int | None = None
        self._obj_sense = highspy.ObjSense.kMinimize
        self._last_primal: np.ndarray | None = None

    # --- structural / bound access ---
    def has_reaction(self, rxn_id: str) -> bool:
        return rxn_id in self.index

    def bounds(self, rxn_id: str) -> tuple[float, float]:
        j = self.index[rxn_id]
        return float(self.cur_lb[j]), float(self.cur_ub[j])

    def set_bound(self, rxn_id, lb=None, ub=None):
        j = self.index[rxn_id]
        if lb is not None:
            self.cur_lb[j] = lb
        if ub is not None:
            self.cur_ub[j] = ub
        self._h.changeColBounds(j, float(self.cur_lb[j]), float(self.cur_ub[j]))

    def knockout(self, rxn_id):
        self.set_bound(rxn_id, 0.0, 0.0)

    def apply_action(self, rxn_id, action, forced_flux):
        _apply_action_bounds(self.bounds, self.set_bound, rxn_id, action, forced_flux)

    # --- objective / solve ---
    def _set_objective(self, rxn_id, direction):
        j = self.index[rxn_id]
        if self._obj_col is not None and self._obj_col != j:
            self._h.changeColCost(self._obj_col, 0.0)
        self._h.changeColCost(j, 1.0)
        sense = (self._hs.ObjSense.kMaximize if direction == "max"
                 else self._hs.ObjSense.kMinimize)
        self._h.changeObjectiveSense(sense)
        self._obj_col = j
        self._obj_sense = sense

    def _status(self):
        st = self._h.getModelStatus()
        MS = self._hs.HighsModelStatus
        if st == MS.kOptimal:
            return "optimal"
        if st == MS.kInfeasible:
            return "infeasible"
        if st == MS.kUnbounded or st == MS.kUnboundedOrInfeasible:
            return "unbounded"
        if st in (MS.kTimeLimit, MS.kIterationLimit):
            return "time_limit"
        return str(st)

    def _run(self, rxn_id, direction):
        self._set_objective(rxn_id, direction)
        self._h.run()
        status = self._status()
        if status != "optimal":
            self._last_primal = None
            return float("nan"), status
        value = float(self._h.getObjectiveValue())
        if not np.isfinite(value):
            return float("nan"), "error"
        return value, "optimal"

    def maximize(self, rxn_id):
        return self._run(rxn_id, "max")

    def minimize(self, rxn_id):
        return self._run(rxn_id, "min")

    def solve_primal(self, rxn_id, direction):
        value, status = self._run(rxn_id, direction)
        if status != "optimal":
            return value, status, {}
        col_value = np.asarray(self._h.getSolution().col_value, dtype=float)
        fluxes = {rid: float(col_value[self.index[rid]]) for rid in self.rxn_ids}
        return value, status, fluxes

    @contextmanager
    def scope(self):
        lb_snap = self.cur_lb.copy()
        ub_snap = self.cur_ub.copy()
        obj_col, obj_sense = self._obj_col, self._obj_sense
        try:
            yield
        finally:
            # Restore every column bound and the objective in one pass so no
            # per-candidate mutation survives into the next evaluation.
            n = len(self.rxn_ids)
            all_idx = np.arange(n, dtype=np.int32)
            self._h.changeColsBounds(n, all_idx, lb_snap, ub_snap)
            self.cur_lb = lb_snap
            self.cur_ub = ub_snap
            if self._obj_col is not None and self._obj_col != obj_col:
                self._h.changeColCost(self._obj_col, 0.0)
            if obj_col is not None:
                self._h.changeColCost(obj_col, 1.0)
                self._h.changeObjectiveSense(obj_sense)
            self._obj_col, self._obj_sense = obj_col, obj_sense


def make_lp(model: cobra.Model, backend: str, timeout: float | None = None):
    """Build the LP handle for ``backend``.

    ``'cobra'`` is the GLPK reference path. Any backend whose name contains
    ``'highs'`` (``'highs'``, ``'design_engine_highs'``, and the legacy
    ``'optforce_highs'`` / ``'optknock_highs'``) builds the persistent HiGHS
    kernel, falling back with a clear error if ``highspy`` is not installed.
    """
    if backend == "cobra":
        return CobraLP(model)
    if "highs" in backend:
        try:
            return LPKernel(model, timeout=timeout)
        except ImportError as exc:  # highspy not installed
            raise ImportError(
                "backend='highs' requires the optional 'highspy' package "
                "(pip install highspy). Use backend='cobra' otherwise."
            ) from exc
    raise ValueError(
        f"Unsupported backend {backend!r}; use 'cobra' or 'highs'.")
