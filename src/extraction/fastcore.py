"""
src/fastcore.py

FASTCORE (Vlassis, Pacheco & Sauter, PLoS Comput Biol 2014) — Stage 1 of the
flux-essential bioproduction pathway extraction pipeline.

Given a GEM and a core reaction set C, FASTCORE returns the minimal (greedy)
flux-consistent subnetwork F ⊇ C: every reaction in F carries |flux| ≥ ε in at
least one feasible steady-state mode, and F contains no blocked reactions.

Implementation notes
--------------------
Two LPs drive the search (paper's LP-7 and LP-9), solved here with
scipy HiGHS on the model's own stoichiometry (signed flux variables with
lb/ub), avoiding cobra's forward/reverse variable splitting:

  LP-7  max Σ_{i∈J} z_i,  0 ≤ z_i ≤ ε,  v_i ≥ z_i,  Sv = 0,  lb ≤ v ≤ ub
        (activate as many core reactions as possible — flux-splitting)
  LP-9  min Σ_{i∈P} |v_i|,  v_i ≥ ε ∀i∈K,  Sv = 0,  lb ≤ v ≤ ub
        (add as few non-core reactions as possible — L1 sparsity on P)

The reversible-reaction orientation flipping follows the reference MATLAB
implementation exactly (singleton / flipped state machine).

References
----------
Vlassis, Pacheco & Sauter. Fast Reconstruction of Compact Context-Specific
Metabolic Network Models. PLoS Comput Biol 10(1):e1003424, 2014.
"""

from __future__ import annotations

import logging

import cobra
import numpy as np
import scipy.sparse as sp
from cobra.flux_analysis import pfba
from scipy.optimize import linprog

logger = logging.getLogger(__name__)

DEFAULT_EPSILON = 1e-5  # COBRA feasibility tolerance × 100 (blueprint default)


# ─────────────────────────────────────────────────────────────────────────────
# Core LPs
# ─────────────────────────────────────────────────────────────────────────────

def _lp7(J: list[int], S: sp.spmatrix, lb: np.ndarray, ub: np.ndarray,
         epsilon: float) -> np.ndarray:
    """max Σ_{i∈J} z_i s.t. v_i ≥ z_i, 0 ≤ z_i ≤ ε, Sv=0, lb≤v≤ub. Returns v."""
    m, n = S.shape
    nj = len(J)
    if nj == 0:
        return np.zeros(n)

    # variables x = [v (n), z (nj)];  minimise -Σ z
    c = np.concatenate([np.zeros(n), -np.ones(nj)])

    A_eq = sp.hstack([S, sp.csr_matrix((m, nj))], format="csr")
    b_eq = np.zeros(m)

    # z_k - v_{J[k]} ≤ 0   (i.e. v_i ≥ z_i)
    rows = np.repeat(np.arange(nj), 2)
    cols = np.empty(2 * nj, dtype=int)
    data = np.empty(2 * nj)
    cols[0::2] = J;              data[0::2] = -1.0   # -v_i
    cols[1::2] = n + np.arange(nj); data[1::2] = 1.0  # +z_k
    A_ub = sp.csr_matrix((data, (rows, cols)), shape=(nj, n + nj))
    b_ub = np.zeros(nj)

    bounds = list(zip(lb, ub)) + [(0.0, epsilon)] * nj
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, method="highs")
    return res.x[:n] if res.success else np.zeros(n)


def _lp9(K: list[int], P: list[int], S: sp.spmatrix, lb: np.ndarray,
         ub: np.ndarray, epsilon: float) -> np.ndarray:
    """min Σ_{i∈P}|v_i| s.t. v_i≥ε ∀i∈K, Sv=0, lb≤v≤ub. Returns v."""
    m, n = S.shape
    npen = len(P)

    # variables x = [v (n), z (npen)];  minimise Σ z   (z_p ≥ |v_p|)
    c = np.concatenate([np.zeros(n), np.ones(npen)])

    A_eq = sp.hstack([S, sp.csr_matrix((m, npen))], format="csr")
    b_eq = np.zeros(m)

    if npen:
        # per penalty reaction p:  v_p - z_p ≤ 0  and  -v_p - z_p ≤ 0
        rows = np.repeat(np.arange(2 * npen), 2)
        cols = np.empty(4 * npen, dtype=int)
        data = np.empty(4 * npen)
        zcol = n + np.arange(npen)
        # row 2k:  +v_p - z_p
        cols[0::4] = P;    data[0::4] = 1.0
        cols[1::4] = zcol; data[1::4] = -1.0
        # row 2k+1: -v_p - z_p
        cols[2::4] = P;    data[2::4] = -1.0
        cols[3::4] = zcol; data[3::4] = -1.0
        A_ub = sp.csr_matrix((data, (rows, cols)), shape=(2 * npen, n + npen))
        b_ub = np.zeros(2 * npen)
    else:
        A_ub, b_ub = None, None

    lb2 = lb.copy()
    for i in K:
        lb2[i] = epsilon  # force core-active reactions above threshold
    bounds = list(zip(lb2, ub)) + [(0.0, None)] * npen
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, method="highs")
    return res.x[:n] if res.success else np.zeros(n)


def _find_sparse_mode(J: set[int], P: set[int], singleton: bool,
                      S: sp.spmatrix, lb: np.ndarray, ub: np.ndarray,
                      epsilon: float) -> set[int]:
    """Find a sparse mode supporting (part of) J while penalising P."""
    if not J:
        return set()
    Jw = [min(J)] if singleton else list(J)

    v = _lp7(Jw, S, lb, ub, epsilon)
    K = [i for i in Jw if v[i] >= 0.99 * epsilon]
    if not K:
        return set()

    v = _lp9(K, list(P), S, lb, ub, epsilon)
    return {i for i in range(S.shape[1]) if abs(v[i]) >= 0.99 * epsilon}


# ─────────────────────────────────────────────────────────────────────────────
# FASTCORE
# ─────────────────────────────────────────────────────────────────────────────

def fastcore(model: cobra.Model, core: list[str],
             epsilon: float = DEFAULT_EPSILON) -> list[str]:
    """
    Compute the minimal flux-consistent subnetwork containing `core`.

    Parameters
    ----------
    model    cobra model (its bounds define the condition/medium)
    core     reaction bigg_ids that must be in the output (e.g. [EX_M, uptakes])
    epsilon  flux threshold defining "carries flux"

    Returns
    -------
    Sorted list of reaction bigg_ids forming the flux-consistent subnetwork F.

    Raises
    ------
    ValueError if a core reaction cannot carry flux in the given model/medium
    (the core set is inconsistent — usually a missing substrate uptake).
    """
    rxn_ids = [r.id for r in model.reactions]
    idx = {r: i for i, r in enumerate(rxn_ids)}
    n = len(rxn_ids)

    missing = [r for r in core if r not in idx]
    if missing:
        raise ValueError(f"Core reactions not in model: {missing}")

    # signed stoichiometry + bounds as mutable copies (columns get flipped)
    S = sp.csc_matrix(cobra.util.create_stoichiometric_matrix(model, array_type="dense"))
    lb = np.array([r.lower_bound for r in model.reactions], dtype=float)
    ub = np.array([r.upper_bound for r in model.reactions], dtype=float)

    def flip(cols: list[int]) -> None:
        for j in cols:
            S.data[S.indptr[j]:S.indptr[j + 1]] *= -1.0
        new_lb = -ub[cols]
        ub[cols] = -lb[cols]
        lb[cols] = new_lb

    N = set(range(n))
    C = {idx[r] for r in core}
    I = {i for i in range(n) if lb[i] >= 0}     # irreversible reactions

    flipped = False
    singleton = False

    # ── initial pass: irreversible core ──────────────────────────────────────
    J = C & I
    P = N - C
    supp = _find_sparse_mode(J, P, singleton, S, lb, ub, epsilon)
    inconsistent = J - supp
    if inconsistent:
        bad = [rxn_ids[i] for i in inconsistent]
        raise ValueError(
            f"Core set inconsistent - these irreversible core reactions cannot "
            f"carry flux in this medium: {bad[:5]}"
            + (f" (+{len(bad) - 5} more)" if len(bad) > 5 else "")
        )
    A = set(supp)
    J = C - A

    # ── main loop ────────────────────────────────────────────────────────────
    while J:
        P = P - A
        supp = _find_sparse_mode(J, P, singleton, S, lb, ub, epsilon)
        A |= supp

        if J & A:
            J = J - A
            flipped = False
        else:
            if singleton:
                j1 = min(J)
                JiRev = {j1} - I
            else:
                JiRev = J - I

            if flipped or not JiRev:
                if singleton:
                    stuck = [rxn_ids[i] for i in J]
                    raise ValueError(
                        f"Network inconsistent - cannot add core reaction(s) "
                        f"{stuck[:5]} to a flux-consistent mode."
                    )
                flipped = False
                singleton = True
            else:
                flip(list(JiRev))
                flipped = True

    return sorted(rxn_ids[i] for i in A)


# ─────────────────────────────────────────────────────────────────────────────
# Bioproduction core seeding (blueprint Stage 1)
# ─────────────────────────────────────────────────────────────────────────────

def bioproduction_core(model: cobra.Model, target_exchange: str,
                       include_uptakes: bool = True,
                       epsilon: float = DEFAULT_EPSILON) -> list[str]:
    """
    Build C_FAST for bioproduction seeding (blueprint Stage 1):

        C_FAST = {EX_M} ∪ {required exchange/uptake reactions for v_EX_M > 0}

    The uptake set is discovered by maximising the target under **pFBA**
    (parsimonious FBA) and taking the exchanges that carry import flux ≥ epsilon
    in that minimal-flux state — the substrates the productive mode actually
    consumes. Two filters keep the seed flux-consistent (so FASTCORE can always
    satisfy it):

      * pFBA (not plain FBA) strips incidental uptake flux from alternative-
        optima degeneracy — e.g. a reversible phosphate exchange drifting to a
        nonzero value — that would otherwise seed a blocked reaction into F.
      * the `epsilon` cutoff drops trace cofactor/vitamin uptakes (nicotinate,
        pantothenate, ...) whose flux sits below FASTCORE's own flux threshold.
        Such reactions cannot carry ≥ epsilon flux consistently, so forcing them
        into the core would make it inconsistent. This is the same threshold
        FASTCORE uses to decide a reaction "carries flux", so the seed is
        consistent-by-construction at that resolution.

    Seeding these anchors F on the biologically real substrate pathways feeding
    EX_M. Set include_uptakes=False to recover the target-only seed ({EX_M}),
    letting FASTCORE discover the support itself.
    """
    if target_exchange not in model.reactions:
        raise ValueError(f"Target exchange '{target_exchange}' not in model.")

    core = {target_exchange}
    if include_uptakes:
        with model:
            model.objective = target_exchange
            try:
                fluxes = pfba(model).fluxes
            except Exception:  # pFBA infeasible/unbounded — fall back to FBA
                fluxes = model.optimize().fluxes
            core |= {ex.id for ex in model.exchanges
                     if ex.id != target_exchange
                     and fluxes.get(ex.id, 0.0) <= -epsilon}
    return sorted(core)
