"""Stage E1 (lean): FASTCC flux-consistency boundary.
Vlassis, Pacheco & Sauter, PLoS Comput Biol 10(1):e1003424, 2014.

Flux-consistent reaction subset of a bounded network: every reaction that
can carry |v| >= eps in some feasible mode (Sv=0, lb<=v<=ub).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog


@dataclass(frozen=True)
class FastccResult:
    consistent: np.ndarray   # int reaction indices carrying flux in some mode
    blocked: np.ndarray      # complement (cannot carry |v| >= eps)
    n_lp: int                # LP-7 solves used (frugality diagnostic)


def _lp7(S, lb, ub, J, eps):
    """LP-7: max sum_{i in J} z_i  s.t.  z in [0,eps], v_i >= z_i (i in J),
    Sv = 0, lb <= v <= ub. Returns v (len n) or None if infeasible."""
    m, n = S.shape
    k = J.size
    if k == 0:
        return np.zeros(n)
    c = np.zeros(n + k); c[n:] = -1.0                      # minimise -sum z
    A_eq = sp.hstack([S, sp.csr_matrix((m, k))], format="csr")
    rows = np.arange(k)
    A_ub = sp.csr_matrix(                                   # z_j - v_{J[j]} <= 0
        (np.r_[np.ones(k), -np.ones(k)],
         (np.r_[rows, rows], np.r_[n + rows, J])),
        shape=(k, n + k))
    bounds = np.vstack([np.column_stack([lb, ub]),
                        np.column_stack([np.zeros(k), np.full(k, eps)])])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(k),
                  A_eq=A_eq, b_eq=np.zeros(m),
                  bounds=bounds, method="highs")
    return res.x[:n] if res.success else None


def fastcc(S, lb, ub, eps=1e-4, tol=None):
    """S: (m,n) sparse stoichiometry; lb, ub: (n,). tol defaults to
    0.99*eps (COBRA convention). Returns FastccResult."""
    S = sp.csr_matrix(S).astype(float)
    lb = np.asarray(lb, float); ub = np.asarray(ub, float)
    m, n = S.shape
    tol = 0.99 * eps if tol is None else tol
    N = set(range(n))
    irr_mask = lb >= 0                                      # irreversible set I
    irr = np.flatnonzero(irr_mask)

    sign = np.ones(n)                                       # column orientations
    lbw, ubw = lb.copy(), ub.copy()
    Seff = lambda: S.multiply(sign.reshape(1, -1)).tocsr()

    def flip(i):
        sign[i] = -sign[i]
        lbw[i], ubw[i] = -ubw[i], -lbw[i]

    n_lp = 0
    v = _lp7(Seff(), lbw, ubw, irr, eps); n_lp += 1        # seed on irreversibles
    A = set(np.flatnonzero(np.abs(v) >= tol).tolist()) if v is not None else set()
    J = N - A - (set(irr.tolist()) - A)                    # drop blocked irreversibles

    flipped = singleton = False
    while J:
        Ji = np.array([min(J)] if singleton else sorted(J))
        v = _lp7(Seff(), lbw, ubw, Ji, eps); n_lp += 1
        if v is not None:
            A |= set(np.flatnonzero(np.abs(v) >= tol).tolist())
        if J & A:
            J -= A
            flipped = False
        else:
            Ji_rev = [i for i in Ji.tolist() if not irr_mask[i]]
            if flipped or not Ji_rev:
                flipped = False
                if singleton:
                    J -= set(Ji.tolist())                  # provably inconsistent
                else:
                    singleton = True
            else:
                for i in Ji_rev:
                    flip(i)
                flipped = True

    consistent = np.array(sorted(A), dtype=int)
    blocked = np.setdiff1d(np.arange(n), consistent, assume_unique=True)
    return FastccResult(consistent, blocked, n_lp)


def fastcc_from_cobra(model, eps=1e-4):
    """Adapter: run fastcc on a cobra model, return (result, consistent_ids)."""
    from cobra.util.array import create_stoichiometric_matrix
    S = sp.csr_matrix(create_stoichiometric_matrix(model, array_type="dok"))
    lb = np.array([r.lower_bound for r in model.reactions], float)
    ub = np.array([r.upper_bound for r in model.reactions], float)
    res = fastcc(S, lb, ub, eps)
    ids = [model.reactions[i].id for i in res.consistent]
    return res, ids
