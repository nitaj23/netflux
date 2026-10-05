# src/extraction/gudie.py
"""
GUDIE subgraph extraction: PPR-VUDIE + FCA-LUDIE.

Consumes the FASTCORE flux-consistent boundary G_F (a NetworkX DiGraph built by
ingest.build_mfg / graph_ops.build_mfg, reaction-id nodes, 'mass_flow' edge attr)
and returns the bioproduction-relevant subgraph H, upstream of the target EX_M.

Pipeline:
    2.1  log-median edge-weight compression        -> w_hat
    2.2  FCA coupling of each reaction to EX_M      -> c
    2.4  target-seeded Personalized PageRank        -> p            (VUDIE)
    2.3  L(i,j) = sqrt(c_i * w_hat_ij)              -> edge interest (LUDIE)
    2.6  interest propagation (gamma_1, h hops)     -> I^h
    2.7  seed expansion (theta_2 decay, threshold k)-> H
"""
from __future__ import annotations

import logging
import math
from collections import deque
from dataclasses import dataclass, field
from typing import Callable, Iterable, Literal, Mapping, Sequence

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

logger = logging.getLogger(__name__)

__all__ = ["GudieConfig", "GudieResult", "run_gudie",
           "log_median_compress", "personalized_pagerank", "flux_coupling_to_target"]

Decay = Literal["theta1", "theta2"]
PPRNorm = Literal["symmetric", "rw"]


 
# Config / result containers
 
@dataclass(slots=True)
class GudieConfig:
    alpha_r: float = 0.3                    # restart probability in (0,1); higher => tighter
                                            # localisation around the target (0.15-0.3 typical)
    ppr_normalization: PPRNorm = "rw"       # proper column-stochastic walk on G^R; converges on
                                            # directed graphs where "symmetric" may not (design note 2)
    ppr_tol: float = 1e-9                   # L1 convergence tolerance
    ppr_max_iter: int = 200
    # compression (Step 2.1)
    compress_eps: float = 1e-9              # MAD denominator floor
    # propagation (Step 2.6)
    hops: int = 8                           # deeper spread reaches multi-step upstream pathways
    # seed expansion (Step 2.7)
    k: float = 0.1                          # inclusion threshold factor delta = k * ref, ref =
                                            # strongest interest in the target's first upstream ring
                                            # (see _threshold_reference); key size knob
    decay: Decay = "theta1"                 # theta_1 (1/|L|) decays slower than theta_2 -> admits
                                            # longer pathways to the product
    # FCA coupling (Step 2.2)
    use_fca: bool = True
    fca_gamma: float = 1.0                  # fraction of optimal target flux to fix
    fca_pure_stoichiometric: bool = False   # relax magnitude bounds -> bound-free ratios
    fca_zero_span_to_zero: bool = True      # 0 in [Rmin,Rmax] -> uncoupled (c=0)
    fca_max_reactions: int | None = None    # guard against pathological F sizes
    # Phase B (advisory FPE embedding blend; lam=0 => exact Phase A)
    lam: float = 0.05                       # mixing coefficient in [0,1]; 0.05 default -- a light
                                            # Flow-Profile-Encoding nudge that lifts flux recall and
                                            # rescues near-collapsed extractions while staying compact
                                            # (validate_gudie_all_organisms sweep). >0.1 grows H fast;
                                            # only used when embedding_sim is supplied to run_gudie.


@dataclass(slots=True)
class GudieResult:
    subgraph: nx.DiGraph                    # H, original MFG edge directions restored
    nodes: list[str]                        # reaction ids retained in H
    node_order: list[str]                   # index order of all F nodes for the arrays below
    ppr: np.ndarray                         # VUDIE p_i, L1-normalised (aligned to node_order)
    coupling: np.ndarray | None             # c_i in [0,1] (None if FCA unavailable -> fallback)
    propagated_interest: np.ndarray         # I^h aligned to node_order
    seed_interest: float                    # I^h at the target (diagnostic; see _threshold_reference)
    threshold: float                        # delta = k * ref (ref = max I^h over seed's producers)
    ppr_converged: bool
    ludie_fallback: bool                    # True if edge-PPR throughput LUDIE was used
    source_ppr: np.ndarray | None = None    # dual-anchor: downstream field seeded at source_ids (None if unused)


 
# Step 2.1 — log-median edge weight compression
 
def log_median_compress(weights: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """
    Map strictly-positive, heavy-tailed mass flows to [0,1] via a log transform
    centred on the median with logistic-MAD squashing (blueprint Step 2.1;
    Almaas et al. 2004 power-law flux distribution).

        f_ij   = log(A_ij / median(A))
        w_hat  = sigmoid( (f_ij - median(f)) / (MAD(f) + eps) )

    median(f) is 0 by construction (log of median/median); computed explicitly
    for numerical safety. Robust (50% breakdown) to the high-flux backbone.
    """
    w = np.asarray(weights, dtype=np.float64)
    if w.size == 0:
        return w.astype(np.float64)
    if np.any(w <= 0):
        # build_mfg thresholds edges at >1e-9; defend anyway.
        raise ValueError("log_median_compress requires strictly positive weights")
    mu = np.median(w)
    f = np.log(w / mu)
    med_f = np.median(f)
    mad = np.median(np.abs(f - med_f))
    z = (f - med_f) / (mad + eps)
    return 1.0 / (1.0 + np.exp(-z))


def _compressed_adjacency(
    G_F: nx.DiGraph, node_index: Mapping[str, int], n: int, eps: float
) -> sp.csr_matrix:
    """Original compressed adjacency A_hat (row=source/producer, col=target/consumer)."""
    if G_F.number_of_edges() == 0:
        return sp.csr_matrix((n, n), dtype=np.float64)
    rows = np.empty(G_F.number_of_edges(), dtype=np.int64)
    cols = np.empty_like(rows)
    data = np.empty(rows.shape[0], dtype=np.float64)
    for e, (u, v, d) in enumerate(G_F.edges(data=True)):
        rows[e] = node_index[u]
        cols[e] = node_index[v]
        data[e] = float(d.get("mass_flow", 0.0))
    keep = data > 0.0
    rows, cols, data = rows[keep], cols[keep], data[keep]
    data = log_median_compress(data, eps)
    A = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    A.sum_duplicates()
    return A


 
# Step 2.4 — target-seeded Personalized PageRank (VUDIE prerequisite)
 
def _ppr_operator(A_hat: sp.csr_matrix, normalization: PPRNorm):
    """
    Build the (matvec) operator that spreads mass EX_M -> upstream producers.

    Design note (2): under the left-multiply update p <- a*p0 + (1-a)*W p, the
    operator that moves mass *against* MFG edges (producer<-EX_M) is the ORIGINAL
    A_hat with D = its column sums (== blueprint's D_ii = sum_j A^R_ij). We do NOT
    transpose. See module docstring / accompanying design notes.
    """
    d_col = np.asarray(A_hat.sum(axis=0)).ravel()          # column sums = in-strength
    if normalization == "symmetric":
        dm = np.zeros_like(d_col)
        nz = d_col > 0
        dm[nz] = 1.0 / np.sqrt(d_col[nz])                  # D^{-1/2}, 0 on dangling
        Dm = sp.diags(dm)
        W = (Dm @ A_hat @ Dm).tocsr()
    elif normalization == "rw":
        di = np.zeros_like(d_col)
        nz = d_col > 0
        di[nz] = 1.0 / d_col[nz]                           # column-stochastic on G^R
        W = (A_hat @ sp.diags(di)).tocsr()
    else:  # pragma: no cover
        raise ValueError(f"unknown normalization {normalization!r}")
    return W


def personalized_pagerank(
    A_hat: sp.csr_matrix,
    seed_idx: int | Sequence[int],
    alpha_r: float = 0.2,
    normalization: PPRNorm = "symmetric",
    tol: float = 1e-9,
    max_iter: int = 200,
) -> tuple[np.ndarray, bool]:
    """
    Iterative PPR (blueprint Step 2.4):  p <- a*p0 + (1-a) W p,  p0 = unit@seed.
    Each step is one sparse matvec, O(|E|). Returns (L1-normalised p, converged).
    seed_idx may be a single node or a set of nodes (restart mass spread uniformly);
    a lone int reproduces the single-seed vector exactly.
    """
    n = A_hat.shape[0]
    seeds = [seed_idx] if isinstance(seed_idx, (int, np.integer)) else list(seed_idx)
    p0 = np.zeros(n, dtype=np.float64)
    p0[seeds] = 1.0 / len(seeds)
    if A_hat.nnz == 0:
        return p0, True
    W = _ppr_operator(A_hat, normalization)
    p = p0.copy()
    converged = False
    for _ in range(max_iter):
        p_new = alpha_r * p0 + (1.0 - alpha_r) * (W @ p)
        if np.abs(p_new - p).sum() < tol:
            p = p_new
            converged = True
            break
        p = p_new
    s = p.sum()
    if s > 0:
        p = p / s
    if not converged:
        logger.warning("PPR did not converge in %d iters (||W||_2 may exceed 1 under "
                       "symmetric norm on a directed graph); consider ppr_normalization='rw'.",
                       max_iter)
    return p, converged


def _katz_fallback(A_hat: sp.csr_matrix, seed_idx: int, beta: float = 1.0) -> np.ndarray:
    """
    Target-personalised Katz fallback (blueprint degradation table; VUDIE cross-check):
        x = alpha * A_hat x + beta * e_seed,   alpha < 1/lambda_max(A_hat).
    Same left-multiply operator direction as PPR (upstream spread).
    """
    n = A_hat.shape[0]
    e = np.zeros(n); e[seed_idx] = beta
    if A_hat.nnz == 0:
        return e / (e.sum() or 1.0)
    try:
        lam = float(np.abs(spla.eigs(A_hat, k=1, which="LM", return_eigenvectors=False)[0]))
    except (spla.ArpackNoConvergence, spla.ArpackError, ValueError):
        lam = float(np.abs(A_hat).sum(axis=1).max())       # crude row-sum bound
    alpha = 0.9 / lam if lam > 1e-12 else 0.0
    x = e.copy()
    for _ in range(200):
        x_new = alpha * (A_hat @ x) + e
        if np.abs(x_new - x).sum() < 1e-9:
            x = x_new
            break
        x = x_new
    s = x.sum()
    return x / s if s > 0 else x


 
# Step 2.2 — FCA coupling of each reaction to the target (LUDIE prerequisite)
 
def flux_coupling_to_target(
    model,                                  # cobra.Model
    target_id: str,
    rxn_subset: Sequence[str],
    *,
    gamma: float = 1.0,
    pure_stoichiometric: bool = False,
    zero_span_to_zero: bool = True,
    cache: dict | None = None,
) -> dict[str, float] | None:
    """
    Continuous FCA coupling strength c_i in [0,1] of each reaction to EX_M
    (blueprint Step 2.2). Fixes target flux, then min/maxes each reaction:

        c_i = min(|Rmin|,|Rmax|) / max(|Rmin|,|Rmax|)   if 0 not in [Rmin,Rmax]
        c_i = 0                                          if 0 in [Rmin,Rmax] (decoupled)
              (fully coupled -> 1; can idle while target active -> 0)

    Note: with finite magnitude bounds this is target-conditional coupling (bounds
    influence the ratio). pure_stoichiometric=True relaxes magnitude bounds (keeps
    only directionality) to recover bound-free stoichiometric coupling ~ Charnes-
    Cooper. c_i is MFGHub's continuous construct, not attributable to Burgard 2004.

    Returns None (caller falls back to edge-PPR LUDIE) if the target is infeasible.
    """
    if cache is not None:
        key = (id(model), target_id, gamma, pure_stoichiometric, tuple(rxn_subset))
        if key in cache:
            return cache[key]

    try:
        import cobra  # noqa: F401  (only needed here)
    except ImportError:
        logger.warning("cobra unavailable; FCA coupling skipped -> LUDIE fallback")
        return None

    with model:  # context manager restores bounds/objective on exit
        target = model.reactions.get_by_id(target_id)
        model.objective = target
        model.objective.direction = "max"
        v_star = model.slim_optimize()
        if v_star is None or math.isnan(v_star) or abs(v_star) < 1e-9:
            logger.warning("Target %s infeasible or carries ~0 flux; FCA skipped.", target_id)
            return None

        if pure_stoichiometric:
            for r in model.reactions:
                lo = -1e6 if r.lower_bound < 0 else 0.0    # keep sign; drop magnitude
                hi = 1e6 if r.upper_bound > 0 else 0.0
                r.bounds = (lo, hi)

        t_val = gamma * v_star
        target.bounds = (t_val, t_val)                      # fix target flux at gamma * v_star

        c: dict[str, float] = {}
        subset = list(rxn_subset)
        for idx, rid in enumerate(subset):
            r = model.reactions.get_by_id(rid)
            model.objective = r
            model.objective.direction = "max"
            r_max = model.slim_optimize()
            model.objective.direction = "min"
            r_min = model.slim_optimize()
            c[rid] = _coupling_from_bounds(r_min, r_max, zero_span_to_zero)
            if logger.isEnabledFor(logging.DEBUG) and idx % 500 == 0 and idx:
                logger.debug("FCA coupling: %d/%d reactions", idx, len(subset))

    if cache is not None:
        cache[key] = c
    return c


def _coupling_from_bounds(r_min, r_max, zero_span_to_zero: bool) -> float:
    if r_min is None or r_max is None or math.isnan(r_min) or math.isnan(r_max):
        return 0.0                                          # unbounded/infeasible -> uncoupled
    lo, hi = min(r_min, r_max), max(r_min, r_max)
    if lo <= 1e-9 and hi >= -1e-9 and zero_span_to_zero:    # 0 in [Rmin,Rmax]
        return 0.0
    a, b = abs(lo), abs(hi)
    denom = max(a, b)
    if denom <= 1e-12 or not math.isfinite(denom):
        return 0.0
    return float(min(a, b) / denom)                         # in [0,1]; 1 == fully coupled


 
# Steps 2.3/2.5/2.6/2.7 — interest functions, propagation, seed expansion
 
def _norm_max(x: np.ndarray) -> np.ndarray:
    """Scale to peak 1 for combining fields (product of normalized scores). Safe on all-zero."""
    m = float(x.max()) if x.size else 0.0
    return x / m if m > 0.0 else x


def _edge_interest_matrix(
    A_rev: sp.csr_matrix,           # reversed compressed adjacency (G^R): A_rev[i,j]=w_hat(orig j->i)
    coupling: np.ndarray | None,
    ppr: np.ndarray,
) -> tuple[sp.csr_matrix, bool]:
    """
    LUDIE on G^R edges (Step 2.3):  L(i,j) = sqrt(c_i * w_hat^R_ij)   [source coupling].
    Fallback (Step 2.3, c unavailable): edge-PPR throughput  L(i,j) = p_i * W^R_ij
    with W^R row-normalised. Geometric fusion => conjunctive zero-forcing.
    """
    A = A_rev.tocoo(copy=True)
    if coupling is not None:
        vals = np.sqrt(np.clip(coupling[A.row], 0.0, None) * A.data)
        fallback = False
    else:
        rowsum = np.asarray(A_rev.sum(axis=1)).ravel()
        inv = np.zeros_like(rowsum); nz = rowsum > 0
        inv[nz] = 1.0 / rowsum[nz]
        vals = ppr[A.row] * (A.data * inv[A.row])           # p_i * normalised transition
        fallback = True
    B = sp.coo_matrix((vals, (A.row, A.col)), shape=A_rev.shape).tocsr()
    B.eliminate_zeros()
    return B, fallback


def _propagate(B: sp.csr_matrix, A_rev: sp.csr_matrix, I0: np.ndarray, hops: int) -> np.ndarray:
    """
    Interest propagation with GUDIE gamma_1 aggregation (Silva Eq. 1; Step 2.6):

        phi(m,n)  = I_m * B[m,n]                              (message m -> n on G^R)
        I_n^(t)   = I_n^(t-1)/2 + mean_{m in N(n)} phi(m,n)/2

    N(n) = G^R-predecessors of n (structural, so decoupled edges still dilute the
    mean -> P4: degree alone cannot inflate interest). Vectorised: B^T @ I.
    """
    BT = B.T.tocsr()
    deg = A_rev.getnnz(axis=0).astype(np.float64)            # |N(n)|: G^R in-degree (structural)
    inv_deg = np.zeros_like(deg); nz = deg > 0
    inv_deg[nz] = 1.0 / deg[nz]
    I = I0.astype(np.float64, copy=True)
    for _ in range(max(0, hops)):
        I = 0.5 * I + 0.5 * (inv_deg * (BT @ I))
    return I


def _retain_factor(decay: Decay) -> Callable[[int], float]:
    """1 - theta(|L|): theta_1 = 1 - 1/|L| -> 1/|L|;  theta_2 = 1 - e^{1-|L|} -> e^{1-|L|}."""
    if decay == "theta2":
        return lambda L: math.exp(1 - L)
    return lambda L: 1.0 / L


def _threshold_reference(A_rev: sp.csr_matrix, I: np.ndarray, seed_idx: int) -> float:
    """
    Reference interest level the inclusion threshold delta = k * ref is anchored to.

    The GUDIE paper anchors delta to the seed's own interest, but that is unsound
    when the target EX_M is a pure sink (no MFG out-edges) — the usual case here.
    Interest only flows *toward* the seed on G^R, so the seed never receives a
    propagated message: I[seed] is just its trapped PPR restart mass decaying as
    0.5^hops (see _propagate), on a scale unrelated to any producer's. Anchoring
    delta to it makes k*I_seed exceed even the strongest direct producer, so seed
    expansion returns only {seed} (see design note / validation on iNF517).

    Instead anchor to the strongest interest in the seed's first upstream ring
    (its direct G^R-producers). This guarantees the primary producer clears the
    bar for any k<=1 (retain(1)=1), so H can never collapse to the target alone,
    and makes k a meaningful "fraction of the primary upstream signal". Falls back
    to the global non-seed max, then to I[seed], for degenerate seeds.
    """
    producers = A_rev.indices[A_rev.indptr[seed_idx]:A_rev.indptr[seed_idx + 1]]
    if producers.size:
        return float(I[producers].max())
    if I.shape[0] > 1:
        masked = I.copy()
        masked[seed_idx] = -np.inf
        m = float(masked.max())
        if math.isfinite(m):
            return m
    return float(I[seed_idx])


def _seed_expansion(
    A_rev: sp.csr_matrix, I: np.ndarray, seed_idx: int, delta: float, decay: Decay
) -> set[int]:
    """
    Single-seed GUDIE seed expansion (Silva Alg. 4; Step 2.7) as a pruned BFS on G^R.

    A node m is included iff, at its shortest reachable path length |L| through
    already-included nodes,  (1 - theta(|L|)) * I_m >= delta.  delta = k * ref is
    computed by the caller from _threshold_reference (anchored to the seed's first
    upstream ring, not I_seed itself — see that function). Because (1 - theta) is
    monotone decreasing in |L|, first (BFS) reach is the best shot; only included
    nodes propagate. O(V + E).
    """
    retain = _retain_factor(decay)
    indptr, indices = A_rev.indptr, A_rev.indices           # row i -> G^R successors (orig preds)
    included = {seed_idx}
    dist = {seed_idx: 0}
    q: deque[int] = deque([seed_idx])
    while q:
        n = q.popleft()
        L = dist[n] + 1
        rf = retain(L)
        for j in indices[indptr[n]:indptr[n + 1]]:
            if j in included:
                continue
            if rf * I[j] >= delta:
                included.add(j)
                dist[j] = L
                q.append(j)
    return included


 
# Orchestration
 
def run_gudie(
    G_F: nx.DiGraph,
    target_id: str,
    *,
    source_ids: Sequence[str] | None = None,    # dual-anchor: substrate/uptake reaction ids
    model=None,                                 # cobra.Model, required for FCA (Step 2.2)
    coupling: Mapping[str, float] | None = None,  # inject precomputed F2C2/FCA result
    embedding_sim: Mapping[str, float] | np.ndarray | None = None,  # Phase B advisory
    config: GudieConfig | None = None,
    fca_cache: dict | None = None,
) -> GudieResult:
    """
    Full GUDIE stage on the FASTCORE boundary. Returns H (induced subgraph of G_F,
    original MFG edge directions) plus diagnostics.

    Parameters
    ----------
    G_F           : FASTCORE-induced MFG subgraph (nx.DiGraph, 'mass_flow' edges).
    target_id     : EX_M reaction id (must be a node of G_F).
    source_ids    : optional substrate/uptake reaction ids for dual-anchor localization.
                    When given, a second PPR seeded here is propagated *downstream* (on
                    A_rev — the operator that moves *along* MFG edges) and multiplied,
                    field-wise after peak-normalisation, with the target's upstream field.
                    This selects the source->target corridor and preserves uptake
                    reactions. None => identical to the single-seed implementation.
    model         : cobra.Model for FCA coupling; omit -> edge-PPR LUDIE fallback.
    coupling      : optional precomputed c_i (overrides in-module FCA).
    embedding_sim : optional cosine sim(z_i, z_EX_M) for Phase B VUDIE blend (lam>0).
    """
    cfg = config or GudieConfig()
    if target_id not in G_F:
        raise KeyError(f"target {target_id!r} not in FASTCORE boundary G_F")

    node_order: list[str] = list(G_F.nodes())
    node_index = {r: i for i, r in enumerate(node_order)}
    n = len(node_order)
    seed_idx = node_index[target_id]

    # Step 2.1 — compressed adjacency (original + reversed for upstream ops)
    A_hat = _compressed_adjacency(G_F, node_index, n, cfg.compress_eps)
    A_rev = A_hat.T.tocsr()                                   # G^R: row=source, entry=w_hat orig->

    # Step 2.4/2.5 — VUDIE = PPR (with Katz fallback on non-convergence)
    p, converged = personalized_pagerank(
        A_hat, seed_idx, cfg.alpha_r, cfg.ppr_normalization, cfg.ppr_tol, cfg.ppr_max_iter)

    # Dual-anchor: second PPR seeded at source_ids, propagated DOWNSTREAM (A_rev), then
    # combined with the target field via product of peak-normalised scores. This selects
    # the source->target corridor; None => single-seed behaviour below is untouched.
    interest = p                                             # VUDIE field feeding I0/edge interest
    src_idx: list[int] = []
    source_ppr: np.ndarray | None = None
    if source_ids:
        src_idx = [node_index[s] for s in source_ids if s in node_index]
        missing = [s for s in source_ids if s not in node_index]
        if missing:
            logger.warning("GUDIE dual-anchor: source_ids not in G_F, ignored: %s", missing)
        if src_idx:
            source_ppr, _ = personalized_pagerank(
                A_rev, src_idx, cfg.alpha_r, cfg.ppr_normalization, cfg.ppr_tol, cfg.ppr_max_iter)
            interest = _norm_max(p) * _norm_max(source_ppr)

    # Step 2.2 — FCA coupling c_i
    c_vec: np.ndarray | None = None
    if coupling is not None:
        c_vec = np.array([float(coupling.get(r, 0.0)) for r in node_order])
    elif cfg.use_fca and model is not None:
        subset = node_order if cfg.fca_max_reactions is None else node_order[:cfg.fca_max_reactions]
        c_map = flux_coupling_to_target(
            model, target_id, subset,
            gamma=cfg.fca_gamma, pure_stoichiometric=cfg.fca_pure_stoichiometric,
            zero_span_to_zero=cfg.fca_zero_span_to_zero, cache=fca_cache)
        if c_map is not None:
            c_vec = np.array([float(c_map.get(r, 0.0)) for r in node_order])

    # Step 2.3 — LUDIE edge interest (geometric fusion, or edge-PPR fallback)
    B, ludie_fallback = _edge_interest_matrix(A_rev, c_vec, interest)

    # VUDIE init, with optional Phase B blend (lam=0 => exact Phase A)
    I0 = interest.copy()
    lam = float(np.clip(cfg.lam, 0.0, 1.0))
    if lam > 0.0 and embedding_sim is not None:
        if isinstance(embedding_sim, np.ndarray):
            sim = embedding_sim.astype(np.float64)
        else:
            sim = np.array([float(embedding_sim.get(r, 0.0)) for r in node_order])
        sim = np.clip(sim, 0.0, None)                        # negative cosine -> 0
        I0 = (1.0 - lam) * interest + lam * sim

    # Step 2.6 — interest propagation
    I_h = _propagate(B, A_rev, I0, cfg.hops)

    # Step 2.7 — seed expansion (restore original edge directions in H)
    ref = _threshold_reference(A_rev, I_h, seed_idx)
    delta = cfg.k * ref
    included_idx = _seed_expansion(A_rev, I_h, seed_idx, delta, cfg.decay)
    if source_ids:                                           # force-include sources + target
        included_idx |= set(src_idx)
        included_idx.add(seed_idx)
    included_nodes = [node_order[i] for i in sorted(included_idx)]
    H = G_F.subgraph(included_nodes).copy()

    if source_ids:                                           # advisory corridor check (warn only)
        for s in source_ids:
            if s in node_index and not nx.has_path(H, s, target_id):
                logger.warning("GUDIE dual-anchor: no path %s -> %s within extracted H",
                               s, target_id)

    logger.info("GUDIE: |F|=%d -> |H|=%d (target=%s, ludie_fallback=%s, ppr_converged=%s)",
                n, H.number_of_nodes(), target_id, ludie_fallback, converged)

    return GudieResult(
        subgraph=H,
        nodes=included_nodes,
        node_order=node_order,
        ppr=p,
        coupling=c_vec,
        propagated_interest=I_h,
        seed_interest=float(I_h[seed_idx]),
        threshold=float(delta),
        ppr_converged=converged,
        ludie_fallback=ludie_fallback,
        source_ppr=source_ppr,
    )