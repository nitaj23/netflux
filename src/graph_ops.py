# src/graph_ops.py
import io
import json
import logging
import warnings
import cobra
import networkx as nx
import numpy as np
import torch
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from torch_geometric.data import Data

logger = logging.getLogger(__name__)

# 2. FBA execution


def run_fba(model: cobra.Model) -> cobra.Solution | None:
    """
    Thin wrapper around model.optimize().
    Returns None on infeasibility so callers raise HTTP 400 without leaking
    solver internals. Orth et al. (2010, Nature Biotechnology).
    """

    solution = model.optimize()
    if solution.status != "optimal":
        logger.warning("FBA returned non-optimal status '%s'", solution.status)
        return None
    return solution


 
# 3. MFG construction
 

def build_mfg(model: cobra.Model, solution: cobra.Solution) -> nx.DiGraph:
    """
    Constructs a Mass Flow Graph via proportional attribution in 2m-dimensional
    space. Edge weight A[i,j] = total mass flux from reaction i to reaction j,
    proportionally attributed by each mediating metabolite's fractional share
    of total production.

    References: Beguerisse-Díaz et al. (NFG/MFG, Eq. 2);
                Islam et al. (2023, Bioinformatics) Section 2.1.
    """

    #  1. Core matrices
    S_dense          = cobra.util.create_stoichiometric_matrix(model)
    S                = sp.csc_matrix(S_dense)  # [n_met, n_rxn]
    rxn_ids          = [r.id for r in model.reactions]
    met_ids          = [m.id for m in model.metabolites]
    num_metabolites, num_reactions = S.shape
    v                = np.array([solution.fluxes.get(rid, 0.0) for rid in rxn_ids])
    reversible_bool  = np.array([1 if r.lower_bound < 0 else 0 for r in model.reactions])

    #  2. 2m-dimensional unfolding
    v_abs = np.abs(v)
    v_2m  = 0.5 * np.concatenate([v_abs + v, v_abs - v])

    I_m    = sp.eye(num_reactions, format='csc')
    diag_r = sp.diags(reversible_bool, format='csc')

    block_diagonal = sp.bmat([
        [I_m,  None],
        [None, diag_r],
    ], format='csc')

    S_hstack = sp.hstack([S, -S], format='csc')

    S_2m     = S_hstack @ block_diagonal
    # Safe sparse thresholding for positive components
    S_2m_pos = S_2m.copy()
    S_2m_pos.data[S_2m_pos.data < 0] = 0
    S_2m_pos.eliminate_zeros()

    # Safe sparse thresholding for negative components
    S_2m_neg = S_2m.copy()
    S_2m_neg.data[S_2m_neg.data > 0] = 0
    S_2m_neg.data = np.abs(S_2m_neg.data)
    S_2m_neg.eliminate_zeros()

    # 3. Adjacency matrix
    V_star   = sp.diags(v_2m)
    j_vector = S_2m_pos @ v_2m

    near_zero = int(np.sum((j_vector > 0) & (j_vector <= 1e-9)))
    if near_zero:
        warnings.warn(
            f"{near_zero} metabolite(s) with 0 < total_production_flow ≤ 1e-9 "
            f"silently excluded from MFG edge set.",
            RuntimeWarning, stacklevel=2,
        )

    j_vector_inv = np.divide(
        1.0, j_vector,
        out=np.zeros_like(j_vector),
        where=j_vector > 1e-9,
    )
    J_pinv     = sp.diags(j_vector_inv)

    PS = S_2m_pos @ V_star   # [n_met, 2*n_rxn] — production weighted
    CS = S_2m_neg @ V_star   # [n_met, 2*n_rxn] — consumption weighted

    A_unfolded = PS.T @ J_pinv @ CS

    A = (A_unfolded[:num_reactions, :num_reactions] +
         A_unfolded[:num_reactions, num_reactions:] +
         A_unfolded[num_reactions:, :num_reactions] +
         A_unfolded[num_reactions:, num_reactions:])

    # 4. NetworkX graph
    # Enumerate model.reactions directly: rxn_ids and model.reactions share
    # the same order, so no .get_by_id() lookup is needed here.
    G = nx.DiGraph()
    for i, rid in enumerate(rxn_ids):
        rxn = model.reactions.get_by_id(rid)
        G.add_node(
            rid,
            reaction_id=rid,
            reaction_name=rxn.name,
            flux=float(v[i]),
            subsystem=rxn.subsystem or "Unknown",
            lower_bound=rxn.lower_bound,
            upper_bound=rxn.upper_bound,
            gene_rule=rxn.gene_reaction_rule or "",
            # sorted(): rxn.genes is a frozenset whose iteration order varies
            # across processes (hash randomization), which would make the
            # serialized GraphML non-reproducible between rebuilds.
            genes=sorted(g.id for g in rxn.genes),
        )

    A_coo = A.tocoo()
    for i, j, mass_flow in zip(A_coo.row, A_coo.col, A_coo.data):
        if mass_flow > 1e-9:
            G.add_edge(rxn_ids[i], rxn_ids[j], mass_flow=float(mass_flow), metabolites=[])

    # 4b. Flow Profile Encoding (FlowGAT Eqs. 4–6) as node attributes.
    # Reuse the adjacency A already computed (rows/cols in rxn_ids order);
    # thresholded the same way as the edge set so attr and edges agree.
    A_fpe = A.tocsr(copy=True)
    A_fpe.data[A_fpe.data <= 1e-9] = 0.0
    A_fpe.eliminate_zeros()
    fpe = _fpe_from_adjacency(A_fpe, FPE_KM, FPE_ALPHA)   # [n_rxn, 2*FPE_KM]
    for i, rid in enumerate(rxn_ids):
        # Plain Python list so it survives GraphML / torch.save round-trips.
        G.nodes[rid]["fpe"] = [float(x) for x in fpe[i]]

    # 5. Metabolite annotation
    P_matrix = PS[:, :num_reactions] + PS[:, num_reactions:]
    C_matrix = CS[:, :num_reactions] + CS[:, num_reactions:]
    P_bool = P_matrix > 1e-9
    C_bool = C_matrix > 1e-9

    for k, met_id in enumerate(met_ids):
        for p in np.where(P_bool[k, :].toarray())[1]:
            for c in np.where(C_bool[k, :].toarray())[1]:
                if G.has_edge(rxn_ids[p], rxn_ids[c]):
                    G[rxn_ids[p]][rxn_ids[c]]["metabolites"].append(met_id)

    return G

 
# Node feature contract — the FPE (Flow Profile Embedding)
# Index 0 is flux — relied on by compute.py Strategy 3 (file-upload):
#     flux_vector = {name: float(pyg_data.x[i, 0]) for i, name in enumerate(node_names)}
# Store NODE_FEATURE_NAMES verbatim in gat_models.node_feature_schema when
# registering a GAT model.


BASE_FEATURE_NAMES: list[str] = [
    "flux",          # 0 — net FBA flux (signed)
    "in_degree",     # 1 — MFG in-degree (unweighted count)
    "out_degree",    # 2 — MFG out-degree (unweighted count)
    "is_reversible", # 3 — 1.0 if lower_bound < 0 else 0.0
]

# Flow Profile Encoding (FPE) defaults.
# FlowGAT (Islam et al. 2023, npj Syst Biol Appl 10:24), Eqs. (4)–(6),
# after Cooper & Barahona (flow profiles for directed graphs).
#   inflow^k_i  = [A^k 1]_i        (Eq. 4)
#   outflow^k_i = [(A^T)^k 1]_i     (Eq. 5)
#   FPE_i = [ β^1·inflow^1, …, β^km·inflow^km,
#            β^1·outflow^1, …, β^km·outflow^km ]   (Eq. 6)
#   β = α / λ1,  λ1 = largest-magnitude eigenvalue of A.
# At k=1, inflow/outflow equal the *weighted* in/out-degree — the higher-order
# generalisation of base features 1–2 (which remain unweighted counts).
FPE_KM: int = 3       # km — maximum path length
FPE_ALPHA: float = 1.0  # α — short/long-path weighting


def _fpe_feature_names(km: int) -> list[str]:
    return (
        [f"inflow_{k}" for k in range(1, km + 1)]
        + [f"outflow_{k}" for k in range(1, km + 1)]
    )

NODE_FEATURE_NAMES: list[str] = BASE_FEATURE_NAMES + _fpe_feature_names(FPE_KM)
NODE_FEATURE_DIM: int = len(NODE_FEATURE_NAMES)  # 4+ 2*FPE_KM


def _fpe_from_adjacency(A: "sp.spmatrix | np.ndarray", km: int, alpha: float) -> np.ndarray:
    """
    Core FPE kernel (FlowGAT Eqs. 4–6) on a dense adjacency A, rows/cols in
    the caller's chosen node order. Returns [n, 2*km]: km inflow columns then
    km outflow columns. Column k-1 holds the length-k flow scaled by β^k.

    β = alpha / λ1, λ1 = spectral radius of A; λ1 ≤ 1e-12 (e.g. a nilpotent
    pure-DAG adjacency) → β = 0, zeroing the FPE block only.
    """

    if not sp.issparse(A):
        A = sp.csr_matrix(A)

    n = A.shape[0]
    if n == 0:
        return np.zeros((0, 2 * km), dtype=np.float32)

    lam1 = 0.0
    try:
        if n < 3:
            # Dense fallback for exceptionally tiny/degenerate graphs
            if A.nnz > 0:
                lam1 = float(np.max(np.abs(np.linalg.eigvals(A.toarray()))))
        elif A.nnz > 0:
            # Sparse extraction of largest magnitude eigenvalue
            evals = spla.eigs(A, k=1, which='LM', return_eigenvectors=False)
            lam1 = float(np.abs(evals[0]))
    except (spla.ArpackNoConvergence, np.linalg.LinAlgError):
        lam1 = 0.0
    beta = (alpha / lam1) if lam1 > 1e-12 else 0.0

    out = np.zeros((n, 2 * km), dtype=np.float64)
    in_vec = np.ones(n, dtype=np.float64)
    out_vec = np.ones(n, dtype=np.float64)
    AT = A.T.tocsr()
    for k in range(1, km + 1):
        in_vec = A @ in_vec
        out_vec = AT @ out_vec
        scale = beta ** k
        out[:, k - 1] = scale * in_vec           # inflow_k
        out[:, km + k - 1] = scale * out_vec     # outflow_k
    return out.astype(np.float32)


def compute_fpe(
    G: nx.DiGraph,
    node_names: list[str],
    km: int = FPE_KM,
    alpha: float = FPE_ALPHA,
) -> np.ndarray:
    """
    Flow Profile Encoding per FlowGAT Eqs. (4)–(6).

    Returns an [n, 2*km] array: km inflow columns followed by km outflow
    columns, row-aligned to `node_names`. Builds the adjacency from the
    'mass_flow' edge weights in `node_names` order, then delegates to
    _fpe_from_adjacency. Used when no prebuilt adjacency is available;
    build_mfg passes its own matrix to _fpe_from_adjacency directly.
    """

    n = len(node_names)
    if n == 0:
        return np.zeros((0, 2 * km), dtype=np.float32)
    idx = {name: i for i, name in enumerate(node_names)}
    row = []
    col = []
    data = []
    for u, v, d in G.edges(data=True):
        row.append(idx[u])
        col.append(idx[v])
        data.append(float(d.get("mass_flow", 0.0)))
        
    A = sp.coo_matrix((data, (row, col)), shape=(n, n)).tocsr()

    return _fpe_from_adjacency(A, km, alpha)


def to_pyg(G: nx.DiGraph) -> Data:
    """
    Convert a Mass Flow Graph (NetworkX DiGraph) to a PyG Data object.

    node_names is stored as a plain list[str] — never converted to a tensor.
    PyG serialises non-tensor attributes as Python objects across torch.save/load.
    Converting to a tensor would silently break GAT seed resolution and
    Strategy 3 flux_vector reconstruction in compute.py.
    """

    node_names: list[str] = list(G.nodes)
    node_index: dict[str, int] = {n: i for i, n in enumerate(node_names)}
    n = len(node_names)

    x = torch.zeros((n, NODE_FEATURE_DIM), dtype=torch.float)
    n_base = len(BASE_FEATURE_NAMES)
    for i, node in enumerate(node_names):
        d  = G.nodes[node]
        lb = float(d.get("lower_bound", 0.0))
        x[i, :n_base] = torch.tensor([
            float(d.get("flux", 0.0)),
            float(G.in_degree(node)),
            float(G.out_degree(node)),
            1.0 if lb < 0 else 0.0,
        ])

    # FPE block (Eqs. 4–6): prefer the attribute stored by build_mfg; only
    # recompute if absent (e.g. a graph loaded without the 'fpe' attr).
    if n:
        first = G.nodes[node_names[0]]
        if "fpe" in first and len(first["fpe"]) == 2 * FPE_KM:
            for i, node in enumerate(node_names):
                x[i, n_base:] = torch.tensor(G.nodes[node]["fpe"], dtype=torch.float)
        else:
            fpe = compute_fpe(G, node_names)            # [n, 2*FPE_KM]
            x[:, n_base:] = torch.from_numpy(fpe)

    edges = list(G.edges(data=True))
    if edges:
        src        = torch.tensor([node_index[u] for u, _, _ in edges], dtype=torch.long)
        dst        = torch.tensor([node_index[v] for _, v, _ in edges], dtype=torch.long)
        edge_index = torch.stack([src, dst], dim=0)                           # [2, E]
        edge_attr  = torch.tensor(
            [[float(d.get("mass_flow", 0.0))] for _, _, d in edges],
            dtype=torch.float,
        )                                                                      # [E, 1]
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_attr  = torch.zeros((0, 1), dtype=torch.float)

    data            = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    data.node_names = node_names   # plain list — do NOT convert to tensor
    data.num_nodes  = n
    return data


def _graphml_ready(G: nx.DiGraph) -> nx.DiGraph:
    """
    Copy of G with list/dict node & edge attributes JSON-encoded to strings.

    GraphML only supports scalar attribute values, so the MFG's list-valued
    attributes (node 'genes', node 'fpe', edge 'metabolites') are stringified.
    The copy is shallow-per-attr-dict: each attr dict is fresh, so reassigning a
    key does not mutate the original G — the in-memory graph keeps its native
    Python lists for to_pyg() and index_gem_and_mfg(). Decode on read with
    json.loads().
    """
    H = G.copy()
    for _, d in H.nodes(data=True):
        for k, val in list(d.items()):
            if isinstance(val, (list, dict)):
                d[k] = json.dumps(val)
    for _, _, d in H.edges(data=True):
        for k, val in list(d.items()):
            if isinstance(val, (list, dict)):
                d[k] = json.dumps(val)
    return H


def write_graphml(G: nx.DiGraph, path: str) -> None:
    """GraphML export to a file `path` (see _graphml_ready for attr encoding)."""
    nx.write_graphml(_graphml_ready(G), path)


def graphml_bytes(G: nx.DiGraph) -> bytes:
    """
    GraphML for `G` as an in-memory UTF-8 bytes object — byte-identical to what
    write_graphml() would put on disk, but with no file involved. Hand straight
    to an HTTP download response (media type application/graphml+xml) or persist
    as a blob. Same list/dict attr encoding as write_graphml (decode with
    json.loads on read).
    """
    buf = io.BytesIO()
    nx.write_graphml(_graphml_ready(G), buf)
    return buf.getvalue()


def graph_json(G: nx.DiGraph, *, as_str: bool = False) -> dict | str:
    """
    Node-link JSON for `G`: a plain dict (as_str=False, default) or a JSON string
    (as_str=True) for an in-memory downloadable object (media type
    application/json). Round-trips via nx.node_link_graph(data, edges="edges").

    The MFG's list attributes (node 'genes'/'fpe', edge 'metabolites') and scalar
    floats/strings are all natively JSON-serialisable, so — unlike GraphML — no
    stringification is needed and the structure survives intact.
    """
    data = nx.node_link_data(G, edges="edges")
    return json.dumps(data) if as_str else data


def fpe_similarity_to_target(
    G: nx.DiGraph,
    target_id: str,
    node_order: list[str] | None = None,
    *,
    zscore: bool = True,
) -> dict[str, float]:
    """
    Cosine similarity of every node's Flow Profile Encoding to the target's FPE
    (the advisory signal for GUDIE's Phase B blend, gudie.run_gudie(embedding_sim=)).

    Reads the per-node ``fpe`` attribute stored by build_mfg (FlowGAT Eqs. 4-6:
    3 inflow + 3 outflow flow-profile channels). A reaction whose local mass-flow
    profile *looks like* the target reaction's is nudged up in interest, letting the
    extractor prefer producers whose flow signature matches the product's demand.

    zscore=True (default) standardises each FPE channel across nodes before the
    cosine, so the raw inflow_k scale (which grows ~beta^k * lambda^k with path
    length) can't let one channel dominate; it also lets similarity go negative
    (opposite profile), which run_gudie clips to 0. zscore=False uses the raw
    non-negative FPE (cosine in [0,1]) -- every node then shares a positive floor
    of similarity, a weaker discriminator.

    Returns {reaction_id: cosine in [-1, 1]}. Missing/degenerate FPE -> 0.0.
    Target absent or zero-norm -> all zeros (blend degenerates to pure Phase A).
    """
    names = list(node_order) if node_order is not None else list(G.nodes())
    dim = 2 * FPE_KM
    X = np.zeros((len(names), dim), dtype=np.float64)
    for i, rid in enumerate(names):
        f = G.nodes[rid].get("fpe")
        if f is not None and len(f) == dim:
            X[i] = f
    if zscore:
        X = fpe_zscored(X).astype(np.float64)

    if target_id not in names:
        return {rid: 0.0 for rid in names}
    t = names.index(target_id)
    tv = X[t]
    tnorm = float(np.linalg.norm(tv))
    if tnorm <= 1e-12:
        return {rid: 0.0 for rid in names}

    row_norms = np.linalg.norm(X, axis=1)
    denom = row_norms * tnorm
    sims = np.zeros(len(names), dtype=np.float64)
    nz = denom > 1e-12
    sims[nz] = (X[nz] @ tv) / denom[nz]
    return {names[i]: float(sims[i]) for i in range(len(names))}


def fpe_zscored(x: "np.ndarray | torch.Tensor") -> np.ndarray:
    
    arr = x.detach().cpu().numpy() if isinstance(x, torch.Tensor) else np.asarray(x)
    arr = arr.astype(np.float64)
    if arr.shape[0] == 0:
        return arr.astype(np.float32)
    mean = arr.mean(axis=0)
    std = arr.std(axis=0)
    # Avoid divide-by-zero on constant dimensions: those map to 0 (mean-centered).
    safe = np.where(std > 1e-12, std, 1.0)
    z = (arr - mean) / safe
    z[:, std <= 1e-12] = 0.0
    return z.astype(np.float32)
