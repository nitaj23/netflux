"""
scripts/validate_gudie.py

DIAGNOSTIC / VALIDATION for the GUDIE localized-subgraph extractor
(src/extraction/gudie.py). Single-run tool that feeds run_gudie a G_F
built from fastcc.py's flux-consistent set and reports:

  1. the SIZE of the extracted subgraph H          (nodes + edges, |H|/|G_F|)
  2. the NUMBER OF REACTIONS in it                  (len(GudieResult.nodes))
  3. a MAP of the subgraph                          (node/edge JSON of H, enriched
                                                      with per-reaction GUDIE
                                                      diagnostics, for the map frontend)

Pipeline (consumes fastcc.py as G_F):
    load GEM
      -> fastcc_from_cobra              : flux-consistent reaction universe
      -> restrict model to it
      -> pick a producible M            : seeded RNG, demand-reaction gated
                                          (--metabolite pins a specific one)
      -> add DM_M, maximise, FBA
      -> build_mfg -> active subgraph   : G_F  (the "graph from fastcc")
      -> run_gudie(seed = DM_M)         : localized subgraph H around M

Run
---
    .venv/Scripts/python.exe -u -m scripts.validate_gudie --model iNF517 --seed 0
    .venv/Scripts/python.exe -u -m scripts.validate_gudie --model iML1515 --metabolite succ_c
    .venv/Scripts/python.exe -u -m scripts.validate_gudie --model iML1515 --alpha 0.3 --k 0.1 --hops 8
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import tempfile
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cobra

from src.extraction.fastcc import (
    fastcc_from_cobra, load_fastcc_model, fastcc_metadata, available_models)
from src.extraction.gudie import run_gudie, GudieConfig
from src.graph_ops import build_mfg, fpe_similarity_to_target
from scripts.fastcc_all_organisms import load_catalogue, max_demand_flux

FASTCC_EPS = 1e-4
EPS = 1e-6


def resolve_path(model_key: str, catalogue: dict) -> Path:
    """Map a metadata key (e.g. 'iML1515') to its GEM file path."""
    return PROJECT_ROOT / catalogue[model_key]["file_path"]


def restrict_to_consistent(model, consistent_ids: set[str]):
    """Delete blocked reactions -> a model whose reactions are exactly the
    FASTCC-consistent set (the 'graph from fastcc')."""
    reduced = model.copy()
    reduced.remove_reactions([r for r in reduced.reactions
                              if r.id not in consistent_ids], remove_orphans=True)
    return reduced


def load_reduced_model(key: str, catalogue: dict, eps: float, use_cache: bool):
    """Return (reduced_model, n_full).

    Default: load the pre-built FASTCC-reduced model from
    extraction/fastcc/models/<key>.json (drop-in reduced GEM, computed once at
    eps=1e-6 with growth preservation) -- no LP solves at run time. n_full comes
    from the cache metadata's num_reactions_full.

    use_cache=False (or no cached model for `key`): fall back to loading the full
    GEM and running fastcc live at `eps`, exactly as before."""
    if use_cache and key in available_models():
        reduced = load_fastcc_model(key)
        n_full = int(fastcc_metadata().get(key, {}).get(
            "num_reactions_full", len(reduced.reactions)))
        return reduced, n_full
    model = cobra.io.load_json_model(str(resolve_path(key, catalogue)))
    n_full = len(model.reactions)
    _, consistent_ids = fastcc_from_cobra(model, eps)
    return restrict_to_consistent(model, set(consistent_ids)), n_full


def pick_producible(model, rng, eps, max_tries=40):
    """Random metabolite whose demand reaction (M ->) can carry flux (producible).
    Replaces the old currency/cofactor-filtered pick_proper; mirrors the sibling
    fastcc+GUDIE scripts."""
    mets = list(model.metabolites)
    rng.shuffle(mets)
    for m in mets[:max_tries]:
        if abs(max_demand_flux(model, m.id, eps)) > eps:
            return m
    return None


def check(label: str, ok: bool, detail: str = "") -> bool:
    """PASS/FAIL sanity line."""
    print(f"   [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
    return ok


def build_map_payload(model, sol, G_F, res, dm_id, meta) -> dict:
    """Node/edge JSON map of H, enriched with per-reaction GUDIE diagnostics.
    Arrays in GudieResult are aligned to res.node_order -> index by reaction id."""
    idx = {rid: i for i, rid in enumerate(res.node_order)}
    H = res.subgraph
    nodes = []
    for rid in res.nodes:
        r = model.reactions.get_by_id(rid)
        i = idx.get(rid)
        nodes.append({
            "id": rid,
            "name": r.name,
            "flux": float(sol.fluxes.get(rid, 0.0)),
            "subsystem": r.subsystem or "",
            "reaction": r.reaction,
            "ppr": float(res.ppr[i]) if i is not None else None,
            "interest": float(res.propagated_interest[i]) if i is not None else None,
            "coupling": (float(res.coupling[i]) if (res.coupling is not None and i is not None)
                         else None),
            "is_target": rid == dm_id,
        })
    edges = [{"source": u, "target": v, "mass_flow": float(d["mass_flow"])}
             for u, v, d in H.edges(data=True)]
    return {
        "model": meta["model"], "compound": meta["compound"], "target_node": dm_id,
        "vmax": meta["vmax"], "lam": meta.get("lam", 0.0),
        "flux_capture": meta.get("flux_capture"),
        "full": meta["full"], "consistent": meta["consistent"],
        "mfg_active": G_F.number_of_nodes(), "H": H.number_of_nodes(),
        "ppr_converged": res.ppr_converged, "ludie_fallback": res.ludie_fallback,
        "threshold": float(res.threshold), "seed_interest": float(res.seed_interest),
        "nodes": nodes, "edges": edges,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Validate GUDIE: size, reaction count, and map of the extracted subgraph H.")
    ap.add_argument("--model", default="iML1515", help="GEM metadata key.")
    ap.add_argument("--metabolite", default=None,
                    help="pin a specific target metabolite bigg_id (default: random proper).")
    ap.add_argument("--source", default=None, nargs="+", metavar="RXN_ID",
                    help="dual-anchor: substrate/uptake reaction id(s) to seed the "
                         "downstream field (e.g. EX_glc__D_e). Omit for single-seed.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epsilon", type=float, default=FASTCC_EPS)
    ap.add_argument("--alpha", type=float, default=0.3, help="GUDIE PPR restart alpha_r.")
    ap.add_argument("--k", type=float, default=0.1, help="GUDIE inclusion threshold factor.")
    ap.add_argument("--hops", type=int, default=8, help="GUDIE interest propagation hops.")
    ap.add_argument("--lam", type=float, default=0.05,
                    help="Phase B FPE blend coefficient in [0,1] (default 0.05). >0 mixes the target "
                         "Flow-Profile-Encoding cosine into VUDIE interest (lam=0 => exact Phase A).")
    ap.add_argument("--fpe-raw", action="store_true",
                    help="use raw non-negative FPE cosine (default z-scores each FPE channel first).")
    ap.add_argument("--no-cache", action="store_true",
                    help="ignore the pre-built FASTCC-reduced model; run fastcc live at --epsilon.")
    ap.add_argument("--out", default=None, help="JSON map output path (default: temp dir).")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    cfg = GudieConfig(alpha_r=args.alpha, k=args.k, hops=args.hops, lam=args.lam)
    catalogue = load_catalogue()
    if args.model not in catalogue:
        raise SystemExit(f"unknown model {args.model!r}; keys: {sorted(catalogue)}")
    info = catalogue[args.model]
    t0 = time.time()

    print("=" * 84)
    print(f"VALIDATE GUDIE  |  {args.model}  {info['organism']}  |  seed={args.seed}  "
          f"alpha_r={args.alpha} k={args.k} hops={args.hops} lam={args.lam}")
    print("=" * 84)

    # 1) FASTCC-consistent universe -> restricted model (cached reduced GEM by
    #    default; --no-cache runs fastcc live). reduced.reactions IS the consistent set.
    reduced, n_full = load_reduced_model(args.model, catalogue, args.epsilon,
                                         use_cache=not args.no_cache)
    consistent = {r.id for r in reduced.reactions}
    src = "fastcc live" if args.no_cache else "cached reduced GEM"
    print(f"   reduced model: {src}  ({len(consistent)} consistent rxns)")

    # 2) pick target metabolite (pinned, or random producible)
    if args.metabolite is not None:
        if args.metabolite not in reduced.metabolites:
            raise SystemExit(f"'{args.metabolite}' not in fastcc-consistent {args.model}")
        met = reduced.metabolites.get_by_id(args.metabolite)
        if abs(max_demand_flux(reduced, met.id, EPS)) <= EPS:
            print(f"   note: {met.id} is not producible via a demand reaction (pinned anyway)")
    else:
        met = pick_producible(reduced, rng, EPS)
        if met is None:
            raise SystemExit(f"no producible metabolite found in {args.model}")

    # 3) add demand for M, maximise it, build the MFG on that flux state
    dm = reduced.add_boundary(reduced.metabolites.get_by_id(met.id), type="demand")
    reduced.objective = dm
    sol = reduced.optimize()
    if sol.status != "optimal":
        raise SystemExit(f"target {met.id} infeasible (status={sol.status})")
    vmax = float(sol.objective_value)

    G = build_mfg(reduced, sol)
    active = [nid for nid in G.nodes if G.degree(nid) > 0]
    G_F = G.subgraph(active).copy()
    if dm.id not in G_F:
        raise SystemExit(f"target {met.id} (DM {dm.id}) not in active MFG G_F")

    # 4) GUDIE localization on the fastcc-derived G_F
    if args.source:
        missing = [s for s in args.source if s not in G_F]
        if missing:
            print(f"   note: source(s) not in active G_F, will be ignored: {missing}")
        print(f"   dual-anchor sources: {args.source}")

    # Phase B advisory signal: cosine of each node's Flow-Profile-Encoding to the
    # target's FPE. Only consulted when lam>0 (else run_gudie ignores it).
    embedding_sim = None
    if args.lam > 0.0:
        embedding_sim = fpe_similarity_to_target(
            G_F, dm.id, list(G_F.nodes()), zscore=not args.fpe_raw)

    res = run_gudie(G_F, dm.id, source_ids=args.source, model=reduced,
                    config=cfg, embedding_sim=embedding_sim)
    H = res.subgraph
    mfg_n = G_F.number_of_nodes()
    h_n, h_e = H.number_of_nodes(), H.number_of_edges()
    C = met.elements.get("C", 0) if met.formula else 0

    # Flux capture: share of the active graph's total |flux| retained by H. A good
    # localization keeps the high-flux backbone leading to the product while staying
    # compact, so flux_capture (recall proxy) should stay high even as |H|/|G_F| drops.
    flux_gf = sum(abs(sol.fluxes.get(n, 0.0)) for n in G_F.nodes)
    flux_h = sum(abs(sol.fluxes.get(n, 0.0)) for n in res.nodes)
    flux_capture = flux_h / flux_gf if flux_gf > 0 else 0.0

    # --- report -----------------------------------------------------------
    print(f"\ntarget compound : {met.id}  ({met.name})  C={C}  formula={met.formula}")
    print(f"target node     : {dm.id}   max production vmax = {vmax:.4f}")
    print(f"reactions       : full={n_full}  fastcc-consistent={len(consistent)}  "
          f"MFG active |G_F|={mfg_n}")
    print("-" * 84)
    print("EXTRACTED SUBGRAPH H")
    print(f"   size          : {h_n} reactions (nodes),  {h_e} edges")
    print(f"   fraction      : |H|/|G_F| = {100*h_n/max(mfg_n,1):.1f}%")
    print(f"   flux capture  : {100*flux_capture:.1f}% of active |flux|  "
          f"(lam={args.lam}, fpe={'raw' if args.fpe_raw else 'zscore' if args.lam>0 else 'off'})")
    print(f"   reactions     : len(res.nodes) = {len(res.nodes)}")
    print(f"   ppr_converged : {res.ppr_converged}   ludie_fallback : {res.ludie_fallback}   "
          f"fca_coupling : {res.coupling is not None}")
    print(f"   threshold     : delta = {res.threshold:.3e}   seed_interest = {res.seed_interest:.3e}")

    print("\nflux funnel - H reactions by |flux| (high -> low):")
    for rid in sorted(res.nodes, key=lambda r: -abs(sol.fluxes.get(r, 0.0))):
        r = reduced.reactions.get_by_id(rid)
        tag = "  <- TARGET" if rid == dm.id else ""
        print(f"   {sol.fluxes.get(rid, 0.0):10.3f}  {rid:16} "
              f"{(r.subsystem or '')[:22]:22} {r.name[:30]:30}{tag}")

    # --- sanity assertions ------------------------------------------------
    print("\nsanity checks:")
    gf_nodes, h_nodes = set(G_F.nodes), set(H.nodes)
    edges_ok = all(G_F.has_edge(u, v) for u, v in H.edges())
    all_ok = all([
        check("target in G_F", dm.id in G_F),
        check("target in H", dm.id in h_nodes),
        check("H nodes subset of G_F", h_nodes <= gf_nodes,
              f"{len(h_nodes - gf_nodes)} stray"),
        check("H edges are G_F edges", edges_ok),
        check("|H| <= |G_F|", h_n <= mfg_n, f"{h_n} <= {mfg_n}"),
        check("res.nodes matches H", len(res.nodes) == h_n),
        check("H non-trivial (more than just target)", h_n > 1, f"|H|={h_n}"),
    ])

    # --- map export -------------------------------------------------------
    meta = {"model": args.model, "compound": met.id, "vmax": vmax,
            "full": n_full, "consistent": len(consistent),
            "lam": args.lam, "flux_capture": flux_capture}
    payload = build_map_payload(reduced, sol, G_F, res, dm.id, meta)
    if args.out:
        out = Path(args.out)
    else:
        out = Path(tempfile.gettempdir()) / "mfghub_gudie" / \
            f"gudie_subgraph_{args.model}_{met.id}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))

    print(f"\nmap: nodes={len(payload['nodes'])} edges={len(payload['edges'])}  ->  {out}")
    print(f"\n{'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED'}  [{time.time()-t0:.1f}s]")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
