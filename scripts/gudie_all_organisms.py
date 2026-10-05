"""
scripts/gudie_all_organisms.py

END-TO-END: for every organism, run the full localization pipeline on the
FASTCC-consistent network and a random (producible) compound of interest:

    load GEM
      -> fastcc                     : flux-consistent universe
      -> restrict model to it       : "the graph from fastcc"
      -> pick random producible M   : demand-reaction gate (max DM_M > eps)
      -> add DM_M, maximise it, FBA
      -> build_mfg                  : mass-flow graph G_F (active subgraph)
      -> run_gudie(seed = DM_M)     : localized subgraph H around M

Reports |consistent|, |MFG active|, |H|, PPR convergence, LUDIE fallback, and --
answering the earlier open question -- how much of central metabolism GUDIE
retains in H (central_in_H / central_in_active).

Run
---
    .venv/Scripts/python.exe -u -m scripts.gudie_all_organisms                 # all organisms
    .venv/Scripts/python.exe -u -m scripts.gudie_all_organisms --model iML1515 # one
    .venv/Scripts/python.exe -u -m scripts.gudie_all_organisms --alpha 0.3 --k 0.1 --seed 2
"""

from __future__ import annotations

import argparse
import random
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cobra

from src.extraction.fastcc import (
    fastcc_from_cobra, load_fastcc_model, available_models)
from src.extraction.gudie import run_gudie, GudieConfig
from src.graph_ops import build_mfg, fpe_similarity_to_target
from scripts.fastcc_all_organisms import load_catalogue, max_demand_flux

FASTCC_EPS = 1e-4


def resolve_path(model_key: str) -> Path:
    """Map a metadata key to its full GEM file path (only needed with --no-cache)."""
    return PROJECT_ROOT / load_catalogue()[model_key]["file_path"]


def central_reactions(model):
    """Central-metabolism reaction ids for the coverage diagnostic. The dedicated
    scripts.central_metab_experiment module was removed from the repo, so fall back
    to a subsystem-keyword heuristic (glycolysis / TCA / PPP / oxidative phos.).
    Returns (set_of_reaction_ids, method_label)."""
    try:  # prefer the curated list if the module is ever restored
        from scripts.central_metab_experiment import central_reactions as _cr
        return _cr(model)
    except Exception:
        keys = ("glycolysis", "gluconeogen", "citric acid", "tca", "tricarbox",
                "pentose phosphate", "oxidative phosphor", "pyruvate metab")
        ids = {r.id for r in model.reactions
               if any(k in (r.subsystem or "").lower() for k in keys)}
        return ids, "subsystem-heuristic"


def restrict_to_consistent(model, consistent_ids: set[str]):
    """Delete blocked reactions -> a model whose reactions are exactly the
    FASTCC-consistent set (the 'graph from fastcc')."""
    reduced = model.copy()
    drop = [r for r in reduced.reactions if r.id not in consistent_ids]
    reduced.remove_reactions(drop, remove_orphans=True)
    return reduced


def load_reduced(key: str, eps: float, use_cache: bool):
    """Cached FASTCC-reduced model by default (extraction/fastcc/models/<key>.json,
    computed once at eps=1e-6 with growth preservation). use_cache=False (or no
    cached model) runs fastcc live at `eps`."""
    if use_cache and key in available_models():
        return load_fastcc_model(key)
    model = cobra.io.load_json_model(str(resolve_path(key)))
    _, consistent_ids = fastcc_from_cobra(model, eps)
    return restrict_to_consistent(model, set(consistent_ids))


def pick_producible(model, rng, eps, max_tries=40):
    """Random metabolite whose demand can carry flux (producible)."""
    mets = list(model.metabolites)
    rng.shuffle(mets)
    for m in mets[:max_tries]:
        v = max_demand_flux(model, m.id, eps)
        if abs(v) > eps:
            return m.id, float(v)
    return None, 0.0


def run_one(key: str, info: dict, rng: random.Random, eps: float,
            cfg: GudieConfig, use_cache: bool) -> None:
    t0 = time.time()
    reduced = load_reduced(key, eps, use_cache)
    consistent = {r.id for r in reduced.reactions}

    mid, _ = pick_producible(reduced, rng, eps)
    if mid is None:
        print(f"{key:9} {info['organism'][:24]:25}  no producible compound found; skipped")
        return

    # add a demand reaction for M, maximise it, build the MFG on that flux state
    dm = reduced.add_boundary(reduced.metabolites.get_by_id(mid), type="demand")
    reduced.objective = dm
    sol = reduced.optimize()
    if sol.status != "optimal":
        print(f"{key:9} {info['organism'][:24]:25}  target={mid} infeasible; skipped")
        return

    G = build_mfg(reduced, sol)
    # the mass-flow graph is the active subgraph: nodes that actually carry mass
    active = [nid for nid in G.nodes if G.degree(nid) > 0]
    G_F = G.subgraph(active).copy()
    if dm.id not in G_F:
        print(f"{key:9} {info['organism'][:24]:25}  target {mid} not in active MFG; skipped")
        return

    # central-metabolism reactions present in the active graph (baseline)
    central, cmethod = central_reactions(reduced)
    central_active = central & set(G_F.nodes)

    embedding_sim = (fpe_similarity_to_target(G_F, dm.id, list(G_F.nodes()))
                     if cfg.lam > 0.0 else None)
    res = run_gudie(G_F, dm.id, model=reduced, config=cfg, embedding_sim=embedding_sim)
    H = set(res.nodes)
    central_in_H = central_active & H

    dt = time.time() - t0
    cov = (f"{len(central_in_H)}/{len(central_active)}"
           if central_active else "0/0")
    print(f"{key:9} {info['organism'][:24]:25} {mid[:14]:14} "
          f"cons={len(consistent):>4} mfg={G_F.number_of_nodes():>4} "
          f"H={len(H):>4} ({100*len(H)/max(G_F.number_of_nodes(),1):>4.1f}%)  "
          f"cen_in_H={cov:>7}  ppr={'Y' if res.ppr_converged else 'n'} "
          f"ludie_fb={'Y' if res.ludie_fallback else 'n'}  [{dt:.1f}s]")


def main() -> None:
    ap = argparse.ArgumentParser(description="End-to-end fastcc -> MFG -> GUDIE over organisms.")
    ap.add_argument("--model", default=None, help="single metadata key (default: all).")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epsilon", type=float, default=FASTCC_EPS)
    ap.add_argument("--alpha", type=float, default=0.3, help="GUDIE PPR restart alpha_r.")
    ap.add_argument("--k", type=float, default=0.1, help="GUDIE inclusion threshold factor.")
    ap.add_argument("--hops", type=int, default=8, help="GUDIE interest propagation hops.")
    ap.add_argument("--lam", type=float, default=0.05,
                    help="Phase B FPE blend coefficient in [0,1] (default 0.05); >0 mixes the target "
                         "Flow-Profile-Encoding cosine into interest (lam=0 => exact Phase A).")
    ap.add_argument("--no-cache", action="store_true",
                    help="ignore the pre-built FASTCC-reduced models; run fastcc live at --epsilon.")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    cfg = GudieConfig(alpha_r=args.alpha, k=args.k, hops=args.hops, lam=args.lam)
    catalogue = load_catalogue()
    if args.model:
        catalogue = {args.model: catalogue[args.model]}

    print("=" * 110)
    print(f"{'fastcc(live)' if args.no_cache else 'fastcc(cached)'} -> MFG -> GUDIE  |  "
          f"organisms={len(catalogue)}  seed={args.seed}  "
          f"alpha_r={args.alpha}  k={args.k}  hops={args.hops}  lam={args.lam}")
    print("=" * 110)
    for key, info in catalogue.items():
        try:
            run_one(key, info, rng, args.epsilon, cfg, use_cache=not args.no_cache)
        except Exception as e:  # keep the sweep going; report the failure
            print(f"{key:9} {info['organism'][:24]:25}  ERROR: {type(e).__name__}: {str(e)[:60]}")


if __name__ == "__main__":
    main()
