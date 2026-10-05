"""Export every (carbon source -> target) subgraph H from the corpus sweep as a
self-describing JSON map, for Escher-style rendering.

This re-runs EXACTLY the sweep's extraction call -- same per-organism seeded target
draw, same auto-detected carbon source, dual-anchor + repair + FVA design
annotation -- and persists H, which the sweep itself computed and discarded. It
does NOT recompute R-FVA / M5 / M5b (the sweep's expensive reference work); a map
needs the subgraph, not the reference metrics, so an export run is much cheaper
than a sweep run.

Output (one file per organism, so a hung organism cannot lose the others):
    maps/<organism>.json = {"organism":..., "pairs":[ {pair}, ... ]}

Each pair carries what an Escher-style map needs:
    reactions  id, name, flux, subsystem, gene_rule, design tag, FCA coupling,
               bounds, and whether repair added it
    edges      source reaction -> target reaction, mass_flow, and the shared
               metabolite ids (the MFG is reaction-to-reaction; the metabolite
               list is what makes the bipartite metabolite/reaction map
               reconstructible)
    metabolites  id, name, formula, carbon count, currency flag
    anchors    target_id, source_ids, repair path

Run (one organism, timeout-guarded -- see run_export_all.sh):
    .venv/Scripts/python.exe -u -m scripts.export_pathway_maps --model iML1515
    .venv/Scripts/python.exe -u -m scripts.export_pathway_maps --all
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
import warnings
from pathlib import Path

import networkx as nx

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.extraction.fastcc import available_models
from src.extraction.pathway import extract_bioproduction_pathway
from scripts.run_corpus_sweep import DEFAULT_EXCLUDE, FVA_F, _detect_sources
from scripts.validate_gudie_all_organisms import (
    load_reduced, pick_proper_producible_many, EPS)

OUT_DIR = PROJECT_ROOT / "maps"

# Currency/cofactor metabolites: present on most edges, so drawing them turns any
# map into a hairball. Kept in the payload but flagged, so the viewer can hide
# them by default and still let you switch them on (Escher's usual treatment).
CURRENCY = {
    "h", "h2o", "atp", "adp", "amp", "pi", "ppi", "nad", "nadh", "nadp", "nadph",
    "co2", "o2", "coa", "nh4", "fad", "fadh2", "q8", "q8h2", "h2o2", "gtp", "gdp",
    "udp", "utp", "ctp", "cdp", "so4", "na1", "k", "mg2", "fe2", "fe3", "cl",
    "thf", "mlthf", "10fthf", "so3", "no3", "no2", "n2", "h2", "hco3", "acp",
}


def _base(mid: str) -> str:
    return mid.rsplit("_", 1)[0] if "_" in mid else mid


def _design_tags(result, f: float) -> dict[str, dict]:
    """{reaction id -> {design, cls, coupling, v_min, v_max}} at operating point f."""
    ann = getattr(result, "annotation", None)
    if ann is None:
        return {}
    out = {}
    for rid, ra in ann.per_reaction.items():
        op = next((o for o in ra.per_op if abs(o.get("frac", -1) - f) < 1e-9),
                  ra.per_op[-1] if ra.per_op else None)
        if op is None:
            continue
        out[rid] = {"design": op.get("design", ""), "cls": op.get("cls", ""),
                    "coupling": ra.coupling, "v_min": op.get("v_min"),
                    "v_max": op.get("v_max"), "v": op.get("v")}
    return out


def export_pair(key, model, met, max_sources) -> dict:
    t0 = time.perf_counter()
    rec = {"organism": key, "target_met": met.id,
           "target_met_name": met.name or met.id}
    source_ids = _detect_sources(model, met.id, max_sources)
    rec["source_ids"] = source_ids
    try:
        result = extract_bioproduction_pathway(
            model, met.id, source_ids=source_ids or None, repair=True,
            fva_operating_points=(FVA_F,))
    except Exception as e:
        rec.update(qualified=False,
                   attrition_reason=f"{type(e).__name__}: {str(e)[:80]}",
                   elapsed_s=round(time.perf_counter() - t0, 2))
        return rec

    H = result.subgraph
    tags = _design_tags(result, FVA_F)
    added = set(result.repair.added) if result.repair is not None else set()

    reactions = []
    for rid in result.reactions:
        d = H.nodes[rid]
        t = tags.get(rid, {})
        reactions.append({
            "id": rid,
            "name": d.get("reaction_name", rid) or rid,
            "flux": round(float(d.get("flux", 0.0)), 6),
            "subsystem": d.get("subsystem", "Unknown") or "Unknown",
            "gene_rule": d.get("gene_rule", "") or "",
            "lower_bound": float(d.get("lower_bound", 0.0)),
            "upper_bound": float(d.get("upper_bound", 0.0)),
            "design": t.get("design", ""),
            "cls": t.get("cls", ""),
            "coupling": (round(float(t["coupling"]), 4)
                         if t.get("coupling") is not None else None),
            "v_min": (round(float(t["v_min"]), 6) if t.get("v_min") is not None else None),
            "v_max": (round(float(t["v_max"]), 6) if t.get("v_max") is not None else None),
            "repair_added": rid in added,
            "is_target": rid == result.target_id,
            "is_source": rid in source_ids,
        })

    met_ids: set[str] = set()
    edges = []
    for u, v, d in H.edges(data=True):
        mets = list(d.get("metabolites", []) or [])
        met_ids.update(mets)
        edges.append({"source": u, "target": v,
                      "mass_flow": round(float(d.get("mass_flow", 0.0)), 8),
                      "metabolites": mets})

    metabolites = []
    for mid in sorted(met_ids):
        try:
            m = model.metabolites.get_by_id(mid)
        except KeyError:
            metabolites.append({"id": mid, "name": mid, "formula": "", "C": 0,
                                "currency": _base(mid) in CURRENCY})
            continue
        metabolites.append({
            "id": mid, "name": m.name or mid, "formula": m.formula or "",
            "C": int(m.elements.get("C", 0)) if m.formula else 0,
            "compartment": m.compartment or "",
            "currency": _base(mid) in CURRENCY,
        })

    connected = any(s in H and result.target_id in H and nx.has_path(H, s, result.target_id)
                    for s in source_ids)
    rec.update(
        qualified=True,
        target_id=result.target_id,
        n_consistent=result.n_consistent,
        n_active=result.n_active,
        n_H=len(result.reactions),
        h_over_active=round(len(result.reactions) / max(1, result.n_active), 4),
        connected_post=connected,
        repair_fired=(result.repair.repaired if result.repair is not None else False),
        repair_reason=(result.repair.reason if result.repair is not None else ""),
        repair_added=sorted(added),
        repair_path=(list(result.repair.path) if result.repair is not None else []),
        connected_pre=connected and not (result.repair is not None and result.repair.repaired),
        ppr_converged=result.gudie.ppr_converged,
        ludie_fallback=result.gudie.ludie_fallback,
        gudie_threshold=float(result.gudie.threshold),
        operating_point=FVA_F,
        v_star_target=(float(result.annotation.v_star_target)
                       if result.annotation is not None else None),
        reactions=reactions, edges=edges, metabolites=metabolites,
        elapsed_s=round(time.perf_counter() - t0, 2),
    )
    return rec


def export_organism(key: str, n_targets: int, max_sources: int, seed: int) -> dict:
    model = load_reduced(key)
    rng = random.Random(f"{seed}:{key}")          # same per-organism draw as the sweep
    targets = pick_proper_producible_many(model, rng, EPS, n_targets)
    pairs = []
    for m in targets:
        try:
            with model:                            # isolate demand/bounds per pair
                rec = export_pair(key, model, m, max_sources)
        except Exception as e:
            rec = {"organism": key, "target_met": m.id, "qualified": False,
                   "attrition_reason": f"pair_fatal: {type(e).__name__}: {str(e)[:60]}"}
        pairs.append(rec)
        if rec.get("qualified"):
            print(f"   {m.id:14} src={','.join(rec['source_ids']) or '-':16} "
                  f"|H|={rec['n_H']:>4}  edges={len(rec['edges']):>5}  "
                  f"mets={len(rec['metabolites']):>4}  {rec['elapsed_s']}s", flush=True)
        else:
            print(f"   {m.id:14} SKIP: {rec['attrition_reason']}", flush=True)
    return {"organism": key, "n_targets": n_targets, "seed": seed,
            "operating_point": FVA_F, "pairs": pairs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--n-targets", type=int, default=8)
    ap.add_argument("--max-sources", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.all:
        keys = [k for k in available_models() if k not in DEFAULT_EXCLUDE]
    elif args.model:
        keys = [args.model]
    else:
        ap.error("pass --model KEY or --all")

    for key in keys:
        print(f"### {key} ###", flush=True)
        try:
            bundle = export_organism(key, args.n_targets, args.max_sources, args.seed)
        except Exception as e:
            print(f"   organism failed: {type(e).__name__}: {e}", flush=True)
            continue
        path = out_dir / f"{key}.json"
        path.write_text(json.dumps(bundle, separators=(",", ":"), sort_keys=True),
                        encoding="utf-8")
        n_ok = sum(1 for p in bundle["pairs"] if p.get("qualified"))
        print(f"   -> {path}  ({n_ok}/{len(bundle['pairs'])} qualified, "
              f"{path.stat().st_size // 1024} KB)", flush=True)


if __name__ == "__main__":
    main()
