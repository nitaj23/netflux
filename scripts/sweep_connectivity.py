"""Source->target connectivity sweep across the corpus.

For many (source, target) pairs per organism, measure whether a DIRECTED
source->target corridor exists inside the extracted subgraph H, before vs after
connectivity repair. Pure connectivity -- no R-FVA -- so it is cheap enough to
sweep many pairs. Reproducible: per-organism seeding (draws independent of corpus
composition); deterministic pipeline.

Per (source, target) pair it records, in BOTH the GUDIE-only H (repair off) and the
repaired H+ (repair on):
    src_present   source reaction retained in H
    connected     directed path source -> target exists in H
    path_len      shortest source->target hop count in H (- if none)
plus the repair record (fired / reason / reactions added / added path).

Run
---
    .venv/Scripts/python.exe -u -m scripts.sweep_connectivity
    .venv/Scripts/python.exe -u -m scripts.sweep_connectivity --n-targets 8 --max-sources 2
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import networkx as nx

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.extraction.fastcc import available_models
from src.extraction.pathway import extract_bioproduction_pathway
from scripts.validate_gudie_all_organisms import (
    load_reduced, pick_proper_producible_many, carbon_uptake_sources, EPS)

DEFAULT_EXCLUDE = {"iLB1027"}   # phototroph: no single carbon-uptake anchor

CSV_COLUMNS = [
    "organism", "target_met", "target_id", "source", "src_C", "src_uptake",
    "n_active", "n_H_gudie", "n_H_repaired",
    "src_present_gudie", "connected_pre", "path_len_pre",
    "src_present_repaired", "connected_post", "path_len_post",
    "repair_fired", "repair_reason", "repair_n_added", "repair_added",
    "note",
]


def _sources(model, met_id, max_sources):
    """Detected carbon uptake(s) in the target's own demand-maximised state:
    [(exchange_id, carbon, |uptake flux|), ...]."""
    with model:
        dm = model.add_boundary(model.metabolites.get_by_id(met_id), type="demand")
        model.objective = dm
        sol = model.optimize()
        if sol.status != "optimal":
            return []
        return carbon_uptake_sources(model, sol, EPS, max_sources)


def _conn(H, source, target):
    """(src_present, connected, path_len) for one source in subgraph H."""
    if source not in H or target not in H:
        return (source in H, False, "")
    if nx.has_path(H, source, target):
        return (True, True, nx.shortest_path_length(H, source, target))
    return (True, False, "")


def sweep_target(key, model, met, max_sources) -> list[dict]:
    srcs = _sources(model, met.id, max_sources)
    src_ids = [rid for rid, _, _ in srcs]
    base = {c: "" for c in CSV_COLUMNS}
    base.update(organism=key, target_met=met.id)
    if not src_ids:
        base["note"] = "no carbon source detected"
        return [base]

    # GUDIE-only H (repair off) and repaired H+ (repair on). No FVA -> fast.
    try:
        r0 = extract_bioproduction_pathway(model, met.id, source_ids=src_ids, repair=False)
        r1 = extract_bioproduction_pathway(model, met.id, source_ids=src_ids, repair=True)
    except Exception as e:
        base["note"] = f"{type(e).__name__}: {str(e)[:60]}"
        return [base]

    H0, H1 = r0.subgraph, r1.subgraph
    tid = r1.target_id
    rep = r1.repair
    rows = []
    for rid, C, upt in srcs:
        sp0, c0, l0 = _conn(H0, rid, tid)
        sp1, c1, l1 = _conn(H1, rid, tid)
        row = {c: "" for c in CSV_COLUMNS}
        row.update(
            organism=key, target_met=met.id, target_id=tid, source=rid,
            src_C=C, src_uptake=round(upt, 3),
            n_active=r1.n_active, n_H_gudie=len(r0.reactions), n_H_repaired=len(r1.reactions),
            src_present_gudie=sp0, connected_pre=c0, path_len_pre=l0,
            src_present_repaired=sp1, connected_post=c1, path_len_post=l1,
            repair_fired=(rep.repaired if rep else False),
            repair_reason=(rep.reason if rep else ""),
            repair_n_added=(rep.n_added if rep else 0),
            repair_added=(";".join(rep.added) if rep and rep.added else ""),
        )
        rows.append(row)
    return rows


def _truth(v):
    return str(v) == "True"


def compile_results(rows):
    pairs = [r for r in rows if r["source"]]           # real (source,target) pairs
    nose = [r for r in rows if not r["source"]]
    print("\n" + "=" * 94)
    print(f"CONNECTIVITY SWEEP   organisms={len({r['organism'] for r in rows})}   "
          f"(source,target) pairs={len(pairs)}   targets w/o source={len(nose)}")
    print("=" * 94)

    # ── exhaustive per-pair listing, grouped by organism ──────────────────────
    byorg = defaultdict(list)
    for r in pairs:
        byorg[r["organism"]].append(r)
    hdr = (f"   {'target':13}{'source':16}{'|H|g':>5}{'|H|+':>5}"
           f"{'pre':>5}{'Lpre':>5}{'post':>6}{'Lpost':>6}{'repair':>8}{'+n':>4}")
    for org in sorted(byorg):
        print(f"\n### {org} ###")
        print(hdr)
        for r in sorted(byorg[org], key=lambda x: x["target_met"]):
            pre = "Y" if _truth(r["connected_pre"]) else ("n" if _truth(r["src_present_gudie"]) else "·")
            post = "Y" if _truth(r["connected_post"]) else ("n" if _truth(r["src_present_repaired"]) else "·")
            rf = "fired" if _truth(r["repair_fired"]) else ""
            print(f"   {r['target_met']:13}{r['source']:16}{r['n_H_gudie']:>5}{r['n_H_repaired']:>5}"
                  f"{pre:>5}{str(r['path_len_pre']):>5}{post:>6}{str(r['path_len_post']):>6}"
                  f"{rf:>8}{str(r['repair_n_added']):>4}")
        for r in [x for x in rows if x["organism"] == org and not x["source"]]:
            print(f"   {r['target_met']:13}{'(no source / '+r['note'][:30]+')':40}")

    # ── aggregate ─────────────────────────────────────────────────────────────
    n = len(pairs)
    if not n:
        print("\nno pairs."); return
    pre = sum(_truth(r["connected_pre"]) for r in pairs)
    post = sum(_truth(r["connected_post"]) for r in pairs)
    fired = sum(_truth(r["repair_fired"]) for r in pairs)
    still = [r for r in pairs if not _truth(r["connected_post"])]
    absent0 = sum(1 for r in pairs if not _truth(r["src_present_gudie"]))
    absent1 = sum(1 for r in pairs if not _truth(r["src_present_repaired"]))
    print("\n" + "-" * 94)
    print("LEGEND  Y=connected  n=source in H but no directed path  ·=source absent from H")
    print(f"\nPRE-REPAIR  (dual-anchor GUDIE only): {pre}/{n} connected ({100*pre/n:.0f}%)   "
          f"source absent from H: {absent0}/{n}")
    print(f"POST-REPAIR (full pipeline)         : {post}/{n} connected ({100*post/n:.0f}%)   "
          f"source absent from H: {absent1}/{n}")
    print(f"repair fired on {fired}/{n} pairs")
    if still:
        print(f"STILL DISCONNECTED post-repair ({len(still)}): "
              + ", ".join(f"{r['organism']}/{r['target_met']}<-{r['source']}" for r in still))
    lens = [int(r["path_len_post"]) for r in pairs if str(r["path_len_post"]).isdigit()]
    if lens:
        import statistics as st
        print(f"corridor length (post): mean={st.mean(lens):.1f}  min={min(lens)}  max={max(lens)} hops")

    # per-organism roll-up
    print("\nPER-ORGANISM:")
    print(f"   {'organism':10}{'pairs':>6}{'conn_pre':>10}{'conn_post':>11}{'repair_fired':>14}")
    for org in sorted(byorg):
        rs = byorg[org]
        cp = sum(_truth(r["connected_pre"]) for r in rs)
        cq = sum(_truth(r["connected_post"]) for r in rs)
        fr = sum(_truth(r["repair_fired"]) for r in rs)
        print(f"   {org:10}{len(rs):>6}{f'{cp}/{len(rs)}':>10}{f'{cq}/{len(rs)}':>11}{fr:>14}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-targets", type=int, default=8)
    ap.add_argument("--max-sources", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exclude", nargs="*", default=sorted(DEFAULT_EXCLUDE))
    ap.add_argument("--out", default=str(PROJECT_ROOT / "eval_connectivity.csv"))
    ap.add_argument("--diagnose-only", action="store_true")
    args = ap.parse_args()

    if args.diagnose_only:
        with open(args.out, newline="", encoding="utf-8") as f:
            compile_results(list(csv.DictReader(f)))
        return

    exclude = set(args.exclude)
    models = [k for k in available_models() if k not in exclude]
    print(f"CONNECTIVITY SWEEP  models={models}  excluded={sorted(exclude)}")
    print(f"n_targets={args.n_targets}  max_sources={args.max_sources}  seed={args.seed} "
          f"(per-organism seed, repair off vs on)")

    rows = []
    for key in models:
        print(f"\n### {key} ###", flush=True)
        try:
            model = load_reduced(key)
            rng = random.Random(f"{args.seed}:{key}")
            targets = pick_proper_producible_many(model, rng, EPS, args.n_targets)
        except Exception as e:
            print(f"   load/pick failed: {type(e).__name__}: {e}")
            continue
        if not targets:
            print("   no proper producible target")
            continue
        for met in targets:
            try:
                with model:
                    rs = sweep_target(key, model, met, args.max_sources)
            except Exception as e:
                rs = [{**{c: "" for c in CSV_COLUMNS}, "organism": key,
                       "target_met": met.id, "note": f"fatal: {type(e).__name__}"}]
            rows.extend(rs)
            for r in rs:
                if r["source"]:
                    print(f"   {met.id:12} <- {r['source']:14} "
                          f"pre={r['connected_pre']} post={r['connected_post']} "
                          f"repair={r['repair_fired']}+{r['repair_n_added']}")
                else:
                    print(f"   {met.id:12} -- {r['note']}")

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    compile_results(rows)
    print(f"\nwrote {len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
