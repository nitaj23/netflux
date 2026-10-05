"""Source x target connectivity GRID per organism.

The plain connectivity sweep gives one carbon source per target because most GEMs
run on a single-carbon medium. To get genuinely MULTIPLE source->target pairs per
organism, this builds a grid: for each candidate carbon source it swaps the medium
so that source is the SOLE carbon input (all other carbon exchanges blocked; non-
carbon nutrients untouched), then for every target records producibility and the
directed source->target corridor inside H, before vs after repair.

A cell is one (source, target) pair. Cells where the product is not makeable from
that substrate are reported as "not producible" -- itself a substrate-scope result.
No R-FVA -> cheap. Reproducible: per-organism target seed; deterministic pipeline.

Run
---
    .venv/Scripts/python.exe -u -m scripts.sweep_grid_connectivity
    .venv/Scripts/python.exe -u -m scripts.sweep_grid_connectivity --n-targets 5 --k-sources 3
"""
from __future__ import annotations

import argparse
import csv
import random
import statistics as st
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
    load_reduced, pick_proper_producible_many, base_id, carbon, EPS,
    CURRENCY, NON_SUBSTRATE_C)
from scripts.fastcc_all_organisms import max_demand_flux

DEFAULT_EXCLUDE = {"iLB1027"}
UPTAKE_LB = -10.0     # forced uptake bound for the chosen sole carbon source

CSV_COLUMNS = [
    "organism", "source", "src_C", "target_met", "target_id",
    "producible", "max_demand_flux",
    "n_active", "n_H_gudie", "n_H_repaired",
    "src_present_gudie", "connected_pre", "path_len_pre",
    "src_present_repaired", "connected_post", "path_len_post",
    "repair_fired", "repair_n_added", "repair_added", "note",
]


def carbon_exchanges(model):
    """[(exch_id, C)] for boundary exchanges of a single carbon substrate
    (>=1 C so C1 substrates like methanol/formate count, not a currency
    metabolite, not CO2/bicarbonate)."""
    out = []
    for r in model.boundary:
        mets = list(r.metabolites)
        if len(mets) != 1:
            continue
        m = mets[0]
        if base_id(m) in CURRENCY or base_id(m) in NON_SUBSTRATE_C:
            continue
        C = carbon(m)
        if C >= 1:
            out.append((r.id, C))
    return out


def swap_sole_source(model, carbon_exch_ids, chosen):
    """Make `chosen` the only carbon uptake: block uptake on every other carbon
    exchange (lb=0), open the chosen one (lb=UPTAKE_LB). Non-carbon nutrient
    exchanges are left exactly as the medium had them.

    Bounds are assigned as a pair (r.bounds) with ub normalised to >=0 first, so a
    forced-uptake exchange (ub<0) never triggers the transient lb>ub check that a
    one-at-a-time assignment would."""
    for cid in carbon_exch_ids:
        r = model.reactions.get_by_id(cid)
        ub = max(float(r.upper_bound), 0.0)              # allow secretion; drop forced uptake
        r.bounds = (UPTAKE_LB, ub) if cid == chosen else (0.0, ub)


def _producible(model, met_id, eps=EPS):
    v = max_demand_flux(model, met_id, eps)
    return (abs(v) > eps, v)


def _conn(H, source, target):
    if source not in H or target not in H:
        return (source in H, False, "")
    if nx.has_path(H, source, target):
        return (True, True, nx.shortest_path_length(H, source, target))
    return (True, False, "")


SCAN_CAP = 45   # bound the alt-source usability scan on exchange-rich models


def pick_sources(model, cexch, targets, k):
    """Up to k carbon sources for this organism: the NATIVE medium carbon input(s)
    first (exchanges open with lb<0 in the loaded model), then alternative carbon
    substrates ranked by how many of `targets` stay producible when that substrate
    is the sole carbon input. Deterministic (fixed scan order + id tie-break)."""
    cids = [c for c, _ in cexch]
    Cof = dict(cexch)
    native = [c for c in cids if model.reactions.get_by_id(c).lower_bound < 0]
    chosen = list(dict.fromkeys(native))                 # native first, de-duped, ordered

    # rank the alternatives (bounded scan: highest-carbon first, then id, for determinism)
    alts = sorted((c for c in cids if c not in chosen),
                  key=lambda c: (-Cof[c], c))[:SCAN_CAP]
    scored = []
    for cid in alts:
        with model:
            swap_sole_source(model, cids, cid)
            cnt = sum(1 for t in targets if _producible(model, t.id)[0])
        if cnt:
            scored.append((cnt, -Cof[cid], cid))
    scored.sort(key=lambda t: (-t[0], t[1], t[2]))       # most targets, then bigger C, then id
    for _, _, cid in scored:
        if len(chosen) >= k:
            break
        if cid not in chosen:
            chosen.append(cid)
    return [(cid, Cof[cid]) for cid in chosen[:k]]


def grid_cell(model, key, source, srcC, met, cids) -> dict:
    row = {c: "" for c in CSV_COLUMNS}
    row.update(organism=key, source=source, src_C=srcC, target_met=met.id)
    with model:
        swap_sole_source(model, cids, source)
        ok, v = _producible(model, met.id)
        row["producible"] = ok
        row["max_demand_flux"] = round(v, 4)
        if not ok:
            row["note"] = "not producible on this source"
            return row
        try:
            r0 = extract_bioproduction_pathway(model, met.id, source_ids=[source], repair=False)
            r1 = extract_bioproduction_pathway(model, met.id, source_ids=[source], repair=True)
        except Exception as e:
            row["note"] = f"{type(e).__name__}: {str(e)[:55]}"
            return row
    H0, H1, tid, rep = r0.subgraph, r1.subgraph, r1.target_id, r1.repair
    sp0, c0, l0 = _conn(H0, source, tid)
    sp1, c1, l1 = _conn(H1, source, tid)
    row.update(
        target_id=tid, n_active=r1.n_active,
        n_H_gudie=len(r0.reactions), n_H_repaired=len(r1.reactions),
        src_present_gudie=sp0, connected_pre=c0, path_len_pre=l0,
        src_present_repaired=sp1, connected_post=c1, path_len_post=l1,
        repair_fired=(rep.repaired if rep else False),
        repair_n_added=(rep.n_added if rep else 0),
        repair_added=(";".join(rep.added) if rep and rep.added else ""))
    return row


def _T(v):
    return str(v) == "True"


def compile_results(rows):
    cells = [r for r in rows if r["source"]]
    prod = [r for r in cells if _T(r["producible"])]
    byorg = defaultdict(list)
    for r in cells:
        byorg[r["organism"]].append(r)

    print("\n" + "=" * 100)
    print(f"SOURCE x TARGET CONNECTIVITY GRID   organisms={len(byorg)}   "
          f"cells={len(cells)}   producible={len(prod)}")
    print("=" * 100)
    print("LEGEND per cell (pre/post):  Y=connected  n=source in H, no path  "
          "X=not producible on source  E=extraction error")

    for org in sorted(byorg):
        rs = byorg[org]
        srcs = sorted({r["source"] for r in rs})
        tgts = sorted({r["target_met"] for r in rs})
        cell = {(r["source"], r["target_met"]): r for r in rs}
        print(f"\n### {org} ###   sources={srcs}")
        head = f"   {'target':14}" + "".join(f"{s[:13]:>15}" for s in srcs)
        print(head)
        for t in tgts:
            line = f"   {t[:14]:14}"
            for s in srcs:
                r = cell.get((s, t))
                if r is None:
                    cellstr = "-"
                elif not _T(r["producible"]):
                    cellstr = "X"
                elif r["note"] or str(r["connected_post"]) == "":
                    cellstr = "E"                        # producible but extraction failed
                else:
                    pre = "Y" if _T(r["connected_pre"]) else "n"
                    post = "Y" if _T(r["connected_post"]) else "n"
                    rf = "+%s" % r["repair_n_added"] if _T(r["repair_fired"]) else ""
                    cellstr = f"{pre}/{post}{rf}"
                line += f"{cellstr:>15}"
            print(line)

    # connectivity is only defined on producible cells that actually extracted
    ext = [r for r in prod if not r["note"] and str(r["connected_post"]) != ""]
    errs = [r for r in prod if r not in ext]
    if not ext:
        print("\nno extracted cells."); return
    pre = sum(_T(r["connected_pre"]) for r in ext)
    post = sum(_T(r["connected_post"]) for r in ext)
    fired = sum(_T(r["repair_fired"]) for r in ext)
    still = [r for r in ext if not _T(r["connected_post"])]
    lens = [int(r["path_len_post"]) for r in ext if str(r["path_len_post"]).isdigit()]
    print("\n" + "-" * 100)
    print(f"CELLS: {len(cells)} total | producible {len(prod)} | extracted {len(ext)} | "
          f"extraction errors {len(errs)} | not-producible {len(cells)-len(prod)}")
    print(f"   (a cell = one source->target pair; X = product not makeable from that substrate)")
    print(f"PRE-REPAIR  connected: {pre}/{len(ext)} ({100*pre/len(ext):.0f}%)")
    print(f"POST-REPAIR connected: {post}/{len(ext)} ({100*post/len(ext):.0f}%)   "
          f"repair fired on {fired} cells")
    if still:
        print("STILL DISCONNECTED: " +
              ", ".join(f"{r['organism']}/{r['target_met']}<-{r['source']}" for r in still))
    if lens:
        print(f"corridor length (post): mean={st.mean(lens):.1f} min={min(lens)} max={max(lens)} hops")

    print("\nPER-ORGANISM:")
    print(f"   {'organism':10}{'srcs':>6}{'tgts':>6}{'cells':>7}{'prod':>6}{'extr':>6}"
          f"{'conn_pre':>10}{'conn_post':>11}{'repaired':>10}")
    for org in sorted(byorg):
        rs = byorg[org]
        p = [r for r in rs if _T(r["producible"])]
        e = [r for r in p if not r["note"] and str(r["connected_post"]) != ""]
        ns = len({r["source"] for r in rs}); nt = len({r["target_met"] for r in rs})
        cp = sum(_T(r["connected_pre"]) for r in e)
        cq = sum(_T(r["connected_post"]) for r in e)
        fr = sum(_T(r["repair_fired"]) for r in e)
        print(f"   {org:10}{ns:>6}{nt:>6}{len(rs):>7}{len(p):>6}{len(e):>6}"
              f"{f'{cp}/{len(e)}':>10}{f'{cq}/{len(e)}':>11}{fr:>10}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-targets", type=int, default=5)
    ap.add_argument("--k-sources", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exclude", nargs="*", default=sorted(DEFAULT_EXCLUDE))
    ap.add_argument("--out", default=str(PROJECT_ROOT / "eval_grid_connectivity.csv"))
    ap.add_argument("--diagnose-only", action="store_true")
    args = ap.parse_args()

    if args.diagnose_only:
        with open(args.out, newline="", encoding="utf-8") as f:
            compile_results(list(csv.DictReader(f)))
        return

    exclude = set(args.exclude)
    models = [k for k in available_models() if k not in exclude]
    print(f"GRID CONNECTIVITY  models={models}  excluded={sorted(exclude)}")
    print(f"n_targets={args.n_targets}  k_sources={args.k_sources}  seed={args.seed} "
          f"(sole-carbon medium swap, repair off vs on)")

    rows = []
    for key in models:
        print(f"\n### {key} ###", flush=True)
        try:
            model = load_reduced(key)
            rng = random.Random(f"{args.seed}:{key}")
            targets = pick_proper_producible_many(model, rng, EPS, args.n_targets)
            cexch = carbon_exchanges(model)
            if not targets:
                print("   no proper producible target"); continue
            if not cexch:
                print("   no carbon-source exchange found"); continue
            cids = [c for c, _ in cexch]
            sources = pick_sources(model, cexch, targets, args.k_sources)
        except Exception as e:
            print(f"   setup failed: {type(e).__name__}: {str(e)[:70]}")
            continue
        print(f"   sources={[s for s,_ in sources]}  targets={[t.id for t in targets]}")
        for source, srcC in sources:
            for met in targets:
                try:
                    row = grid_cell(model, key, source, srcC, met, cids)
                except Exception as e:
                    row = {c: "" for c in CSV_COLUMNS}
                    row.update(organism=key, source=source, target_met=met.id,
                               note=f"fatal: {type(e).__name__}")
                rows.append(row)
                tag = ("X" if not _T(row["producible"]) else
                       f"pre={row['connected_pre']} post={row['connected_post']} "
                       f"rep={row['repair_fired']}+{row['repair_n_added']}")
                print(f"   {source:16} -> {met.id:12} {tag}")

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader(); w.writerows(rows)
    compile_results(rows)
    print(f"\nwrote {len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
