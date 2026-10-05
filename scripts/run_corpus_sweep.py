"""Reproducible full-corpus Tier-1 + compactness sweep of the FULL extraction
pipeline (dual-anchor source seeding + connectivity repair + FVA design
annotation), over every cached FASTCC-reduced organism.

Why this script exists (reproducibility)
----------------------------------------
The extraction pipeline itself is deterministic: for a FIXED (organism, target,
source) it is byte-identical across processes and PYTHONHASHSEED values (FASTCC ->
GLPK FBA -> MFG -> GUDIE -> repair -> FVA carry no RNG; GLPK is single-threaded).
The run-to-run drift in the older `run_dual_eval` sweep came from the HARNESS, not
the pipeline:

  1. ONE shared RNG consumed across organisms in sequence -> excluding an organism
     (iLB1027), a model failing to load, or changing --n-targets reshuffles every
     DOWNSTREAM organism's target draw. Fixed here by seeding PER ORGANISM
     (`random.Random(f"{seed}:{key}")`), so each organism's targets are independent
     of corpus composition and order.
  2. R-KEGG pulled over the network into the OS temp dir (which clears) -> KEGG
     columns drift. Tier-1 needs only R-FVA (intrinsic, offline, circularity-safe),
     so this sweep makes NO network calls and is fully offline-reproducible.

Diagnostics (per spec section 2 validity matrix)
------------------------------------------------
TIER 1 (load-bearing for a strain-design candidate scope):
  * PRODUCIBILITY / CONNECTIVITY   source->target reachable inside H (post-repair);
                                   repair_fired + n_added record the closed gap
                                   (spec Connectivity_Repair A1: 68% -> 100%).
  * M2 BACKBONE GATE vs R-FVA      did H retain every FVA-essential reaction at the
                                   operating flux f? (the core invariant; pass>=0.99)
  * M5 DESIGN-TAG CONSISTENCY      do ko_candidate/alt_route KOs leave production
    vs R-FVA growth-essentiality   alive, and do protected/nontarget_essential KOs
                                   actually kill growth? (the design payload)
TIER 2 (scope shape):
  * COMPACTNESS                    |H|, |active MFG|, H/active, |consistent|.

Run
---
    .venv/Scripts/python.exe -u -m scripts.run_corpus_sweep
    .venv/Scripts/python.exe -u -m scripts.run_corpus_sweep --n-targets 6 --seed 1
    .venv/Scripts/python.exe -u -m scripts.run_corpus_sweep --diagnose-only --out eval_corpus_sweep.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import statistics as st
import sys
import time
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.extraction.fastcc import available_models
from src.extraction.pathway import extract_bioproduction_pathway
from src.evaluation import gem_reference as gr
from src.evaluation import metrics as M
from scripts.validate_gudie_all_organisms import (
    load_reduced, pick_proper_producible_many, carbon_uptake_sources, EPS)

# Excluded from the corpus by request:
#   iLB1027  phototroph (Phaeodactylum tricornutum): light-driven, so the
#            demand-maximised flux state has no single dominant carbon UPTAKE to
#            serve as the dual anchor -- source detection is ill-posed there.
#   iYli21   Yarrowia lipolytica: excluded by request from this sweep.
DEFAULT_EXCLUDE = {"iLB1027", "iYli21"}

CSV_COLUMNS = [
    "organism", "target_met", "carbon_source", "n_sources",
    # compactness (Tier 2)
    "n_consistent", "n_active", "n_H_pre", "n_H", "h_over_active",
    # producibility / connectivity (Tier 1)
    "connected_post", "connected_pre", "repair_fired", "repair_reason", "repair_n_added",
    # M2 backbone gate vs R-FVA (Tier 1)
    "n_fva_essential", "backbone_recall", "backbone_gate_pass",
    # M2 attribution: was a missed essential even SELECTABLE (present in the active
    # MFG F that GUDIE carves H from)?  in_F  => GUDIE dropped it (a real M2 fault)
    #                                  out_F => it never entered the candidate pool
    "n_ess_in_active", "backbone_recall_active",
    "n_missed", "n_missed_in_F", "n_missed_out_F", "missed_json",
    # M5 design-tag consistency vs R-FVA (Tier 1)
    "m5_ko_consistent", "m5_ko_total", "m5_ess_consistent", "m5_ess_total",
    "m5_confusion_json",
    # design payload size (the deliverable): tag counts over H
    "n_protected", "n_nontarget_essential", "n_alt_route", "n_ko_candidate",
    # M5b supplementary: tag-matched single-KO test at the SAME operating point the
    # tag was defined at (see _m5b_tag_matched)
    "m5b_flex_ok", "m5b_flex_tot", "m5b_bb_ok", "m5b_bb_tot", "m5b_json",
    # pipeline health / cost
    "ppr_converged", "ludie_fallback", "elapsed_s",
    # bookkeeping
    "qualified", "attrition_reason",
]

FVA_F = 0.5   # operating point: target pinned at f * v* for R-FVA + annotation


def _detect_sources(model, met_id, max_sources):
    """Carbon source(s) in the target's own demand-maximised flux state (the active
    carbon uptake), rolled back via the model context. Deterministic under GLPK."""
    with model:
        dm = model.add_boundary(model.metabolites.get_by_id(met_id), type="demand")
        model.objective = dm
        sol = model.optimize()
        if sol.status != "optimal":
            return []
        srcs = carbon_uptake_sources(model, sol, EPS, max_sources)
    return [rid for rid, _, _ in srcs]


def _m5_scalars(confusion: dict) -> tuple[int, int, int, int]:
    """Collapse the M5 confusion into two consistency ratios:
      ko  : ko_candidate/alt_route single-KO leaves production alive (production_survives)
      ess : protected/nontarget_essential single-KO is growth-essential (growth_essential)"""
    ko_ok = ko_tot = ess_ok = ess_tot = 0
    for tag in ("ko_candidate", "alt_route"):
        oc = confusion.get(tag, {})
        ko_ok += oc.get("production_survives", 0)
        ko_tot += sum(oc.values())
    for tag in ("protected", "nontarget_essential"):
        oc = confusion.get(tag, {})
        ess_ok += oc.get("growth_essential", 0)
        ess_tot += sum(oc.values())
    return ko_ok, ko_tot, ess_ok, ess_tot


def _m5b_tag_matched(result, model, target, f: float) -> dict | None:
    """M5b -- the tag-matched supplement to M5.

    M5 (metrics.py:134) scores the FVA design tags against GROWTH-essentiality: it
    asks whether knocking out a `protected` / `nontarget_essential` reaction kills
    BIOMASS. But those tags are defined by target-production backbone-ness at the
    operating point f (0 not in [v_min, v_max] when the target is pinned at f*v*),
    not by growth. Judging them on growth is a category mismatch, so a low M5
    essential-side number is expected and does not by itself indicate a tagging
    fault.

    M5b closes that gap by testing each tag against the exact predicate that
    defined it -- single-KO, then re-maximise the TARGET:

        backbone  (protected / nontarget_essential)  =>  v_max_target < f * v*
                  i.e. the KO must cost the operating point
        flexible  (alt_route / ko_candidate)         =>  v_max_target >= f * v*
                  i.e. the KO must preserve the operating point

    This is stricter than M5's flexible-side check (which only asks for any
    production at all, > EPS). Both are reported; neither replaces the other.
    """
    ann = getattr(result, "annotation", None)
    if ann is None:
        return None
    v_star = ann.v_star_target
    level = f * abs(v_star)
    tol = max(1e-9, 1e-6 * max(1.0, level))
    m = model.copy()
    tid = gr._seed_target(m, target)
    m.objective = tid
    tags = {}
    for rid, ra in ann.per_reaction.items():
        op = next((o for o in ra.per_op if abs(o.get("f", o.get("frac", -1)) - f) < 1e-9),
                  ra.per_op[-1] if ra.per_op else None)
        if op:
            tags[rid] = op.get("design", "")
    BACKBONE = ("protected", "nontarget_essential")
    FLEX = ("alt_route", "ko_candidate")
    per_tag = {t: {"ok": 0, "tot": 0} for t in BACKBONE + FLEX}
    for rid, tag in tags.items():
        if tag not in per_tag or rid == tid or not m.reactions.has_id(rid):
            continue
        r = m.reactions.get_by_id(rid)
        lo, hi = r.bounds
        r.bounds = (0.0, 0.0)
        v = m.slim_optimize()
        r.bounds = (lo, hi)
        keeps_op = v is not None and abs(v) >= level - tol
        ok = (not keeps_op) if tag in BACKBONE else keeps_op
        per_tag[tag]["tot"] += 1
        per_tag[tag]["ok"] += int(ok)
    flex_ok = sum(per_tag[t]["ok"] for t in FLEX)
    flex_tot = sum(per_tag[t]["tot"] for t in FLEX)
    bb_ok = sum(per_tag[t]["ok"] for t in BACKBONE)
    bb_tot = sum(per_tag[t]["tot"] for t in BACKBONE)
    return {"flex_ok": flex_ok, "flex_tot": flex_tot, "bb_ok": bb_ok,
            "bb_tot": bb_tot, "per_tag": per_tag}


def _base_row(key, target):
    row = {c: "" for c in CSV_COLUMNS}
    row.update(organism=key, target_met=target, qualified=False)
    return row


def evaluate_pair(key, model, met, max_sources) -> dict:
    t0 = time.perf_counter()
    row = _base_row(key, met.id)
    source_ids = _detect_sources(model, met.id, max_sources)
    row["n_sources"] = len(source_ids)
    row["carbon_source"] = ";".join(source_ids)

    # FULL pipeline: dual-anchor (source_ids) + connectivity repair + FVA annotation.
    try:
        result = extract_bioproduction_pathway(
            model, met.id, source_ids=source_ids or None, repair=True,
            fva_operating_points=(FVA_F,))
    except Exception as e:
        row["attrition_reason"] = f"{type(e).__name__}: {str(e)[:70]}"
        row["elapsed_s"] = round(time.perf_counter() - t0, 2)
        return row

    row["qualified"] = True
    # compactness (Tier 2)
    n_added = result.repair.n_added if result.repair is not None else 0
    row.update(n_consistent=result.n_consistent, n_active=result.n_active,
               n_H=len(result.reactions), n_H_pre=len(result.reactions) - n_added,
               h_over_active=round(len(result.reactions) / max(1, result.n_active), 4),
               ppr_converged=result.gudie.ppr_converged,
               ludie_fallback=result.gudie.ludie_fallback)

    # producibility / connectivity in H (Tier 1)
    connected = any(
        s in result.subgraph and result.target_id in result.subgraph
        and nx.has_path(result.subgraph, s, result.target_id)
        for s in source_ids)
    row["connected_post"] = connected
    if result.repair is not None:
        row.update(repair_fired=result.repair.repaired,
                   repair_reason=result.repair.reason,
                   repair_n_added=result.repair.n_added)
        # repair fires ONLY on a disconnected H, so pre-repair connectivity is
        # exactly (connected_post AND repair did NOT fire).
        row["connected_pre"] = connected and not result.repair.repaired
    else:
        row["connected_pre"] = connected

    # metrics (extraction succeeded; a metric hiccup must not sink the row)
    try:
        R_fva = gr.r_fva(model, met.id, f=FVA_F)
        if R_fva.available:
            row["n_fva_essential"] = len(R_fva.reactions)
            m2 = M.m2_backbone_gate(result, R_fva)
            row.update(backbone_recall=m2["backbone_recall"],
                       backbone_gate_pass=m2["gate_pass"])
            # Attribute every miss. F = the active MFG node set GUDIE selected from
            # (gudie.node_order); an essential reaction absent from F was never a
            # candidate (killed upstream by the pipeline's own FASTCC/mass-flow
            # restriction), so blaming GUDIE's localization for it is wrong.
            F = set(result.gudie.node_order)
            ess = R_fva.reactions
            missed = list(m2["missed"])
            in_F = sorted(r for r in missed if r in F)
            out_F = sorted(r for r in missed if r not in F)
            ess_in_F = ess & F
            row.update(
                n_ess_in_active=len(ess_in_F),
                backbone_recall_active=(round(len(ess_in_F & set(result.reactions))
                                              / len(ess_in_F), 4) if ess_in_F else ""),
                n_missed=len(missed), n_missed_in_F=len(in_F), n_missed_out_F=len(out_F),
                missed_json=json.dumps({"in_F": in_F[:40], "out_F": out_F[:40]},
                                       separators=(",", ":"), sort_keys=True))
            m5 = M.m5_design_agreement(result, R_fva, model, met.id)
            if m5.get("available"):
                conf = m5["confusion"]
                # sort_keys => byte-stable column (dict/Counter key order is otherwise
                # hash-randomized across processes; the counts are already deterministic).
                row["m5_confusion_json"] = json.dumps(conf, separators=(",", ":"),
                                                      sort_keys=True)
                ko_ok, ko_tot, ess_ok, ess_tot = _m5_scalars(conf)
                row.update(m5_ko_consistent=ko_ok, m5_ko_total=ko_tot,
                           m5_ess_consistent=ess_ok, m5_ess_total=ess_tot)
                for tag in ("protected", "nontarget_essential", "alt_route", "ko_candidate"):
                    row[f"n_{tag}"] = sum(conf.get(tag, {}).values())
            m5b = _m5b_tag_matched(result, model, met.id, FVA_F)
            if m5b:
                row.update(m5b_flex_ok=m5b["flex_ok"], m5b_flex_tot=m5b["flex_tot"],
                           m5b_bb_ok=m5b["bb_ok"], m5b_bb_tot=m5b["bb_tot"],
                           m5b_json=json.dumps(m5b["per_tag"], separators=(",", ":"),
                                               sort_keys=True))
    except Exception as e:
        row["attrition_reason"] = f"metrics_partial: {type(e).__name__}: {str(e)[:50]}"
    row["elapsed_s"] = round(time.perf_counter() - t0, 2)
    return row


# ── analysis ──────────────────────────────────────────────────────────────────
def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _truth(v):
    return str(v) == "True"


def _stats(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "n/a"
    return (f"mean={st.mean(xs):.3f}  sd={(st.pstdev(xs) if len(xs) > 1 else 0):.3f}  "
            f"min={min(xs):.3f}  max={max(xs):.3f}  (n={len(xs)})")


def diagnose(rows):
    q = [r for r in rows if _truth(r["qualified"])]
    att = [r for r in rows if not _truth(r["qualified"])]
    print("\n" + "=" * 80)
    print(f"CORPUS SWEEP DIAGNOSIS   pairs={len(rows)}  qualified={len(q)}  "
          f"organisms={len({r['organism'] for r in rows})}")
    print("=" * 80)

    if att:
        print("\nATTRITION:")
        for reason, n in Counter(r["attrition_reason"].split(":")[0] for r in att).most_common():
            print(f"   {n:>3}  {reason}")
    if not q:
        print("\nno qualified pairs."); return

    # ── TIER 1 ────────────────────────────────────────────────────────────────
    print("\n--- TIER 1 ---")

    # producibility / connectivity + repair
    pre = sum(1 for r in q if _truth(r["connected_pre"]))
    post = sum(1 for r in q if _truth(r["connected_post"]))
    fired = [r for r in q if _truth(r["repair_fired"])]
    added = [_f(r["repair_n_added"]) for r in fired]
    print(f"\nPRODUCIBILITY / SOURCE->TARGET CONNECTIVITY IN H:")
    print(f"   pre-repair  (dual-anchor GUDIE only): {pre}/{len(q)} ({100*pre/len(q):.0f}%)")
    print(f"   post-repair (full pipeline)         : {post}/{len(q)} ({100*post/len(q):.0f}%)")
    print(f"   repair fired: {len(fired)}/{len(q)} pairs" +
          (f"   |added| {_stats(added)}" if added else ""))
    rsn = Counter(r["repair_reason"] for r in q if r["repair_reason"])
    if rsn:
        print("   repair reasons: " + ", ".join(f"{k}={v}" for k, v in rsn.most_common()))
    disc = [(r["organism"], r["target_met"]) for r in q if not _truth(r["connected_post"])]
    if disc:
        print(f"   STILL DISCONNECTED post-repair: {disc}")

    # M2 backbone gate
    gated = [r for r in q if r["backbone_gate_pass"] != ""]
    bp = sum(1 for r in gated if _truth(r["backbone_gate_pass"]))
    print(f"\nM2 BACKBONE GATE vs R-FVA (core invariant, pass>=0.99):  "
          f"{bp}/{len(gated)} ({100*bp/max(1,len(gated)):.0f}%) pass")
    print(f"   backbone_recall: {_stats([_f(r['backbone_recall']) for r in gated])}")
    print(f"   backbone_recall vs SELECTABLE essentials (ess in active MFG F): "
          f"{_stats([_f(r['backbone_recall_active']) for r in gated])}")
    mi = sum(_f(r["n_missed_in_F"]) or 0 for r in gated)
    mo = sum(_f(r["n_missed_out_F"]) or 0 for r in gated)
    if mi + mo:
        print(f"   missed-essential attribution: {int(mi)} in F (GUDIE dropped a "
              f"selectable reaction) / {int(mo)} outside F (never a candidate: the "
              f"pipeline's FASTCC+mass-flow restriction removed it) "
              f"= {100*mo/(mi+mo):.0f}% not GUDIE's doing")
    fails = [(r["organism"], r["target_met"], r["backbone_recall"],
              r["backbone_recall_active"], r["n_missed_in_F"], r["n_missed_out_F"])
             for r in gated if not _truth(r["backbone_gate_pass"])]
    for org, tgt, rec, reca, mif, mof in fails[:25]:
        print(f"      FAIL {org:9} {tgt:12} recall={rec} recall_active={reca} "
              f"missed_in_F={mif} missed_out_F={mof}")
    if len(fails) > 25:
        print(f"      ... and {len(fails)-25} more")

    # M5 design-tag consistency
    ko_ok = sum(_f(r["m5_ko_consistent"]) or 0 for r in q if r["m5_ko_consistent"] != "")
    ko_tot = sum(_f(r["m5_ko_total"]) or 0 for r in q if r["m5_ko_total"] != "")
    ess_ok = sum(_f(r["m5_ess_consistent"]) or 0 for r in q if r["m5_ess_consistent"] != "")
    ess_tot = sum(_f(r["m5_ess_total"]) or 0 for r in q if r["m5_ess_total"] != "")
    print("\nM5 DESIGN-TAG CONSISTENCY vs R-FVA:")
    if ko_tot:
        print(f"   ko_candidate/alt_route KO leaves production alive: "
              f"{int(ko_ok)}/{int(ko_tot)} ({100*ko_ok/ko_tot:.0f}%)")
    if ess_tot:
        print(f"   protected/nontarget_essential KO kills growth    : "
              f"{int(ess_ok)}/{int(ess_tot)} ({100*ess_ok/ess_tot:.0f}%)")
        print("      (growth is NOT what these tags predict -- see M5b below)")

    fo = sum(_f(r["m5b_flex_ok"]) or 0 for r in q if r["m5b_flex_ok"] != "")
    ft = sum(_f(r["m5b_flex_tot"]) or 0 for r in q if r["m5b_flex_tot"] != "")
    bo = sum(_f(r["m5b_bb_ok"]) or 0 for r in q if r["m5b_bb_ok"] != "")
    bt = sum(_f(r["m5b_bb_tot"]) or 0 for r in q if r["m5b_bb_tot"] != "")
    print("")
    print(f"M5b TAG-MATCHED SINGLE-KO TEST (same operating point f={FVA_F}):")
    if ft:
        print(f"   flexible (ko_candidate/alt_route) KO KEEPS f*v* production: "
              f"{int(fo)}/{int(ft)} ({100*fo/ft:.0f}%)")
    if bt:
        print(f"   backbone (protected/nontarget_essential) KO COSTS f*v*    : "
              f"{int(bo)}/{int(bt)} ({100*bo/bt:.0f}%)")

    # ── TIER 2 ────────────────────────────────────────────────────────────────
    print("\n--- TIER 2: COMPACTNESS ---")
    print(f"   |H|            : {_stats([_f(r['n_H']) for r in q])}")
    print(f"   |active MFG|   : {_stats([_f(r['n_active']) for r in q])}")
    print(f"   H / active     : {_stats([_f(r['h_over_active']) for r in q])}")
    print(f"   |consistent|   : {_stats([_f(r['n_consistent']) for r in q])}")
    print(f"   |H| pre-repair : {_stats([_f(r['n_H_pre']) for r in q])}")
    print("")
    print("   DESIGN PAYLOAD over H (KO-enumerable set):")
    for tag in ("protected", "nontarget_essential", "alt_route", "ko_candidate"):
        print(f"      {tag:20}: {_stats([_f(r[f'n_{tag}']) for r in q])}")
    print("")
    print(f"   cost/pair (s) : {_stats([_f(r['elapsed_s']) for r in q])}")

    # ── per-organism ──────────────────────────────────────────────────────────
    print("\nPER-ORGANISM:")
    byorg = defaultdict(list)
    for r in q:
        byorg[r["organism"]].append(r)
    print(f"   {'organism':10}{'pairs':>6}{'conn_pre':>9}{'conn_post':>10}"
          f"{'backbone%':>10}{'|H|':>7}{'H/act':>7}")
    for org, rs in sorted(byorg.items()):
        cpre = 100 * sum(1 for r in rs if _truth(r["connected_pre"])) / len(rs)
        cpost = 100 * sum(1 for r in rs if _truth(r["connected_post"])) / len(rs)
        g = [r for r in rs if r["backbone_gate_pass"] != ""]
        bpr = 100 * sum(1 for r in g if _truth(r["backbone_gate_pass"])) / max(1, len(g))
        nH = st.mean([_f(r["n_H"]) for r in rs])
        hoa = st.mean([_f(r["h_over_active"]) for r in rs])
        print(f"   {org:10}{len(rs):>6}{cpre:>8.0f}%{cpost:>9.0f}%"
              f"{bpr:>9.0f}%{nH:>7.0f}{hoa:>7.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-targets", type=int, default=5)
    ap.add_argument("--max-sources", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exclude", nargs="*", default=sorted(DEFAULT_EXCLUDE))
    ap.add_argument("--out", default=str(PROJECT_ROOT / "eval_corpus_sweep.csv"))
    ap.add_argument("--diagnose-only", action="store_true",
                    help="re-run analysis on an existing --out CSV (no extraction)")
    args = ap.parse_args()

    if args.diagnose_only:
        with open(args.out, newline="", encoding="utf-8") as f:
            diagnose(list(csv.DictReader(f)))
        return

    exclude = set(args.exclude)
    models = [k for k in available_models() if k not in exclude]
    print(f"CORPUS SWEEP  models={models}  excluded={sorted(exclude)}")
    print(f"n_targets={args.n_targets}  seed={args.seed}  (PER-ORGANISM seeding, "
          f"repair=ON, dual-anchor=ON, FVA f={FVA_F})")

    rows = []
    for key in models:
        print(f"\n### {key} ###", flush=True)
        try:
            model = load_reduced(key)
            # PER-ORGANISM seed: target draw is independent of corpus composition/order,
            # so excluding iLB1027 (or any model) leaves every other organism identical.
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
                with model:                              # isolate demand/bounds per pair
                    row = evaluate_pair(key, model, met, args.max_sources)
            except Exception as e:
                row = _base_row(key, met.id)
                row["attrition_reason"] = f"pair_fatal: {type(e).__name__}: {str(e)[:50]}"
            rows.append(row)
            if row["qualified"]:
                print(f"   {met.id:12} src={row['carbon_source'] or '-':14} "
                      f"|H|={row['n_H']:>3}  conn_pre={row['connected_pre']} "
                      f"conn_post={row['connected_post']} repair+{row['repair_n_added']}  "
                      f"backbone={row['backbone_gate_pass']}")
            else:
                print(f"   {met.id:12} SKIP: {row['attrition_reason']}")

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    diagnose(rows)
    print(f"\nwrote {len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
