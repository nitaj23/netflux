"""Dual-reference evaluation of GUDIE's OWN output (as produced by
src.extraction.pathway.extract_bioproduction_pathway), under
validate_gudie_all_organisms-style conditions:

  * per organism, N random PROPER producible targets (not currency/cofactors)
  * carbon source(s) AUTO-DETECTED from each target's own demand-maximised flux
    state (the active carbon uptake), passed as GUDIE's dual anchor
  * cached FASTCC-reduced models, GUDIE default hyper-parameters (GudieConfig())
  * each pair in an isolated model context

No baselines: this focuses solely on GUDIE, scores every dual-reference metric on
its H, and then diagnoses the distribution (backbone-gate failures, recovery,
KEGG reliability, fragmentation, design-tag consistency) across all pairs.

Emits eval_dual_latest.csv (one row per organism,target) + a printed diagnosis.

Run:  .venv/Scripts/python.exe -u -m scripts.run_dual_eval --n-targets 4 --seed 1
      .venv/Scripts/python.exe -u -m scripts.run_dual_eval --model iML1515
"""
from __future__ import annotations

import argparse, csv, json, random, statistics as st, sys, warnings
from collections import Counter, defaultdict
from pathlib import Path
import networkx as nx
warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.extraction.pathway import extract_bioproduction_pathway
from src.evaluation import Reference
from src.evaluation import gem_reference as gr
from src.evaluation import kegg_reference as kr
from src.evaluation import metrics as M
from src.evaluation.divergence import divergence
from scripts.validate_gudie_all_organisms import (
    load_reduced, pick_proper_producible_many, carbon_uptake_sources, EPS)

KEGG_ORG = {"iML1515": "eco", "e_coli_core": "eco", "iMM904": "sce",
            "iAF692": "mba", "iNF517": "lla", "iJN1463": "ppu", "iJN1026": "ppu",
            "iCN1361": "", "iCTH669": "", "iAF987": "gme", "iLB1027": "pti"}

DEFAULT_MODELS = ["iML1515", "iMM904"]

CSV_COLUMNS = [
    "organism", "kegg_org", "target_met", "condition_id", "method", "k",
    "n_sources", "carbon_source",
    "auprc_fva", "auprc_kegg", "auprc_efm",
    "backbone_recall_f05", "backbone_gate_pass",
    "sn_kegg", "ppv_kegg_strict", "ppv_kegg_charitable",
    "accg_kegg_strict", "accg_kegg_charitable", "kegg_coverage",
    "sn_sub", "ppv_sub_strict", "accg_sub_strict",
    "sn_efm", "ppv_efm_strict", "accg_efm_strict", "efm_tractable",
    "h_over_active", "n_H",
    "m5_confusion_json", "m6_oddsratio_mcs", "m6_p_mcs", "m6_precision_at_k_lit",
    "n_subsystems_touched", "mean_subsystem_capture",
    "namespace_loss", "real_miss", "model_variant",
    "ppr_converged", "ludie_fallback", "qualified", "attrition_reason",
    "repair_enabled", "connected", "repair_fired", "repair_reason", "repair_n_added",
]


def _kegg_ref(model, target, key):
    """R-KEGG, but a network/fetch failure degrades to an unavailable reference
    (KEGG is the only network dependency; a blip must not abort a corpus sweep).
    Matches the spec's treatment of KEGG as organism/coverage-gated."""
    try:
        return kr.r_kegg(model, target, KEGG_ORG.get(key))
    except Exception as e:
        return Reference("R-KEGG", "external", "kegg", frozenset(), "route",
                         available=False, reason=f"fetch failed: {type(e).__name__}")


def _build_refs(model, target, key):
    return {
        "sub": gr.r_subsystem(model, target),
        "fva": gr.r_fva(model, target, f=0.5),
        "lit": gr.r_literature(target),
        "efm": gr.r_efm(model, target),
        "mcs": gr.r_mcs(model, target),
        "kegg": _kegg_ref(model, target, key),
    }


def _fill_metrics(row, result, refs, seeds, model, target):
    R_sub, R_fva, R_lit = refs["sub"], refs["fva"], refs["lit"]
    R_efm, R_mcs, R_kegg = refs["efm"], refs["mcs"], refs["kegg"]
    row["efm_tractable"] = R_efm.available
    row["n_H"] = len(result.reactions)
    row["h_over_active"] = round(len(result.reactions) / max(1, result.n_active), 4)

    # M1 AUPRC (dual)
    if R_fva.available:
        row["auprc_fva"] = M.m1_auprc(result, R_fva)["auprc"]
    if R_kegg.available:
        row["auprc_kegg"] = M.m1_auprc(result, R_kegg)["auprc"]
    if R_efm.available:
        row["auprc_efm"] = M.m1_auprc(result, R_efm)["auprc"]

    # M2 backbone gate
    if R_fva.available:
        m2 = M.m2_backbone_gate(result, R_fva)
        row["backbone_recall_f05"] = m2["backbone_recall"]
        row["backbone_gate_pass"] = m2["gate_pass"]

    # M4 topological (dual)
    if R_kegg.available:
        m4 = M.m4_topological(result, R_kegg, seeds=seeds)
        row.update(sn_kegg=m4["Sn"], ppv_kegg_strict=m4["PPV_strict"],
                   ppv_kegg_charitable=m4["PPV_charitable"],
                   accg_kegg_strict=m4["accg_strict"],
                   accg_kegg_charitable=m4["accg_charitable"],
                   kegg_coverage=R_kegg.coverage)
    if R_sub.available and R_sub.reactions:
        m4s = M.m4_topological(result, R_sub, seeds=seeds)
        row.update(sn_sub=m4s["Sn"], ppv_sub_strict=m4s["PPV_strict"],
                   accg_sub_strict=m4s["accg_strict"])
    if R_efm.available:
        m4e = M.m4_topological(result, R_efm, seeds=seeds)
        row.update(sn_efm=m4e["Sn"], ppv_efm_strict=m4e["PPV_strict"],
                   accg_efm_strict=m4e["accg_strict"])

    # M8 fragmentation
    if R_sub.available:
        m8 = M.m8_fragmentation(result, R_sub)
        row["n_subsystems_touched"] = m8["n_subsystems_touched"]
        row["mean_subsystem_capture"] = m8["mean_subsystem_capture"]

    # M5 / M6 design + intervention
    if R_fva.available:
        m5 = M.m5_design_agreement(result, R_fva, model, target)
        if m5.get("available"):
            row["m5_confusion_json"] = json.dumps(m5["confusion"], separators=(",", ":"))
    if R_lit.available:
        row["m6_precision_at_k_lit"] = M.m6_enrichment(result, R_lit).get("precision_at_k", "")
    if R_mcs.available:
        m6m = M.m6_enrichment(result, R_mcs)
        row["m6_oddsratio_mcs"] = m6m.get("odds_ratio", "")
        row["m6_p_mcs"] = m6m.get("hypergeom_p", "")

    # divergence
    if R_kegg.available and R_fva.available:
        d = divergence(result, R_kegg, R_fva)
        row.update(namespace_loss=d["namespace_loss"], real_miss=d["real_miss"],
                   model_variant=d["model_variant"])
    return row


def _base_row(key, target):
    row = {c: "" for c in CSV_COLUMNS}
    row.update(organism=key, kegg_org=KEGG_ORG.get(key, ""), target_met=target,
               condition_id="default_medium", method="GUDIE", qualified=False)
    return row


def _detect_sources(model, met_id, max_sources):
    """Carbon source(s) in the target's own demand-maximised flux state, rolled back."""
    with model:
        dm = model.add_boundary(model.metabolites.get_by_id(met_id), type="demand")
        model.objective = dm
        sol = model.optimize()
        if sol.status != "optimal":
            return []
        srcs = carbon_uptake_sources(model, sol, EPS, max_sources)
    return [rid for rid, _, _ in srcs]


def evaluate_pair(key, model, met, max_sources, repair=False) -> dict:
    row = _base_row(key, met.id)
    row["repair_enabled"] = repair
    source_ids = _detect_sources(model, met.id, max_sources)
    row["n_sources"] = len(source_ids)
    row["carbon_source"] = ";".join(source_ids)
    try:
        result = extract_bioproduction_pathway(
            model, met.id, source_ids=source_ids or None, repair=repair,
            fva_operating_points=(0.5,))
    except Exception as e:
        row["attrition_reason"] = f"{type(e).__name__}: {str(e)[:60]}"
        return row
    seeds = set(source_ids)
    row.update(qualified=True, k=result.gudie.threshold,
               ppr_converged=result.gudie.ppr_converged,
               ludie_fallback=result.gudie.ludie_fallback)
    # Direct test of the repair claim: is the target reachable from a source in H?
    row["connected"] = any(s in result.subgraph and result.target_id in result.subgraph
                           and nx.has_path(result.subgraph, s, result.target_id)
                           for s in source_ids)
    if result.repair is not None:
        row.update(repair_fired=result.repair.repaired,
                   repair_reason=result.repair.reason,
                   repair_n_added=result.repair.n_added)
    try:
        refs = _build_refs(model, met.id, key)
        _fill_metrics(row, result, refs, seeds, model, met.id)
    except Exception as e:  # extraction succeeded; a metric hiccup must not sink the sweep
        row["attrition_reason"] = f"metrics_partial: {type(e).__name__}: {str(e)[:50]}"
    return row


# ── diagnosis ────────────────────────────────────────────────────────────────
def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _truth(v):
    """backbone_gate_pass may be an in-memory bool or a CSV string -> normalize."""
    return str(v) == "True"


def _stats(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "n/a"
    return f"mean={st.mean(xs):.3f}  sd={(st.pstdev(xs) if len(xs) > 1 else 0):.3f}  " \
           f"min={min(xs):.3f}  max={max(xs):.3f}  (n={len(xs)})"


def diagnose(rows):
    q = [r for r in rows if _truth(r["qualified"])]
    att = [r for r in rows if not _truth(r["qualified"])]
    print("\n" + "=" * 78)
    print(f"DIAGNOSIS  GUDIE dual-reference eval   pairs={len(rows)}  qualified={len(q)}")
    print("=" * 78)

    if att:
        print("\nATTRITION (extraction failures):")
        for reason, n in Counter(r["attrition_reason"].split(":")[0] for r in att).most_common():
            print(f"   {n:>3}  {reason}")
    if not q:
        print("\nno qualified pairs."); return

    # 1. backbone gate -- the core invariant
    bp = [_truth(r["backbone_gate_pass"]) for r in q]
    fails = [(r["organism"], r["target_met"], r["backbone_recall_f05"])
             for r in q if not _truth(r["backbone_gate_pass"])]
    print(f"\nM2 BACKBONE GATE (core invariant):  pass {sum(bp)}/{len(bp)} "
          f"({100*sum(bp)/len(bp):.0f}%)")
    if fails:
        print("   FAILURES (GUDIE dropped an FVA-essential reaction):")
        for org, tgt, rec in fails:
            print(f"      {org:9} {tgt:10} backbone_recall={rec}")

    # 1b. connectivity + repair (the feature under test)
    conn = sum(1 for r in q if _truth(r["connected"]))
    print(f"\nSOURCE->TARGET CONNECTIVITY IN H:  {conn}/{len(q)} "
          f"({100*conn/len(q):.0f}%) connected")
    if any(_truth(r["repair_enabled"]) for r in q):
        fired = [r for r in q if _truth(r["repair_fired"])]
        added = [_f(r["repair_n_added"]) for r in fired]
        print(f"   repair fired: {len(fired)}/{len(q)} pairs"
              + (f"  |added| {_stats(added)}" if added else ""))
        rsn = Counter(r["repair_reason"] for r in q if r["repair_reason"])
        if rsn:
            print("   repair reasons: " + ", ".join(f"{k}={v}" for k, v in rsn.most_common()))
        disc = [(r["organism"], r["target_met"]) for r in fired if not _truth(r["connected"])]
        if disc:
            print(f"   WARNING: {len(disc)} repaired pairs still not connected: {disc}")

    # 2. topological recovery
    print("\nM4 RECOVERY (absolute -- interpret only relative to a baseline):")
    print(f"   accg_kegg strict : {_stats([_f(r['accg_kegg_strict']) for r in q])}")
    print(f"   accg_kegg charit.: {_stats([_f(r['accg_kegg_charitable']) for r in q])}")
    print(f"   accg_sub  strict : {_stats([_f(r['accg_sub_strict']) for r in q])}")
    print(f"   AUPRC vs KEGG    : {_stats([_f(r['auprc_kegg']) for r in q])}")
    print(f"   AUPRC vs FVA     : {_stats([_f(r['auprc_fva']) for r in q])}")

    # 3. KEGG reliability (divergence)
    cov = [_f(r["kegg_coverage"]) for r in q if _f(r["kegg_coverage"]) is not None]
    unrel = sum(1 for r in q if _f(r["namespace_loss"]) is not None
                and _f(r["real_miss"]) is not None
                and _f(r["namespace_loss"]) >= 2 * max(_f(r["real_miss"]), 1))
    print(f"\nKEGG RELIABILITY:  coverage {_stats(cov)}")
    print(f"   'KEGG-axis-unreliable' (namespace_loss >= 2*real_miss): {unrel}/{len(q)} pairs")
    mv = [_f(r["model_variant"]) for r in q if _f(r["model_variant"]) is not None]
    print(f"   model_variant reactions (GEM route KEGG omits): {_stats(mv)}")

    # 4. compactness + fragmentation
    print(f"\nSCOPE:  |H| {_stats([_f(r['n_H']) for r in q])}")
    print(f"   H / active MFG   : {_stats([_f(r['h_over_active']) for r in q])}")
    print(f"   subsystems touched: {_stats([_f(r['n_subsystems_touched']) for r in q])}")
    print(f"   mean subsys capture: {_stats([_f(r['mean_subsystem_capture']) for r in q])}")

    # 5. design-tag consistency (M5) -- do KO interventions match the tag?
    good = bad = 0
    for r in q:
        if not r["m5_confusion_json"]:
            continue
        conf = json.loads(r["m5_confusion_json"])
        for tag in ("ko_candidate", "alt_route"):
            oc = conf.get(tag, {})
            good += oc.get("production_survives", 0)
            bad += sum(v for k, v in oc.items() if k != "production_survives")
    tot = good + bad
    if tot:
        print(f"\nM5 DESIGN-TAG CONSISTENCY (ko_candidate/alt_route KO should NOT kill "
              f"production): {good}/{tot} consistent ({100*good/tot:.0f}%)")

    # 6. per-organism
    print("\nPER-ORGANISM:")
    byorg = defaultdict(list)
    for r in q:
        byorg[r["organism"]].append(r)
    print(f"   {'organism':10}{'pairs':>6}{'backbone%':>10}{'accg_kegg':>11}{'accg_sub':>10}")
    for org, rs in sorted(byorg.items()):
        bpr = 100 * sum(1 for r in rs if _truth(r["backbone_gate_pass"])) / len(rs)
        ak = [_f(r["accg_kegg_strict"]) for r in rs if _f(r["accg_kegg_strict"]) is not None]
        asub = [_f(r["accg_sub_strict"]) for r in rs if _f(r["accg_sub_strict"]) is not None]
        print(f"   {org:10}{len(rs):>6}{bpr:>9.0f}%"
              f"{(st.mean(ak) if ak else 0):>11.3f}{(st.mean(asub) if asub else 0):>10.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", help="single organism key (default: DEFAULT_MODELS)")
    ap.add_argument("--all", action="store_true",
                    help="evaluate every cached FASTCC-reduced model (the full corpus)")
    ap.add_argument("--n-targets", type=int, default=4)
    ap.add_argument("--max-sources", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--repair", action="store_true",
                    help="enable connectivity repair (extract_bioproduction_pathway repair=True)")
    ap.add_argument("--out", default=str(PROJECT_ROOT / "eval_dual_latest.csv"))
    ap.add_argument("--diagnose-only", action="store_true",
                    help="re-run the diagnosis on an existing --out CSV (no extraction)")
    args = ap.parse_args()
    if args.diagnose_only:
        with open(args.out, newline="", encoding="utf-8") as f:
            diagnose(list(csv.DictReader(f)))
        return
    if args.all:
        from src.extraction.fastcc import available_models
        models = sorted(available_models())
    elif args.model:
        models = [args.model]
    else:
        models = DEFAULT_MODELS
    rng = random.Random(args.seed)

    print(f"GUDIE dual-reference eval   models={models}  n_targets={args.n_targets}  "
          f"seed={args.seed}  (default GudieConfig, auto carbon source)")
    rows = []
    for key in models:
        print(f"\n### {key} ###", flush=True)
        try:
            model = load_reduced(key)
            targets = pick_proper_producible_many(model, rng, EPS, args.n_targets)
        except Exception as e:
            print(f"   load/pick failed: {type(e).__name__}: {e}"); continue
        if not targets:
            print("   no proper producible target"); continue
        for met in targets:
            try:
                with model:                               # isolate demand/bounds per pair
                    row = evaluate_pair(key, model, met, args.max_sources, args.repair)
            except Exception as e:                         # never let one pair sink the corpus
                row = _base_row(key, met.id)
                row["attrition_reason"] = f"pair_fatal: {type(e).__name__}: {str(e)[:50]}"
            rows.append(row)
            if row["qualified"]:
                print(f"   {met.id:12} src={row['carbon_source'] or '-':14} |H|={row['n_H']:>3}  "
                      f"backbone={row['backbone_gate_pass']}  accg_kegg={row['accg_kegg_strict']}  "
                      f"accg_sub={row['accg_sub_strict']}  kegg_cov={row['kegg_coverage']}")
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
