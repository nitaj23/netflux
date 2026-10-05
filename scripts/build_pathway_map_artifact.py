"""Assemble maps/*.json (from scripts.export_pathway_maps) plus designs/*.json (from
scripts.export_design_results) into a single self-contained HTML atlas: one
Escher-style map per (carbon source -> target) pair, each with the strain-design
engine's verified KO front attached.

The data is inlined into scripts/pathway_map_template.html at the
<script id="map-data"> placeholder, because a published artifact's CSP blocks
fetch() to any external host -- the page has to carry its own data. That is also why
the design results are precomputed and inlined rather than called live: the tool layer
is a Python service, so the atlas carries the tool layer's *output* (the
``design_table`` / ``design_csv`` view, verbatim) instead of calling it.

designs/ is OPTIONAL. Without it the atlas builds exactly as before, minus the design
panel -- so a map-only rebuild never blocks on a design run.

Run: .venv/Scripts/python.exe -u -m scripts.build_pathway_map_artifact
     .venv/Scripts/python.exe -u -m scripts.build_pathway_map_artifact --out atlas.html
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

TEMPLATE = PROJECT_ROOT / "scripts" / "pathway_map_template.html"
MAPS_DIR = PROJECT_ROOT / "maps"
DESIGNS_DIR = PROJECT_ROOT / "designs"
RESULTS_DIR = PROJECT_ROOT / "results"
OPTKNOCK_CSVS = ("optknock_k1.csv", "optknock_k2.csv", "optknock_k3.csv")
DEFAULT_OUT = PROJECT_ROOT / "pathway_atlas.html"

# Provenance keys worth showing in the design panel. The full provenance dict carries
# hashes and null placeholders that only add payload; these are the fields that tell
# you how the front was reached.
PROV_KEEP = ("tier", "mode", "F_source", "f_reuse_note", "biomass_id", "target_id",
             "n_reactions", "n_F", "n_candidates", "n_masked", "n_evaluated",
             "n_verified", "n_front", "search_iterations", "threshold",
             "v_star_target", "floor_fraction", "time_limited", "workers",
             "solver", "scip_version", "milp_status", "timings_s",
             "tiers_run", "tier_front_sizes", "tier_failures")

# Fields the viewer actually reads. Anything else is dropped so the inlined
# payload stays as small as the maps allow.
RXN_KEEP = ("id", "name", "flux", "subsystem", "gene_rule", "lower_bound",
            "upper_bound", "design", "cls", "coupling", "v_min", "v_max",
            "repair_added", "is_target", "is_source")
MET_KEEP = ("id", "name", "formula", "C", "currency")
PAIR_KEEP = ("organism", "target_met", "target_met_name", "source_ids", "target_id",
             "n_consistent", "n_active", "n_H", "h_over_active", "connected_pre",
             "connected_post", "repair_fired", "repair_reason", "repair_added",
             "operating_point", "v_star_target", "qualified", "attrition_reason")


def _round(v, nd):
    return None if v is None else round(float(v), nd)


def trim_pair(p: dict, pid: str) -> dict:
    out = {k: p[k] for k in PAIR_KEEP if k in p}
    out["_id"] = pid
    if not p.get("qualified"):
        out.setdefault("attrition_reason", "")
        return out
    out["v_star_target"] = _round(p.get("v_star_target"), 6)
    out["reactions"] = [
        {**{k: r[k] for k in RXN_KEEP if k in r},
         "flux": _round(r.get("flux"), 5),
         "v_min": _round(r.get("v_min"), 5),
         "v_max": _round(r.get("v_max"), 5)}
        for r in p["reactions"]
    ]
    out["metabolites"] = [{k: m[k] for k in MET_KEEP if k in m} for m in p["metabolites"]]
    out["edges"] = [{"source": e["source"], "target": e["target"],
                     "mass_flow": _round(e["mass_flow"], 6),
                     "metabolites": e["metabolites"]}
                    for e in p["edges"]]
    return out


def trim_design(rec: dict) -> dict:
    """One pair's design-engine result, trimmed to what the panel renders. The table's
    ``columns`` are carried verbatim -- that self-describing schema is exactly what
    lets the viewer build the grid without hardcoding the engine's columns."""
    out = {"status": rec["status"], "target_id": rec.get("target_id", ""),
           "elapsed_s": rec.get("elapsed_s")}
    # trivial-KO pruning stats ride along on both outcomes: on a "none" they are the
    # evidence that the search still had a real candidate set to work with.
    if rec.get("prune"):
        out["prune"] = rec["prune"]
    if rec["status"] != "ok":
        out["reason"] = rec.get("reason", "")
        out["hint"] = rec.get("hint", "")
        # Per-tier detail on a pair that produced nothing: it separates "the MILP was
        # capped at the greedy incumbent size" from "no design exists" from a solver
        # fault, which all otherwise read as the same empty panel.
        if rec.get("tier_failures"):
            out["tier_failures"] = rec["tier_failures"]
        return out
    t = rec["table"]
    out["columns"] = t["columns"]
    out["rows"] = t["rows"]
    out["n_designs"] = t["n_designs"]
    out["guarantee_label"] = rec["guarantee_label"]
    out["csv"] = rec["csv"]
    out["params"] = rec["params"]
    prov = rec.get("provenance", {})
    out["provenance"] = {k: prov[k] for k in PROV_KEEP if k in prov and prov[k] is not None}
    return out


def load_designs(designs_dir: Path) -> dict[str, dict[str, dict]]:
    """{organism key -> {target_met -> trimmed design record}}. Missing dir -> empty."""
    idx: dict[str, dict[str, dict]] = {}
    if not designs_dir.is_dir():
        return idx
    for f in sorted(designs_dir.glob("*.json")):
        bundle = json.loads(f.read_text(encoding="utf-8"))
        idx[bundle["organism"]] = {r["target_met"]: trim_design(r)
                                   for r in bundle.get("pairs", [])}
    return idx


# ---------------------------------------------------------------------------
# OptKnock results -> design-panel records.
#
# The optknock engine (src/optknock) targets SECRETED products (acetate, ethanol,
# ...), which do not line up with the intracellular corridor targets the maps are
# drawn for. So its designs cannot attach to an existing map pair; instead each
# (organism -> secreted target) with a buildable design is added as a "design-only"
# pair: no corridor geometry, just the verified KO front in the design drawer.
#
# The record uses the same self-describing {columns, rows} schema the viewer
# already renders (report.py), so the grid is built from the data with no viewer
# change; `engine: "optknock"` only swaps the wording of the notes/provenance.
# ---------------------------------------------------------------------------

# Self-describing grid schema for an optknock design. `knockouts` is the array key
# the viewer chips; the three `axis` columns are the growth-coupling trade-off.
_OK_COLUMNS = [
    {"key": "rank", "label": "#", "numeric": True, "type": "int", "unit": None, "axis": False},
    {"key": "engine_rank", "label": "Engine rank", "numeric": True, "type": "int", "unit": None, "axis": False},
    {"key": "knockouts", "label": "Reaction knockouts", "numeric": False, "type": "str", "unit": None, "axis": False},
    {"key": "ko_count", "label": "k", "numeric": True, "type": "int", "unit": None, "axis": True},
    {"key": "guaranteed_floor", "label": "Min product @ max growth", "numeric": True, "type": "float", "unit": "mmol/gDW/h", "axis": True},
    {"key": "yield_or_titer", "label": "Max product", "numeric": True, "type": "float", "unit": "mmol/gDW/h", "axis": True},
    {"key": "verified_growth", "label": "Growth at design", "numeric": True, "type": "float", "unit": "1/h", "axis": False},
    {"key": "retained_mu", "label": "Retained growth", "numeric": True, "type": "float", "unit": None, "axis": False},
    {"key": "growth_robustness_ratio", "label": "Robustness", "numeric": True, "type": "float", "unit": None, "axis": False},
    {"key": "coupling_class", "label": "Coupling", "numeric": False, "type": "str", "unit": None, "axis": False},
    {"key": "gene_ko", "label": "Gene KO", "numeric": False, "type": "str", "unit": None, "axis": False},
    {"key": "source", "label": "Search", "numeric": False, "type": "str", "unit": None, "axis": False},
]


def _ok_float(s):
    try:
        return round(float(s), 5)
    except (TypeError, ValueError):
        return None


def _parse_kos(interventions: str) -> list[str]:
    """"ECHH_10:KO, POR2:KO" -> ["ECHH_10", "POR2"] (reaction ids, action stripped)."""
    out = []
    for tok in interventions.split(","):
        tok = tok.strip()
        if tok:
            out.append(tok.rsplit(":", 1)[0].strip())
    return out


def load_optknock(results_dir: Path) -> dict[str, list[dict]]:
    """{organism -> [design-only pair dict, ...]}.

    Reads results/optknock_k{1,2,3}.csv, keeps only rows with a BUILDABLE design
    (a gene-KO set realizes the reaction knockouts), and groups them into one
    design-only pair per (organism, secreted target). Missing files -> {}.
    """
    # (organism, target_met) -> {"meta": {...}, "rows": [raw csv rows], "n_total": int}
    groups: dict[tuple[str, str], dict] = {}
    for name in OPTKNOCK_CSVS:
        path = results_dir / name
        if not path.exists():
            continue
        with path.open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh):
                key = (r["organism"], r["target_met"])
                g = groups.setdefault(key, {"meta": r, "rows": [], "n_total": 0})
                g["n_total"] += 1
                if r.get("buildable", "").strip().lower() == "true":
                    g["rows"].append(r)

    out: dict[str, list[dict]] = {}
    for (org, target_met), g in groups.items():
        if not g["rows"]:                       # no buildable design for this target
            continue
        meta = g["meta"]
        # Sort by the engine's Pareto rank (each CSV's `rank` folds in pareto_front
        # + tie-breaks), cardinality only as a tie-break when ranks collide across
        # the single-cardinality runs. The displayed `#` is a fresh global index
        # over that order (and the unique row id the viewer uses).
        def _int(x, default=0):
            try:
                return int(x)
            except (TypeError, ValueError):
                return default

        raw = sorted(g["rows"], key=lambda r: (_int(r["rank"], 10 ** 9),
                                               _int(r["cardinality"])))
        rows, by_k = [], {}
        for i, r in enumerate(raw, 1):
            kos = _parse_kos(r["interventions"])
            k = int(r["cardinality"])
            by_k[str(k)] = by_k.get(str(k), 0) + 1
            rows.append({
                "rank": i,
                "engine_rank": _int(r["rank"], None) if str(r.get("rank", "")).strip() else None,
                "pareto_front": _int(r["pareto_front"], None) if str(r.get("pareto_front", "")).strip() else None,
                "knockouts": kos,
                "ko_count": k,
                "guaranteed_floor": _ok_float(r["verified_floor"]),
                "yield_or_titer": _ok_float(r["verified_ceiling"]),
                "verified_growth": _ok_float(r["verified_growth"]),
                "retained_mu": _ok_float(r["retained_mu"]),
                "growth_robustness_ratio": _ok_float(r["growth_robustness_ratio"]),
                "coupling_class": r.get("coupling_class", ""),
                "gene_ko": r.get("gene_ko", ""),
                "source": r.get("source", ""),
                "coupling_provenance": {rid: r.get("coupling_class", "") for rid in kos},
                # kept for the shared KO-detail view; optknock has no floor threshold
                "threshold_used": None, "floor_margin": None, "tier": k,
            })
        max_k = max((row["ko_count"] for row in rows), default=1)
        design = {
            "status": "ok",
            "engine": "optknock",
            "target_id": meta.get("target_rxn", ""),
            "n_designs": len(rows),
            "guarantee_label": (
                "OptKnock (verified): production floor > 0 at maximum growth, "
                "k ≤ 3 reaction knockouts. “Buildable” means a gene-deletion "
                "set realizes the reaction knockouts. Rows are sorted by the engine’s "
                "Pareto rank (Engine rank column); there is no composite score."),
            "elapsed_s": None,
            "columns": _OK_COLUMNS,
            "rows": rows,
            "csv": _ok_csv(rows),
            "params": {"engine": "optknock", "max_knockouts": max_k,
                       "substrate": meta.get("substrate", ""),
                       "solver": "GLPK (k=1) / HiGHS (k≥2)"},
            "provenance": {"engine": "optknock", "substrate": meta.get("substrate", ""),
                           "target_rxn": meta.get("target_rxn", ""),
                           "target_met": target_met,
                           "wild_type_growth": _ok_float(meta.get("mu_wt")),
                           "n_buildable": len(rows), "n_total": g["n_total"],
                           "by_k": by_k, "solver": "GLPK / HiGHS"},
        }
        pair = {
            "organism": org, "target_met": target_met,
            "target_met_name": meta.get("target_name", target_met),
            "source_ids": [meta.get("substrate", "")] if meta.get("substrate") else [],
            "target_id": meta.get("target_rxn", ""),
            "qualified": True, "design_only": True,
            "n_H": 0, "reactions": [], "edges": [], "metabolites": [],
            "design_result": design,
        }
        out.setdefault(org, []).append(pair)
    # stable target order within an organism
    for org in out:
        out[org].sort(key=lambda p: p["target_met"])
    return out


def _ok_csv(rows: list[dict]) -> str:
    """Downloadable CSV of an optknock design front (mirrors the grid columns)."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow([c["label"] for c in _OK_COLUMNS])
    for row in rows:
        w.writerow(["+".join(row[c["key"]]) if c["key"] == "knockouts"
                    else ("" if row.get(c["key"]) is None else row.get(c["key"]))
                    for c in _OK_COLUMNS])
    return buf.getvalue()


def build(maps_dir: Path, template: Path, designs_dir: Path,
          results_dir: Path | None = None) -> tuple[str, dict]:
    meta_path = PROJECT_ROOT / "gem_models" / "metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.exists() else {}
    designs = load_designs(designs_dir)
    optknock = load_optknock(results_dir) if results_dir else {}

    organisms, n_pairs, n_ok, n_des, n_designs, n_ok_pairs = [], 0, 0, 0, 0, 0
    for f in sorted(maps_dir.glob("*.json")):
        bundle = json.loads(f.read_text(encoding="utf-8"))
        key = bundle["organism"]
        dmap = designs.get(key, {})
        pairs = []
        for i, p in enumerate(bundle["pairs"]):
            tp = trim_pair(p, f"{key}:{i}")
            d = dmap.get(p.get("target_met"))
            if d is not None:
                tp["design_result"] = d
                if d["status"] == "ok" and d["n_designs"]:
                    n_des += 1
                    n_designs += d["n_designs"]
            pairs.append(tp)
            n_pairs += 1
            n_ok += bool(p.get("qualified"))
        # OptKnock design-only pairs for this organism (secreted targets, no corridor).
        for j, op in enumerate(optknock.pop(key, [])):
            op["_id"] = f"{key}:ok{j}"
            pairs.append(op)
            n_pairs += 1
            n_ok += 1
            n_des += 1
            n_ok_pairs += 1
            n_designs += op["design_result"]["n_designs"]
        organisms.append({"key": key,
                          "organism": meta.get(key, {}).get("organism", ""),
                          "pairs": pairs})
    # Any optknock organism with no map bundle still gets its own block.
    for key, ops in optknock.items():
        pairs = []
        for j, op in enumerate(ops):
            op["_id"] = f"{key}:ok{j}"
            pairs.append(op)
            n_pairs += 1
            n_ok += 1
            n_des += 1
            n_ok_pairs += 1
            n_designs += op["design_result"]["n_designs"]
        organisms.append({"key": key,
                          "organism": meta.get(key, {}).get("organism", ""),
                          "pairs": pairs})
    data = {"organisms": organisms, "has_designs": bool(designs or n_ok_pairs)}
    payload = json.dumps(data, separators=(",", ":"))
    html = template.read_text(encoding="utf-8")
    if "__DATA__" not in html:
        raise SystemExit("template is missing the __DATA__ placeholder")
    # </script> inside the JSON would close the host <script> tag early.
    html = html.replace("__DATA__", payload.replace("</", "<\\/"))
    return html, {"organisms": len(organisms), "pairs": n_pairs, "qualified": n_ok,
                  "with_designs": n_des, "designs": n_designs,
                  "optknock_pairs": n_ok_pairs,
                  "payload_kb": len(payload) // 1024}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maps-dir", default=str(MAPS_DIR))
    ap.add_argument("--designs-dir", default=str(DESIGNS_DIR))
    ap.add_argument("--results-dir", default=str(RESULTS_DIR),
                    help="optknock_k{1,2,3}.csv location; '' to skip optknock designs")
    ap.add_argument("--template", default=str(TEMPLATE))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()
    results_dir = Path(args.results_dir) if args.results_dir else None
    html, stats = build(Path(args.maps_dir), Path(args.template),
                        Path(args.designs_dir), results_dir)
    out = Path(args.out)
    out.write_text(html, encoding="utf-8")
    print(f"organisms={stats['organisms']}  pairs={stats['pairs']}  "
          f"qualified={stats['qualified']}  data={stats['payload_kb']} KB")
    print(f"pairs with KO designs={stats['with_designs']}  "
          f"verified designs={stats['designs']}  "
          f"(optknock design-only pairs={stats['optknock_pairs']})")
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
