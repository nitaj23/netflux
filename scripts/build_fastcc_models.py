"""
scripts/build_fastcc_models.py

Generate the FASTCC-reduced GEMs into
    src/extraction/fastcc/models/<key>.json

For each organism in gem_models/metadata.json: take the full GEM, remove every
flux-BLOCKED reaction (keeping the FASTCC-consistent set from the cache), and
save the result as a standalone cobra JSON model -- loadable exactly like a
genome-scale model via src.extraction.fastcc.load_fastcc_model(key).

Also writes models/metadata.json cataloguing each reduced model.

Run
---
    .venv/Scripts/python.exe -u -m scripts.build_fastcc_models
    .venv/Scripts/python.exe -u -m scripts.build_fastcc_models --model iML1515
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cobra

from scripts.fastcc_all_organisms import load_catalogue
from scripts.fastcc_cache import get_consistent_ids, FASTCC_EPS

MODELS_DIR = PROJECT_ROOT / "src" / "extraction" / "fastcc" / "models"


def resolve_path(model_key: str) -> Path:
    return PROJECT_ROOT / load_catalogue()[model_key]["file_path"]


def build_one(key: str, info: dict, eps: float) -> dict:
    src = resolve_path(key)
    model = cobra.io.load_json_model(str(src))
    n_full = len(model.reactions)
    growth_full = float(model.slim_optimize() or 0.0)

    keep = set(get_consistent_ids(key, model, eps))     # cached; computes on miss
    n_consistent = len(keep)

    # Growth-preserving union: also keep every reaction active at max growth.
    # FASTCC at any numerically safe eps prunes trace-flux (<eps) but biologically
    # essential biosynthesis (cofactors, cell wall), leaving most reduced GEMs
    # unable to grow. Unioning the growth-mode support restores viability while
    # keeping the network flux-consistent (each kept reaction carries flux in a
    # fastcc mode OR the growth mode).
    sol = model.optimize()
    n_growth_added = 0
    if sol.status == "optimal" and growth_full > 1e-9:
        gact = {r.id for r in model.reactions if sol.fluxes[r.id] != 0.0}
        n_growth_added = len(gact - keep)
        keep |= gact

    reduced = model.copy()
    reduced.remove_reactions([r for r in reduced.reactions if r.id not in keep],
                             remove_orphans=True)
    reduced.id = f"{model.id or key}_fastcc"

    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    out = MODELS_DIR / f"{key}.json"
    cobra.io.save_json_model(reduced, str(out))
    growth_reduced = float(reduced.slim_optimize() or 0.0)

    return {
        "model_id": key,
        "file_path": f"src/extraction/fastcc/models/{key}.json",
        "source_gem": info["file_path"],
        "organism": info.get("organism", ""),
        "epsilon": eps,
        "growth_preserving": True,
        "num_reactions_full": n_full,
        "num_reactions": len(reduced.reactions),
        "num_consistent": n_consistent,
        "num_growth_added": n_growth_added,
        "num_metabolites": len(reduced.metabolites),
        "num_removed": n_full - len(reduced.reactions),
        "growth_full": round(growth_full, 5),
        "growth_reduced": round(growth_reduced, 5),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Build FASTCC-reduced GEMs into the fastcc package.")
    ap.add_argument("--model", default=None, help="single metadata key (default: all).")
    ap.add_argument("--epsilon", type=float, default=FASTCC_EPS)
    args = ap.parse_args()

    catalogue = load_catalogue()
    keys = [args.model] if args.model else list(catalogue)

    print("=" * 80)
    print(f"Building FASTCC-reduced GEMs -> {MODELS_DIR}")
    print("=" * 80)
    print(f"{'model':9} {'full':>5} {'consist':>7} {'+grow':>5} {'reduced':>7} "
          f"{'growth':>7} {'':>5}  time")
    print("-" * 80)

    meta = {}
    for key in keys:
        t0 = time.time()
        entry = build_one(key, catalogue[key], args.epsilon)
        meta[key] = entry
        grows = "GROWS" if entry["growth_reduced"] > 1e-6 else "DEAD"
        print(f"{key:9} {entry['num_reactions_full']:>5} {entry['num_consistent']:>7} "
              f"{entry['num_growth_added']:>5} {entry['num_reactions']:>7} "
              f"{entry['growth_reduced']:>7.4f} {grows:>5}  [{time.time()-t0:.1f}s]")

    # merge into metadata.json (so single-model rebuilds don't wipe the rest)
    meta_path = MODELS_DIR / "metadata.json"
    existing = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    existing.update(meta)
    meta_path.write_text(json.dumps(existing, indent=1))

    print("-" * 80)
    print(f"wrote {len(meta)} model(s) + metadata.json  ({len(existing)} total catalogued)")


if __name__ == "__main__":
    main()
