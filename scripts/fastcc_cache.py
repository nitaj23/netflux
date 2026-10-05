"""
scripts/fastcc_cache.py

Cache FASTCC flux-consistent reaction sets per organism so downstream scripts
(GUDIE, pFBA, OptForce, validation sweeps) don't recompute the expensive
consistency boundary every run.

Cache layout: fastcc_cache/<model_key>.json  (gitignored, regenerable)
    { "model", "file_hash", "epsilon", "n_reactions", "n_consistent", "consistent": [ids] }

The cache is keyed by the GEM file's SHA1 + epsilon, so editing a model or
changing epsilon transparently invalidates the stale entry.

Use in code
-----------
    from scripts.fastcc_cache import get_consistent_ids
    consistent = get_consistent_ids("iML1515", model)   # set[str]; computes+caches on miss

Build the whole cache
---------------------
    .venv/Scripts/python.exe -u -m scripts.fastcc_cache            # all organisms
    .venv/Scripts/python.exe -u -m scripts.fastcc_cache --rebuild  # force recompute
    .venv/Scripts/python.exe -u -m scripts.fastcc_cache --model iML1515
"""

from __future__ import annotations

import argparse
import hashlib
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

from src.extraction.fastcc import fastcc_from_cobra
from scripts.fastcc_all_organisms import load_catalogue

# 1e-6 (not 1e-4): the coarser 1e-4 threshold prunes trace-flux (<1e-4) but
# biologically essential reactions (cofactor/cell-wall biosynthesis), leaving the
# reduced model unable to grow in 9/12 GEMs. 1e-6 sits at the solver feasibility
# floor -- it removes only genuinely blocked reactions and preserves full growth.
FASTCC_EPS = 1e-6
CACHE_DIR = PROJECT_ROOT / "fastcc_cache"


def resolve_path(model_key: str) -> Path:
    return PROJECT_ROOT / load_catalogue()[model_key]["file_path"]


def _file_hash(path: Path) -> str:
    h = hashlib.sha1()
    h.update(path.read_bytes())
    return h.hexdigest()


def _cache_file(model_key: str) -> Path:
    return CACHE_DIR / f"{model_key}.json"


def get_consistent_ids(model_key: str, model=None, eps: float = FASTCC_EPS,
                       rebuild: bool = False) -> set[str]:
    """Return the FASTCC-consistent reaction ids for `model_key`, using the cache
    when the GEM file hash + epsilon match. Loads the model only on a cache miss
    (or if `model` is supplied, reuses it). Writes/refreshes the cache entry."""
    path = resolve_path(model_key)
    fhash = _file_hash(path)
    cf = _cache_file(model_key)

    if not rebuild and cf.exists():
        try:
            entry = json.loads(cf.read_text())
            if entry.get("file_hash") == fhash and entry.get("epsilon") == eps:
                return set(entry["consistent"])
        except (json.JSONDecodeError, KeyError):
            pass  # corrupt/old cache -> recompute

    if model is None:
        model = cobra.io.load_json_model(str(path))
    _, consistent = fastcc_from_cobra(model, eps)
    consistent = list(consistent)

    CACHE_DIR.mkdir(exist_ok=True)
    cf.write_text(json.dumps({
        "model": model_key, "file_hash": fhash, "epsilon": eps,
        "n_reactions": len(model.reactions), "n_consistent": len(consistent),
        "consistent": sorted(consistent),
    }, indent=1))
    return set(consistent)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build/refresh the FASTCC consistent-set cache.")
    ap.add_argument("--model", default=None, help="single metadata key (default: all).")
    ap.add_argument("--epsilon", type=float, default=FASTCC_EPS)
    ap.add_argument("--rebuild", action="store_true", help="ignore existing cache and recompute.")
    args = ap.parse_args()

    catalogue = load_catalogue()
    keys = [args.model] if args.model else list(catalogue)

    print("=" * 78)
    print(f"FASTCC cache  ->  {CACHE_DIR}   (eps={args.epsilon}, rebuild={args.rebuild})")
    print("=" * 78)
    print(f"{'model':9} {'reactions':>9} {'consistent':>10} {'%':>4}  {'status':>8}  time")
    print("-" * 78)

    for key in keys:
        cf = _cache_file(key)
        hit = (not args.rebuild and cf.exists()
               and json.loads(cf.read_text()).get("file_hash") == _file_hash(resolve_path(key))
               and json.loads(cf.read_text()).get("epsilon") == args.epsilon)
        t0 = time.time()
        ids = get_consistent_ids(key, eps=args.epsilon, rebuild=args.rebuild)
        entry = json.loads(cf.read_text())
        n = entry["n_reactions"]
        status = "cached" if hit else "computed"
        print(f"{key:9} {n:>9} {len(ids):>10} {100*len(ids)//max(n,1):>3}%  "
              f"{status:>8}  [{time.time()-t0:.1f}s]")

    print("-" * 78)
    print(f"cache dir: {CACHE_DIR}  ({len(list(CACHE_DIR.glob('*.json')))} entries)")


if __name__ == "__main__":
    main()
