"""
scripts/fastcc_all_organisms.py

Run FASTCC (flux-consistency boundary) for EVERY organism in metadata.json and,
per organism, choose a random compound of interest and report whether it is
actually producible.

Why the producibility column matters
------------------------------------
FASTCC returns the flux-consistent universe (removes only blocked/dead-end
reactions). But a random metabolite may still be a dead-end with no drain, or a
genuine synthesis gap -- in which case NO extraction method can build a flux-
consistent graph around it. We test producibility the validated way: add a
demand reaction (M ->) and maximise it. max flux > epsilon => producible.

Columns
-------
  model / organism / reactions / consistent (fastcc) / %cons
  compound (random) / compartment / producible? / max demand flux
  #consistent reactions touching the compound

Run
---
    .venv/Scripts/python.exe -m scripts.fastcc_all_organisms
    .venv/Scripts/python.exe -m scripts.fastcc_all_organisms --seed 7
    .venv/Scripts/python.exe -m scripts.fastcc_all_organisms --producible-only
    .venv/Scripts/python.exe -m scripts.fastcc_all_organisms --compartment c --seed 1
"""

from __future__ import annotations

import argparse
import json
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

from src.extraction.fastcc import fastcc_from_cobra

# fastcc.py exposes fastcc_from_cobra(model, eps=1e-4); no DEFAULT_EPSILON there.
FASTCC_EPS = 1e-4

METADATA_PATH = PROJECT_ROOT / "gem_models" / "metadata.json"


def load_catalogue() -> dict:
    return json.loads(METADATA_PATH.read_text())


def max_demand_flux(model, mid: str, eps: float) -> float:
    """Add a demand reaction (M ->), maximise it, return max flux. Producible if
    result > eps. Restores the model via the context manager."""
    with model:
        m = model.metabolites.get_by_id(mid)
        try:
            dm = model.add_boundary(m, type="demand")
        except ValueError:
            # a demand/boundary for M already exists -> reuse it
            dm = next((r for r in model.reactions
                       if len(r.metabolites) == 1 and m in r.metabolites), None)
            if dm is None:
                return 0.0
        model.objective = dm
        v = model.slim_optimize()
    return float(v) if v is not None and v == v else 0.0  # v==v filters NaN


def pick_compound(model, rng: random.Random, compartment: str | None,
                  eps: float, producible_only: bool, max_tries: int):
    """Return (mid, compartment, producible, max_flux). If producible_only,
    resample up to max_tries until a producible compound is found."""
    pool = [m for m in model.metabolites
            if compartment is None or m.compartment == compartment]
    if not pool:
        pool = list(model.metabolites)
    rng.shuffle(pool)

    tries = max_tries if producible_only else 1
    last = None
    for m in pool[:tries]:
        v = max_demand_flux(model, m.id, eps)
        prod = abs(v) > eps
        last = (m.id, m.compartment or "?", prod, v)
        if prod or not producible_only:
            return last
    return last  # exhausted tries; return the last (non-producible) one


def main() -> None:
    ap = argparse.ArgumentParser(description="FASTCC over all organisms + a random compound each.")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed (reproducible compound picks).")
    ap.add_argument("--epsilon", type=float, default=FASTCC_EPS, help="FASTCC/producibility threshold.")
    ap.add_argument("--compartment", default=None,
                    help="restrict the random compound to this compartment (e.g. c). Default: any.")
    ap.add_argument("--producible-only", action="store_true",
                    help="resample until a producible compound is found (up to --max-tries).")
    ap.add_argument("--max-tries", type=int, default=25, help="resample budget for --producible-only.")
    args = ap.parse_args()

    eps = args.epsilon
    rng = random.Random(args.seed)
    catalogue = load_catalogue()

    print("=" * 100)
    print(f"FASTCC over {len(catalogue)} organisms  |  eps={eps}  seed={args.seed}"
          f"  compartment={args.compartment or 'any'}  producible_only={args.producible_only}")
    print("=" * 100)
    hdr = (f"{'model':9} {'organism':26} {'rxns':>5} {'consist':>7} {'%':>4}  "
           f"{'compound':16} {'cmp':>3} {'prod?':>5} {'maxflux':>10}  {'#cons_rxns':>9}")
    print(hdr)
    print("-" * len(hdr))

    n_prod = 0
    for key, info in catalogue.items():
        path = PROJECT_ROOT / info["file_path"]
        model = cobra.io.load_json_model(str(path))
        n = len(model.reactions)

        t0 = time.time()
        res, consistent_ids = fastcc_from_cobra(model, eps)
        consistent = set(consistent_ids)
        dt = time.time() - t0

        mid, cmp, prod, v = pick_compound(
            model, rng, args.compartment, eps, args.producible_only, args.max_tries)
        n_prod += int(prod)

        m_cons = len({r.id for r in model.metabolites.get_by_id(mid).reactions} & consistent)

        print(f"{key:9} {info['organism'][:25]:26} {n:>5} {len(consistent):>7} "
              f"{100*len(consistent)//max(n,1):>3}%  "
              f"{mid[:16]:16} {cmp[:3]:>3} {('yes' if prod else 'NO'):>5} {v:>10.4f}  "
              f"{m_cons:>9}   [{dt:.1f}s]")

    print("-" * len(hdr))
    print(f"producible compounds: {n_prod}/{len(catalogue)}")


if __name__ == "__main__":
    main()
