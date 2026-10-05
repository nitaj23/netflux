"""
scripts/validate_gudie_all_organisms.py

For EVERY organism, pin the target to a PROPER metabolite (not a cofactor/
currency/ion) and detect the organism's DEFAULT CARBON SOURCE(s) (the exchange
reaction(s) actively taking up carbon in the demand-maximised flux state), then
compare GUDIE extraction with and without those sources as a second anchor:

    single   target-seeded VUDIE only          (upstream field of the product)
    dual     bidirectional interest scoring     (product upstream field  x
             source downstream field), i.e. run_gudie(source_ids=<carbon source>)

Everything runs under GUDIE's DEFAULT hyper-parameters (GudieConfig()), so this
driver measures what the dual-anchor / bidirectional signal buys you out of the
box rather than sweeping a knob.

Does the dual bloat BUY quality?  (the question this driver answers)
-------------------------------------------------------------------
Dual-anchor extraction produces a larger H than single-anchor. Size alone is not
quality: |H|, edges, flux_cap, genes and subsys ALL grow monotonically with node
count, so a bigger graph trivially "wins" on them regardless of whether the extra
nodes are relevant. To decide whether the bloat is *useful*, this driver reports
precision metrics that FALL when dual adds junk:

    reach%   fraction of H's nodes that can actually reach the target in H.
             Off-corridor filler that doesn't feed the product drags this down.
    corr%    fraction of H on a genuine source->target corridor (reachable FROM a
             source AND able to reach the target). The core precision metric.
    ccap%    share of the active graph's total |flux| captured by CORRIDOR nodes
             only (not all of H) -- useful flux on real substrate->product routes.
    ndp      node-disjoint source->target routes in H (design redundancy: more
             independent routes => more robust knockout/overexpression scope).
    +prec    of the nodes dual ADDS over single, the fraction that land on the
             dual corridor. This is the direct "is the added bloat relevant?"
             number: high => dual spends its extra size on corridor; low => filler.

Good bloat raises corr%/ccap%/ndp/src>tgt and keeps reach% high; bad bloat inflates
|H| while reach%/corr% collapse.

Generalizability / reproducibility
-----------------------------------
One random (source, target) pair per organism cannot separate a real effect from a
lucky draw. This driver evaluates N targets per organism (--n-targets), each
demand-maximised in an isolated model context, with the carbon source auto-detected
per target's own flux state -- so BOTH endpoints vary across pairs. The summary
reports, per metric, mean(single), mean(dual), the mean delta, and a WIN RATE
(fraction of source-bearing pairs where dual >= single), so the verdict rests on the
distribution over pairs rather than any single extraction.

Speed/quality choices
---------------------
  * loads the pre-built FASTCC-reduced models from extraction/fastcc/models/
    (load_fastcc_model) instead of running fastcc live per organism -- the cached
    reductions were computed once at eps=1e-6 with growth preservation.
  * builds the MFG + FCA coupling ONCE per (organism, target), then reuses them
    across the single/dual comparison (coupling is source-anchor-independent, but
    target-specific) -- so each pair costs one FCA solve.

Run
---
    .venv/Scripts/python.exe -u -m scripts.validate_gudie_all_organisms --seed 1
    .venv/Scripts/python.exe -u -m scripts.validate_gudie_all_organisms --n-targets 5
    .venv/Scripts/python.exe -u -m scripts.validate_gudie_all_organisms --model iML1515
"""

from __future__ import annotations

import argparse
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

from src.extraction.fastcc import load_fastcc_model, available_models
from src.extraction.gudie import (
    run_gudie, GudieConfig, flux_coupling_to_target)
from src.graph_ops import build_mfg, fpe_similarity_to_target
from scripts.fastcc_all_organisms import load_catalogue, max_demand_flux

EPS = 1e-6

# currency/cofactor metabolites to exclude -- keep only real bioproduction-relevant
# compounds (mirrors the filter previously in optforce_all_organisms.py, which has
# since been removed from the repo; reproduced here so this driver is self-contained).
CURRENCY = {
    "h", "h2o", "atp", "adp", "amp", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp",
    "utp", "udp", "ump", "itp", "idp", "nad", "nadh", "nadp", "nadph", "coa",
    "accoa", "fad", "fadh2", "fmn", "fmnh2", "pi", "ppi", "pppi", "co2", "o2",
    "nh4", "so4", "so3", "h2s", "hco3", "q8", "q8h2", "mqn8", "mql8", "2dmmq8",
    "thf", "5mthf", "10fthf", "methf", "mlthf", "h2", "no", "no2", "no3", "amet",
    "ahcys", "glu__L", "akg",
}

# inorganic carbon carriers never count as a "carbon source" for design purposes.
NON_SUBSTRATE_C = {"co2", "hco3", "co", "cyan", "hcn"}


def base_id(met) -> str:
    bid = met.id
    c = met.compartment
    if c and bid.endswith("_" + c):
        bid = bid[: -(len(c) + 1)]
    return bid


def carbon(met) -> int:
    try:
        return met.elements.get("C", 0) if met.formula else 0
    except Exception:
        return 0


def is_proper(met) -> bool:
    """Real carbon metabolite (2-12 C), not a currency/cofactor/ion, not a
    macromolecule (ACP-bound / tRNA)."""
    C = carbon(met)
    if C < 2 or C > 12:
        return False
    if base_id(met) in CURRENCY:
        return False
    idl, nm = met.id.lower(), (met.name or "").lower()
    if "acp" in idl or "trna" in idl or "trna" in nm:
        return False
    return True


def load_reduced(key: str):
    """Cached FASTCC-reduced model (drop-in reduced GEM). Raises if not built."""
    if key not in available_models():
        raise FileNotFoundError(f"no cached FASTCC-reduced model for {key!r}")
    return load_fastcc_model(key)


def pick_proper_producible(model, rng, eps, max_tries=60):
    cyto = [m for m in model.metabolites if m.compartment == "c" and is_proper(m)]
    other = [m for m in model.metabolites if m.compartment != "c" and is_proper(m)]
    rng.shuffle(cyto); rng.shuffle(other)
    for m in (cyto + other)[:max_tries]:
        if abs(max_demand_flux(model, m.id, eps)) > eps:
            return m
    return None


def pick_proper_producible_many(model, rng, eps, n: int, max_tries=120) -> list:
    """Up to `n` DISTINCT proper, producible metabolites (for source/target-pair
    generalizability). Shares pick_proper_producible's candidate ordering (cytosol
    first, then other compartments) but keeps scanning past the first hit."""
    cyto = [m for m in model.metabolites if m.compartment == "c" and is_proper(m)]
    other = [m for m in model.metabolites if m.compartment != "c" and is_proper(m)]
    rng.shuffle(cyto); rng.shuffle(other)
    found: list = []
    for m in (cyto + other)[:max_tries]:
        if abs(max_demand_flux(model, m.id, eps)) > eps:
            found.append(m)
            if len(found) >= n:
                break
    return found


def carbon_uptake_sources(model, sol, eps, top: int) -> list[tuple[str, int, float]]:
    """The organism's default carbon source(s) in *this* flux state: boundary
    exchange reactions running as uptake (flux < 0) whose metabolite carries
    carbon (excluding CO2/bicarbonate), ranked by carbon uptake rate |flux|*C.
    Returns [(rxn_id, C, |flux|), ...] for the top `top` substrates."""
    cands: list[tuple[float, str, int, float]] = []
    for r in model.boundary:
        f = sol.fluxes.get(r.id, 0.0)
        if f >= -eps:                                   # keep uptake only
            continue
        C, skip = 0, False
        for m in r.metabolites:
            if base_id(m) in NON_SUBSTRATE_C:
                skip = True
            C = max(C, carbon(m))
        if skip or C <= 0:
            continue
        cands.append((abs(f) * C, r.id, C, abs(f)))
    cands.sort(reverse=True)
    return [(rid, C, flux) for _, rid, C, flux in cands[:top]]



# ---- corridor / reachability quality primitives -------------------------------
# A subgraph's *quality* for strain design is about the source->product route it
# contains, not its size. These operate on H with MFG edge directions restored
# (produced-by -> consumes), so a directed path source ... -> target is a genuine
# mass-flow corridor from substrate to product.

def _reaches_target(H, target_id) -> set:
    """Nodes of H that can reach the target along directed edges (incl. target)."""
    if target_id not in H:
        return set()
    return nx.ancestors(H, target_id) | {target_id}


def _reachable_from_sources(H, source_ids) -> set:
    """Nodes of H reachable FROM any source along directed edges (incl. sources)."""
    seen: set = set()
    for s in source_ids:
        if s in H and s not in seen:
            seen |= nx.descendants(H, s) | {s}
    return seen


def corridor_nodes(H, source_ids, target_id) -> set:
    """Nodes on a genuine source->target corridor: reachable from a source AND able
    to reach the target. Empty if no source is present or no corridor exists."""
    present = [s for s in source_ids if s in H]
    if not present or target_id not in H:
        return set()
    return _reachable_from_sources(H, present) & _reaches_target(H, target_id)


def disjoint_routes(H, source_ids, target_id) -> int:
    """Number of NODE-disjoint source->target routes in H (a super-source feeds all
    present sources so distinct-source routes count too). Corridor robustness: more
    independent routes => a wider, more fault-tolerant design scope. 0 if no path."""
    present = [s for s in source_ids if s in H and s != target_id]
    if not present or target_id not in H:
        return 0
    G2 = H.copy()
    ss = "__SUPER_SRC__"
    G2.add_node(ss)
    for s in present:
        G2.add_edge(ss, s)
    try:
        return len(list(nx.node_disjoint_paths(G2, ss, target_id)))
    except (nx.NetworkXNoPath, nx.NetworkXError):
        return 0


def subgraph_metrics(H, target_id, source_ids, sol, flux_gf, gf_n) -> dict:
    """Size, pathway-quality, corridor-precision and annotation metrics for H."""
    nodes = list(H.nodes())
    n = len(nodes)
    frac = n / max(gf_n, 1)
    cap = sum(abs(sol.fluxes.get(x, 0.0)) for x in nodes) / flux_gf
    dens = cap / frac if frac > 0 else 0.0

    genes: set[str] = set()
    subsys: set[str] = set()
    n_gpr = n_rev = 0
    for x in nodes:
        d = H.nodes[x]
        genes.update(d.get("genes") or [])
        s = d.get("subsystem") or "Unknown"
        if s and s != "Unknown":
            subsys.add(s)
        if d.get("gene_rule"):
            n_gpr += 1
        if float(d.get("lower_bound", 0.0)) < 0:
            n_rev += 1

    # substrate->product corridor completeness (mass-flow direction restored in H)
    conn, plen = 0, None
    for s in source_ids:
        if s in H and nx.has_path(H, s, target_id):
            conn += 1
            L = nx.shortest_path_length(H, s, target_id)
            plen = L if plen is None else min(plen, L)

    # ---- precision: does H's mass actually serve the product / corridor? ----
    reach_set = _reaches_target(H, target_id)
    corr_set = corridor_nodes(H, source_ids, target_id)
    reach_frac = len(reach_set) / max(n, 1)
    corr_frac = len(corr_set) / max(n, 1)
    corr_cap = sum(abs(sol.fluxes.get(x, 0.0)) for x in corr_set) / flux_gf
    ndp = disjoint_routes(H, source_ids, target_id)

    return {
        "n": n, "frac": frac, "edges": H.number_of_edges(),
        "cap": cap, "dens": dens,
        "reach": reach_frac, "corr": corr_frac, "corrcap": corr_cap, "ndp": ndp,
        "genes": len(genes), "gpr": n_gpr, "gpr_frac": n_gpr / max(n, 1),
        "subsys": len(subsys), "rev": n_rev,
        "conn": conn, "nsrc": len(source_ids), "plen": plen,
        "node_set": frozenset(nodes),   # for cross-mode added-node precision
    }


ROW_HDR = (f"   {'mode':>6} {'|H|':>4} {'H/GF':>6} {'flux_cap':>8} {'cap/sz':>6} "
           f"{'reach':>6} {'corr':>6} {'ccap':>6} {'ndp':>4} {'genes':>5} "
           f"{'subsys':>6} {'src>tgt':>7} {'plen':>4} {'+prec':>6}")


def print_row(mode: str, m: dict, add_prec: float | None, add_n: int) -> None:
    src = f"{m['conn']}/{m['nsrc']}" if m["nsrc"] else "-"
    plen = str(m["plen"]) if m["plen"] is not None else "-"
    if add_prec is None:
        ap = "-"
    else:
        ap = f"{100*add_prec:.0f}%(+{add_n})"
    print(f"   {mode:>6} {m['n']:>4} {100*m['frac']:>5.1f}% {100*m['cap']:>7.1f}% "
          f"{m['dens']:>6.2f} {100*m['reach']:>5.1f}% {100*m['corr']:>5.1f}% "
          f"{100*m['corrcap']:>5.1f}% {m['ndp']:>4} {m['genes']:>5} {m['subsys']:>6} "
          f"{src:>7} {plen:>4} {ap:>6}")


def eval_pair(reduced, met, args, cfg) -> dict | None:
    """One (source, target) pair: build MFG + coupling once, run single vs dual,
    return a record with per-mode metrics. Isolated in a model context by caller."""
    C = carbon(met)
    dm = reduced.add_boundary(reduced.metabolites.get_by_id(met.id), type="demand")
    reduced.objective = dm
    sol = reduced.optimize()
    if sol.status != "optimal":
        print(f"   target {met.id}: infeasible (status={sol.status}); skipped")
        return None

    G = build_mfg(reduced, sol)
    active = [n for n in G.nodes if G.degree(n) > 0]
    G_F = G.subgraph(active).copy()
    if dm.id not in G_F:
        print(f"   target {met.id} (DM {dm.id}) not in active MFG; skipped")
        return None

    node_order = list(G_F.nodes())
    gf_n = G_F.number_of_nodes()
    flux_gf = sum(abs(sol.fluxes.get(n, 0.0)) for n in node_order) or 1.0

    # default carbon source(s): substrate exchanges actively taken up in this state
    if args.substrate:
        srcs = []
        for s in args.substrate:
            if s not in reduced.reactions:
                continue
            mets = list(reduced.reactions.get_by_id(s).metabolites)
            srcs.append((s, max((carbon(m) for m in mets), default=0),
                         abs(sol.fluxes.get(s, 0.0))))
    else:
        srcs = carbon_uptake_sources(reduced, sol, EPS, args.max_sources)
    source_ids = [rid for rid, _, _ in srcs if rid in G_F]

    src_str = (", ".join(f"{rid}(C{Cs},|v|={fx:.2f})"
                         + ("" if rid in G_F else "[not in G_F]")
                         for rid, Cs, fx in srcs)
               if srcs else "none detected (dual == single)")
    print(f"   target: {met.id}  ({met.name})  C={C}   |G_F|={gf_n}")
    print(f"   carbon source(s): {src_str}")

    coupling = flux_coupling_to_target(reduced, dm.id, node_order, gamma=1.0)
    embedding_sim = fpe_similarity_to_target(
        G_F, dm.id, node_order, zscore=not args.fpe_raw)
    print(f"   FCA coupling {'ok' if coupling else 'FALLBACK'}")
    print(ROW_HDR)

    modes: dict[str, dict] = {}
    single_nodes: frozenset | None = None
    for mode, srcs_arg in (("single", None), ("dual", source_ids or None)):
        res = run_gudie(G_F, dm.id, model=reduced, coupling=coupling,
                        embedding_sim=(embedding_sim if cfg.lam > 0 else None),
                        source_ids=srcs_arg, config=cfg)
        m = subgraph_metrics(res.subgraph, dm.id, source_ids, sol, flux_gf, gf_n)
        modes[mode] = m

        # added-node precision: of the nodes dual adds over single, how many land on
        # the (dual) corridor -- the direct "is the extra bloat relevant?" number.
        add_prec = add_n = None
        if mode == "single":
            single_nodes = m["node_set"]
        elif source_ids and single_nodes is not None:
            added = m["node_set"] - single_nodes
            add_n = len(added)
            if add_n:
                corr = corridor_nodes(res.subgraph, source_ids, dm.id)
                add_prec = len(added & corr) / add_n
            else:
                add_prec = 1.0    # nothing added -> vacuously "clean"
        print_row(mode, m, add_prec, add_n or 0)
        if add_prec is not None:
            m["add_prec"] = add_prec
            m["add_n"] = add_n

        if not source_ids:          # dual identical to single -> skip the reprint
            modes["dual"] = m
            break

    return {"target": met.id, "modes": modes, "nsrc": len(source_ids)}


def eval_organism(key: str, rng, args, cfg) -> list[dict]:
    """Evaluate up to --n-targets (source, target) pairs for one organism, each in an
    isolated model context so demand reactions / bound changes never leak across
    pairs. Returns the list of per-pair records."""
    reduced = load_reduced(key)
    targets = pick_proper_producible_many(reduced, rng, EPS, args.n_targets)
    if not targets:
        print(f"   no proper producible metabolite found for {key}; skipped")
        return []

    records: list[dict] = []
    for i, met in enumerate(targets):
        print(f"   -- pair {i + 1}/{len(targets)} " + "-" * 40)
        t0 = time.time()
        with reduced:                                   # roll back demand rxn + bounds
            try:
                rec = eval_pair(reduced, met, args, cfg)
            except Exception as e:
                print(f"   ERROR on {met.id}: {type(e).__name__}: {str(e)[:80]}")
                rec = None
        if rec is not None:
            rec["key"] = key
            records.append(rec)
        print(f"   [{time.time() - t0:.1f}s]")
    return records



# ---- cross-pair summary helpers ----------------------------------------------
def _mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def _std(xs):
    if len(xs) < 2:
        return 0.0
    mu = _mean(xs)
    return (sum((x - mu) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5


def _summarise(records: list[dict]) -> None:
    print("\n" + "=" * 100)
    print("SUMMARY  (per (source,target) pair)   dual = bidirectional source-anchored "
          "interest scoring vs single-anchor")
    print("=" * 100)
    n_pairs = len(records)
    src_pairs = [r for r in records if r["nsrc"]]
    print(f"   pairs={n_pairs}   with a detected carbon source={len(src_pairs)}   "
          f"organisms={len({r['key'] for r in records})}")

    # size/mass metrics (grow with |H|; NOT evidence of quality on their own)
    # precision metrics (fall when dual adds off-corridor junk; the real verdict)
    print("\n   size / mass  (rise trivially with node count -- context, not verdict)")
    print(f"   {'metric':>10} {'single':>10} {'dual':>10} {'delta':>10} {'win%':>7}")
    for label, keyf, scale, pct in (
            ("H/GF", "frac", 100, True),
            ("flux_cap", "cap", 100, True),
            ("genes", "genes", 1, False),
            ("subsys", "subsys", 1, False)):
        _print_metric(label, keyf, src_pairs or records, scale, pct)

    print("\n   PRECISION  (fall when bloat is off-corridor filler -- the verdict)")
    print(f"   {'metric':>10} {'single':>10} {'dual':>10} {'delta':>10} {'win%':>7}")
    for label, keyf, scale, pct in (
            ("cap/size", "dens", 1, False),
            ("reach%", "reach", 100, True),
            ("corr%", "corr", 100, True),
            ("corrcap%", "corrcap", 100, True),
            ("ndp", "ndp", 1, False)):
        _print_metric(label, keyf, src_pairs or records, scale, pct)

    # corridor completeness + pathway depth: only over source-bearing pairs
    if src_pairs:
        for mode in ("single", "dual"):
            crate = [r["modes"][mode]["conn"] / r["modes"][mode]["nsrc"] for r in src_pairs]
            plens = [r["modes"][mode]["plen"] for r in src_pairs
                     if r["modes"][mode]["plen"] is not None]
            print(f"   {mode:>6} src>tgt corridor rate={100 * _mean(crate):>5.1f}%   "
                  f"mean plen={_mean(plens):>4.2f} (over {len(plens)} connected)")

        # the headline: how relevant is the extra size dual spends?
        aps = [r["modes"]["dual"].get("add_prec") for r in src_pairs
               if r["modes"]["dual"].get("add_prec") is not None]
        adds = [r["modes"]["dual"].get("add_n", 0) for r in src_pairs]
        if aps:
            print(f"\n   dual added-node corridor precision: mean={100 * _mean(aps):.1f}% "
                  f"(+/-{100 * _std(aps):.1f}%)   mean nodes added/pair={_mean(adds):.1f}")
            print("      (fraction of the nodes dual adds over single that lie on the "
                  "source->target corridor)")


def _print_metric(label, keyf, pairs, scale, pct) -> None:
    s = [r["modes"]["single"][keyf] for r in pairs]
    d = [r["modes"]["dual"][keyf] for r in pairs]
    delta = [dd - ss for ss, dd in zip(s, d)]
    wins = sum(1 for x in delta if x > 1e-9)
    ties = sum(1 for x in delta if abs(x) <= 1e-9)
    win_pct = 100 * (wins + 0.5 * ties) / len(delta) if delta else 0.0
    suffix = "%" if pct else " "
    fmt = (lambda x: f"{scale * x:>9.1f}{suffix}") if pct else (lambda x: f"{scale * x:>9.2f}{suffix}")
    print(f"   {label:>10} {fmt(_mean(s))} {fmt(_mean(d))} {fmt(_mean(delta))} "
          f"{win_pct:>6.0f}%")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="validate_gudie: bidirectional (dual-anchor) interest scoring vs "
                    "single-anchor, over all organisms x N targets, measuring whether "
                    "the dual bloat buys corridor quality (not just size).")
    ap.add_argument("--model", default=None, help="single metadata key (default: all).")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-targets", type=int, default=3,
                    help="distinct producible targets per organism (source/target-pair "
                         "generalizability). More pairs => stronger reproducibility signal.")
    ap.add_argument("--substrate", nargs="+", default=None,
                    help="force these exchange reaction id(s) as the carbon source "
                         "(default: auto-detect the active carbon uptake).")
    ap.add_argument("--max-sources", type=int, default=1,
                    help="how many carbon sources to anchor on, ranked by uptake rate.")
    # GUDIE default hyper-parameters (GudieConfig()); exposed for override only.
    # slots=True dataclass -> read defaults off an instance, not the class.
    dflt = GudieConfig()
    ap.add_argument("--alpha", type=float, default=dflt.alpha_r)
    ap.add_argument("--k", type=float, default=dflt.k)
    ap.add_argument("--hops", type=int, default=dflt.hops)
    ap.add_argument("--lam", type=float, default=dflt.lam,
                    help="Phase B FPE blend coefficient (default: GudieConfig default).")
    ap.add_argument("--fpe-raw", action="store_true",
                    help="use raw non-negative FPE cosine (default z-scores each channel).")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    catalogue = load_catalogue()
    if args.model:
        catalogue = {args.model: catalogue[args.model]}

    print("=" * 100)
    print(f"validate_gudie  single vs dual-anchor (bidirectional interest scoring)   "
          f"organisms={len(catalogue)}  n_targets={args.n_targets}  seed={args.seed}")
    print(f"default hyper-parameters: alpha={args.alpha} k={args.k} hops={args.hops} "
          f"lam={args.lam}   fpe={'raw' if args.fpe_raw else 'zscore'}   "
          f"max_sources={args.max_sources}")
    print("=" * 100)

    cfg = GudieConfig(alpha_r=args.alpha, k=args.k, hops=args.hops, lam=args.lam)
    records: list[dict] = []
    for key, info in catalogue.items():
        print(f"\n{'#' * 3} {key}  {info['organism']} {'#' * 3}")
        try:
            records.extend(eval_organism(key, rng, args, cfg))
        except Exception as e:
            print(f"   ERROR: {type(e).__name__}: {str(e)[:80]}")

    if records:
        _summarise(records)
    else:
        print("\nno pairs evaluated.")


if __name__ == "__main__":
    main()
