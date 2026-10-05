# src/extraction/repair.py
"""
Connectivity Repair: post-GUDIE source->target corridor completion.

Output-level post-processor that runs AFTER GUDIE localization and BEFORE FVA
design-annotation (see Connectivity_Repair_spec_v1.md). It takes GUDIE's subgraph
H and the active mass-flow graph G_F it was carved from and -- only when H has no
directed source->target route -- adds back the minimum-cost flux-carrying path
that reconnects them. It never removes a GUDIE-selected reaction and never fires
when H is already connected.

The reconnection is a two-tier flux-weighted shortest path (Dijkstra):
    w(u,v) = 0                        if both u,v already in H  (hug GUDIE's pick)
    w(u,v) = log(mf_max / mf(u,v))    otherwise                (dominant-flux route
                                                                is cheapest; Croes
                                                                2006, MFG weighting)
Zeroing in-H edges makes Dijkstra detour out of H only across the actual hole, so
the number of NEW reactions is minimized -- not merely the geometric path length.
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Sequence

import networkx as nx

from src.extraction.gudie import GudieConfig, GudieResult

logger = logging.getLogger(__name__)

__all__ = ["RepairResult", "repair_connectivity"]


@dataclass(slots=True)
class RepairResult:
    nodes: frozenset[str]                    # H+ = H union V(path) (superset of H)
    added: list[str] = field(default_factory=list)   # reactions added by repair, in path order
    n_added: int = 0
    repaired: bool = False
    path: list[str] = field(default_factory=list)     # the source->target corridor (empty if no repair)
    reason: str = ""                         # "", "already_connected", or "no_corridor_in_G_F"


def repair_connectivity(
    res: GudieResult,
    G_F: nx.DiGraph,
    source_ids: Sequence[str],
    target_id: str,
    *,
    config: GudieConfig | None = None,
) -> RepairResult:
    """
    Reconnect `target_id` to a source inside GUDIE's H when the localized subgraph
    dropped the corridor. Additive only: H subset of the returned nodes.

    Returns a RepairResult with repaired=False (nodes == H) when H already carries
    a directed source->target route ("already_connected") or when no such route
    exists anywhere in G_F ("no_corridor_in_G_F", an unreachable guard for a
    well-formed producing G_F). Otherwise repaired=True with the added reactions.
    """
    H = set(res.nodes)
    sources = [s for s in source_ids if s in G_F]

    # Step 1 -- connectivity gate (no-op guard).
    H_sub = G_F.subgraph(H)
    if target_id in H_sub and any(
        s in H_sub and nx.has_path(H_sub, s, target_id) for s in sources
    ):
        return RepairResult(nodes=frozenset(H), reason="already_connected")

    # Step 2 -- two-tier flux-weighted edge cost, evaluated on demand (never written
    # to G_F). In-H edges are free so the path reuses GUDIE's selection and pays only
    # for the detour across the hole.
    eps = (config if config is not None else GudieConfig()).compress_eps
    mf_max = max(
        (mf for _, _, mf in G_F.edges(data="mass_flow") if mf is not None), default=1.0)
    if mf_max <= 0.0:
        mf_max = 1.0

    def weight(u: str, v: str, d: dict) -> float:
        if u in H and v in H:
            return 0.0
        mf = d.get("mass_flow", 0.0)
        if mf <= 0.0:
            mf = eps                                    # guard degenerate/missing flux
        return math.log(mf_max / mf)

    # Step 3 -- cheapest source->target corridor over the full G_F.
    best_path: list[str] | None = None
    best_cost = math.inf
    for s in sources:
        try:
            cost, path = nx.single_source_dijkstra(G_F, s, target_id, weight=weight)
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            continue
        if cost < best_cost:
            best_cost, best_path = cost, path

    if best_path is None:
        return RepairResult(nodes=frozenset(H), reason="no_corridor_in_G_F")

    # Step 4 -- union path into H (additive only).
    added = [n for n in best_path if n not in H]
    logger.info("repair: reconnected %s (n_added=%d, cost=%.4f, path_len=%d)",
                target_id, len(added), best_cost, len(best_path))
    return RepairResult(
        nodes=frozenset(H | set(best_path)),
        added=added,
        n_added=len(added),
        path=list(best_path),
        repaired=True,
    )
