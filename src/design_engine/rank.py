"""
rank.py - Pareto ranking and finalization, shared by the unified pipeline.

Ranking uses full-model *verified* values only; discovery-model numbers never
enter here. Candidates from every discovery strategy are pooled before this
runs, so a single Pareto front is computed across the whole batch regardless of
which engine proposed each design.
"""

from __future__ import annotations

import collections
from dataclasses import replace

from .schema import InterventionCandidate


# Objectives of the Pareto vector, in order (all maximized).
PARETO_OBJECTIVES = ("verified_floor", "retained_mu",
                     "growth_robustness_ratio", "verified_ceiling")


def _q(x: float, tol: float) -> float:
    return round(x / tol) * tol


def _dominates(a: dict, b: dict, objs: list[str], tol: float) -> bool:
    """A dominates B: no worse on every objective, materially better on one."""
    no_worse = all(a[o] >= b[o] - tol for o in objs)
    strictly = any(a[o] > b[o] + tol for o in objs)
    return no_worse and strictly


def _assign_fronts(cands: list[InterventionCandidate],
                   class_tol: float) -> tuple[list[int], tuple[str, ...]]:
    # Choose objectives: drop growth_robustness_ratio if undefined for any in batch.
    objs = list(PARETO_OBJECTIVES)
    if any(c.growth_robustness_ratio is None for c in cands):
        objs = [o for o in objs if o != "growth_robustness_ratio"]
    vecs = []
    for c in cands:
        vecs.append({
            "verified_floor": _q(c.verified_floor, class_tol),
            "retained_mu": _q(c.retained_mu, class_tol),
            "growth_robustness_ratio": _q(c.growth_robustness_ratio or 0.0, class_tol),
            "verified_ceiling": _q(c.verified_ceiling, class_tol),
        })
    n = len(cands)
    fronts = [0] * n
    remaining = set(range(n))
    front_num = 1
    while remaining:
        current = []
        for i in remaining:
            if not any(_dominates(vecs[j], vecs[i], objs, class_tol)
                       for j in remaining if j != i):
                current.append(i)
        if not current:  # tolerance safety: should not happen with antisymmetric rule
            current = list(remaining)
        for i in current:
            fronts[i] = front_num
        remaining -= set(current)
        front_num += 1
    return fronts, tuple(objs)


def batch_objectives(cands: list[InterventionCandidate]) -> tuple[str, ...]:
    if any(c.growth_robustness_ratio is None for c in cands):
        return tuple(o for o in PARETO_OBJECTIVES if o != "growth_robustness_ratio")
    return PARETO_OBJECTIVES


def finalize(accepted: list[InterventionCandidate], best_k1_floor,
             mu_wt: float, class_tol: float, top_k) -> list[InterventionCandidate]:
    """Flag, Pareto-rank, sort, and number the accepted (buildable) designs."""
    if not accepted:
        return []

    # adds_nothing: a k>=2 set whose verified floor <= the best verified k=1 floor.
    # ties: >= 5 accepted candidates share an identical verified floor.
    floor_counts = collections.Counter(_q(c.verified_floor, class_tol) for c in accepted)

    rebuilt: list[InterventionCandidate] = []
    for c in accepted:
        flags = list(c.artifact_flags)
        if (c.cardinality >= 2 and best_k1_floor is not None
                and c.verified_floor <= best_k1_floor + class_tol):
            flags.append("adds_nothing")
        if floor_counts[_q(c.verified_floor, class_tol)] >= 5:
            flags.append("ties")
        ties = "ties" in flags
        collapse = "collapse" in flags
        starved = "starved" in flags
        adds_nothing = "adds_nothing" in flags
        passes = not (adds_nothing or (ties and (collapse or starved)))
        rebuilt.append(replace(c, artifact_flags=tuple(flags), passes_checks=passes))

    fronts, _objs = _assign_fronts(rebuilt, class_tol)
    for i, f in enumerate(fronts):
        rebuilt[i] = replace(rebuilt[i], pareto_front=f)

    def sort_key(c: InterventionCandidate):
        return (
            c.pareto_front,
            -_q(c.verified_floor, class_tol),
            -_q(c.retained_mu, class_tol),
            c.cardinality,
            tuple(sorted(a.identifier() for a in c.reaction_actions)),
        )

    rebuilt.sort(key=sort_key)
    ranked = [replace(c, rank=i) for i, c in enumerate(rebuilt, 1)]
    if top_k is not None:
        ranked = ranked[:top_k]
    return ranked
