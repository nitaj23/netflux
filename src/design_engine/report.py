"""
report.py - Display/export formatting for the unified design-engine pipeline.

:func:`rank_intervention_candidates` returns a rich :class:`RankingResult` of
dataclasses. The tool layer (and any frontend) wants two flat, self-describing
artifacts instead:

    ranking_table(result)  -> {"columns": [...], "rows": [...]}   render grid
    ranking_csv(result)    -> str                                 download blob
    candidate_to_dict(c)   -> dict                                one JSON row

These are pure functions of a :class:`RankingResult`; they never touch a model
or a solver, so they are safe to call anywhere the result has already been
computed. The accepted (ranked) candidates are what a user acts on, so they are
what the table and CSV describe; the rejected candidates remain available on the
result for diagnostics and via :func:`candidate_to_dict`.
"""

from __future__ import annotations

import csv
import io
from typing import Optional

from .schema import InterventionCandidate, RankingResult, ReactionAction

# Compact, ASCII-safe action labels (CSV / terminal friendly).
_ACTION_LABEL = {
    "upregulation": "up",
    "downregulation": "down",
    "knockout": "KO",
}


def _fmt(x: Optional[float], places: int = 4) -> Optional[float]:
    """Round a float for display; pass through None unchanged."""
    return None if x is None else round(float(x), places)


def _action_dict(a: ReactionAction) -> dict:
    return {
        "reaction_id": a.reaction_id,
        "action": a.action,
        "forced_flux": _fmt(a.forced_flux),
    }


def intervention_label(actions: tuple[ReactionAction, ...]) -> str:
    """A one-line human-readable summary, e.g. ``"PGI:up@8.2, FUM:KO"``."""
    parts = []
    for a in actions:
        tag = _ACTION_LABEL.get(a.action, a.action)
        if a.action == "knockout":
            parts.append(f"{a.reaction_id}:{tag}")
        else:
            parts.append(f"{a.reaction_id}:{tag}@{a.forced_flux:g}")
    return ", ".join(parts)


def _gene_options(c: InterventionCandidate) -> list[dict]:
    return [
        {"genes": list(o.genes), "kind": o.kind, "complete": o.complete}
        for o in c.gene_options
    ]


def _gene_constructs(c: InterventionCandidate) -> list[dict]:
    return [
        {
            "genes": list(gc.genes),
            "silenced_reactions": list(gc.silenced_reactions),
            "verified_growth": _fmt(gc.verified_growth),
            "verified_floor": _fmt(gc.verified_floor),
            "verified_ceiling": _fmt(gc.verified_ceiling),
            "retained_mu": _fmt(gc.retained_mu),
            "import_loop_status": gc.import_loop_status,
            "buildable": gc.buildable,
        }
        for gc in c.gene_constructs
    ]


def _best_gene_ko(c: InterventionCandidate) -> str:
    """A short 'genes' cell for the table: the buildable gene-deletion construct
    if one exists (that is the design the user would actually build), else the
    first realizable KO gene option."""
    for gc in c.gene_constructs:
        if gc.buildable and gc.genes:
            return " + ".join(gc.genes)
    for o in c.gene_options:
        if o.kind == "knockout" and o.genes:
            suffix = "" if o.complete else " (partial)"
            return " + ".join(o.genes) + suffix
    return ""


def candidate_to_dict(c: InterventionCandidate) -> dict:
    """Flatten one :class:`InterventionCandidate` to a JSON-serializable dict."""
    return {
        "rank": c.rank,
        "pareto_front": c.pareto_front,
        "source": c.source,
        "cardinality": c.cardinality,
        "interventions": [_action_dict(a) for a in c.reaction_actions],
        "intervention_label": intervention_label(c.reaction_actions),
        # Full-model verified metrics (the authoritative ranking values).
        "verified_floor": _fmt(c.verified_floor),
        "verified_ceiling": _fmt(c.verified_ceiling),
        "verified_growth": _fmt(c.verified_growth),
        "retained_mu": _fmt(c.retained_mu),
        "floor_delta": _fmt(c.floor_delta),
        "ceiling_delta": _fmt(c.ceiling_delta),
        "coupling_class": c.coupling_class,
        "coupling_potential": _fmt(c.coupling_potential),
        "mandatory_production": c.mandatory_production,
        "guaranteed_growth": _fmt(c.guaranteed_growth),
        "pinned_growth": _fmt(c.pinned_growth),
        "growth_robustness_ratio": _fmt(c.growth_robustness_ratio),
        # Gene translation.
        "translation_status": c.translation_status,
        "all_knockout": c.all_knockout,
        "buildable": c.buildable,
        "gene_options": _gene_options(c),
        "gene_constructs": _gene_constructs(c),
        # Audit / safety.
        "import_loop_status": c.import_loop_status,
        "target_import_reactions": list(c.target_import_reactions),
        "closed_import_floor": _fmt(c.closed_import_floor),
        "artifact_flags": list(c.artifact_flags),
        "passes_checks": c.passes_checks,
        "rejection_reason": c.rejection_reason,
        # Discovery-model diagnostics (never used for ranking).
        "discovery_floor": _fmt(c.discovery_floor),
        "discovery_growth": _fmt(c.discovery_growth),
    }


# The self-describing table the frontend renders: (key, header) in display order.
_COLUMNS: tuple[tuple[str, str], ...] = (
    ("rank", "Rank"),
    ("intervention_label", "Interventions"),
    ("cardinality", "k"),
    ("verified_floor", "Floor"),
    ("verified_ceiling", "Ceiling"),
    ("verified_growth", "Growth"),
    ("retained_mu", "Retained mu"),
    ("growth_robustness_ratio", "Robustness"),
    ("coupling_class", "Coupling"),
    ("translation_status", "Translation"),
    ("genes", "Gene KO"),
    ("flags", "Flags"),
)


def _row(c: InterventionCandidate) -> dict:
    return {
        "rank": c.rank,
        "intervention_label": intervention_label(c.reaction_actions),
        "cardinality": c.cardinality,
        "verified_floor": _fmt(c.verified_floor),
        "verified_ceiling": _fmt(c.verified_ceiling),
        "verified_growth": _fmt(c.verified_growth),
        "retained_mu": _fmt(c.retained_mu),
        "growth_robustness_ratio": _fmt(c.growth_robustness_ratio),
        "coupling_class": c.coupling_class,
        "translation_status": c.translation_status,
        "genes": _best_gene_ko(c),
        "flags": ", ".join(c.artifact_flags),
    }


def ranking_table(result: RankingResult) -> dict:
    """A self-describing grid of the ranked (accepted) candidates.

    ``{"columns": [{"key", "label"}...], "rows": [{...}...]}`` -- the frontend
    renders it directly; ``rows`` are ordered by rank (the result already is).
    """
    return {
        "kind": "designs",   # decisive, buildable, ranked designs -- never hypotheses
        "columns": [{"key": k, "label": label} for k, label in _COLUMNS],
        "rows": [_row(c) for c in result.accepted],
    }


def ranking_csv(result: RankingResult) -> str:
    """The ranked candidates as a downloadable CSV string (header + one row each)."""
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([label for _, label in _COLUMNS])
    for c in result.accepted:
        row = _row(c)
        writer.writerow([row[k] if row[k] is not None else "" for k, _ in _COLUMNS])
    return buf.getvalue()


# --- Hypotheses (verified regulatory / mixed designs, never ranked) ---------

_HYP_COLUMNS: tuple[tuple[str, str], ...] = (
    ("intervention_label", "Interventions"),
    ("cardinality", "k"),
    ("verified_floor", "Floor"),
    ("verified_ceiling", "Ceiling"),
    ("verified_growth", "Growth"),
    ("retained_mu", "Retained mu"),
    ("coupling_class", "Coupling"),
    ("translation_status", "Translation"),
    ("genes", "Putative targets"),
)


def _hyp_row(c: InterventionCandidate) -> dict:
    putative = ""
    for o in c.gene_options:
        if o.genes:
            putative = " + ".join(o.genes)
            break
    return {
        "intervention_label": intervention_label(c.reaction_actions),
        "cardinality": c.cardinality,
        "verified_floor": _fmt(c.verified_floor),
        "verified_ceiling": _fmt(c.verified_ceiling),
        "verified_growth": _fmt(c.verified_growth),
        "retained_mu": _fmt(c.retained_mu),
        "coupling_class": c.coupling_class,
        "translation_status": c.translation_status,
        "genes": putative,
    }


def hypotheses_table(result: RankingResult) -> dict:
    """A self-describing grid of the verified regulatory/mixed designs.

    These passed full-model verification and the import audit but are not
    buildable as a single gene-deletion construct, so they are offered as
    hypotheses (e.g. up/down-regulation targets) rather than ranked designs.
    """
    return {
        "kind": "hypotheses",   # regulatory leads -- separate grid, never ranked designs
        "columns": [{"key": k, "label": label} for k, label in _HYP_COLUMNS],
        "rows": [_hyp_row(c) for c in result.hypotheses],
    }


def hypotheses_csv(result: RankingResult) -> str:
    """The verified regulatory/mixed hypotheses as a downloadable CSV string."""
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([label for _, label in _HYP_COLUMNS])
    for c in result.hypotheses:
        row = _hyp_row(c)
        writer.writerow([row[k] if row[k] is not None else "" for k, _ in _HYP_COLUMNS])
    return buf.getvalue()
