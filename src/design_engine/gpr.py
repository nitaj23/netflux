"""
gpr.py - GPR annotation and reaction-knockout -> gene-construct translation.

This is the OptForce-local genetic-interpretation layer (build prompt sect. 7-8).
Its authority is COBRA's Boolean GPR semantics: annotation and knockout
translation are driven by ``reaction.gpr`` (a ``cobra.core.gene.GPR``, whose
body is a Python ``ast`` node) and validated with ``GPR.eval(knockouts)``. No
GPR string is re-parsed by hand and no directionality is read from a GPR --
genetic logic is independent of reaction flux direction.

Two operations are exposed:

* :func:`annotate` -- report the raw rule, referenced genes, and the normalized
  AND-complex / OR-isozyme structure for any discovered reaction.
* :func:`minimal_gene_knockouts` -- for a knockout proposal, the minimal
  gene-deletion sets that make the GPR evaluate False (the minimal transversals
  of the isozyme complexes), deduplicated, capped, and each verified with
  ``GPR.eval``.

Malformed non-empty GPRs fail closed (no guess); empty GPRs are labelled
``no_gpr`` and remain valid reaction-level proposals.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from itertools import product

import cobra


# Default cap on inverse-GPR enumeration (build prompt sect. 8, 24). Bounds both
# the DNF term count and the number of returned gene-deletion alternatives.
DEFAULT_GPR_CAP = 64


@dataclass(frozen=True)
class GeneAnnotation:
    """GPR-derived context for one reaction (build prompt sect. 7)."""
    reaction_id: str
    reaction_name: str
    reaction_equation: str
    gpr_string: str
    gene_ids: tuple[str, ...]
    # Normalized Boolean structure where derivable: each isozyme is one
    # AND-complex (a tuple of gene ids). Empty when the GPR is empty or
    # unparseable.
    isozyme_complexes: tuple[tuple[str, ...], ...] = ()
    parse_ok: bool = True


def _dnf_terms(node, cap: int) -> list[frozenset[str]] | None:
    """
    Expand a cobra GPR ast body into disjunctive normal form.

    Returns a list of AND-terms (each a frozenset of gene ids); the reaction is
    active iff at least one term is fully present. Returns ``None`` for a
    malformed / unsupported node (fail closed), or when expansion would exceed
    ``cap`` terms.
    """
    if node is None:
        return []
    if isinstance(node, ast.Name):
        return [frozenset({node.id})]
    if isinstance(node, ast.BoolOp):
        child_dnfs = []
        for value in node.values:
            sub = _dnf_terms(value, cap)
            if sub is None:
                return None
            child_dnfs.append(sub)
        if isinstance(node.op, ast.Or):
            terms: list[frozenset[str]] = []
            for sub in child_dnfs:
                terms.extend(sub)
                if len(terms) > cap:
                    return None
            return terms
        if isinstance(node.op, ast.And):
            # Cartesian product of the child term lists, unioned per combination.
            combos = 1
            for sub in child_dnfs:
                combos *= max(1, len(sub))
                if combos > cap:
                    return None
            terms = []
            for combo in product(*child_dnfs):
                merged = frozenset().union(*combo) if combo else frozenset()
                terms.append(merged)
            return terms
        return None  # unknown boolean operator
    return None  # Call, Constant, etc. -> unsupported, fail closed


def _minimal_hitting_sets(
    terms: list[frozenset[str]], cap: int
) -> tuple[list[frozenset[str]], bool]:
    """
    Minimal transversals of ``terms``: minimal gene sets intersecting every
    AND-complex, i.e. every isozyme is disabled.

    Returns (sets, complete). ``complete`` is False when the branch search hit
    the node cap and the returned list may be a subset of all minimal solutions.
    """
    # Drop duplicate and superset terms: hitting the smaller term is necessary
    # and a superset adds no new coverage obligation beyond it only if it is a
    # superset -- but a superset term is *easier* to hit, so any transversal of
    # the subset also hits the superset. Removing supersets is therefore safe.
    uniq = sorted({t for t in terms if t}, key=lambda t: (len(t), sorted(t)))
    minimal_terms: list[frozenset[str]] = []
    for t in uniq:
        if not any(m <= t for m in minimal_terms):
            minimal_terms.append(t)

    solutions: list[frozenset[str]] = []
    complete = True
    node_budget = [cap * 8]

    def recurse(chosen: frozenset[str]):
        nonlocal complete
        if node_budget[0] <= 0:
            complete = False
            return
        node_budget[0] -= 1
        # First term not yet hit.
        for t in minimal_terms:
            if chosen.isdisjoint(t):
                uncovered = t
                break
        else:
            solutions.append(chosen)
            return
        for gene in sorted(uncovered):
            recurse(chosen | {gene})

    recurse(frozenset())

    # Reduce to minimal (drop any solution that is a superset of another) and dedup.
    solutions = sorted(set(solutions), key=lambda s: (len(s), sorted(s)))
    minimal: list[frozenset[str]] = []
    for s in solutions:
        if not any(m <= s for m in minimal):
            minimal.append(s)
        if len(minimal) >= cap:
            complete = False
            break
    return minimal, complete


def annotate(reaction: cobra.Reaction, cap: int = DEFAULT_GPR_CAP) -> GeneAnnotation:
    """Report GPR-derived context for one reaction (build prompt sect. 7)."""
    rule = (reaction.gene_reaction_rule or "").strip()
    genes = tuple(sorted(g.id for g in reaction.genes))
    try:
        equation = reaction.build_reaction_string(use_metabolite_names=False)
    except Exception:
        equation = ""
    if not rule:
        return GeneAnnotation(
            reaction.id, reaction.name or "", equation, "", genes,
            isozyme_complexes=(), parse_ok=True,
        )
    terms = _dnf_terms(reaction.gpr.body, cap)
    if terms is None:
        # Malformed / too complex: fail closed for structure, keep the raw rule.
        return GeneAnnotation(
            reaction.id, reaction.name or "", equation, rule, genes,
            isozyme_complexes=(), parse_ok=False,
        )
    complexes = tuple(sorted((tuple(sorted(t)) for t in terms), key=lambda c: (len(c), c)))
    return GeneAnnotation(
        reaction.id, reaction.name or "", equation, rule, genes,
        isozyme_complexes=complexes, parse_ok=True,
    )


def minimal_gene_knockouts(
    reaction: cobra.Reaction, cap: int = DEFAULT_GPR_CAP
) -> tuple[list[frozenset[str]], bool, bool]:
    """
    Minimal gene-deletion sets that make ``reaction`` inactive.

    Returns ``(sets, complete, parse_ok)``. Each set is verified with
    ``GPR.eval`` so a returned set provably disables the reaction. ``complete``
    is False if the cap was hit (translation incomplete -- keep reaction-level
    output, sect. 8). ``parse_ok`` is False for a malformed GPR (fail closed:
    no gene sets returned).
    """
    rule = (reaction.gene_reaction_rule or "").strip()
    if not rule:
        return [], True, True  # no_gpr: reaction-level proposal only
    terms = _dnf_terms(reaction.gpr.body, cap)
    if terms is None:
        return [], True, False  # parse failure -> fail closed
    if any(len(t) == 0 for t in terms):
        # A constantly-true term (shouldn't happen for a real gene rule) means
        # the reaction cannot be disabled by gene deletion. Fail closed.
        return [], True, True
    sets, complete = _minimal_hitting_sets(terms, cap)
    gpr = reaction.gpr
    verified = [s for s in sets if not gpr.eval(set(s))]
    return verified, complete, True
