"""
design_engine - Unified metabolic strain-design engine.

A single verification-first pipeline that runs one or more discovery strategies
(OptForce and/or OptKnock), pools their reaction-level proposals, and co-ranks
them on identical full-model verified metrics. Every ranked result is buildable
by default: a knockout design with a concrete gene-deletion construct that both
grows and holds its production floor on the full model. Regulatory and mixed
designs are returned as verified hypotheses, never ranked.

The shared core (GPR translation, full-model verification, import-loop audit,
Pareto ranking, reporting) is defined once here; discovery strategies are
pluggable and live under :mod:`design_engine.discovery`.

Usage
-----
>>> import cobra
>>> from design_engine import rank_intervention_candidates
>>> model = cobra.io.load_json_model("iML1515.json")
>>> res = rank_intervention_candidates(
...     model, "EX_succ_e", "BIOMASS_Ec_iML1515_core_75p37M",
...     methods=("optforce", "optknock"), max_knockouts=2)
>>> for c in res.accepted:
...     print(c.rank, c.source, c.verified_floor, c.buildable)

References
----------
Ranganathan S, Suthers PF, Maranas CD (2010) OptForce. PLoS Comput Biol 6(4).
Burgard AP, Pharkya P, Maranas CD (2003) OptKnock. Biotechnol Bioeng 84(6).
"""

from .pipeline import rank_intervention_candidates, PARETO_OBJECTIVES, DEFAULT_METHODS
from .schema import (
    ReactionAction, Proposal, GeneOption, GeneConstructResult,
    InterventionCandidate, RankingResult, SUPPORTED_ACTIONS,
)
from .gpr import annotate, minimal_gene_knockouts, GeneAnnotation
from .lpkernel import make_lp, LPKernel, CobraLP
from .report import (
    ranking_table, ranking_csv, candidate_to_dict, intervention_label,
    hypotheses_table, hypotheses_csv,
)
from .discovery import STRATEGIES, get_strategy
from . import errors

__version__ = "0.1.0"
__all__ = [
    "rank_intervention_candidates", "PARETO_OBJECTIVES", "DEFAULT_METHODS",
    "ReactionAction", "Proposal", "GeneOption", "GeneConstructResult",
    "InterventionCandidate", "RankingResult", "SUPPORTED_ACTIONS",
    "annotate", "minimal_gene_knockouts", "GeneAnnotation",
    "make_lp", "LPKernel", "CobraLP",
    "ranking_table", "ranking_csv", "candidate_to_dict", "intervention_label",
    "hypotheses_table", "hypotheses_csv",
    "STRATEGIES", "get_strategy",
    "errors",
]
