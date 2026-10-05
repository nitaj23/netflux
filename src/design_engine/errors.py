"""
errors.py - Typed failures for the unified design-engine ranking pipeline.

Two kinds of failure exist in the pipeline and they are handled differently
(build prompt sect. 20):

* **Base-model contract failures** are raised as exceptions and abort the whole
  operation -- there is no point ranking candidates against a model that cannot
  grow, cannot produce the target, or was called with an invalid gamma.

* **Candidate-level failures** (a single design turning out infeasible, losing
  its floor on the full model, requiring a direction flip, or failing the
  import audit) are *not* exceptions. They become a ``rejection_reason`` string
  on the candidate so one bad design never aborts an otherwise valid batch. The
  reason strings live in :data:`REJECTION_REASONS` and the audit-status strings
  in :data:`IMPORT_AUDIT_STATUSES`; both are plain constants, not classes,
  because they are data attached to a result rather than control flow.

Every exception carries an actionable ``hint`` so the message never bottoms out
in a bare solver status.
"""

from __future__ import annotations


class DesignEngineError(Exception):
    """Base class for design-engine pipeline errors. Carries an actionable hint."""

    def __init__(self, message: str, hint: str | None = None):
        self.hint = hint
        super().__init__(message if hint is None else f"{message}  Hint: {hint}")


# Backward-compatible aliases: the former per-engine base classes. Code that
# does ``except OptForceError`` / ``except OptKnockError`` keeps working, since
# both now resolve to the unified base.
OptForceError = DesignEngineError
OptKnockError = DesignEngineError


# --- Base-model contract failures (abort the whole operation) ---------------

class BaseModelInfeasible(OptForceError):
    """The full model is infeasible under its fixed medium."""


class BiomassNotGrowing(OptForceError):
    """Wild-type biomass optimum is at or below tolerance."""


class TargetNotProducible(OptForceError):
    """The target reaction cannot carry finite flux above tolerance."""


class InvalidGamma(OptForceError):
    """gamma is outside the open interval (0, 1)."""


class DiscoveryModelInvalid(OptForceError):
    """A supplied discovery model is missing biomass, target, or medium exchanges."""


class MissingFullModelReaction(OptForceError):
    """A discovered reaction id does not resolve on the full model."""


class AmbiguousTargetIdentity(OptForceError):
    """The produced target metabolite could not be resolved unambiguously."""


class SolverFailure(OptForceError):
    """An optimization returned a non-optimal, non-finite, or errored result."""


class SolverTimeout(OptForceError):
    """An optimization exceeded the configured solver timeout."""


# --- Candidate-level rejection reasons (attached to results, never raised) ---
# These are the values `InterventionCandidate.rejection_reason` may take.
CANDIDATE_INFEASIBLE = "full_model_infeasible"
MISSING_FULL_MODEL_REACTION = "missing_full_model_reaction"
FULL_MODEL_BYPASS = "full_model_bypass"
BELOW_GROWTH_FLOOR = "below_growth_floor"
BELOW_PRODUCTION_THRESHOLD = "below_production_threshold"
UNSUPPORTED_DIRECTION_FLIP = "unsupported_direction_flip"
IMPORT_LOOP = "import_loop"
GPR_PARSE_FAILURE = "gpr_parse_failure"
# A knockout design that verified at the reaction level but has no gene-deletion
# construct that both grows and holds its production floor on the full model.
# Only applied when the caller requires buildable results (require_buildable).
NOT_BUILDABLE = "not_buildable"

REJECTION_REASONS = frozenset({
    CANDIDATE_INFEASIBLE,
    MISSING_FULL_MODEL_REACTION,
    FULL_MODEL_BYPASS,
    BELOW_GROWTH_FLOOR,
    BELOW_PRODUCTION_THRESHOLD,
    UNSUPPORTED_DIRECTION_FLIP,
    IMPORT_LOOP,
    GPR_PARSE_FAILURE,
    NOT_BUILDABLE,
})

# --- Import-loop audit outcomes (sect. 13.5) --------------------------------
IMPORT_PASSES = "passes"
IMPORT_LOOP_STATUS = "import_loop"
IMPORT_AMBIGUOUS = "ambiguous_target_identity"
IMPORT_UNAVAILABLE = "audit_unavailable"

IMPORT_AUDIT_STATUSES = frozenset({
    IMPORT_PASSES, IMPORT_LOOP_STATUS, IMPORT_AMBIGUOUS, IMPORT_UNAVAILABLE,
})

# --- GPR translation statuses (sect. 7) -------------------------------------
GPR_EXACT_KO = "exact_gene_knockout"
GPR_ALTERNATIVE_KO = "alternative_gene_knockouts"
GPR_PUTATIVE = "putative_regulation"
GPR_NONE = "no_gpr"
GPR_FAILURE = "gpr_parse_failure"
GPR_REACTION_ONLY = "reaction_only"
GPR_TRANSLATION_INCOMPLETE = "translation_incomplete"
