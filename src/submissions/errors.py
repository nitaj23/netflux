"""
errors.py - Typed failures for the community-submission pipeline.

Two kinds of failure exist here and they are reported differently:

* **Pipeline failures** (bad upload, illegal transition, lost race, duplicate
  model) are exceptions. Each carries an actionable ``hint`` so the caller
  never has to read a stack trace to know what to do next.

* **Validation failures** are *not* exceptions. A GEM that fails a check is a
  legitimate outcome of validation -- the report records every check and the
  submission moves to ``validation_failed`` with its reasons intact. Only a
  crash inside validation raises.

Mirrors the convention in ``optforce/errors.py``: one base class, a hint on
every instance, no error that bottoms out in a bare library message.
"""

from __future__ import annotations


class SubmissionError(Exception):
    """Base class for submission-pipeline errors. Carries an actionable hint."""

    def __init__(self, message: str, hint: str | None = None):
        super().__init__(message)
        self.hint = hint


class UnsupportedFormat(SubmissionError):
    """Upload is neither cobra-JSON nor SBML (S-1, S-4)."""


class UploadTooLarge(SubmissionError):
    """Upload exceeds the per-format size cap (V14)."""


class DuplicateSubmission(SubmissionError):
    """An identical upload or canonical model is already live (V13, S-3b/f)."""


class DuplicateModel(SubmissionError):
    """The corpus already holds this model (V13 / V13b)."""


class NameCollision(SubmissionError):
    """model_name is taken by an unrelated model (V13c)."""


class IllegalTransition(SubmissionError):
    """The submission is not in a state where this action is legal (§5.1)."""


class PublicationInProgress(SubmissionError):
    """Another publication holds the single-flight slot (U-7), or this row is live."""


class ClaimConflict(SubmissionError):
    """The file_hash claim belongs to another submission (§5.2.1)."""


class PublicationIncomplete(SubmissionError):
    """A promotion assertion failed -- the model stays invisible (§5.5)."""
