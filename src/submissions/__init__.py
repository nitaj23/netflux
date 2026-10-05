"""
submissions - Community GEM submissions for MFGHub.

A quarantined upload becomes a curated repository entry through one path:

    create -> validate -> review -> approve -> publish -> ingest_gem -> visible

The package is **additive**: ``src`` is not modified. Publication
delegates GEM -> MFG/FBA/GraphML/PyG to ``src.ingest.ingest_gem`` and
keeps its own state in two new tables, so the catalogue schema, the tool layer
and ``build_db`` are untouched.

Usage
-----
>>> from src.db import SessionLocal, engine
>>> from submissions import init_submissions_schema, create_submission
>>> init_submissions_schema(engine)
>>> with SessionLocal() as db:
...     sub = create_submission(open("iYO844.xml", "rb").read(),
...                             {"model_name": "iYO844",
...                              "organism": "Bacillus subtilis 168",
...                              "slug": "b_subtilis_168",
...                              "conditions": ["default_medium"]}, db=db)
...     report = validate_submission(sub.submission_id, db=db)
...     if report["passed"]:
...         approve(sub.submission_id, "BiGG-derived, DOI checked", db=db)
...         publish_submission(sub.submission_id, db=db)

Or from the command line:

    python -m submissions submit model.xml --name iYO844 --organism "B. subtilis" \
        --slug b_subtilis_168
    python -m submissions list
    python -m submissions validate <id>
    python -m submissions approve <id> --notes "checked"
    python -m submissions publish <id>
    python -m submissions recover
"""

from .errors import (ClaimConflict, DuplicateModel, DuplicateSubmission,
                     IllegalTransition, NameCollision, PublicationInProgress,
                     PublicationIncomplete, SubmissionError, UnsupportedFormat,
                     UploadTooLarge)
from .identity import (canon_gpr, canonical_bytes, config_fingerprint,
                       db_structure_fingerprint, file_hash,
                       structure_fingerprint)
from .schema import (Submission, SubmissionModel, init_submissions_schema,
                     schema_ready, set_sqlite_busy_timeout)
from .service import (Paths, approve, create_submission, publish_submission,
                      recover_publications, reject, release_model_claim,
                      submission_entries, validate_submission, visible_models)
from .validation import validate_gem

__all__ = [
    # schema
    "Submission", "SubmissionModel", "init_submissions_schema", "schema_ready",
    "set_sqlite_busy_timeout",
    # pipeline
    "create_submission", "validate_submission", "approve", "reject",
    "publish_submission", "recover_publications", "release_model_claim",
    "visible_models", "submission_entries", "Paths",
    # identity
    "canonical_bytes", "file_hash", "structure_fingerprint",
    "config_fingerprint", "db_structure_fingerprint", "canon_gpr",
    # validation
    "validate_gem",
    # errors
    "SubmissionError", "UnsupportedFormat", "UploadTooLarge",
    "DuplicateSubmission", "DuplicateModel", "NameCollision",
    "IllegalTransition", "PublicationInProgress", "PublicationIncomplete",
    "ClaimConflict",
]
