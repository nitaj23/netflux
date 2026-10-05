"""
schema.py - Submission tables and the DB-enforced uniqueness rules.

**Additive by construction.** ``src/schema.py`` is not edited. The
per-model publication state the spec puts on ``models`` lives here instead, in
a side table keyed by ``model_id``:

    models (untouched)  1 --- 0..1  submission_models  1 --- 1  submissions

A model with no ``submission_models`` row is an ordinary curated entry
registered by ``build_db`` -- visible, current, published. That is what makes
the default "published" semantics free: the 11 curated rows need no backfill.

Everything the spec asks the database to arbitrate is an index here, not an
application check (§3.1):

* one model row per submission        -> retry cannot create a second (§5.2.1)
* one live submission per upload/file -> resubmission after rejection stays legal
* one publication at a time           -> serializes the catalogue-insert race (U-7)
* one current version per base_name   -> "two current versions" unrepresentable
* (base_name, version) unique         -> version allocation cannot race (§5.3.1)

Partial and expression indexes are used by both backends the repo targets
(SQLite >= 3.8, PostgreSQL), so the same DDL serves ``mfghub.db`` and a
Supabase ``DATABASE_URL``.
"""

from __future__ import annotations

from datetime import date

from sqlalchemy import (Boolean, Date, Float, ForeignKey, Integer, Text,
                        inspect, text)
from sqlalchemy.orm import Mapped, mapped_column

from src.schema import Base

# ── FSM states (§5.1) ────────────────────────────────────────────────────────
PENDING = "pending_upload"
VALIDATING = "validating"
VALIDATION_FAILED = "validation_failed"
AWAITING_REVIEW = "awaiting_review"
APPROVED = "approved"
PUBLISHING = "publishing"
PUBLISHED = "published"
REJECTED = "rejected"

#: Statuses that still reserve their upload/canonical hash. A rejected or
#: failed submission releases both, so an appeal can resubmit the same file.
LIVE_STATUSES = (PENDING, VALIDATING, AWAITING_REVIEW, APPROVED, PUBLISHING, PUBLISHED)

CLAIMED = "claimed"


class Submission(Base):
    """One community upload and everything known about it."""

    __tablename__ = "submissions"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    submission_id: Mapped[str] = mapped_column(Text, unique=True, nullable=False)

    # identity
    model_name: Mapped[str] = mapped_column(Text, nullable=False)
    organism: Mapped[str] = mapped_column(Text, nullable=False)
    slug: Mapped[str] = mapped_column(Text, nullable=False)

    # provenance (unverified -- there are no user accounts)
    source: Mapped[str | None] = mapped_column(Text)
    publication_doi: Mapped[str | None] = mapped_column(Text)
    license: Mapped[str | None] = mapped_column(Text)
    submitter_email: Mapped[str | None] = mapped_column(Text)
    trophic_mode: Mapped[str | None] = mapped_column(Text)

    # payload
    upload_format: Mapped[str] = mapped_column(Text, nullable=False)
    upload_path: Mapped[str] = mapped_column(Text, nullable=False)
    upload_hash: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    model_path: Mapped[str] = mapped_column(Text, nullable=False)
    file_hash: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    structure_fp: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    config_fp: Mapped[str] = mapped_column(Text, nullable=False)
    fingerprint_version: Mapped[str] = mapped_column(Text, nullable=False)
    conditions_json: Mapped[str] = mapped_column(Text, nullable=False, default="[]")

    # state
    status: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    validation_report_json: Mapped[str | None] = mapped_column(Text)
    review_notes: Mapped[str | None] = mapped_column(Text)
    allow_version: Mapped[bool] = mapped_column(Boolean, default=False)
    # Saga marker: a HINT that lets a resume skip re-checking. Never the truth --
    # recovery rebuilds progress from the DB and filesystem (§5.2.2).
    publish_progress_json: Mapped[str | None] = mapped_column(Text)
    publish_started_at: Mapped[float | None] = mapped_column(Float)
    publish_heartbeat_at: Mapped[float | None] = mapped_column(Float)

    submitted_at: Mapped[date] = mapped_column(Date, default=date.today)
    reviewed_at: Mapped[date | None] = mapped_column(Date)
    published_model_id: Mapped[int | None] = mapped_column(ForeignKey("models.id"))
    # No published_entry_id: one submission with N conditions publishes N
    # entries. The relation is mfg_entries(model_id, condition_id_fk) (§5.4).


class SubmissionModel(Base):
    """Publication state for a model created by a submission.

    Exists so ``models`` stays untouched. A row here means "this model row was
    created by the submission pipeline"; its absence means "curated entry".
    """

    __tablename__ = "submission_models"

    model_id: Mapped[int] = mapped_column(ForeignKey("models.id"), primary_key=True)
    submission_id: Mapped[str] = mapped_column(
        ForeignKey("submissions.submission_id"), nullable=False)
    file_hash: Mapped[str] = mapped_column(Text, nullable=False)
    name: Mapped[str] = mapped_column(Text, nullable=False)
    base_name: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    version: Mapped[int] = mapped_column(Integer, nullable=False, default=1)
    is_current: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    #: 'claimed' until every promotion assertion holds, then 'published' (§5.5).
    publication_state: Mapped[str] = mapped_column(Text, nullable=False, default=CLAIMED)
    supersedes: Mapped[int | None] = mapped_column(ForeignKey("models.id"))


_LIVE = ", ".join(f"'{s}'" for s in LIVE_STATUSES)

#: Raw DDL rather than ``Index(...)`` objects: partial and expression indexes
#: need per-dialect kwargs in SQLAlchemy, and this form is identical on both
#: backends. Every statement here was exercised against SQLite 3.45.3.
_INDEXES = (
    # One model row per submission -- the retry guard (§5.2.1).
    "CREATE UNIQUE INDEX IF NOT EXISTS uq_submodels_submission "
    "ON submission_models (submission_id)",
    # A canonical artifact is claimed once.
    "CREATE UNIQUE INDEX IF NOT EXISTS uq_submodels_file_hash "
    "ON submission_models (file_hash)",
    # Case-insensitive: gem_models/ is an NTFS directory, so iYO844 and iyo844
    # are the same file and must not be two catalogue rows.
    "CREATE UNIQUE INDEX IF NOT EXISTS uq_submodels_name_nocase "
    "ON submission_models (lower(name))",
    # Version allocation cannot race (§5.3.1).
    "CREATE UNIQUE INDEX IF NOT EXISTS uq_submodels_base_version "
    "ON submission_models (base_name, version)",
    # Exactly one current version per lineage.
    "CREATE UNIQUE INDEX IF NOT EXISTS uq_submodels_one_current "
    "ON submission_models (base_name) WHERE is_current = 1",
    # Live-only hash reservations: rejected/failed submissions are resubmittable.
    f"CREATE UNIQUE INDEX IF NOT EXISTS uq_submissions_live_upload "
    f"ON submissions (upload_hash) WHERE status IN ({_LIVE})",
    f"CREATE UNIQUE INDEX IF NOT EXISTS uq_submissions_live_file "
    f"ON submissions (file_hash) WHERE status IN ({_LIVE})",
    # Single-flight publication (U-7). Every row in this index has the same
    # status value, so the second concurrent publisher collides.
    "CREATE UNIQUE INDEX IF NOT EXISTS uq_one_publication "
    "ON submissions (status) WHERE status = 'publishing'",
)


def init_submissions_schema(engine) -> None:
    """Create the submission tables and their indexes. Safe to call repeatedly.

    Only ``submissions`` and ``submission_models`` are created -- passing
    ``tables=`` keeps ``create_all`` from touching the catalogue tables.
    """
    Base.metadata.create_all(
        engine, tables=[Submission.__table__, SubmissionModel.__table__])
    with engine.begin() as conn:
        for ddl in _INDEXES:
            conn.exec_driver_sql(ddl)


def schema_ready(engine) -> bool:
    """True when the submission tables exist (used by the CLI to fail helpfully)."""
    return inspect(engine).has_table(Submission.__tablename__)


def set_sqlite_busy_timeout(engine, ms: int = 5000) -> None:
    """Wait rather than fail on a concurrent writer (U-4).

    ``db.py`` sets WAL and foreign_keys but no busy timeout, so a second writer
    raises SQLITE_BUSY immediately. Publication is minutes long; waiting a few
    seconds for a lock is always better than failing the saga.
    """
    if engine.dialect.name != "sqlite":
        return
    with engine.begin() as conn:
        conn.exec_driver_sql(f"PRAGMA busy_timeout={int(ms)}")
        conn.execute(text("SELECT 1"))
