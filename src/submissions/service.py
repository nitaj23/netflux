"""
service.py - The submission FSM and the publication saga.

    pending_upload -> validating -> validation_failed
                                 -> awaiting_review -> rejected
                                                    -> approved -> publishing -> published

Two properties do the heavy lifting here.

**The database arbitrates, not a prior SELECT.** Every uniqueness rule is an
index (``schema._INDEXES``) and every state change is a compare-and-swap
``UPDATE ... WHERE status = <expected>``. A duplicate query still runs, but
only to produce a friendly message; the insert is the gate.

**Publication is a saga, not a transaction.** ``ingest_gem`` commits internally
and writes files, so a publication cannot be rolled back as a unit. Instead
every step is idempotent and ordered so the cheapest-to-undo comes first:

    claim model row -> copy canonical JSON -> merge metadata.json
                    -> ingest per condition -> promote to visible

A crash anywhere is resolved by *rolling forward* (``recover_publications``),
because unwinding a committed ingest would mean deleting entry rows and graph
files. The claim is owned (``submission_models.submission_id``, unique), so a
retry recognizes its own work instead of mistaking it for someone else's
duplicate -- the failure that made the naive version fail permanently.

Nothing in ``src`` is modified: publication delegates to ``ingest_gem``
and records its own state in the side table.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import time
import uuid
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import cobra
from sqlalchemy import func, or_, update
from sqlalchemy.exc import IntegrityError

from src.ingest import ingest_gem, make_entry_id
from src.schema import (Condition, EntryMetabolite, EntryReaction,
                                MFGEntry, Metabolite, Model, Reaction)

from . import identity
from .errors import (ClaimConflict, DuplicateModel, DuplicateSubmission,
                     IllegalTransition, NameCollision, PublicationInProgress,
                     PublicationIncomplete, UnsupportedFormat, UploadTooLarge)
from .schema import (APPROVED, AWAITING_REVIEW, CLAIMED, LIVE_STATUSES,
                     PENDING, PUBLISHED, PUBLISHING, REJECTED, Submission,
                     SubmissionModel, VALIDATING, VALIDATION_FAILED)
from .validation import validate_gem

PROJECT_ROOT = Path(__file__).resolve().parent.parent

MAX_JSON_BYTES = 50 * 1024 * 1024     # matches tool.MAX_UPLOAD_BYTES
MAX_SBML_BYTES = 200 * 1024 * 1024    # SBML runs ~4x the JSON for one model
STALE_AFTER = 3600                    # well above a genome-scale ingest
SOLVER_KEY = "glpk"                   # ingest_gem's default; part of the entry id

_SLUG_RE = re.compile(r"^[a-z0-9_]+$")
_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class Paths:
    """Filesystem layout. Injectable so tests never touch the real corpus."""
    root: Path = PROJECT_ROOT

    @property
    def uploads(self) -> Path: return self.root / "uploads"
    @property
    def gem_models(self) -> Path: return self.root / "gem_models"
    @property
    def graphs(self) -> Path: return self.root / "graphs"
    @property
    def metadata(self) -> Path: return self.root / "gem_models" / "metadata.json"


DEFAULT_PATHS = Paths()


# ── helpers ──────────────────────────────────────────────────────────────────
def _cas(db, submission_id: str, expected: str, **values) -> bool:
    """Compare-and-swap a status. False means someone else got there first."""
    result = db.execute(
        update(Submission)
        .where(Submission.submission_id == submission_id, Submission.status == expected)
        .values(**values))
    db.commit()
    return result.rowcount == 1


def _get(db, submission_id: str) -> Submission:
    row = db.query(Submission).filter_by(submission_id=submission_id).one_or_none()
    if row is None:
        raise IllegalTransition(f"no submission {submission_id!r}",
                                hint="list submissions to find the id")
    return row


def _atomic_write(path: Path, data: bytes) -> None:
    """Write via a sibling temp file + os.replace: a reader sees old or new."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _raw_model_hash(row: Submission) -> str:
    """sha256 of the canonical ``model.json`` **bytes**.

    ``models.file_hash`` keeps the meaning ``build_db`` and ``ingest_gem`` give
    it -- the hash of the file on disk -- because ``ingest_gem`` finds its model
    row with exactly this value (``ingest.py:327``). The *canonical* hash, which
    ignores key order and float spelling, stays in the submission tables where
    it belongs. Publication copies these bytes verbatim, so the value computed
    here before the copy is the value the published file will hash to.
    """
    return identity.sha256_bytes(Path(row.model_path).read_bytes())


def _progress(row: Submission) -> dict:
    return json.loads(row.publish_progress_json or "{}")


def _save_progress(db, row: Submission, progress: dict) -> None:
    row.publish_progress_json = json.dumps(progress)
    db.commit()


def _heartbeat(db, row: Submission) -> None:
    """Liveness, so a slow multi-condition publish is not 'recovered' mid-run."""
    row.publish_heartbeat_at = time.time()
    db.commit()


# ── visibility (§5.5) ────────────────────────────────────────────────────────
def visible_models(db):
    """Catalogue query filtered to fully published models.

    A model with no ``submission_models`` row is a curated entry registered by
    ``build_db`` -- visible. A submission-created model is visible only once
    every promotion assertion held. Read paths inside publication
    (``ingest_gem``, ``build_db.register_models``) deliberately do *not* use
    this: they must be able to find a claim.
    """
    return (db.query(Model)
              .outerjoin(SubmissionModel, SubmissionModel.model_id == Model.id)
              .filter(or_(SubmissionModel.model_id.is_(None),
                          SubmissionModel.publication_state == PUBLISHED)))


def is_visible(db, model_id: int) -> bool:
    return visible_models(db).filter(Model.id == model_id).count() == 1


def submission_entries(db, submission_id: str) -> list[MFGEntry]:
    """Every MFG entry this submission published, one per condition (§5.4)."""
    return (db.query(MFGEntry)
              .join(SubmissionModel, SubmissionModel.model_id == MFGEntry.model_id)
              .filter(SubmissionModel.submission_id == submission_id)
              .order_by(MFGEntry.condition_id_fk).all())


# ── create ───────────────────────────────────────────────────────────────────
def create_submission(raw: bytes, meta: dict, *, db, paths: Paths = DEFAULT_PATHS
                      ) -> Submission:
    """Quarantine an upload, canonicalize it, and open a submission.

    The upload is stored verbatim as provenance; ``model.json`` beside it is
    the canonical artifact that gets validated, hashed and published.
    """
    for field in ("model_name", "organism", "slug"):
        if not meta.get(field):
            raise UnsupportedFormat(f"{field} is required",
                                    hint="model_name, organism and slug identify the entry")
    if not _NAME_RE.match(meta["model_name"]):
        raise UnsupportedFormat("model_name must match [A-Za-z0-9_.-]+",
                                hint="it becomes a gem_models/ filename")
    if not _SLUG_RE.match(meta["slug"]):
        raise UnsupportedFormat("slug must match [a-z0-9_]+",
                                hint="it becomes a graphs/<slug>/ directory")
    if not raw:
        raise UnsupportedFormat("empty upload", hint="POST the GEM file bytes")

    fmt = identity.sniff_format(raw)
    if fmt == "gzip":
        raw = identity.decompress(raw, MAX_SBML_BYTES)
        fmt = identity.sniff_format(raw)
    cap = MAX_SBML_BYTES if fmt == "sbml" else MAX_JSON_BYTES
    if len(raw) > cap:
        raise UploadTooLarge(f"upload is {len(raw)} bytes (> {cap})",
                             hint="reduce the model or submit cobra-JSON")

    submission_id = uuid.uuid4().hex
    folder = paths.uploads / submission_id
    folder.mkdir(parents=True, exist_ok=True)
    source = folder / ("source.xml" if fmt == "sbml" else "source.json")
    source.write_bytes(raw)

    try:
        model = identity.load_model(source, fmt)
    except Exception as exc:
        shutil.rmtree(folder, ignore_errors=True)
        raise UnsupportedFormat(f"could not parse the upload: {exc}",
                                hint="accepted formats are cobra-JSON and SBML") from exc

    canonical = folder / "model.json"
    cobra.io.save_json_model(model, str(canonical))
    if canonical.stat().st_size > MAX_JSON_BYTES:
        shutil.rmtree(folder, ignore_errors=True)
        raise UploadTooLarge("canonical model.json exceeds the JSON cap",
                             hint="the model is too large for the repository")

    row = Submission(
        submission_id=submission_id,
        model_name=meta["model_name"], organism=meta["organism"], slug=meta["slug"],
        source=meta.get("source"), publication_doi=meta.get("publication_doi"),
        license=meta.get("license"), submitter_email=meta.get("submitter_email"),
        trophic_mode=meta.get("trophic_mode"),
        upload_format=fmt, upload_path=str(source),
        upload_hash=identity.sha256_bytes(raw),
        model_path=str(canonical), file_hash=identity.file_hash(model),
        structure_fp=identity.structure_fingerprint(model),
        config_fp=identity.config_fingerprint(model),
        fingerprint_version=identity.FINGERPRINT_VERSION,
        conditions_json=json.dumps(meta.get("conditions") or ["default_medium"]),
        status=PENDING, allow_version=bool(meta.get("allow_version")),
        submitted_at=date.today())
    db.add(row)
    try:
        db.commit()
    except IntegrityError as exc:
        db.rollback()
        shutil.rmtree(folder, ignore_errors=True)
        raise _explain_live_duplicate(db, row) from exc
    return row


def _explain_live_duplicate(db, row: Submission) -> DuplicateSubmission:
    """Name the live submission that already holds this upload or model."""
    clash = (db.query(Submission)
               .filter(Submission.status.in_(LIVE_STATUSES),
                       or_(Submission.upload_hash == row.upload_hash,
                           Submission.file_hash == row.file_hash))
               .first())
    if clash is None:
        return DuplicateSubmission("submission rejected by a uniqueness constraint",
                                   hint="retry; another submission may have just landed")
    same = "file" if clash.upload_hash == row.upload_hash else "model"
    return DuplicateSubmission(
        f"the same {same} is already live as submission {clash.submission_id} "
        f"({clash.status})",
        hint="a rejected or validation-failed submission may be resubmitted; "
             "a live one should be corrected in place")


# ── validate ─────────────────────────────────────────────────────────────────
def validate_submission(submission_id: str, *, db, paths: Paths = DEFAULT_PATHS) -> dict:
    """Run the gate and move to awaiting_review or validation_failed."""
    row = _get(db, submission_id)
    if row.status not in (PENDING, VALIDATION_FAILED):
        raise IllegalTransition(
            f"cannot validate from {row.status!r}",
            hint="validation runs from pending_upload or validation_failed")
    if not _cas(db, submission_id, row.status, status=VALIDATING):
        raise IllegalTransition("the submission changed state concurrently",
                                hint="re-read it and retry")

    model = cobra.io.load_json_model(row.model_path)
    conditions, unknown = _resolve_conditions(db, json.loads(row.conditions_json))
    report = validate_gem(model, conditions=conditions,
                          trophic_mode=row.trophic_mode,
                          upload_format=row.upload_format)
    report.setdefault("canonicalizer", identity.CANONICALIZER_VERSION)
    if unknown:
        report["checks"].append({"id": "S-5", "name": "unknown_conditions",
                                 "severity": "warn", "passed": False,
                                 "detail": {"skipped": unknown}})
        report["warnings"].append("unknown_conditions")

    _catalogue_checks(db, row, model, report, paths)

    row = _get(db, submission_id)
    row.validation_report_json = json.dumps(report)
    row.status = AWAITING_REVIEW if report["passed"] else VALIDATION_FAILED
    db.commit()
    return report


def _resolve_conditions(db, condition_ids: list[str]) -> tuple[list[dict], list[str]]:
    resolved, unknown = [], []
    for condition_id in condition_ids:
        cond = db.query(Condition).filter_by(condition_id=condition_id).one_or_none()
        if cond is None:
            unknown.append(condition_id)
            continue
        resolved.append({"condition_id": condition_id,
                         "exchange_bounds": json.loads(cond.exchange_bounds_json or "{}")})
    return resolved, unknown


def _catalogue_checks(db, row: Submission, model: cobra.Model, report: dict,
                      paths: Paths) -> None:
    """V12 namespace overlap, V13 exact duplicate, V13b semantic, V13c name."""
    rxn_ids = [r.id for r in model.reactions]
    met_ids = [m.id for m in model.metabolites]
    known_rxn = (db.query(func.count(Reaction.id))
                   .filter(Reaction.bigg_id.in_(rxn_ids)).scalar() or 0)
    known_met = (db.query(func.count(Metabolite.id))
                   .filter(Metabolite.bigg_id.in_(met_ids)).scalar() or 0)
    report["checks"].append({
        "id": "V12", "name": "namespace_overlap", "severity": "warn", "passed": True,
        "detail": {"reactions": round(known_rxn / max(1, len(rxn_ids)), 3),
                   "metabolites": round(known_met / max(1, len(met_ids)), 3)}})

    exact = _exact_duplicate(db, row)
    _record(report, "V13", "exact_duplicate", exact is None,
            **({"model_id": exact} if exact else {}))

    twin = _structural_twin(db, row)
    if twin is not None:
        model_id, same_config = twin
        passed = bool(row.allow_version) and not same_config
        _record(report, "V13b", "semantic_duplicate", passed, model_id=model_id,
                same_configuration=same_config,
                hint="identical structure; " + (
                    "identical configuration too -- there is nothing new to publish"
                    if same_config else
                    "re-parameterized: set allow_version to publish it as a version"))

    clash = _name_clash(db, row, paths)
    _record(report, "V13c", "name_collision", clash is None,
            **({"held_by": clash} if clash else {}))


def _record(report: dict, cid: str, name: str, passed: bool, **detail) -> None:
    report["checks"].append({"id": cid, "name": name, "severity": "error",
                             "passed": passed, "detail": detail})
    if not passed:
        report["passed"] = False
        report["errors"].append(name)


def _exact_duplicate(db, row: Submission) -> int | None:
    """A model whose canonical bytes are identical (V13)."""
    claimed = (db.query(SubmissionModel)
                 .filter(SubmissionModel.file_hash == row.file_hash,
                         SubmissionModel.submission_id != row.submission_id)
                 .first())
    if claimed:
        return claimed.model_id
    curated = db.query(Model).filter(Model.file_hash == _raw_model_hash(row)).first()
    return curated.id if curated else None


def _structural_twin(db, row: Submission) -> tuple[int, bool] | None:
    """A registered model with the same semantics (V13b).

    Returns (model_id, same_config). ``config_fp`` for a registered model is
    rebuilt from ``model_reactions`` bounds only -- the objective is not in the
    catalogue schema -- which is enough for its one job: telling a
    re-parameterized copy from an identical one.
    """
    for model in visible_models(db).all():
        if identity.db_structure_fingerprint(db, model.id) != row.structure_fp:
            continue
        if model.file_hash == row.file_hash:
            return model.id, True
        same_config = _db_config_fingerprint(db, model.id) == row.config_fp
        return model.id, same_config
    return None


def _db_config_fingerprint(db, model_id: int) -> str:
    """Bounds-only config fingerprint, marked by its own prefix so it can never
    be confused with the full in-memory one."""
    from src.schema import ModelReaction
    rows = (db.query(Reaction.bigg_id, ModelReaction.lower_bound, ModelReaction.upper_bound)
              .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
              .filter(ModelReaction.model_id == model_id).all())
    parts = sorted(f"{identity._norm(b)}\x1f{identity._num(lb)}\x1f{identity._num(ub)}"
                   for b, lb, ub in rows)
    return identity.sha256_bytes(("partial\n" + "\n".join(parts)).encode())


def _name_clash(db, row: Submission, paths: Paths) -> str | None:
    """model_name held by an unrelated model (V13c), case-insensitively."""
    held = (db.query(SubmissionModel)
              .filter(func.lower(SubmissionModel.base_name) == row.model_name.lower(),
                      SubmissionModel.submission_id != row.submission_id)
              .first())
    if row.allow_version:
        return None                      # a deliberate new version of that lineage
    if held is not None:
        return held.name
    if row.model_name in _metadata(paths):
        return row.model_name            # a curated entry holds the name
    return None


# ── review ───────────────────────────────────────────────────────────────────
def approve(submission_id: str, notes: str | None = None, *, db) -> Submission:
    if not _cas(db, submission_id, AWAITING_REVIEW, status=APPROVED,
                review_notes=notes, reviewed_at=date.today()):
        raise IllegalTransition(f"cannot approve from {_get(db, submission_id).status!r}",
                                hint="only an awaiting_review submission can be approved")
    return _get(db, submission_id)


def reject(submission_id: str, notes: str | None = None, *, db) -> Submission:
    if not _cas(db, submission_id, AWAITING_REVIEW, status=REJECTED,
                review_notes=notes, reviewed_at=date.today()):
        raise IllegalTransition(f"cannot reject from {_get(db, submission_id).status!r}",
                                hint="only an awaiting_review submission can be rejected")
    return _get(db, submission_id)


# ── publication saga ─────────────────────────────────────────────────────────
def publish_submission(submission_id: str, *, db, paths: Paths = DEFAULT_PATHS) -> dict:
    """Claim -> copy -> metadata -> ingest -> promote. Idempotent throughout."""
    row = _get(db, submission_id)
    try:
        claimed = _cas(db, submission_id, APPROVED, status=PUBLISHING,
                       publish_started_at=time.time(), publish_heartbeat_at=time.time())
    except IntegrityError:
        db.rollback()
        raise PublicationInProgress(
            "another publication is already running",
            hint="publication is single-flight; retry when it finishes") from None
    if not claimed:
        raise IllegalTransition(f"cannot publish from {row.status!r}",
                                hint="approve the submission first, or it is already publishing")
    return _run_saga(db, _get(db, submission_id), paths)


def _run_saga(db, row: Submission, paths: Paths) -> dict:
    """The publication steps. Safe to re-enter after any crash."""
    try:
        progress = _progress(row)
        model, claim = _acquire_claim(db, row, paths)
        progress["model_id"] = model.id
        _save_progress(db, row, progress)

        dest = Path(model.model_path)
        _publish_file(row, dest)
        _merge_metadata(paths, claim, row, dest)

        entries = progress.get("entries", {})
        for condition_id in json.loads(row.conditions_json):
            if condition_id in entries:
                continue
            _heartbeat(db, row)
            entry_id = _ingest_condition(db, model, row, dest, condition_id, paths)
            if entry_id:
                entries[condition_id] = entry_id
                progress["entries"] = entries
                _save_progress(db, row, progress)

        _promote(db, row, model, claim, dest, paths)
    except Exception:
        # Compensation: return to approved. The claim is kept -- it is owned, so
        # the retry adopts it rather than seeing a duplicate (§5.2.1).
        db.rollback()
        _cas(db, row.submission_id, PUBLISHING, status=APPROVED)
        raise

    _cas(db, row.submission_id, PUBLISHING, status=PUBLISHED,
         published_model_id=_progress(_get(db, row.submission_id)).get("model_id"))
    fresh = _get(db, row.submission_id)
    return {"submission_id": fresh.submission_id, "status": fresh.status,
            "model_id": fresh.published_model_id,
            "entries": _progress(fresh).get("entries", {})}


def _acquire_claim(db, row: Submission, paths: Paths) -> tuple[Model, SubmissionModel]:
    """Get this submission's model row, creating it if it does not exist yet.

    Never trusts ``publish_progress_json``: a crash can land between the claim
    commit and the progress write, and the DB is the truth.
    """
    mine = (db.query(SubmissionModel)
              .filter_by(submission_id=row.submission_id).one_or_none())
    if mine is not None:
        if mine.file_hash != row.file_hash:
            raise ClaimConflict(
                f"model {mine.model_id} claims a different file_hash",
                hint="the upload changed under an existing claim; release it first")
        return db.get(Model, mine.model_id), mine

    holder = (db.query(SubmissionModel)
                .filter_by(file_hash=row.file_hash).one_or_none())
    if holder is not None:
        raise _classify_holder(db, holder)
    curated = db.query(Model).filter_by(file_hash=_raw_model_hash(row)).first()
    if curated is not None:
        raise DuplicateModel(f"already published as model_id={curated.id}",
                             hint="the corpus holds this exact model")

    for _ in range(5):
        base, version, name, supersedes = _allocate_version(db, row, paths)
        dest = paths.gem_models / f"{name}.json"
        model = Model(organism=row.organism, source_db=row.source,
                      model_path=str(dest), file_hash=_raw_model_hash(row))
        db.add(model)
        try:
            db.flush()
            claim = SubmissionModel(
                model_id=model.id, submission_id=row.submission_id,
                file_hash=row.file_hash, name=name, base_name=base, version=version,
                is_current=False, publication_state=CLAIMED, supersedes=supersedes)
            db.add(claim)
            db.commit()
            return model, claim
        except IntegrityError:
            db.rollback()                      # someone took this version; recompute
    raise PublicationInProgress("could not allocate a model version",
                                hint="concurrent publications kept winning; retry")


def _classify_holder(db, holder: SubmissionModel) -> Exception:
    """Whose claim is this, and is it live? (§5.2.1)"""
    if holder.publication_state == PUBLISHED:
        return DuplicateModel(f"already published as model_id={holder.model_id}",
                              hint="the corpus holds this exact model")
    owner = (db.query(Submission)
               .filter_by(submission_id=holder.submission_id).one_or_none())
    if owner is None or owner.status == REJECTED:
        return ClaimConflict(
            f"model_id={holder.model_id} is an orphaned claim",
            hint="release it with release_model_claim, then retry")
    if owner.status == PUBLISHING:
        return PublicationInProgress(
            f"submission {owner.submission_id} is publishing this model",
            hint="wait for it to finish or be recovered")
    return ClaimConflict(
        f"model_id={holder.model_id} is claimed by submission {owner.submission_id}",
        hint="let recovery resolve that submission first -- a live claim is never stolen")


def _allocate_version(db, row: Submission, paths: Paths
                      ) -> tuple[str, int, str, int | None]:
    """Next version for this lineage. The DB enforces it; this only proposes."""
    base = row.model_name
    siblings = (db.query(SubmissionModel)
                  .filter(func.lower(SubmissionModel.base_name) == base.lower())
                  .order_by(SubmissionModel.version.desc()).all())
    if not siblings:
        if base in _metadata(paths):
            raise NameCollision(f"model_name {base!r} is already in metadata.json",
                                hint="pick another name; first writer keeps it")
        return base, 1, base, None
    if not row.allow_version:
        raise NameCollision(
            f"model_name {base!r} is held by model_id={siblings[0].model_id}",
            hint="set allow_version to publish this as a new version")
    version = siblings[0].version + 1
    return base, version, f"{base}_v{version}", siblings[0].model_id


def _publish_file(row: Submission, dest: Path) -> None:
    """Copy the canonical JSON into gem_models/ and verify what landed.

    Both hashes are checked: the raw one because ``ingest_gem`` resolves its
    model row by it, and the canonical one because that is what the submission
    was validated and de-duplicated against.
    """
    if not dest.exists():
        _atomic_write(dest, Path(row.model_path).read_bytes())
    if identity.sha256_bytes(dest.read_bytes()) != _raw_model_hash(row) or \
            identity.file_hash(cobra.io.load_json_model(str(dest))) != row.file_hash:
        raise PublicationIncomplete(
            f"{dest.name} does not match the submitted model",
            hint="delete the file and retry; the quarantine copy is authoritative")


def _metadata(paths: Paths) -> dict:
    return json.loads(paths.metadata.read_text()) if paths.metadata.exists() else {}


def _merge_metadata(paths: Paths, claim: SubmissionModel, row: Submission,
                    dest: Path) -> None:
    """Read-merge-replace, keyed by name so a redo is a no-op (U-6)."""
    model = cobra.io.load_json_model(str(dest))
    entry = {"file_path": f"gem_models/{dest.name}", "model_id": claim.base_name,
             "slug": row.slug, "organism": row.organism, "version": str(claim.version),
             "num_metabolites": len(model.metabolites),
             "num_reactions": len(model.reactions), "num_genes": len(model.genes)}
    data = _metadata(paths)
    if data.get(claim.name) == entry:
        return
    data[claim.name] = entry
    _atomic_write(paths.metadata, (json.dumps(data, indent=1) + "\n").encode())


def _ingest_condition(db, model: Model, row: Submission, dest: Path,
                      condition_id: str, paths: Paths) -> str | None:
    """One condition through ingest_gem, with the two compensations it needs."""
    cond = db.query(Condition).filter_by(condition_id=condition_id).one_or_none()
    if cond is None:
        return None                                    # S-5: warned at validation
    # A-2: index_gem_and_mfg always writes is_latest=True and never demotes the
    # prior entry for the same (model, condition). Demote every *other* entry
    # here -- the target is named up front because make_entry_id is deterministic,
    # so this is correct whether ingest_gem creates the entry or skips it.
    target_id = make_entry_id(identity.sha256_bytes(dest.read_bytes()),
                              condition_id, SOLVER_KEY)
    db.execute(update(MFGEntry)
               .where(MFGEntry.model_id == model.id,
                      MFGEntry.condition_id_fk == cond.id,
                      MFGEntry.entry_id != target_id)
               .values(is_latest=False))
    db.commit()

    result = ingest_gem(str(dest), db, condition_id, solver_key=SOLVER_KEY,
                        output_root=str(paths.graphs), output_subdir=row.slug)
    entry_id = result["entry_id"]
    if result.get("skipped") and not _artifacts_ok(db, entry_id):
        # A-1: the entry_id guard is DB-only. An entry whose artifacts were lost
        # would otherwise be "skipped" forever and promoted with no outputs.
        _purge_entry(db, entry_id)
        result = ingest_gem(str(dest), db, condition_id,
                            output_root=str(paths.graphs), output_subdir=row.slug)
        entry_id = result["entry_id"]
    if not _artifacts_ok(db, entry_id):
        raise PublicationIncomplete(
            f"artifacts missing for entry {entry_id}",
            hint="graph files could not be produced; check disk space and retry")
    # A skipped ingest writes nothing, so the entry keeps whatever is_latest it
    # already had. Assert it here rather than assuming ingest_gem set it.
    db.execute(update(MFGEntry).where(MFGEntry.entry_id == entry_id)
               .values(is_latest=True))
    db.commit()
    return entry_id


def _artifacts_ok(db, entry_id: str) -> bool:
    entry = db.query(MFGEntry).filter_by(entry_id=entry_id).one_or_none()
    if entry is None:
        return False
    paths = (entry.graphml_path, entry.pyg_path, entry.fba_solution_path)
    return all(p and Path(p).exists() and Path(p).stat().st_size > 0 for p in paths)


def _purge_entry(db, entry_id: str) -> None:
    """Delete one broken entry and its membership rows, in one transaction.

    The only place this pipeline deletes committed rows, and it is confined to
    an entry already proven to have no artifacts.
    """
    db.query(EntryReaction).filter_by(entry_id=entry_id).delete()
    db.query(EntryMetabolite).filter_by(entry_id=entry_id).delete()
    db.query(MFGEntry).filter_by(entry_id=entry_id).delete()
    db.commit()


def _promote(db, row: Submission, model: Model, claim: SubmissionModel,
             dest: Path, paths: Paths) -> None:
    """Make the model visible -- only once every assertion holds (§5.5)."""
    if not dest.exists():
        raise PublicationIncomplete(f"{dest} is missing", hint="retry publication")
    if claim.name not in _metadata(paths):
        raise PublicationIncomplete("metadata.json entry is missing",
                                    hint="retry publication")
    for condition_id in json.loads(row.conditions_json):
        cond = db.query(Condition).filter_by(condition_id=condition_id).one_or_none()
        if cond is None:
            continue
        entry = (db.query(MFGEntry)
                   .filter_by(model_id=model.id, condition_id_fk=cond.id)
                   .order_by(MFGEntry.id.desc()).first())
        if entry is None or not _artifacts_ok(db, entry.entry_id):
            raise PublicationIncomplete(
                f"condition {condition_id!r} has no complete entry",
                hint="retry publication; recovery will resume from here")

    # One transaction: demote the previous current version, then promote this
    # one. uq_submodels_one_current makes an overlap unrepresentable.
    db.execute(update(SubmissionModel)
               .where(SubmissionModel.base_name == claim.base_name,
                      SubmissionModel.is_current.is_(True))
               .values(is_current=False))
    db.execute(update(SubmissionModel)
               .where(SubmissionModel.model_id == claim.model_id)
               .values(is_current=True, publication_state=PUBLISHED))
    db.commit()


# ── recovery ─────────────────────────────────────────────────────────────────
def recover_publications(*, db, paths: Paths = DEFAULT_PATHS,
                         stale_after: int = STALE_AFTER) -> list[dict]:
    """Roll stranded publications forward; roll back only if nothing was claimed.

    Evidence comes from the database and filesystem, never from the progress
    hint, so a crash between a commit and its progress write cannot mislead us.
    """
    outcomes = []
    now = time.time()
    for row in db.query(Submission).filter_by(status=PUBLISHING).all():
        # Explicit None checks: 0.0 is a legitimate (epoch, very stale)
        # timestamp, and `or` would silently treat it as "no heartbeat".
        beat = row.publish_heartbeat_at
        if beat is None:
            beat = row.publish_started_at
        if beat is not None and now - beat < stale_after:
            outcomes.append({"submission_id": row.submission_id, "action": "left_running"})
            continue
        claim = (db.query(SubmissionModel)
                   .filter_by(submission_id=row.submission_id).one_or_none())
        if claim is None:
            _cas(db, row.submission_id, PUBLISHING, status=APPROVED,
                 publish_progress_json=None)
            outcomes.append({"submission_id": row.submission_id, "action": "rolled_back"})
            continue
        try:
            _run_saga(db, row, paths)
            outcomes.append({"submission_id": row.submission_id, "action": "rolled_forward"})
        except Exception as exc:
            outcomes.append({"submission_id": row.submission_id,
                             "action": "failed", "error": str(exc)})
    return outcomes


def release_model_claim(model_id: int, *, db) -> None:
    """Free an abandoned claim so its file_hash stops being reserved.

    Refuses a published model, a claim with entries, and a claim whose owner is
    still publishing -- a live claim is never stolen.
    """
    claim = db.query(SubmissionModel).filter_by(model_id=model_id).one_or_none()
    if claim is None:
        raise ClaimConflict(f"model {model_id} was not created by a submission",
                            hint="curated models are not releasable")
    if claim.publication_state == PUBLISHED:
        raise ClaimConflict(f"model {model_id} is published",
                            hint="published models are never released")
    owner = (db.query(Submission)
               .filter_by(submission_id=claim.submission_id).one_or_none())
    if owner is not None and owner.status == PUBLISHING:
        raise PublicationInProgress(f"submission {owner.submission_id} is publishing",
                                    hint="wait for it to finish or be recovered")
    if db.query(MFGEntry).filter_by(model_id=model_id).count():
        raise ClaimConflict(f"model {model_id} has MFG entries",
                            hint="a partially ingested model must be recovered, not released")
    db.delete(claim)
    db.query(Model).filter_by(id=model_id).delete()
    db.commit()
