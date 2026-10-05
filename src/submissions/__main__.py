"""
__main__.py - Operator CLI for community submissions.

There are no user accounts, so review is an operator action taken here (or on a
deployment that is not publicly exposed) rather than an authenticated route.

    python -m submissions init
    python -m submissions submit <file> --name iYO844 --organism "B. subtilis 168" \
                                        --slug b_subtilis_168 [--condition default_medium]
    python -m submissions list [--status awaiting_review]
    python -m submissions show <id>
    python -m submissions validate <id> [--allow-version]
    python -m submissions approve <id> [--notes "..."]
    python -m submissions reject  <id> [--notes "..."]
    python -m submissions publish <id>
    python -m submissions recover [--stale-after 3600]
    python -m submissions release <model_id>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from src.db import SessionLocal, engine

from . import service
from .errors import SubmissionError
from .schema import (Submission, init_submissions_schema, schema_ready,
                     set_sqlite_busy_timeout)


def _report_lines(report: dict) -> list[str]:
    lines = []
    for check in report.get("checks", []):
        mark = "ok  " if check["passed"] else ("FAIL" if check["severity"] == "error"
                                               else "warn")
        lines.append(f"  {mark} {check['id']:<5} {check['name']:<28} "
                     f"{json.dumps(check['detail'], default=str)[:110]}")
    return lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m submissions",
                                     description="MFGHub community submissions")
    parser.add_argument("--root", type=Path,
                        help="repository root to publish into (default: the mfghub "
                             "checkout). Point it at a scratch tree to rehearse a "
                             "publication without touching the corpus.")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("init", help="create the submission tables and indexes")

    p = sub.add_parser("submit", help="quarantine a GEM (SBML or cobra-JSON)")
    p.add_argument("file")
    p.add_argument("--name", required=True)
    p.add_argument("--organism", required=True)
    p.add_argument("--slug", required=True)
    p.add_argument("--condition", action="append", dest="conditions")
    p.add_argument("--source"), p.add_argument("--doi"), p.add_argument("--license")
    p.add_argument("--email"), p.add_argument("--trophic-mode")
    p.add_argument("--allow-version", action="store_true")

    p = sub.add_parser("list"); p.add_argument("--status")
    p = sub.add_parser("show"); p.add_argument("submission_id")
    p = sub.add_parser("validate")
    p.add_argument("submission_id"); p.add_argument("--allow-version", action="store_true")
    p = sub.add_parser("approve")
    p.add_argument("submission_id"); p.add_argument("--notes")
    p = sub.add_parser("reject")
    p.add_argument("submission_id"); p.add_argument("--notes")
    p = sub.add_parser("publish"); p.add_argument("submission_id")
    p = sub.add_parser("recover"); p.add_argument("--stale-after", type=int, default=3600)
    p = sub.add_parser("release"); p.add_argument("model_id", type=int)

    args = parser.parse_args(argv)

    if args.command == "init":
        init_submissions_schema(engine)
        print(f"submission tables ready in {engine.url}")
        return 0
    if not schema_ready(engine):
        print("submission tables are missing -- run: python -m submissions init",
              file=sys.stderr)
        return 2
    set_sqlite_busy_timeout(engine)

    paths = service.Paths(root=args.root) if args.root else service.DEFAULT_PATHS
    with SessionLocal() as db:
        try:
            return _dispatch(args, db, paths)
        except SubmissionError as exc:
            print(f"error: {exc}", file=sys.stderr)
            if exc.hint:
                print(f"hint:  {exc.hint}", file=sys.stderr)
            return 1


def _dispatch(args, db, paths) -> int:
    if args.command == "submit":
        for folder in (paths.uploads, paths.gem_models, paths.graphs):
            folder.mkdir(parents=True, exist_ok=True)
        row = service.create_submission(
            Path(args.file).read_bytes(),
            {"model_name": args.name, "organism": args.organism, "slug": args.slug,
             "conditions": args.conditions, "source": args.source,
             "publication_doi": args.doi, "license": args.license,
             "submitter_email": args.email, "trophic_mode": args.trophic_mode,
             "allow_version": args.allow_version},
            db=db, paths=paths)
        print(f"{row.submission_id}  {row.status}  format={row.upload_format}")
        return 0

    if args.command == "list":
        query = db.query(Submission).order_by(Submission.id)
        if args.status:
            query = query.filter_by(status=args.status)
        for row in query.all():
            age = ""
            if row.status == "publishing" and row.publish_heartbeat_at:
                import time
                age = f"  heartbeat {int(time.time() - row.publish_heartbeat_at)}s ago"
            print(f"{row.submission_id}  {row.status:<18} {row.model_name:<16} "
                  f"{row.organism}{age}")
        return 0

    if args.command == "show":
        row = db.query(Submission).filter_by(submission_id=args.submission_id).one()
        print(f"{row.submission_id}  {row.status}")
        print(f"  model      {row.model_name} ({row.organism}, slug={row.slug})")
        print(f"  upload     {row.upload_format}  {row.upload_path}")
        print(f"  file_hash  {row.file_hash[:16]}  structure_fp {row.structure_fp[:16]}")
        print(f"  conditions {row.conditions_json}")
        if row.review_notes:
            print(f"  notes      {row.review_notes}")
        if row.validation_report_json:
            report = json.loads(row.validation_report_json)
            print(f"  validation passed={report['passed']} "
                  f"errors={report['errors']} warnings={report['warnings']}")
            print("\n".join(_report_lines(report)))
        for entry in service.submission_entries(db, row.submission_id):
            print(f"  entry      {entry.entry_id} condition_fk={entry.condition_id_fk} "
                  f"obj={entry.fba_objective_value}")
        return 0

    if args.command == "validate":
        if args.allow_version:
            row = db.query(Submission).filter_by(submission_id=args.submission_id).one()
            row.allow_version = True
            db.commit()
        report = service.validate_submission(args.submission_id, db=db, paths=paths)
        print(f"passed={report['passed']}  errors={report['errors']}  "
              f"warnings={report['warnings']}")
        print("\n".join(_report_lines(report)))
        return 0 if report["passed"] else 1

    if args.command in ("approve", "reject"):
        action = service.approve if args.command == "approve" else service.reject
        row = action(args.submission_id, args.notes, db=db)
        print(f"{row.submission_id}  {row.status}")
        return 0

    if args.command == "publish":
        result = service.publish_submission(args.submission_id, db=db, paths=paths)
        print(f"{result['submission_id']}  {result['status']}  "
              f"model_id={result['model_id']}")
        for condition_id, entry_id in (result.get("entries") or {}).items():
            print(f"  {condition_id:<18} {entry_id}")
        return 0

    if args.command == "recover":
        for outcome in service.recover_publications(db=db, paths=paths, stale_after=args.stale_after):
            print(f"{outcome['submission_id']}  {outcome['action']}"
                  f"{'  ' + outcome.get('error', '') if outcome.get('error') else ''}")
        return 0

    if args.command == "release":
        service.release_model_claim(args.model_id, db=db)
        print(f"released claim on model_id={args.model_id}")
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
