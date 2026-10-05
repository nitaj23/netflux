"""
scripts/populate_test_db.py

End-to-end ingestion test over EVERY GEM in gem_models/metadata.json.

Builds a throwaway ("fake") SQLite database from scratch, registers each model
whose JSON exists on disk, ingests it under a single default-medium condition
(empty exchange bounds -> each model uses its own native medium, so FBA is
feasible regardless of organism), then verifies the database was filled
correctly.

Everything lands in scripts/_test_run/ (DB + graphs), which is wiped at the
start of each run — it never touches the canonical mfghub.db or graphs/.

Usage
-----
    python -m scripts.populate_test_db
    python -m scripts.populate_test_db --keep      # don't wipe previous run first
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker

from src.schema import (
    Base, Model, Reaction, Metabolite, Gene,
    ModelReaction, ModelMetabolite, ModelGene, ReactionGene, ReactionMetabolite,
    Condition, MFGEntry, EntryReaction, EntryMetabolite,
)
from src.ingest import ingest_gem

METADATA_PATH = PROJECT_ROOT / "gem_models" / "metadata.json"
TEST_DIR = PROJECT_ROOT / "scripts" / "_test_run"
DB_PATH = TEST_DIR / "test_mfghub.db"
GRAPHS_DIR = TEST_DIR / "graphs"
CONDITION_ID = "default_medium"


def setup(keep: bool):
    if TEST_DIR.exists() and not keep:
        shutil.rmtree(TEST_DIR)
    TEST_DIR.mkdir(parents=True, exist_ok=True)
    engine = create_engine(f"sqlite:///{DB_PATH}")
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)


def seed_condition(db):
    if not db.query(Condition).filter_by(condition_id=CONDITION_ID).first():
        db.add(Condition(
            condition_id=CONDITION_ID,
            carbon_source=None,
            oxygen_status="default",
            exchange_bounds_json=json.dumps({}),  # use each model's native medium
        ))
        db.commit()


def register_and_ingest(db) -> list[dict]:
    catalogue = json.loads(METADATA_PATH.read_text())
    results = []

    for slug, entry in catalogue.items():
        model_path = PROJECT_ROOT / entry["file_path"]
        rec = {"slug": slug, "organism": entry["organism"], "status": None,
               "seconds": 0.0, "detail": ""}

        if not model_path.exists():
            rec["status"] = "MISSING_FILE"
            rec["detail"] = str(model_path)
            results.append(rec)
            print(f"[skip] {slug:14} file not found")
            continue

        # Register model row (idempotent on file_hash)
        file_hash = hashlib.sha256(model_path.read_bytes()).hexdigest()
        if not db.query(Model).filter_by(file_hash=file_hash).first():
            db.add(Model(
                organism=entry["organism"],
                source_db="BiGG",
                model_path=str(model_path),
                file_hash=file_hash,
            ))
            db.commit()

        print(f"[ingest] {slug:14} ({entry['organism']}) ...", end="", flush=True)
        t0 = time.time()
        try:
            res = ingest_gem(str(model_path), db, CONDITION_ID, output_root=str(GRAPHS_DIR))
            rec["seconds"] = round(time.time() - t0, 1)
            rec["status"] = "SKIPPED" if res.get("skipped") else "OK"
            rec["detail"] = (f"nodes={res.get('n_nodes')} edges={res.get('n_edges')} "
                             f"active={res.get('n_active_reactions')} "
                             f"obj={res.get('fba_objective_value'):.4f}"
                             if not res.get("skipped") else "already present")
            print(f" {rec['status']} ({rec['seconds']}s)  {rec['detail']}")
        except Exception as exc:  # keep going; one bad model shouldn't stop the run
            db.rollback()
            rec["seconds"] = round(time.time() - t0, 1)
            rec["status"] = "FAILED"
            rec["detail"] = f"{type(exc).__name__}: {exc}"
            print(f" FAILED ({rec['seconds']}s)  {rec['detail']}")
            traceback.print_exc()
        results.append(rec)

    return results


def verify(db, results: list[dict]) -> bool:
    print("\n" + "=" * 70)
    print("VERIFICATION")
    print("=" * 70)

    ok = True

    # Global table counts
    counts = {
        "models": db.query(func.count()).select_from(Model).scalar(),
        "reactions": db.query(func.count()).select_from(Reaction).scalar(),
        "metabolites": db.query(func.count()).select_from(Metabolite).scalar(),
        "genes": db.query(func.count()).select_from(Gene).scalar(),
        "model_reactions": db.query(func.count()).select_from(ModelReaction).scalar(),
        "model_metabolites": db.query(func.count()).select_from(ModelMetabolite).scalar(),
        "model_genes": db.query(func.count()).select_from(ModelGene).scalar(),
        "reaction_genes": db.query(func.count()).select_from(ReactionGene).scalar(),
        "reaction_metabolites": db.query(func.count()).select_from(ReactionMetabolite).scalar(),
        "mfg_entries": db.query(func.count()).select_from(MFGEntry).scalar(),
        "entry_reactions": db.query(func.count()).select_from(EntryReaction).scalar(),
        "entry_metabolites": db.query(func.count()).select_from(EntryMetabolite).scalar(),
    }
    print("\nTable row counts:")
    for k, v in counts.items():
        print(f"  {k:20} {v}")

    n_ingested = sum(1 for r in results if r["status"] == "OK")
    n_failed = sum(1 for r in results if r["status"] == "FAILED")
    n_missing = sum(1 for r in results if r["status"] == "MISSING_FILE")

    # Per-model integrity: every ingested model must have an MFG entry plus
    # populated structural + functional membership rows.
    print("\nPer-model integrity:")
    for m in db.query(Model).all():
        n_mr = db.query(ModelReaction).filter_by(model_id=m.id).count()
        n_mm = db.query(ModelMetabolite).filter_by(model_id=m.id).count()
        n_mg = db.query(ModelGene).filter_by(model_id=m.id).count()
        entry = db.query(MFGEntry).filter_by(model_id=m.id).first()
        n_er = db.query(EntryReaction).filter_by(entry_id=entry.entry_id).count() if entry else 0

        problems = []
        if not entry:               problems.append("no MFGEntry")
        if n_mr == 0:               problems.append("no model_reactions")
        if n_mm == 0:               problems.append("no model_metabolites")
        if n_er == 0:               problems.append("no entry_reactions")
        # model_genes may legitimately be 0 (e.g. iRP911 has 0 genes in metadata)
        status = "OK" if not problems else "FAIL: " + ", ".join(problems)
        if problems:
            ok = False
        print(f"  id={m.id:2} {m.organism[:34]:34} "
              f"mr={n_mr:5} mm={n_mm:5} mg={n_mg:5} er={n_er:5}  {status}")

    # Referential sanity: no reaction_genes / entry rows pointing nowhere
    orphan_rg = (db.query(ReactionGene)
                   .outerjoin(Gene, Gene.id == ReactionGene.gene_id)
                   .filter(Gene.id.is_(None)).count())
    orphan_er = (db.query(EntryReaction)
                   .outerjoin(MFGEntry, MFGEntry.entry_id == EntryReaction.entry_id)
                   .filter(MFGEntry.entry_id.is_(None)).count())
    print(f"\nOrphan checks:  reaction_genes->genes: {orphan_rg}   "
          f"entry_reactions->mfg_entries: {orphan_er}")
    if orphan_rg or orphan_er:
        ok = False

    # Reaction global de-dup: bigg_id should be unique (shared across models)
    distinct_rxn = db.query(func.count(func.distinct(Reaction.bigg_id))).scalar()
    if distinct_rxn != counts["reactions"]:
        ok = False
        print(f"  ! reactions.bigg_id not unique: {distinct_rxn} distinct of {counts['reactions']}")

    print("\n" + "-" * 70)
    print(f"Models: {n_ingested} ingested, {n_failed} failed, {n_missing} missing-file")
    print(f"RESULT: {'PASS' if ok and n_failed == 0 else 'CHECK'}")
    print("-" * 70)
    return ok and n_failed == 0


def main():
    parser = argparse.ArgumentParser(description="Ingest every GEM into a throwaway test DB and verify.")
    parser.add_argument("--keep", action="store_true", help="Don't wipe scripts/_test_run/ first.")
    args = parser.parse_args()

    Session = setup(args.keep)
    with Session() as db:
        seed_condition(db)
        results = register_and_ingest(db)
        verify(db, results)

    print(f"\nTest DB: {DB_PATH}")


if __name__ == "__main__":
    main()
