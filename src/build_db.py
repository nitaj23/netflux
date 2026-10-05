"""
src/build_db.py

Build the canonical MFGHub database from the curated GEMs in gem_models/.

This is the production counterpart of scripts/populate_test_db.py: instead of a
throwaway SQLite file, it writes to the real database configured in db.py
(mfghub.db by default, or DATABASE_URL for PostgreSQL/Supabase) and stores the
generated MFGs under graphs/.

Pipeline per model:
    register (models row) -> ingest_gem (FBA -> MFG -> GraphML/PyG files)
        -> populate every catalogue + membership + entry table.

Canonical condition
-------------------
default_medium = each model's own curated medium (empty exchange override). For
a multi-organism repository this is the right baseline: a fixed M9-glucose-
aerobic condition is infeasible for the phototrophs/methanogens in the set.
Additional canonical conditions can be appended to CANONICAL_CONDITIONS.

Run
---
    python -m src.build_db                 # build / top up mfghub.db
    python -m src.build_db --reset         # drop all tables, rebuild
    python -m src.build_db --condition default_medium
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import time
from pathlib import Path

from sqlalchemy import func

from src.db import engine, Session
from src.schema import (
    Base, Model, Reaction, Metabolite, Gene,
    ModelReaction, ModelMetabolite, ModelGene, ReactionGene, ReactionMetabolite,
    Condition, MFGEntry, EntryReaction, EntryMetabolite,
)
from src.ingest import ingest_gem

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("build_db")
logging.getLogger("cobra").setLevel(logging.ERROR)  # silence solver chatter

PROJECT_ROOT = Path(__file__).resolve().parent.parent
METADATA_PATH = PROJECT_ROOT / "gem_models" / "metadata.json"
GRAPHS_ROOT = str(PROJECT_ROOT / "graphs")

# Canonical conditions seeded directly into the DB (no conditions/ JSON files).
CANONICAL_CONDITIONS = [
    {
        "condition_id": "default_medium",
        "carbon_source": None,
        "oxygen_status": "model_default",
        "exchange_bounds": {},  # {} -> keep each model's curated native medium
    },
]
DEFAULT_CONDITION_ID = "default_medium"

# metadata key -> source database (extend when non-BiGG models are added)
SOURCE_DB_OVERRIDES: dict[str, str] = {}


def init_schema(reset: bool) -> None:
    if reset:
        logger.warning("--reset: dropping all tables")
        Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)


def seed_conditions(db) -> None:
    for cond in CANONICAL_CONDITIONS:
        if db.query(Condition).filter_by(condition_id=cond["condition_id"]).first():
            continue
        db.add(Condition(
            condition_id=cond["condition_id"],
            carbon_source=cond["carbon_source"],
            oxygen_status=cond["oxygen_status"],
            exchange_bounds_json=json.dumps(cond["exchange_bounds"]),
        ))
        logger.info("Seeded condition '%s'", cond["condition_id"])
    db.commit()


def register_models(db) -> list[str]:
    """Insert a models row per GEM in metadata.json whose file exists on disk.

    Idempotent on file_hash. Returns the metadata keys that are available to
    ingest (registered now or already present).
    """
    catalogue = json.loads(METADATA_PATH.read_text())
    available: list[str] = []

    for key, entry in catalogue.items():
        model_path = PROJECT_ROOT / entry["file_path"]
        if not model_path.exists():
            logger.warning("Skipping '%s' — file not found: %s", key, model_path)
            continue

        file_hash = hashlib.sha256(model_path.read_bytes()).hexdigest()
        if not db.query(Model).filter_by(file_hash=file_hash).first():
            db.add(Model(
                organism=entry["organism"],
                ncbi_taxon_id=None,
                source_db=SOURCE_DB_OVERRIDES.get(key, "BiGG"),
                model_path=str(model_path),
                file_hash=file_hash,
            ))
            logger.info("Registered '%s' (%s)", key, entry["organism"])
        available.append(key)

    db.commit()
    return available


def ingest_all(db, keys: list[str], condition_ids: list[str]) -> dict:
    catalogue = json.loads(METADATA_PATH.read_text())
    summary = {"ok": 0, "skipped": 0, "failed": 0}

    for key in keys:
        model_path = str(PROJECT_ROOT / catalogue[key]["file_path"])
        # route graph files to the organism slug folder; fall back to model id
        slug = catalogue[key].get("slug")
        for cond_id in condition_ids:
            t0 = time.time()
            try:
                res = ingest_gem(model_path, db, cond_id,
                                 output_root=GRAPHS_ROOT, output_subdir=slug)
            except Exception as exc:
                db.rollback()
                summary["failed"] += 1
                logger.error("FAILED %s / %s: %s", key, cond_id, exc)
                continue
            dt = time.time() - t0
            if res.get("skipped"):
                summary["skipped"] += 1
                logger.info("skip  %-14s %-14s already present", key, cond_id)
            else:
                summary["ok"] += 1
                logger.info("ok    %-14s %-14s nodes=%d edges=%d active=%d obj=%.4f (%.1fs)",
                            key, cond_id, res["n_nodes"], res["n_edges"],
                            res["n_active_reactions"], res["fba_objective_value"], dt)
    return summary


def report(db) -> None:
    tables = [
        ("models", Model), ("conditions", Condition), ("reactions", Reaction),
        ("metabolites", Metabolite), ("genes", Gene),
        ("model_reactions", ModelReaction), ("model_metabolites", ModelMetabolite),
        ("model_genes", ModelGene), ("reaction_genes", ReactionGene),
        ("reaction_metabolites", ReactionMetabolite), ("mfg_entries", MFGEntry),
        ("entry_reactions", EntryReaction), ("entry_metabolites", EntryMetabolite),
    ]
    print("\n" + "=" * 50 + "\nCANONICAL DB TABLE COUNTS\n" + "=" * 50)
    for label, model in tables:
        print(f"  {label:22} {db.query(func.count()).select_from(model).scalar()}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the canonical MFGHub database.")
    parser.add_argument("--reset", action="store_true", help="Drop all tables before building.")
    parser.add_argument("--condition", action="append", dest="conditions",
                         help="Restrict ingestion to this condition_id (repeatable). "
                              f"Default: {DEFAULT_CONDITION_ID}.")
    args = parser.parse_args()

    logger.info("Target DB: %s", engine.url)
    init_schema(args.reset)

    with Session() as db:
        seed_conditions(db)
        keys = register_models(db)
        conds = args.conditions or [DEFAULT_CONDITION_ID]
        logger.info("Ingesting %d model(s) under condition(s): %s", len(keys), conds)
        summary = ingest_all(db, keys, conds)
        report(db)

    print(f"\nSummary: {summary['ok']} ingested, {summary['skipped']} skipped, "
          f"{summary['failed']} failed")
    print(f"Database: {engine.url}")
    print(f"Graphs:   {GRAPHS_ROOT}")


if __name__ == "__main__":
    main()
