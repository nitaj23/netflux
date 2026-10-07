# src/ingest.py
"""
Ingest pipeline.

Public API
----------
make_entry_id             deterministic 16-char primary key
load_model_with_condition apply DB-resolved condition bounds to cached model copy
index_gem_and_mfg         populate all 4 membership tables + MFGEntry record
ingest_gem                top-level orchestrator (condition → files → DB)

build_mfg, run_fba live in graph_ops.py and are re-exported
here for the compute.py import line:
    from src.ingest import run_fba, build_mfg
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import cobra
import networkx as nx
import numpy as np
from sqlalchemy.orm import Session

from src.graph_ops import build_mfg, run_fba, write_graphml
from src.schema import (
    Condition,
    EntryMetabolite,
    EntryReaction,
    Gene,
    MFGEntry,
    Metabolite,
    Model as DBModel,
    ModelGene,
    ModelMetabolite,
    ModelReaction,
    Reaction,
    ReactionGene,
    ReactionMetabolite,
)

logger = logging.getLogger(__name__)

# Re-export so compute.py can keep its single import line unchanged
__all__ = [
    "make_entry_id", "run_fba", "build_mfg",
    "load_model_with_condition", "index_gem_and_mfg", "ingest_gem",
]


 
# 1. Primary key
 

def make_entry_id(model_hash: str, condition_id: str, solver_key: str = "glpk") -> str:
    """
    Deterministic 16-char hex ID for a (file_hash, condition_id, solver_key) triple.
    Uses file_hash rather than slug to guard against slug reuse across Zenodo versions.
    """
    raw = json.dumps([model_hash, condition_id, solver_key], sort_keys=True)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


 
# 4. Condition application
 

def load_model_with_condition(model_path: str, condition: dict) -> cobra.Model:
    """
    Retrieve a pristine model copy from the LRU cache and apply exchange bound
    overrides. Uses get_cached_model() which calls .copy() — the cached base
    model is never mutated. Heirendt et al. (2019, Nature Protocols).
    """
    model = cobra.io.load_json_model(model_path)
    for rxn_id, bounds in condition.get("exchange_bounds", {}).items():
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            lb, ub = bounds
            rxn.lower_bound, rxn.upper_bound = lb, ub
    return model


 
# 5. DB population
 

def index_gem_and_mfg(
    model:            cobra.Model,
    solution:         cobra.Solution,
    G:                nx.DiGraph,
    model_db_id:      int,
    entry_id:         str,
    condition_id_fk:  int,
    db:               Session,
    *,
    graphml_path:     str | None = None,
    pyg_path:         str | None = None,
    fba_solution_path: str | None = None,
) -> None:
    """
    Populate MFGEntry + membership tables in a single transaction:
    model_reactions, model_metabolites, model_genes, reaction_genes
    (GEM-structural, flux-independent) and entry_reactions, entry_metabolites
    (MFG-functional, flux-derived).

    n_reactions uses the flux-threshold definition (|flux| > 1e-9).

    """

    # ── 0. MFGEntry ──────────────────────────────────────────────────────────
    active_node_count = sum(
        1 for n in G.nodes if abs(G.nodes[n].get("flux", 0.0)) > 1e-9
    )
    # Atomically demote any prior latest entry for this (model, condition): the
    # new row below becomes the sole is_latest=True, in the same transaction so a
    # reader never sees two latests or none. A NULL condition is matched as NULL.
    cond_filter = (MFGEntry.condition_id_fk == condition_id_fk if condition_id_fk is not None
                   else MFGEntry.condition_id_fk.is_(None))
    (db.query(MFGEntry)
       .filter(MFGEntry.model_id == model_db_id, cond_filter,
               MFGEntry.is_latest == True)
       .update({MFGEntry.is_latest: False}, synchronize_session=False))
    mfg_entry = MFGEntry(
        entry_id=entry_id,
        model_id=model_db_id,
        condition_id_fk=condition_id_fk,
        fba_objective_value=solution.objective_value,
        n_reactions=active_node_count,
        n_edges=G.number_of_edges(),
        graphml_path=graphml_path,
        pyg_path=pyg_path,
        fba_solution_path=fba_solution_path,
        is_latest=True,
    )
    db.add(mfg_entry)
    db.flush()

    # ── 1. Reactions — bulk pre-fetch to avoid N+1 ───────────────────────────
    model_rxn_ids   = [r.id for r in model.reactions]
    existing_rxns   = {
        r.bigg_id: r.id
        for r in db.query(Reaction).filter(Reaction.bigg_id.in_(model_rxn_ids)).all()
    }
    existing_model_rxns = {
        r.reaction_id
        for r in db.query(ModelReaction.reaction_id)
                    .filter_by(model_id=model_db_id).all()
    }

    # 1a. Genes 
    model_gene_ids = [g.id for g in model.genes]
    existing_genes = {
        g.bigg_id: g.id
        for g in db.query(Gene).filter(Gene.bigg_id.in_(model_gene_ids)).all()
    } if model_gene_ids else {}
    existing_model_genes = {
        g.gene_id
        for g in db.query(ModelGene.gene_id).filter_by(model_id=model_db_id).all()
    }
    # reaction_genes is global (model-agnostic). Only reactions already in the
    # DB can have pre-existing links; new reactions never do.
    existing_rxn_gene_pairs = {
        (rg.reaction_id, rg.gene_id)
        for rg in db.query(ReactionGene.reaction_id, ReactionGene.gene_id)
                    .filter(ReactionGene.reaction_id.in_(list(existing_rxns.values()))).all()
    } if existing_rxns else set()
    # reaction_metabolites is model-specific (stoichiometry/direction varies per
    # GEM), so the dedup guard is scoped to THIS model — not global.
    existing_rxn_met_pairs = {
        (rm.reaction_id, rm.metabolite_id)
        for rm in db.query(ReactionMetabolite.reaction_id, ReactionMetabolite.metabolite_id)
                    .filter_by(model_id=model_db_id).all()
    }

    for gene in model.genes:
        gene_db_id = existing_genes.get(gene.id)
        if not gene_db_id:
            # Safe annotation parsing — values can be list or str (or absent)
            annot = gene.annotation or {}
            ncbi = annot.get("ncbigene")
            ncbi = ncbi[0] if isinstance(ncbi, list) else ncbi
            uniprot = annot.get("uniprot")
            uniprot = uniprot[0] if isinstance(uniprot, list) else uniprot

            new_gene = Gene(
                bigg_id=gene.id,
                name=gene.name or None,
                ncbigene_id=ncbi,
                uniprot_id=uniprot,
            )
            db.add(new_gene)
            db.flush()
            gene_db_id = existing_genes[gene.id] = new_gene.id

        if gene_db_id not in existing_model_genes:
            db.add(ModelGene(model_id=model_db_id, gene_id=gene_db_id))
            existing_model_genes.add(gene_db_id)

    for rxn in model.reactions:
        rxn_db_id = existing_rxns.get(rxn.id)
        if not rxn_db_id:
            new_rxn = Reaction(
                bigg_id=rxn.id,
                name=rxn.name,
                subsystem=rxn.subsystem or "Unknown",
                reaction_formula=rxn.reaction,
                gene_rule=rxn.gene_reaction_rule,   
            )
            db.add(new_rxn)
            db.flush()
            rxn_db_id = existing_rxns[rxn.id] = new_rxn.id

        if rxn_db_id not in existing_model_rxns:
            db.add(ModelReaction(
                model_id=model_db_id, reaction_id=rxn_db_id,
                lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound,
            ))

        # reaction_genes — flatten this reaction's GPR into (reaction, gene)
        # links. Boolean AND/OR logic stays in reactions.gene_rule; this table
        # is the fast many-to-many join.
        for gene in rxn.genes:
            gene_db_id = existing_genes.get(gene.id)
            if gene_db_id and (rxn_db_id, gene_db_id) not in existing_rxn_gene_pairs:
                db.add(ReactionGene(reaction_id=rxn_db_id, gene_id=gene_db_id))
                existing_rxn_gene_pairs.add((rxn_db_id, gene_db_id))

        flux = float(solution.fluxes.get(rxn.id, 0.0))
        db.add(EntryReaction(
            entry_id=entry_id, reaction_id=rxn_db_id, flux=flux,
            in_degree=G.in_degree(rxn.id) if rxn.id in G else 0,
            out_degree=G.out_degree(rxn.id) if rxn.id in G else 0,
            is_active=abs(flux) > 1e-9,
        ))

    # ── 2. Metabolites — bulk pre-fetch to avoid N+1 ─────────────────────────
    model_met_ids   = [m.id for m in model.metabolites]
    existing_mets   = {
        m.bigg_id: m.id
        for m in db.query(Metabolite).filter(Metabolite.bigg_id.in_(model_met_ids)).all()
    }
    existing_model_mets = {
        m.metabolite_id
        for m in db.query(ModelMetabolite.metabolite_id)
                    .filter_by(model_id=model_db_id).all()
    }

    met_flow:      dict[str, float] = {}
    met_rxns:      dict[str, set[str]] = {}

    for u, v, data in G.edges(data=True):
        for mid in data.get("metabolites", []):
            met_flow[mid] = met_flow.get(mid, 0.0) + data.get("mass_flow", 0.0)
            
            if mid not in met_rxns:
                met_rxns[mid] = set()
            met_rxns[mid].update([u, v])  

    for met in model.metabolites:
        met_db_id = existing_mets.get(met.id)
        if not met_db_id:
            # Safe KEGG parsing — annotation values can be list or str
            kegg_annot = (met.annotation or {}).get("kegg.compound")
            kegg = kegg_annot[0] if isinstance(kegg_annot, list) else kegg_annot

            new_met = Metabolite(
                bigg_id=met.id,
                name=met.name,
                formula=met.formula,
                compartment=met.compartment,   
                charge=met.charge,             
                kegg_id=kegg,
            )
            db.add(new_met)
            db.flush()
            met_db_id = existing_mets[met.id] = new_met.id

        if met_db_id not in existing_model_mets:
            db.add(ModelMetabolite(model_id=model_db_id, metabolite_id=met_db_id))

        flow = met_flow.get(met.id, 0.0)
        if flow > 0:
            db.add(EntryMetabolite(
                entry_id=entry_id,
                metabolite_id=met_db_id,
                total_flow=flow,
                n_reactions=len(met_rxns.get(met.id, set())),
            ))

    # ── 3. reaction_metabolites — GEM stoichiometry (structural, flux-free) ───
    # Runs after both loops so reaction and metabolite DB ids are all resolved.
    # cobrapy rxn.metabolites maps Metabolite -> coefficient (neg=substrate,
    # pos=product). Model-specific: each model stores its own stoichiometry,
    # since direction/magnitude for a shared reaction bigg_id varies across GEMs.
    for rxn in model.reactions:
        rxn_db_id = existing_rxns.get(rxn.id)
        for met, coeff in rxn.metabolites.items():
            met_db_id = existing_mets.get(met.id)
            if met_db_id and (rxn_db_id, met_db_id) not in existing_rxn_met_pairs:
                db.add(ReactionMetabolite(
                    model_id=model_db_id,
                    reaction_id=rxn_db_id, metabolite_id=met_db_id,
                    stoichiometry=float(coeff),
                ))
                existing_rxn_met_pairs.add((rxn_db_id, met_db_id))

    db.commit()


 
# 6. Top-level orchestrator
 

def ingest_gem(
    model_path:    str,
    db:            Session,
    condition_id:  str,
    output_root:   str = "graphs",
    solver_key:    str = "glpk",
    output_subdir: str | None = None,
) -> dict:
    """
    Full ingest pipeline for a single GEM under a canonical DB condition.
    Call only for default-FBA entries qualifying for repo storage.
    Ephemeral tool-layer graphs use build_mfg / run_fba directly
    via compute.py and must never call ingest_gem().
    """

    # ── 1. File hash + model DB record ───────────────────────────────────────
    file_hash = hashlib.sha256(Path(model_path).read_bytes()).hexdigest()
    model_row = db.query(DBModel).filter_by(file_hash=file_hash).first()
    if not model_row:
        raise ValueError(
            f"No model record in DB for file_hash={file_hash[:12]}... "
            f"Register the model via populate_db.py --register before ingestion."
        )

    # ── 2. Idempotency guard — before any expensive work ─────────────────────
    entry_id = make_entry_id(file_hash, condition_id, solver_key)
    if db.query(MFGEntry).filter_by(entry_id=entry_id).first():
        logger.info("Entry %s already in DB — skipping.", entry_id)
        return {
            "entry_id":     entry_id,
            "model_id":     model_row.id,
            "condition_id": condition_id,
            "skipped":      True,
        }

    # ── 3. Condition resolution ───────────────────────────────────────────────
    cond_row = db.query(Condition).filter_by(condition_id=condition_id).first()
    if not cond_row:
        raise ValueError(
            f"Condition '{condition_id}' not found in DB. "
            f"Seed it in the conditions table via populate_db.py before ingestion."
        )
    condition = {
        "exchange_bounds": json.loads(cond_row.exchange_bounds_json or "{}"),
        "carbon_source":   cond_row.carbon_source,
    }

    # ── 4. Model load + FBA ──────────────────────────────────────────────────
    cobra_model = load_model_with_condition(model_path, condition)
    solution    = run_fba(cobra_model)
    if solution is None:
        raise ValueError(
            f"FBA infeasible for model_id={model_row.id}, condition='{condition_id}'. "
            f"Check exchange_bounds_json in the conditions table."
        )

    # ── 5. Graph construction ────────────────────────────────────────────────
    G        = build_mfg(cobra_model, solution)

    # Graph stats for the return payload / logs only (never persisted — the
    # MFGEntry row computes its own counts). Derived directly from G:
    n_active_reactions = sum(
        1 for nd in G.nodes if abs(G.nodes[nd].get("flux", 0.0)) > 1e-9
    )
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()

    # ── 6. File I/O ──────────────────────────────────────────────────────────
    # output_subdir (e.g. the organism slug from metadata) routes files to a
    # readable folder; falls back to the model DB id when not supplied.
    out_dir = Path(output_root) / (output_subdir or str(model_row.id))
    out_dir.mkdir(parents=True, exist_ok=True)

    graphml_path      = str(out_dir / f"{entry_id}.graphml")
    fba_solution_path = str(out_dir / f"{entry_id}_fba.json")

    write_graphml(G, graphml_path)
    Path(fba_solution_path).write_text(
        json.dumps({
            "objective_value": float(solution.objective_value),  # key preserved
            "status":          solution.status,                  # preserved for validation
            "fluxes":          solution.fluxes.to_dict(),
        }, indent=2)
    )
    logger.info("Files written → %s", out_dir)

    # ── 7. DB population ─────────────────────────────────────────────────────
    index_gem_and_mfg(
        model=cobra_model,
        solution=solution,
        G=G,
        model_db_id=model_row.id,
        entry_id=entry_id,
        condition_id_fk=cond_row.id,
        db=db,
        graphml_path=graphml_path,
        fba_solution_path=fba_solution_path,
    )

    logger.info(
        "Ingested entry_id=%s | active_nodes=%d | edges=%d | obj=%.6f",
        entry_id, n_active_reactions, n_edges,
        solution.objective_value,
    )

    return {
        "entry_id":            entry_id,
        "model_id":            model_row.id,
        "condition_id":        condition_id,
        "graphml_path":        graphml_path,
        "fba_solution_path":   fba_solution_path,
        "n_nodes":             n_nodes,
        "n_edges":             n_edges,
        "n_active_reactions":  n_active_reactions,
        "fba_objective_value": float(solution.objective_value),
        "skipped":             False,
    }