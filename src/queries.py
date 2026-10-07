from sqlalchemy import or_
from sqlalchemy.orm import Session
from src.schema import Model, Reaction, Metabolite, ModelReaction, ModelMetabolite, MFGEntry, EntryReaction, EntryMetabolite


def _base_compound_id(met_id: str) -> str:
    """Strip a trailing compartment suffix (``_c``/``_e``/``_p``/...) from a met id.

    Shares the single definition in src.tool (imported lazily to keep this module's
    import graph free of the heavy solver stack)."""
    from src.tool import base_compound_id
    return base_compound_id(met_id)


# ── Compound search (base-id grouped) ──────────────────────────────────────

def search_compounds(db: Session, q: str, limit: int = 25) -> list[dict]:
    """Search metabolites by bigg_id or name, grouped by base compound.

    Normalizes compartment suffixes (``_c``/``_e``/``_p``/...) so ``succ_c`` and
    ``succ_e`` fold into one ``succ`` group. Each group carries base_id, name,
    formula, KEGG id and the set of compartments seen. Returns at most `limit`
    groups, ordered by base_id."""
    rows = (db.query(Metabolite)
            .filter(or_(Metabolite.bigg_id.ilike(f"%{q}%"),
                        Metabolite.name.ilike(f"%{q}%")))
            .order_by(Metabolite.bigg_id).all())
    groups: dict[str, dict] = {}
    for m in rows:
        base = _base_compound_id(m.bigg_id)
        g = groups.get(base)
        if g is None:
            g = groups[base] = {
                "base_id": base, "name": m.name, "formula": m.formula,
                "kegg_id": m.kegg_id, "compartments": [],
            }
        # Prefer the first non-empty name/formula/kegg seen across variants.
        g["name"] = g["name"] or m.name
        g["formula"] = g["formula"] or m.formula
        g["kegg_id"] = g["kegg_id"] or m.kegg_id
        if m.compartment and m.compartment not in g["compartments"]:
            g["compartments"].append(m.compartment)
    return list(groups.values())[:limit]


def compound_organisms(db: Session, base_id: str) -> list[dict]:
    """Models carrying a compound (any compartment variant of `base_id`).

    For each model: its catalogue row, ``has_exchange`` (an ``EX_<met>`` reaction
    exists for a variant), ``carries_flow`` (the metabolite appears as edge mass
    flow in the model's latest MFG), and ``suggested_target`` (the preferred
    exchange reaction to overproduce it, extracellular variant first). Models whose
    compound carries flow are returned first."""
    base = _base_compound_id(base_id)
    mets = (db.query(Metabolite)
            .filter(or_(Metabolite.bigg_id == base,
                        Metabolite.bigg_id.like(f"{base}\\_%", escape="\\")))
            .all())
    met_ids = [m.id for m in mets]
    if not met_ids:
        return []
    met_biggs = [m.bigg_id for m in mets]
    ex_biggs = [f"EX_{b}" for b in met_biggs]

    models = (db.query(Model)
              .join(ModelMetabolite, Model.id == ModelMetabolite.model_id)
              .filter(ModelMetabolite.metabolite_id.in_(met_ids))
              .distinct().all())

    flow_model_ids = {row.model_id for row in
                      db.query(MFGEntry.model_id)
                      .join(EntryMetabolite, MFGEntry.entry_id == EntryMetabolite.entry_id)
                      .filter(EntryMetabolite.metabolite_id.in_(met_ids),
                              MFGEntry.is_latest == True).all()}

    ex_by_model: dict[int, list[str]] = {}
    for mid, bigg in (db.query(ModelReaction.model_id, Reaction.bigg_id)
                      .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                      .filter(Reaction.bigg_id.in_(ex_biggs)).all()):
        ex_by_model.setdefault(mid, []).append(bigg)

    out = []
    for model in models:
        exchanges = ex_by_model.get(model.id, [])
        suggested = next((e for e in exchanges if e.endswith("_e")),
                         exchanges[0] if exchanges else None)
        out.append({
            "model": model,
            "has_exchange": bool(exchanges),
            "carries_flow": model.id in flow_model_ids,
            "suggested_target": suggested,
        })
    out.sort(key=lambda r: (not r["carries_flow"], r["model"].organism))
    return out

# ── Reaction queries ───────────────────────────────────────────────────────

def gems_containing_reaction(db: Session, bigg_id: str) -> list[Model]:
    """GEMs that structurally contain this reaction (regardless of flux)."""
    return (db.query(Model)
            .join(ModelReaction, Model.id == ModelReaction.model_id)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.bigg_id == bigg_id).all())

def gems_containing_reaction_by_name(db: Session, name: str) -> list[Model]:
    """GEMs that structurally contain a reaction matching name (case-insensitive substring)."""
    return (db.query(Model)
            .join(ModelReaction, Model.id == ModelReaction.model_id)
            .join(Reaction, Reaction.id == ModelReaction.reaction_id)
            .filter(Reaction.name.ilike(f"%{name}%")).all())

def mfgs_with_active_reaction(db: Session, bigg_id: str) -> list[MFGEntry]:
    """MFGs where this reaction carries non-zero flux under default conditions."""
    return (db.query(MFGEntry)
            .join(EntryReaction, MFGEntry.entry_id == EntryReaction.entry_id)
            .join(Reaction, Reaction.id == EntryReaction.reaction_id)
            .filter(Reaction.bigg_id == bigg_id,
                    EntryReaction.is_active == True,
                    MFGEntry.is_latest == True).all())

def mfgs_with_active_reaction_by_name(db: Session, name: str) -> list[MFGEntry]:
    """MFGs where a reaction matching name carries non-zero flux."""
    return (db.query(MFGEntry)
            .join(EntryReaction, MFGEntry.entry_id == EntryReaction.entry_id)
            .join(Reaction, Reaction.id == EntryReaction.reaction_id)
            .filter(Reaction.name.ilike(f"%{name}%"),
                    EntryReaction.is_active == True,
                    MFGEntry.is_latest == True).all())

# ── Metabolite queries ────────────────────────────────────────────────────

def gems_containing_metabolite(db: Session, bigg_id: str) -> list[Model]:
    """GEMs that structurally contain this metabolite."""
    return (db.query(Model)
            .join(ModelMetabolite, Model.id == ModelMetabolite.model_id)
            .join(Metabolite, Metabolite.id == ModelMetabolite.metabolite_id)
            .filter(Metabolite.bigg_id == bigg_id).all())

def gems_containing_metabolite_by_name(db: Session, name: str) -> list[Model]:
    """GEMs that structurally contain a metabolite matching name (case-insensitive substring)."""
    return (db.query(Model)
            .join(ModelMetabolite, Model.id == ModelMetabolite.model_id)
            .join(Metabolite, Metabolite.id == ModelMetabolite.metabolite_id)
            .filter(Metabolite.name.ilike(f"%{name}%")).all())

def mfgs_carrying_metabolite(db: Session, bigg_id: str) -> list[MFGEntry]:
    """MFGs where this metabolite appears as edge mass flow."""
    return (db.query(MFGEntry)
            .join(EntryMetabolite, MFGEntry.entry_id == EntryMetabolite.entry_id)
            .join(Metabolite, Metabolite.id == EntryMetabolite.metabolite_id)
            .filter(Metabolite.bigg_id == bigg_id,
                    MFGEntry.is_latest == True).all())

def mfgs_carrying_metabolite_by_name(db: Session, name: str) -> list[MFGEntry]:
    """MFGs where a metabolite matching name appears as edge mass flow."""
    return (db.query(MFGEntry)
            .join(EntryMetabolite, MFGEntry.entry_id == EntryMetabolite.entry_id)
            .join(Metabolite, Metabolite.id == EntryMetabolite.metabolite_id)
            .filter(Metabolite.name.ilike(f"%{name}%"),
                    MFGEntry.is_latest == True).all())

# ── Organism queries ──────────────────────────────────────────────────────

def gems_by_organism(db: Session, organism: str) -> list[Model]:
    return db.query(Model).filter(Model.organism.ilike(f"%{organism}%")).all()

def mfgs_by_organism(db: Session, organism: str) -> list[MFGEntry]:
    return (db.query(MFGEntry).join(Model)
            .filter(Model.organism.ilike(f"%{organism}%"),
                    MFGEntry.is_latest == True).all())

# ── Cross-dimension (composable joins) ────────────────────────────────────

def mfgs_by_organism_and_reaction(db: Session, organism: str, bigg_id: str) -> list[MFGEntry]:
    """E. coli MFGs where PFK is active — representative cross query."""
    return (db.query(MFGEntry)
            .join(Model, Model.id == MFGEntry.model_id)
            .join(EntryReaction, MFGEntry.entry_id == EntryReaction.entry_id)
            .join(Reaction, Reaction.id == EntryReaction.reaction_id)
            .filter(Model.organism.ilike(f"%{organism}%"), 
                    Reaction.bigg_id == bigg_id,
                    EntryReaction.is_active == True, 
                    MFGEntry.is_latest == True).all())

def mfgs_by_organism_and_metabolite(db: Session, organism: str, bigg_id: str) -> list[MFGEntry]:
    return (db.query(MFGEntry)
            .join(Model, Model.id == MFGEntry.model_id)
            .join(EntryMetabolite, MFGEntry.entry_id == EntryMetabolite.entry_id)
            .join(Metabolite, Metabolite.id == EntryMetabolite.metabolite_id)
            .filter(Model.organism.ilike(f"%{organism}%"), 
                    Metabolite.bigg_id == bigg_id,
                    MFGEntry.is_latest == True).all()) 
