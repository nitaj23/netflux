from sqlalchemy.orm import Session
from src.schema import Model, Reaction, Metabolite, ModelReaction, ModelMetabolite, MFGEntry, EntryReaction, EntryMetabolite

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
