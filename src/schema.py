from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column
from sqlalchemy import ForeignKey, Integer, Text, Float, Boolean, Date, Index, UniqueConstraint
from datetime import date

class Base(DeclarativeBase): pass

class Model(Base): 
    __tablename__ = "models"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    organism: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    ncbi_taxon_id: Mapped[int|None] = mapped_column(Integer, index=True)
    source_db: Mapped[str|None] = mapped_column(Text)
    model_path: Mapped[str] = mapped_column(Text, nullable=False)
    file_hash: Mapped[str] = mapped_column(Text, nullable=False)

class Reaction(Base):
    __tablename__ = "reactions"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    bigg_id: Mapped[str] = mapped_column(Text, unique=True, nullable=False)
    name: Mapped[str|None] = mapped_column(Text)
    reaction_formula: Mapped[str|None] = mapped_column(Text)
    subsystem: Mapped[str|None] = mapped_column(Text, index=True)
    gene_rule: Mapped[str|None] = mapped_column(Text)

class Metabolite(Base):
    __tablename__ = "metabolites"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    bigg_id: Mapped[str] = mapped_column(Text, unique=True, nullable=False)
    name: Mapped[str|None] = mapped_column(Text)
    formula: Mapped[str|None] = mapped_column(Text)
    compartment: Mapped[str|None] = mapped_column(Text, index=True)
    charge: Mapped[int|None] = mapped_column(Integer)
    kegg_id: Mapped[str|None] = mapped_column(Text)

class ModelReaction(Base):
    __tablename__ = "model_reactions"
    __table_args__ = (UniqueConstraint("model_id", "reaction_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("models.id"), index=True)
    reaction_id: Mapped[int] = mapped_column(ForeignKey("reactions.id"), index=True)
    lower_bound: Mapped[float|None] = mapped_column(Float)
    upper_bound: Mapped[float|None] = mapped_column(Float)

class ModelMetabolite(Base):
    __tablename__ = "model_metabolites"
    __table_args__ = (UniqueConstraint("model_id", "metabolite_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("models.id"), index=True)
    metabolite_id: Mapped[int] = mapped_column(ForeignKey("metabolites.id"), index=True)

class Gene(Base):
    __tablename__ = "genes"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    bigg_id: Mapped[str] = mapped_column(Text, unique=True, nullable=False)
    name: Mapped[str|None] = mapped_column(Text)
    ncbigene_id: Mapped[str|None] = mapped_column(Text)
    uniprot_id: Mapped[str|None] = mapped_column(Text)

class ModelGene(Base):
    __tablename__ = "model_genes"
    __table_args__ = (UniqueConstraint("model_id", "gene_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("models.id"), index=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id"), index=True)

class ReactionGene(Base):
    __tablename__ = "reaction_genes"
    __table_args__ = (UniqueConstraint("reaction_id", "gene_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    reaction_id: Mapped[int] = mapped_column(ForeignKey("reactions.id"), index=True)
    gene_id: Mapped[int] = mapped_column(ForeignKey("genes.id"), index=True)

class ReactionMetabolite(Base):
    # Stoichiometry is model-specific: the same reaction bigg_id is written with
    # different directionality (sign) and occasionally different magnitude across
    # GEMs, so this is keyed by (model_id, reaction_id, metabolite_id), NOT global.
    __tablename__ = "reaction_metabolites"
    __table_args__ = (UniqueConstraint("model_id", "reaction_id", "metabolite_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("models.id"), index=True)
    reaction_id: Mapped[int] = mapped_column(ForeignKey("reactions.id"), index=True)
    metabolite_id: Mapped[int] = mapped_column(ForeignKey("metabolites.id"), index=True)
    # cobrapy convention: negative coefficient = substrate, positive = product
    stoichiometry: Mapped[float] = mapped_column(Float)

class Condition(Base):
    __tablename__ = "conditions"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    condition_id: Mapped[str] = mapped_column(Text, unique=True, nullable=False)
    carbon_source: Mapped[str|None] = mapped_column(Text)
    oxygen_status: Mapped[str|None] = mapped_column(Text)
    exchange_bounds_json: Mapped[str|None] = mapped_column(Text)
    
class MFGEntry(Base): 
    __tablename__ = "mfg_entries"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    entry_id: Mapped[str] = mapped_column(Text, unique=True, nullable=False)
    model_id: Mapped[int] = mapped_column(ForeignKey("models.id"), index=True)
    condition_id_fk: Mapped[int|None] = mapped_column(ForeignKey("conditions.id"))
    fba_objective_value: Mapped[float|None] = mapped_column(Float)
    n_reactions: Mapped[int|None] = mapped_column(Integer)
    n_edges: Mapped[int|None] = mapped_column(Integer)
    graphml_path: Mapped[str|None] = mapped_column(Text)
    pyg_path: Mapped[str|None] = mapped_column(Text)
    fba_solution_path: Mapped[str|None] = mapped_column(Text)
    is_latest: Mapped[bool] = mapped_column(Boolean, default=True, index=True)
    
class EntryReaction(Base): 
    __tablename__ = "entry_reactions"
    __table_args__ = (UniqueConstraint("entry_id", "reaction_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    entry_id: Mapped[str] = mapped_column(ForeignKey("mfg_entries.entry_id"), index=True)
    reaction_id: Mapped[int] = mapped_column(ForeignKey("reactions.id"), index=True)
    flux: Mapped[float|None] = mapped_column(Float)
    in_degree: Mapped[int|None] = mapped_column(Integer)
    out_degree: Mapped[int|None] = mapped_column(Integer)
    is_active: Mapped[bool] = mapped_column(Boolean, default=False, index=True)

class EntryMetabolite(Base): 
    __tablename__ = "entry_metabolites"
    __table_args__ = (UniqueConstraint("entry_id", "metabolite_id"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    entry_id: Mapped[str] = mapped_column(ForeignKey("mfg_entries.entry_id"), index=True)
    metabolite_id: Mapped[int] = mapped_column(ForeignKey("metabolites.id"), index=True)
    total_flow: Mapped[float|None] = mapped_column(Float)
    n_reactions: Mapped[int|None] = mapped_column(Integer)

class GATModel(Base):
    __tablename__ = "gat_models"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    model_slug: Mapped[str] = mapped_column(Text, unique=True)
    architecture_json: Mapped[str] = mapped_column(Text)
    weights_path: Mapped[str] = mapped_column(Text)
    config_path: Mapped[str] = mapped_column(Text)
    node_feature_schema: Mapped[str] = mapped_column(Text)
    val_metrics_json: Mapped[str|None] = mapped_column(Text)
  
