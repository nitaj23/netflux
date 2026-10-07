"""Pydantic request/response models for the NetFlux API.

Boundary types only: they validate and document what crosses the HTTP edge, then
hand plain dicts to the tool layer, which owns all business logic. ``SourceSpec``
/ ``OverrideParams`` / ``GudieConfigSpec`` serialise to exactly the dicts the
tool functions already accept; ``DesignParams`` forbids unknown keys so a typo is
a clean 422 at the edge, while emitting only the keys the client set so every
omitted field keeps its engine default.
"""
from __future__ import annotations

from typing import Any, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator


# ── Model resolution + overrides ─────────────────────────────────────────────
class SourceSpec(BaseModel):
    """One of the three tool-layer resolution strategies; exactly one id is
    required, matching ``strategy``."""
    model_config = ConfigDict(extra="forbid")

    strategy: Literal["repo", "custom", "session"]
    model_id: Optional[int] = None
    gem_token: Optional[str] = None
    session_token: Optional[str] = None
    overrides: Optional[dict[str, Any]] = None  # session-replay override layer

    _REQUIRED = {"repo": "model_id", "custom": "gem_token", "session": "session_token"}

    @model_validator(mode="after")
    def _require_matching_id(self) -> "SourceSpec":
        field = self._REQUIRED[self.strategy]
        if getattr(self, field) is None:
            raise ValueError(f"strategy '{self.strategy}' requires '{field}'")
        return self

    def to_dict(self) -> dict:
        out: dict[str, Any] = {"strategy": self.strategy}
        for key in ("model_id", "gem_token", "session_token", "overrides"):
            if (value := getattr(self, key)) is not None:
                out[key] = value
        return out


class CarbonSource(BaseModel):
    model_config = ConfigDict(extra="forbid")
    exchange_id: str
    uptake: float  # applied as lb = -abs(uptake)


class OverrideParams(BaseModel):
    """The four override axes (applied in order) plus the pFBA toggle; a knockout
    is ``reaction_bounds`` with ``[0, 0]``. Every axis is optional."""
    model_config = ConfigDict(extra="forbid")

    condition_id: Optional[str] = None
    custom_exchange_bounds: Optional[dict[str, tuple[float, float]]] = None
    reaction_bounds: Optional[dict[str, tuple[float, float]]] = None
    carbon_source: Optional[CarbonSource] = None
    objective_rxn: Optional[str] = None
    direction: Optional[Literal["min", "max"]] = None
    use_pfba: bool = False

    def to_params(self) -> dict:
        return self.model_dump(exclude_none=True)


class GudieConfigSpec(BaseModel):
    """Overrides for GUDIE's validated-default hyper-parameters; only the set
    fields are applied, the rest keep their swept defaults."""
    model_config = ConfigDict(extra="forbid")

    alpha_r: Optional[float] = None
    hops: Optional[int] = None
    k: Optional[float] = None
    lam: Optional[float] = None
    use_fca: Optional[bool] = None

    def to_config(self):
        fields = self.model_dump(exclude_none=True)
        if not fields:
            return None
        from src.extraction.gudie import GudieConfig
        return GudieConfig(**fields)


# ── Compute / extraction / design requests ───────────────────────────────────
class MFGRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    source: SourceSpec
    params: OverrideParams = Field(default_factory=OverrideParams)
    output_format: Literal["graphml", "pyg", "fba_json", "stats_json"] = "stats_json"
    keep_in_session: bool = False


class PathwayRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    source: SourceSpec
    target: str  # product reaction id, or a metabolite id to drain
    source_ids: Optional[list[str]] = None        # GUDIE dual-anchor substrates
    config: Optional[GudieConfigSpec] = None
    repair: bool = False                          # opt-in connectivity repair (extractor)
    fva_operating_points: Optional[list[float]] = None  # fractions of v*; gates FVA annotation
    fva_loopless: bool = False
    include_graphml: bool = False


class DesignParams(BaseModel):
    """Tuning knobs for the unified design engine (mirrors the tool-layer
    allow-list); unknown keys are rejected at the edge."""
    model_config = ConfigDict(extra="forbid")

    max_knockouts: Optional[int] = None
    enable_k2: Optional[bool] = None
    min_growth_frac: Optional[float] = None
    min_growth_abs: Optional[float] = None
    gamma: Optional[float] = None
    eps_prod: Optional[float] = None
    class_tol: Optional[float] = None
    backend: Optional[str] = None
    top_k: Optional[int] = None
    solver_timeout: Optional[float] = None
    gpr_cap: Optional[int] = None

    def to_params(self) -> dict:
        return self.model_dump(exclude_unset=True)


class DesignRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    source: SourceSpec
    target: str     # product reaction id whose overproduction is coupled to growth
    biomass: str    # growth/biomass reaction id
    substrate: Optional[str] = None            # auto-detected when omitted
    methods: list[Literal["optforce", "optknock"]] = Field(
        default_factory=lambda: ["optforce", "optknock"])
    require_buildable: bool = True
    parallel: bool = False
    params: Optional[DesignParams] = None
    discovery_source: Optional[SourceSpec] = None


# ── Catalogue / misc response models ─────────────────────────────────────────
class ModelOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)
    id: int
    organism: str
    ncbi_taxon_id: Optional[int] = None
    source_db: Optional[str] = None
    model_path: str
    file_hash: str


class ConditionOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)
    id: int
    condition_id: str
    carbon_source: Optional[str] = None
    oxygen_status: Optional[str] = None
    exchange_bounds_json: Optional[str] = None


class GemRegistrationResponse(BaseModel):
    gem_token: str
    summary: dict[str, Any]
    warnings: list[str]


class HealthResponse(BaseModel):
    status: Literal["ok", "degraded"]
    database: bool


# ── Compound search ───────────────────────────────────────────────────────────
class CompoundOut(BaseModel):
    """One base compound (compartment variants folded together)."""
    base_id: str
    name: Optional[str] = None
    formula: Optional[str] = None
    kegg_id: Optional[str] = None
    compartments: list[str] = Field(default_factory=list)
    suggested_target: Optional[str] = None  # set for custom-model listings


class CompoundOrganismOut(BaseModel):
    """A model carrying a compound, with its production readiness signals."""
    model: ModelOut
    has_exchange: bool
    carries_flow: bool
    suggested_target: Optional[str] = None


# ── Model defaults ──────────────────────────────────────────────────────────
class ModelDefaultsOut(BaseModel):
    objective_rxn: Optional[str] = None
    objective_rxns: list[str] = Field(default_factory=list)
    direction: str
    carbon_source: Optional[str] = None
    carbon_source_options: list[str] = Field(default_factory=list)
    exchange_bounds: dict[str, list[float]] = Field(default_factory=dict)
    condition_id: Optional[str] = None
    condition_options: list[str] = Field(default_factory=list)
    target: Optional[str] = None


# ── Async jobs ────────────────────────────────────────────────────────────────
class JobSubmitResponse(BaseModel):
    job_id: str
    status: Literal["queued", "running", "completed", "failed"]


class JobStatusResponse(BaseModel):
    job_id: str
    status: Literal["queued", "running", "completed", "failed"]
    result: Optional[dict[str, Any]] = None
    error: Optional[dict[str, Any]] = None
