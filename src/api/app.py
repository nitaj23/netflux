"""
NetFlux API - a thin FastAPI edge over the existing backend.

Every endpoint validates its request (Pydantic), injects a request-scoped DB
session, delegates to the tool/query layer (no business logic is re-implemented
here), serialises the in-memory result, and lets one exception handler map the
backend's typed errors onto HTTP status codes.

Run with ``uvicorn src.api.app:app --reload`` or ``python -m src.api``; docs at
``/docs``.
"""
from __future__ import annotations

import base64
import logging
from contextlib import asynccontextmanager
from typing import Iterator

from fastapi import APIRouter, Depends, FastAPI, HTTPException, Path, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from sqlalchemy import text
from sqlalchemy.orm import Session

from src import queries, tool
from src.db import SessionLocal
from src.schema import Condition, Model
from src.api import schemas as S
from src.extraction.pathway import (
    AnnotationBudgetExceeded, InfeasibleModel, InvalidOperatingPoint,
    TargetBlocked, TargetNoMassFlow, TargetNotFound, TargetNotProducible,
)

logger = logging.getLogger(__name__)
API_PREFIX = "/api/v1"


# ── Dependency: request-scoped DB session ────────────────────────────────────
def get_db() -> Iterator[Session]:
    """Yield a DB session and always close it. Tests override this provider."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# ── Error mapping: typed domain errors -> one JSON envelope ──────────────────
# Keyed by exception class NAME so the design-engine errors (whose packages pull
# in heavy solver deps, imported lazily by the tool layer) need no import here.
# A domain error not in the table but carrying a ``.hint`` is a client-actionable
# 400; anything else is an opaque, logged 500. The cheap tool/pathway error bases
# are registered directly in create_app so they convert without re-raising; the
# design-engine errors ride the ``Exception`` fallback (correct status in prod).
_STATUS: dict[str, int] = {
    # custom-GEM upload
    "GemUnparseable": 422, "GemEmpty": 422, "GemNoGenes": 422, "GemTooLarge": 413,
    # compute / feasibility
    "Infeasible": 422, "InfeasibleMaintenance": 422,
    # session cache
    "TokenNotFound": 404, "TokenExpired": 410, "CacheFull": 429,
    "AnnotationUnavailable": 409,
    # pathway extraction
    "TargetNotFound": 404, "TargetBlocked": 422, "TargetNotProducible": 422,
    "TargetNoMassFlow": 422, "InfeasibleModel": 422, "InvalidOperatingPoint": 422,
    "AnnotationBudgetExceeded": 422,
    # design-engine base-model contract failures
    "BaseModelInfeasible": 422, "BiomassNotGrowing": 422, "InvalidGamma": 422,
    "DiscoveryModelInvalid": 422, "MissingFullModelReaction": 422,
    "AmbiguousTargetIdentity": 422, "SolverFailure": 500, "SolverTimeout": 504,
}


async def _handle_domain_error(request: Request, exc: Exception) -> JSONResponse:
    name = type(exc).__name__
    status = _STATUS.get(name)
    if status is None:
        if hasattr(exc, "hint"):          # an unlisted but typed domain error
            status = 400
        else:                             # a real bug: don't leak internals
            logger.exception("unhandled error serving %s", request.url.path)
            return JSONResponse(status_code=500, content={
                "error": "InternalServerError",
                "detail": "an internal error occurred", "hint": None})
    return JSONResponse(status_code=status, content={
        "error": name, "detail": str(exc) or name, "hint": getattr(exc, "hint", None)})


# ── Serialisation of in-memory artifacts ─────────────────────────────────────
def _mfg_payload(out: dict) -> dict:
    """Base64 the binary GraphML/PyG blobs; everything else passes through."""
    result = dict(out)
    for key in ("graphml", "pyg"):
        if isinstance(result.get(key), (bytes, bytearray)):
            result[key] = base64.b64encode(result[key]).decode("ascii")
            result["content_encoding"] = "base64"
    return result


def _pathway_payload(result, *, include_graphml: bool) -> dict:
    g = result.gudie
    ann = result.annotation
    payload = {
        "target_id": result.target_id,
        "reactions": result.reactions,
        "n_pathway": len(result.reactions),
        "n_consistent": result.n_consistent,
        "n_active": result.n_active,
        "gudie": {
            "ppr_converged": g.ppr_converged,
            "ludie_fallback": g.ludie_fallback,
            "threshold": g.threshold,
            "coupling_available": g.coupling is not None,
        },
        "subgraph": result.to_json(),
        "annotation": None if ann is None else {
            "operating_points": ann.operating_points,
            "v_star_target": ann.v_star_target,
            "coupling_tau": ann.coupling_tau,
            "per_reaction": {rid: {"reaction_id": ra.reaction_id, "gene_id": ra.gene_id,
                                   "coupling": ra.coupling, "per_op": ra.per_op}
                             for rid, ra in ann.per_reaction.items()},
        },
    }
    if include_graphml:
        payload["graphml_base64"] = base64.b64encode(result.to_graphml()).decode("ascii")
    return payload


# ── Routes ───────────────────────────────────────────────────────────────────
router = APIRouter(prefix=API_PREFIX)


@router.get("/health", response_model=S.HealthResponse, tags=["health"])
def health(db: Session = Depends(get_db)) -> S.HealthResponse:
    """Service status + a trivial DB probe (never raises)."""
    try:
        db.execute(text("SELECT 1"))
        ok = True
    except Exception:
        ok = False
    return S.HealthResponse(status="ok" if ok else "degraded", database=ok)


@router.get("/models", response_model=list[S.ModelOut], tags=["catalog"])
def list_models(organism: str | None = Query(None), limit: int = Query(100, ge=1, le=1000),
                offset: int = Query(0, ge=0), db: Session = Depends(get_db)):
    """List repository GEMs, optionally filtered by organism substring."""
    if organism:
        return queries.gems_by_organism(db, organism)[offset: offset + limit]
    return db.query(Model).order_by(Model.id).offset(offset).limit(limit).all()


@router.get("/models/{model_id}", response_model=S.ModelOut, tags=["catalog"])
def get_model(model_id: int = Path(..., ge=1), db: Session = Depends(get_db)):
    if (row := db.get(Model, model_id)) is None:
        raise HTTPException(404, f"no model with id={model_id}")
    return row


@router.get("/conditions", response_model=list[S.ConditionOut], tags=["catalog"])
def list_conditions(db: Session = Depends(get_db)):
    """Medium/condition presets usable as an override ``condition_id``."""
    return db.query(Condition).order_by(Condition.condition_id).all()


@router.get("/conditions/{condition_id}", response_model=S.ConditionOut, tags=["catalog"])
def get_condition(condition_id: str, db: Session = Depends(get_db)):
    row = db.query(Condition).filter(Condition.condition_id == condition_id).one_or_none()
    if row is None:
        raise HTTPException(404, f"no condition_id={condition_id!r}")
    return row


@router.get("/reactions/{bigg_id}/gems", response_model=list[S.ModelOut], tags=["catalog"])
def gems_with_reaction(bigg_id: str, db: Session = Depends(get_db)):
    """Repository GEMs whose network contains the reaction (by BiGG id)."""
    return queries.gems_containing_reaction(db, bigg_id)


@router.get("/metabolites/{bigg_id}/gems", response_model=list[S.ModelOut], tags=["catalog"])
def gems_with_metabolite(bigg_id: str, db: Session = Depends(get_db)):
    """Repository GEMs whose network contains the metabolite (by BiGG id)."""
    return queries.gems_containing_metabolite(db, bigg_id)


@router.post("/gems", response_model=S.GemRegistrationResponse, status_code=201, tags=["gems"])
async def register_gem(request: Request) -> S.GemRegistrationResponse:
    """Register a custom GEM (raw cobra-JSON body) -> ephemeral ``gem_token``."""
    token, summary, warnings = tool.register_custom_gem(await request.body())
    return S.GemRegistrationResponse(gem_token=token, summary=summary, warnings=warnings)


@router.post("/mfg", tags=["compute"])
def compute_mfg(req: S.MFGRequest) -> dict:
    """Resolve + override a model, run FBA/pFBA, return the MFG artifact."""
    return _mfg_payload(tool.compute_mfg(
        req.source.to_dict(), req.params.to_params(),
        output_format=req.output_format, keep_in_session=req.keep_in_session))


@router.post("/pathways", tags=["pathway"])
def extract_pathway(req: S.PathwayRequest, db: Session = Depends(get_db)) -> dict:
    """Extract the bioproduction pathway H for ``target`` (+ optional FVA annotation)."""
    result = tool.extract_pathway(
        req.source.to_dict(), req.target, source_ids=req.source_ids,
        config=req.config.to_config() if req.config else None,
        fva_operating_points=req.fva_operating_points, fva_loopless=req.fva_loopless, db=db)
    return _pathway_payload(result, include_graphml=req.include_graphml)


@router.post("/designs", tags=["design"])
def rank_designs(req: S.DesignRequest, db: Session = Depends(get_db)) -> dict:
    """Rank buildable KO designs (unified OptForce + OptKnock engine); regulatory
    leads come back separately as ``hypotheses``."""
    return tool.rank_design_interventions(
        req.source.to_dict(), req.target, req.biomass, substrate=req.substrate,
        methods=tuple(req.methods), require_buildable=req.require_buildable,
        parallel=req.parallel, params=req.params.to_params() if req.params else None,
        discovery_source=req.discovery_source.to_dict() if req.discovery_source else None,
        db=db)


# ── App factory ───────────────────────────────────────────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    try:
        from src.db import init_db
        init_db()  # idempotent; creates only missing tables
    except Exception:
        logger.warning("init_db() failed at startup", exc_info=True)
    yield


def create_app() -> FastAPI:
    app = FastAPI(title="NetFlux API", version="1.0.0", lifespan=lifespan,
                  description="Ephemeral FBA/MFG compute, bioproduction-pathway "
                              "extraction and strain-design ranking over repository "
                              "and user-uploaded GEMs. All compute is RAM-only.")
    app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"],
                       allow_headers=["*"])
    # Register the cheap tool/pathway error bases so they convert cleanly; the
    # Exception fallback catches design-engine errors and unexpected bugs.
    for base in (tool.ToolError, TargetNotFound, TargetBlocked, TargetNotProducible,
                 TargetNoMassFlow, InfeasibleModel, InvalidOperatingPoint,
                 AnnotationBudgetExceeded, Exception):
        app.add_exception_handler(base, _handle_domain_error)
    app.include_router(router)

    @app.get("/", include_in_schema=False)
    def root() -> dict:
        return {"service": "NetFlux API", "version": app.version, "docs": "/docs"}

    return app


app = create_app()
