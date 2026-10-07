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
import os
from contextlib import asynccontextmanager
from typing import Iterator

from fastapi import (APIRouter, Depends, FastAPI, File, HTTPException, Path,
                     Query, Request, UploadFile)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from sqlalchemy import text
from sqlalchemy.orm import Session

from src import queries, tool
from src.db import SessionLocal
from src.schema import Condition, Model
from src.api import schemas as S
from src.api.jobs import JobManager, MAX_CONCURRENT_DESIGN_JOBS
from src.extraction.pathway import (
    AnnotationBudgetExceeded, InfeasibleModel, InvalidOperatingPoint,
    TargetBlocked, TargetNoMassFlow, TargetNotFound, TargetNotProducible,
)

logger = logging.getLogger(__name__)
API_PREFIX = "/api/v1"

# Default allowed origins for local frontend dev. Production overrides via the
# NETFLUX_CORS_ORIGINS env var (comma-separated). A permissive "*" is never the
# built-in default; it is only used if an operator explicitly sets the var to "*".
_DEFAULT_CORS_ORIGINS = ("http://localhost:3000", "http://127.0.0.1:3000",
                         "http://localhost:5173", "http://127.0.0.1:5173")


def allowed_origins() -> list[str]:
    """Explicit, configurable CORS allow-list (env NETFLUX_CORS_ORIGINS or dev defaults)."""
    raw = os.getenv("NETFLUX_CORS_ORIGINS")
    if raw is None:
        return list(_DEFAULT_CORS_ORIGINS)
    return [o.strip() for o in raw.split(",") if o.strip()]


# One process-wide bounded runner for async design jobs (RAM-only, ephemeral).
job_manager = JobManager(max_workers=MAX_CONCURRENT_DESIGN_JOBS,
                         error_serializer=lambda exc: {
                             "error": type(exc).__name__,
                             "detail": str(exc) or type(exc).__name__,
                             "hint": getattr(exc, "hint", None)})


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


@router.get("/compounds", response_model=list[S.CompoundOut], tags=["catalog"])
def search_compounds(q: str = Query(..., min_length=1), limit: int = Query(25, ge=1, le=500),
                     db: Session = Depends(get_db)):
    """Search compounds by id/name, grouped by base compound (compartments folded)."""
    return queries.search_compounds(db, q, limit)


@router.get("/compounds/{base_id}/organisms", response_model=list[S.CompoundOrganismOut],
            tags=["catalog"])
def compound_organisms(base_id: str, db: Session = Depends(get_db)):
    """Models carrying `base_id` with has_exchange / carries_flow / suggested_target;
    flow-carrying models first."""
    return queries.compound_organisms(db, base_id)


@router.get("/models/{model_id}/defaults", response_model=S.ModelDefaultsOut, tags=["catalog"])
def model_defaults(model_id: int = Path(..., ge=1), compound: str | None = Query(None),
                   db: Session = Depends(get_db)):
    """Default operating point (objective, carbon source, exchange bounds, conditions,
    suggested target) the frontend needs to pre-fill a run form."""
    model = tool.resolve_model({"strategy": "repo", "model_id": model_id}, db=db)
    conditions = db.query(Condition).order_by(Condition.condition_id).all()
    return tool.model_defaults(model, conditions=conditions, compound=compound)


@router.get("/gems/{gem_token}/compounds", response_model=list[S.CompoundOut], tags=["gems"])
def gem_compounds(gem_token: str, q: str | None = Query(None),
                  limit: int = Query(25, ge=1, le=500)):
    """Compounds in an uploaded GEM, grouped by base compound with a suggested
    target -- the Route-1 compound picker for custom models (optional `q` filter)."""
    model = tool.resolve_model({"strategy": "custom", "gem_token": gem_token})
    return tool.model_compounds(model, q=q, limit=limit)


@router.get("/gems/{gem_token}/defaults", response_model=S.ModelDefaultsOut, tags=["gems"])
def gem_defaults(gem_token: str, compound: str | None = Query(None),
                 db: Session = Depends(get_db)):
    """Default operating point for an uploaded GEM (parity with the repo
    /models/{id}/defaults): objective, carbon source, exchange bounds, conditions,
    and the suggested target for `compound`."""
    model = tool.resolve_model({"strategy": "custom", "gem_token": gem_token})
    conditions = db.query(Condition).order_by(Condition.condition_id).all()
    return tool.model_defaults(model, conditions=conditions, compound=compound)


@router.post("/gems", response_model=S.GemRegistrationResponse, status_code=201, tags=["gems"])
async def register_gem(request: Request) -> S.GemRegistrationResponse:
    """Register a custom GEM (raw COBRA-JSON or SBML/XML body) -> ephemeral ``gem_token``."""
    token, summary, warnings = tool.register_custom_gem(await request.body())
    return S.GemRegistrationResponse(gem_token=token, summary=summary, warnings=warnings)


@router.post("/gems/upload", response_model=S.GemRegistrationResponse, status_code=201,
             tags=["gems"])
async def upload_gem(file: UploadFile = File(...)) -> S.GemRegistrationResponse:
    """Browser-friendly multipart upload of a custom GEM (COBRA-JSON or SBML/XML).

    Parallels the raw-body ``POST /gems`` (which is preserved); same validation,
    same 50 MB cap, same ephemeral token/summary/warning/error structure. The
    upload is never written to disk or a database."""
    token, summary, warnings = tool.register_custom_gem(await file.read())
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
        config=req.config.to_config() if req.config else None, repair=req.repair,
        fva_operating_points=req.fva_operating_points, fva_loopless=req.fva_loopless, db=db)
    return _pathway_payload(result, include_graphml=req.include_graphml)


def _run_design_job(req: S.DesignRequest) -> dict:
    """Execute a design ranking with its OWN DB session (runs off-request-thread)."""
    db = SessionLocal()
    try:
        return tool.rank_design_interventions(
            req.source.to_dict(), req.target, req.biomass, substrate=req.substrate,
            methods=tuple(req.methods), require_buildable=req.require_buildable,
            parallel=req.parallel, params=req.params.to_params() if req.params else None,
            discovery_source=req.discovery_source.to_dict() if req.discovery_source else None,
            db=db)
    finally:
        db.close()


@router.post("/designs", response_model=S.JobSubmitResponse, status_code=202, tags=["design"])
def rank_designs(req: S.DesignRequest) -> S.JobSubmitResponse:
    """Submit a strain-design ranking as an async job; returns a ``job_id``.

    The unified OptForce + OptKnock engine runs on a bounded in-process executor
    (at most 4 concurrent jobs); poll ``GET /jobs/{job_id}`` for the result. A
    completed job's ``result`` is exactly the payload the engine produces."""
    job_id = job_manager.submit(_run_design_job, req)
    job = job_manager.get(job_id)
    return S.JobSubmitResponse(job_id=job_id, status=job.status if job else "queued")


@router.get("/jobs/{job_id}", response_model=S.JobStatusResponse, tags=["design"])
def get_job(job_id: str) -> S.JobStatusResponse:
    """Poll an async job: queued | running | completed | failed.

    A completed job carries the full design ``result``; a failed job retains the
    typed ``error`` / ``detail`` / ``hint`` from the engine."""
    job = job_manager.get(job_id)
    if job is None:
        raise HTTPException(404, f"no job with id={job_id!r}")
    return S.JobStatusResponse(job_id=job.id, status=job.status,
                               result=job.result if job.status == "completed" else None,
                               error=job.error)


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
    app.add_middleware(CORSMiddleware, allow_origins=allowed_origins(),
                       allow_methods=["*"], allow_headers=["*"])
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
