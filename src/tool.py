# src/tool.py
"""
Ephemeral TOOL-LAYER SERVICE

On-demand FBA + MFG compute and bioproduction-pathway extraction over repo GEMs
*and* user-uploaded custom GEMs, plus an FVA design-annotation layer — with one
hard invariant: **nothing is ever written to disk or the DB**. Every artifact is
returned in memory; any cross-call state lives only in a RAM-only session cache
with a 30-minute TTL. This module never calls ingest_gem, never writes graphs/,
never touches the membership tables.

Public surface
--------------
    register_custom_gem(raw)                 -> (gem_token, summary, warnings)
    resolve_model(source, *, db=None)        -> fresh mutable cobra.Model
    apply_overrides(model, params)           -> cobra.Model (mutated in place)
    compute_mfg(source, params, *, ...)      -> in-memory artifact dict
    extract_pathway(source, target, *, ...)  -> PathwayResult (+ annotation)
    design_strains_table(source, target, *)  -> KO-design table/CSV dict
    rank_design_interventions(source, target, biomass, *) -> unified design payload
    rank_optforce_interventions(source, target, biomass, *) -> OptForce table/CSV dict
    rank_optknock_interventions(source, target, biomass, *) -> OptKnock table/CSV dict

`rank_design_interventions` is the current entry point: it runs the unified
design engine (OptForce + OptKnock pooled under one verification/ranking core)
and returns only the fields that define a buildable design. The two single-engine
functions are retained for backward compatibility.

`source` is a small dict selecting a model-resolution strategy:
    {"strategy": "repo",    "model_id": int}
    {"strategy": "custom",  "gem_token": str}
    {"strategy": "session", "session_token": str}
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import threading
import time
import uuid
from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, ContextManager, Protocol, Sequence, runtime_checkable

import cobra
import torch

from src import db as _db
from src.graph_ops import (
    build_mfg, graph_json, graphml_bytes, run_fba, to_pyg)
from src.schema import Condition, Model
from src.extraction.gudie import GudieConfig
from src.extraction.pathway import (
    PathwayResult, extract_bioproduction_pathway,
    AnnotationBudgetExceeded, InfeasibleModel, InvalidOperatingPoint,
    TargetBlocked, TargetNoMassFlow, TargetNotFound, TargetNotProducible,
)

# ── Tunables (all RAM-only) ──────────────────────────────────────────────────
SESSION_TTL = 1800          # 30-minute session TTL, per spec
MAX_UPLOAD_BYTES = 50 * 1024 * 1024   # custom-GEM upload cap (~50 MB)
MAX_ITEMS = 64              # session cache: LRU cap by item count
MAX_BYTES = 512 * 1024 * 1024         # session cache: LRU cap by total bytes
MAX_CONCURRENT_JOBS = 4    # bound on concurrent extract/annotate jobs


# ── Typed errors (service domain) ────────────────────────────────────────────
class ToolError(Exception):
    """Base for all tool-layer errors. Carries an optional user-actionable hint."""
    def __init__(self, message: str, *, hint: str | None = None):
        super().__init__(message)
        self.hint = hint


class GemUnparseable(ToolError): ...      # E-1  upload is not valid cobra-JSON
class GemTooLarge(ToolError): ...         # E-2  upload exceeds MAX_UPLOAD_BYTES
class GemEmpty(ToolError): ...            # E-2  upload is empty / no reactions/metabolites
class GemNoGenes(ToolError): ...          # E-3  no genes / gene-reaction rules (GPR)
class Infeasible(ToolError): ...          # E-4  FBA/pFBA non-optimal
class InfeasibleMaintenance(ToolError): ...  # E-5  infeasible w/ hard ATPM floor
class TokenNotFound(ToolError): ...       # E-13 unknown session/gem token
class TokenExpired(ToolError): ...        # E-13 token past its TTL
class CacheFull(ToolError): ...           # E-14 too many concurrent jobs
class AnnotationUnavailable(ToolError): ...  # E-15 annotation asked in graph-only mode

# Re-export extract-domain errors for a single import surface.
__all__ = [
    "ToolError", "GemUnparseable", "GemTooLarge", "GemEmpty",
    "GemNoGenes", "Infeasible",
    "InfeasibleMaintenance", "TokenNotFound", "TokenExpired", "CacheFull",
    "AnnotationUnavailable", "TargetNotFound", "TargetBlocked",
    "TargetNotProducible", "TargetNoMassFlow", "InvalidOperatingPoint",
    "AnnotationBudgetExceeded", "CacheBackend", "SessionCache",
    "register_custom_gem", "resolve_model", "apply_overrides",
    "compute_mfg", "extract_pathway", "model_defaults",
    "rank_design_interventions",
]


# ── Cache backend interface ──────────────────────────────────────────────────
# The tool layer talks to its token store only through this small surface, so a
# shared/out-of-process backend (e.g. Redis) can be dropped in later without
# touching any call site. The default implementation below stays RAM-only.
@runtime_checkable
class CacheBackend(Protocol):
    """Minimal token-store contract: TTL'd put/get + a bounded job slot."""
    def put(self, payload: Any, kind: str, *, size_bytes: int = 0) -> str: ...
    def get(self, token: str, *, kind: str | None = None) -> Any: ...
    def job_slot(self) -> ContextManager[None]: ...


# ── RAM-only session cache ───────────────────────────────────────────────────
@dataclass(slots=True)
class _Entry:
    payload: Any
    kind: str          # "gem" | "recipe" | "pathway"
    size_bytes: int
    expires_at: float


class SessionCache:
    """Thread-safe, RAM-only token store: TTL + bounded LRU (by count and bytes).

    Nothing here is ever persisted. Expired entries are dropped lazily on every
    get/put; inserts evict the LRU end until both caps are satisfied. A bounded
    semaphore caps concurrent extract/annotate jobs (model copies dominate RAM).
    """

    def __init__(self, *, ttl: int = SESSION_TTL, max_items: int = MAX_ITEMS,
                 max_bytes: int = MAX_BYTES, max_jobs: int = MAX_CONCURRENT_JOBS):
        self.ttl = ttl
        self.max_items = max_items
        self.max_bytes = max_bytes
        self._data: "OrderedDict[str, _Entry]" = OrderedDict()
        self._bytes = 0
        self._lock = threading.Lock()
        self._sem = threading.BoundedSemaphore(max_jobs)

    def _drop_expired(self, now: float) -> None:
        for tok in [t for t, e in self._data.items() if e.expires_at <= now]:
            self._bytes -= self._data.pop(tok).size_bytes

    def put(self, payload: Any, kind: str, *, size_bytes: int = 0) -> str:
        now = time.monotonic()
        token = uuid.uuid4().hex
        with self._lock:
            self._drop_expired(now)
            self._data[token] = _Entry(payload, kind, size_bytes, now + self.ttl)
            self._bytes += size_bytes
            while (len(self._data) > self.max_items
                   or self._bytes > self.max_bytes) and len(self._data) > 1:
                _, ev = self._data.popitem(last=False)
                self._bytes -= ev.size_bytes
        return token

    def get(self, token: str, *, kind: str | None = None) -> Any:
        now = time.monotonic()
        with self._lock:
            entry = self._data.get(token)
            if entry is None:
                raise TokenNotFound(f"unknown token {token!r}",
                                    hint="the token was never issued or already evicted")
            if entry.expires_at <= now:
                self._bytes -= self._data.pop(token).size_bytes
                raise TokenExpired(f"token {token!r} expired",
                                   hint=f"sessions live {self.ttl}s; recompute to get a fresh token")
            if kind is not None and entry.kind != kind:
                raise TokenNotFound(f"token {token!r} is a {entry.kind}, expected {kind}")
            self._data.move_to_end(token)
            return entry.payload

    @contextmanager
    def job_slot(self):
        """Non-blocking concurrency guard around an extract/annotate job."""
        if not self._sem.acquire(blocking=False):
            raise CacheFull("too many concurrent compute jobs",
                            hint=f"at most {MAX_CONCURRENT_JOBS} extract/annotate jobs at once; retry shortly")
        try:
            yield
        finally:
            self._sem.release()


# Module-level default cache instance.
_CACHE = SessionCache()


# ── Custom-GEM registration ──────────────────────────────────────────────────
# Inorganic-carbon exchanges (CO2, bicarbonate) carry carbon but are not a
# feedstock; excluding them keeps the carbon-source check honest. Matched by
# molecular formula so it is independent of id conventions (mirrors the same
# guard in optforce.pipeline._detect_substrate).
_INORGANIC_CARBON_FORMULAS = {"CO2", "CHO3", "CO3", "HCO3"}

# Compartment suffix on a BiGG-style id (``succ_c`` / ``glc__D_e`` / ``atp_im``):
# a trailing ``_<compartment>`` where the compartment is 1-2 short tokens. Used to
# group compound variants across compartments by their base id.
_COMPARTMENT_SUFFIX = re.compile(r"_([a-z][a-z0-9]?)$")


def base_compound_id(met_id: str) -> str:
    """Strip a trailing compartment suffix (``_c``/``_e``/``_p``/...) from a met id."""
    return _COMPARTMENT_SUFFIX.sub("", met_id)


def _carbon_sources(model: cobra.Model) -> list[str]:
    """Boundary reactions supplying an *organic*-carbon uptake route (lb < 0).

    A carbon source is a boundary/exchange reaction that can carry inward flux
    (negative lower bound) for a carbon-bearing metabolite, excluding pure
    inorganic carbon (CO2 / bicarbonate). Returns the reaction ids, in model order.
    """
    sources: list[str] = []
    for rxn in model.boundary:
        if rxn.lower_bound >= 0:                       # not an open uptake route
            continue
        carbons = sum(met.elements.get("C", 0) for met in rxn.metabolites)
        if carbons <= 0:                               # carries no carbon
            continue
        formulas = {(met.formula or "").replace(" ", "") for met in rxn.metabolites}
        if formulas <= _INORGANIC_CARBON_FORMULAS:     # CO2 / bicarbonate only
            continue
        sources.append(rxn.id)
    return sources


def _looks_like_sbml(text: str) -> bool:
    """Cheap sniff: SBML/XML starts with an XML prolog or an ``<sbml``/``<?xml`` tag."""
    head = text.lstrip()[:512].lower()
    return head.startswith("<?xml") or head.startswith("<sbml") or "<sbml" in head


def _parse_gem_bytes(raw: bytes) -> tuple[cobra.Model, str, str]:
    """Parse an uploaded GEM from either COBRA-JSON or SBML/XML bytes.

    Returns (model, fmt, text) where ``fmt`` is ``"json"`` or ``"sbml"`` and
    ``text`` is the decoded source (cached so resolve_model can re-parse a fresh,
    never-mutated copy). Nothing touches disk: SBML is read from an in-memory
    buffer. Raises GemUnparseable on malformed input.
    """
    try:
        text = raw.decode("utf-8")
    except Exception as exc:
        raise GemUnparseable(f"upload is not UTF-8 text: {exc}",
                             hint="upload a COBRA-JSON or SBML/XML model") from exc
    if _looks_like_sbml(text):
        try:
            model = cobra.io.read_sbml_model(io.StringIO(text))
        except Exception as exc:
            raise GemUnparseable(f"not valid SBML: {exc}",
                                 hint="export with cobra.io.write_sbml_model") from exc
        return model, "sbml", text
    try:
        model = cobra.io.from_json(text)
    except Exception as exc:
        raise GemUnparseable(f"not valid cobra-JSON or SBML: {exc}",
                             hint="export with cobra.io.save_json_model / write_sbml_model") from exc
    return model, "json", text


def _model_from_cached(cached: Any) -> cobra.Model:
    """Re-parse a cached GEM payload into a fresh cobra.Model (format-aware)."""
    if isinstance(cached, dict):          # new {"fmt", "text"} envelope
        if cached.get("fmt") == "sbml":
            return cobra.io.read_sbml_model(io.StringIO(cached["text"]))
        return cobra.io.from_json(cached["text"])
    return cobra.io.from_json(cached)     # legacy: bare JSON string


def register_custom_gem(raw: bytes) -> tuple[str, dict, list[str]]:
    """Validate an uploaded COBRA-JSON or SBML/XML GEM, cache it in RAM, return a token.

    Never persisted. Returns (gem_token, summary, warnings). Metabolites, reactions,
    and genes (GPR) are hard requirements -- each a typed rejection when absent,
    because no downstream engine can run without them. A missing carbon source or
    objective is a warning ("no_carbon_source" / "no_objective"), not a rejection:
    the medium/objective may be set later via apply_overrides.
    """
    if not raw:
        raise GemEmpty("empty upload", hint="POST the GEM's JSON or SBML bytes")
    if len(raw) > MAX_UPLOAD_BYTES:
        raise GemTooLarge(f"upload is {len(raw)} bytes (> {MAX_UPLOAD_BYTES})",
                          hint="reduce the model or raise MAX_UPLOAD_BYTES")
    model, fmt, text = _parse_gem_bytes(raw)

    # 1-2. Metabolites and reactions: a model without either is not a GEM.
    if not model.metabolites or not model.reactions:
        raise GemEmpty("model has no metabolites or reactions",
                       hint="upload a complete GEM with metabolites and reactions")

    # 3. Genes (GPR): strain design (OptForce / KO engines) can only translate a
    #    reaction intervention into a construct through gene-reaction rules.
    if not model.genes or not any((r.gene_reaction_rule or "").strip()
                                  for r in model.reactions):
        raise GemNoGenes(
            "model has no genes / gene-reaction rules (GPR)",
            hint="upload a GEM whose reactions carry gene_reaction_rule (GPR) annotations")

    warnings: list[str] = []

    # 4. Carbon source: a model with no open organic-carbon uptake cannot grow or
    #    produce as-uploaded, but the medium may legitimately be set later via a
    #    carbon_source / condition override -- so this is a warning, not a rejection.
    carbon_sources = _carbon_sources(model)
    if not carbon_sources:
        warnings.append("no_carbon_source")

    try:
        objective_id = str(model.objective.expression)
    except Exception:
        objective_id = ""
    if not objective_id:
        warnings.append("no_objective")

    summary = {
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
        "objective_id": objective_id,
        "carbon_sources": carbon_sources,
        "exchanges": [r.id for r in model.exchanges],
        "format": fmt,
        "content_hash": hashlib.sha256(raw).hexdigest(),
    }
    token = _CACHE.put({"fmt": fmt, "text": text}, "gem", size_bytes=len(raw))
    return token, summary, warnings


# ── Model resolution (three strategies) ──────────────────────────────────────
def resolve_model(source: dict, *, db=None) -> cobra.Model:
    """Return a FRESH mutable cobra.Model for `source`; cached bases never mutated.

    A repo   {"model_id": int}       -> DB path -> load_json_model (fresh from disk)
    C custom {"gem_token": str}      -> from_json of the cached RAM string (fresh parse)
    B session{"session_token": str}  -> replay the stored recipe (resolve + overrides)
    """
    strategy = source.get("strategy")
    if strategy == "repo":
        own = db is None
        session = db or _db.SessionLocal()
        try:
            row = session.get(Model, source["model_id"])
            if row is None:
                raise TokenNotFound(f"no repo model_id={source['model_id']!r}")
            return cobra.io.load_json_model(row.model_path)
        finally:
            if own:
                session.close()

    if strategy == "custom":
        cached = _CACHE.get(source["gem_token"], kind="gem")
        return _model_from_cached(cached)

    if strategy == "session":
        recipe = _CACHE.get(source["session_token"], kind="recipe")
        # Replay faithfully: the original source's own overrides (folded into
        # source_ref) then the recipe's compute params, in the canonical order.
        return resolve_overridden(recipe["source_ref"], recipe["overrides"], db=db)

    raise ToolError(f"unknown resolution strategy {strategy!r}",
                    hint="strategy must be 'repo', 'custom', or 'session'")


# ── Overrides (four axes, applied in order) ──────────────────────────────────
def _apply_bounds(model: cobra.Model, rxn_id: str, lb, ub) -> None:
    """Mirror load_model_with_condition's guarded idiom: skip unknown reactions."""
    if rxn_id in model.reactions:
        model.reactions.get_by_id(rxn_id).bounds = (lb, ub)


def apply_overrides(model: cobra.Model, params: dict, *, db=None) -> cobra.Model:
    """Apply the four override axes in order (condition, reaction_bounds,
    carbon_source, objective) to `model` in place, then return it."""
    params = params or {}

    # 1. condition: DB Condition.exchange_bounds_json, then custom overrides on top.
    cond_id = params.get("condition_id")
    if cond_id:
        own = db is None
        session = db or _db.SessionLocal()
        try:
            row = session.query(Condition).filter(
                Condition.condition_id == cond_id).one_or_none()
            if row is None:
                raise TokenNotFound(f"no condition_id={cond_id!r}")
            bounds = json.loads(row.exchange_bounds_json or "{}")
        finally:
            if own:
                session.close()
        for rxn_id, (lb, ub) in bounds.items():
            _apply_bounds(model, rxn_id, lb, ub)
    for rxn_id, (lb, ub) in (params.get("custom_exchange_bounds") or {}).items():
        _apply_bounds(model, rxn_id, lb, ub)

    # 2. reaction_bounds (knockout = lb=ub=0).
    for rxn_id, (lb, ub) in (params.get("reaction_bounds") or {}).items():
        _apply_bounds(model, rxn_id, lb, ub)

    # 3. carbon-source uptake lower bound.
    cs = params.get("carbon_source")
    if cs and cs.get("exchange_id") in model.reactions:
        rxn = model.reactions.get_by_id(cs["exchange_id"])
        rxn.lower_bound = -abs(float(cs["uptake"]))

    # 4. objective (+ direction).
    obj = params.get("objective_rxn")
    if obj and obj in model.reactions:
        model.objective = obj
        if params.get("direction") in ("min", "max"):
            model.objective_direction = params["direction"]

    return model


# ── Compute ──────────────────────────────────────────────────────────────────
def _recipe_hash(source_ref: dict, overrides: dict) -> str:
    raw = json.dumps([source_ref, overrides], sort_keys=True, default=str)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def _run_fba(model: cobra.Model, use_pfba: bool) -> cobra.Solution:
    """FBA or pFBA, raising a typed Infeasible error on non-optimal status."""
    if use_pfba:
        try:
            return cobra.flux_analysis.pfba(model)
        except Exception as exc:
            raise _infeasible(model, str(exc)) from exc
    sol = run_fba(model)
    if sol is None:
        raise _infeasible(model, "FBA returned non-optimal status")
    return sol


def _infeasible(model: cobra.Model, detail: str) -> ToolError:
    """Distinguish a hard maintenance floor (ATPM) from a plain infeasibility."""
    atpm = next((r for r in model.reactions if r.id.upper() == "ATPM"), None)
    if atpm is not None and atpm.lower_bound > 0:
        return InfeasibleMaintenance(
            f"infeasible with ATPM floor {atpm.lower_bound}: {detail}",
            hint="the medium cannot meet the maintenance (ATPM) demand; loosen uptake bounds")
    return Infeasible(detail, hint="no feasible flux state under these bounds; check the medium/knockouts")


def resolve_overridden(source: dict, extra: dict | None = None, *, db=None) -> cobra.Model:
    """Resolve `source` and apply its override layers, in the one canonical order.

    Applies ``source["overrides"]`` for EVERY strategy (repo / custom / session) --
    not just session -- so ``source.overrides`` behaves consistently across /mfg,
    /pathways and /designs. For a session source the stored recipe is first replayed
    by resolve_model, then ``source.overrides`` is layered on top, then any `extra`
    (e.g. /mfg's explicit params). Within each apply_overrides call the axes keep
    their fixed order: condition -> custom exchange bounds -> reaction bounds/KO ->
    carbon source -> objective/direction.
    """
    model = resolve_model(source, db=db)
    if (ov := source.get("overrides")):
        apply_overrides(model, ov, db=db)
    if extra:
        apply_overrides(model, extra, db=db)
    return model


def _objective_reaction_ids(model: cobra.Model) -> list[str]:
    """Reaction ids carrying a non-zero objective coefficient (model order)."""
    return [r.id for r in model.reactions if r.objective_coefficient != 0]


def _suggested_target(model: cobra.Model, compound: str | None) -> str | None:
    """Best exchange reaction to overproduce `compound` (base id), extracellular
    preferred; None when the compound has no exchange in this model."""
    if not compound:
        return None
    base = base_compound_id(compound)
    candidates: list[tuple[int, str]] = []
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            if base_compound_id(met.id) == base:
                # Prefer an extracellular ('_e') drain, then any exchange.
                candidates.append((0 if met.compartment == "e"
                                   or met.id.endswith("_e") else 1, rxn.id))
                break
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][1]


def model_defaults(model: cobra.Model, *, conditions: Sequence[Condition] | None = None,
                   compound: str | None = None) -> dict:
    """Default operating point the frontend needs to pre-fill a run form.

    Reuses the same primitives as custom-GEM registration (``_carbon_sources``,
    the objective/exchange readout) rather than a second defaults system. Returns
    the objective reaction + direction, the current carbon source and the full set
    of organic-carbon uptake options, every exchange's current bounds, the DB
    condition catalogue (id + options), and -- when `compound` is given -- the
    suggested target exchange to overproduce it.
    """
    objective = _objective_reaction_ids(model)
    carbon_options = _carbon_sources(model)
    try:
        direction = model.objective_direction
    except Exception:
        direction = "max"
    cond_rows = list(conditions or [])
    return {
        "objective_rxn": objective[0] if objective else None,
        "objective_rxns": objective,
        "direction": direction,
        "carbon_source": carbon_options[0] if carbon_options else None,
        "carbon_source_options": carbon_options,
        "exchange_bounds": {r.id: [r.lower_bound, r.upper_bound] for r in model.exchanges},
        "condition_id": cond_rows[0].condition_id if cond_rows else None,
        "condition_options": [c.condition_id for c in cond_rows],
        "target": _suggested_target(model, compound),
    }


def compute_mfg(source: dict, params: dict, *, output_format: str = "stats_json",
                keep_in_session: bool = False) -> dict:
    """Resolve -> apply_overrides -> FBA/pFBA -> build_mfg -> in-memory artifact.

    output_format in {"graphml", "pyg", "fba_json", "stats_json"}. When
    keep_in_session=True the recipe (source + overrides) is memoized under a
    session_token returned in the result for later replay. Nothing is persisted.
    `source.overrides` and `params` are both applied (source.overrides first).
    """
    params = params or {}
    model = resolve_overridden(source, params)
    sol = _run_fba(model, bool(params.get("use_pfba")))
    G = build_mfg(model, sol)

    out: dict[str, Any] = {
        "objective_value": float(sol.objective_value),
        "n_reactions": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
    }
    if output_format == "graphml":
        out["graphml"] = graphml_bytes(G)
    elif output_format == "pyg":
        buf = io.BytesIO()
        torch.save(to_pyg(G), buf)
        out["pyg"] = buf.getvalue()
    elif output_format == "fba_json":
        out["fba_json"] = {r: float(f) for r, f in sol.fluxes.items()}
    elif output_format == "stats_json":
        out["graph_json"] = graph_json(G)
    else:
        raise ToolError(f"unknown output_format {output_format!r}",
                        hint="one of graphml, pyg, fba_json, stats_json")

    if keep_in_session:
        recipe = {"source_ref": source, "overrides": params}
        out["session_token"] = _CACHE.put(recipe, "recipe", size_bytes=1024)
    return out


# ── Pathway extraction (+ FVA design-annotation) ─────────────────────────────
def extract_pathway(source: dict, target: str, *, source_ids: Sequence[str] | None = None,
                    config: GudieConfig | None = None, repair: bool = False,
                    fva_operating_points: Sequence[float] | None = None,
                    fva_loopless: bool = False, db=None) -> PathwayResult:
    """Extract the bioproduction pathway for `target` (+ optional FVA annotation).

    Full mode (repo/custom/session): resolve+override a cobra model, then run
    extract_bioproduction_pathway inside a bounded job slot. Annotation is gated
    behind fva_operating_points (None => today's behaviour, unchanged). `repair`
    (default False, byte-identical output) is forwarded to the existing extractor's
    opt-in connectivity repair -- no new repair algorithm is introduced here.
    `source.overrides` is applied for every strategy (see resolve_overridden).

    Graph-only mode is not available here: without a cobra model neither FASTCC
    nor FVA can run, so a graph-only source that requests annotation raises
    AnnotationUnavailable rather than silently returning an empty set (E-15).
    """
    with _CACHE.job_slot():
        model = resolve_overridden(source, db=db)
        try:
            return extract_bioproduction_pathway(
                model, target, source_ids=source_ids, config=config, repair=repair,
                fva_operating_points=fva_operating_points, fva_loopless=fva_loopless)
        except (TargetNotFound, TargetBlocked, TargetNotProducible, TargetNoMassFlow,
                InfeasibleModel, InvalidOperatingPoint, AnnotationBudgetExceeded):
            raise
        except ValueError as exc:
            raise _infeasible(model, str(exc)) from exc


# ── Unified design engine: tool-layer entry (design-focused payload) ──────────
# Engine kwargs forwarded to design_engine.rank_intervention_candidates. The
# discovery `methods`, the `require_buildable` gate, and `parallel` are first-class
# features (explicit kwargs); the rest are opt-in tuning knobs sent via `params`.
_DESIGN_PARAM_KEYS = frozenset({
    "max_knockouts", "enable_k2", "min_growth_frac", "min_growth_abs", "gamma",
    "eps_prod", "class_tol", "backend", "top_k", "solver_timeout", "gpr_cap",
})
_DESIGN_METHODS = ("optforce", "optknock")


def _design_row(c, intervention_label) -> dict:
    """The fields that DEFINE a buildable design -- what to knock out, what to
    build it as (the verified gene deletion), and what it guarantees. Discovery
    bookkeeping, audit internals, and non-buildable gene alternatives are dropped."""
    build = next((gc for gc in c.gene_constructs if gc.buildable), None)
    return {
        "rank": c.rank,                           # Pareto-ordered position (1 = best)
        "k": c.cardinality,
        # Intervention type is implicit: a ranked design is always a knockout
        # (reaction_knockouts / gene_knockouts, label ":KO"), so engine provenance
        # (source) is not reported; nor is `buildable` (every ranked design is, by
        # construction) or `pareto_front` (rank already encodes the ordering).
        "label": intervention_label(c.reaction_actions),
        "reaction_knockouts": [a.reaction_id for a in c.reaction_actions],
        "gene_knockouts": list(build.genes) if build else [],
        "silenced_reactions": list(build.silenced_reactions) if build else [],
        "guaranteed_floor": c.verified_floor,     # min product flux at the growth optimum
        "max_ceiling": c.verified_ceiling,        # best-case product flux
        "growth": c.verified_growth,              # biomass at the design optimum
        "retained_mu": c.retained_mu,             # fraction of wild-type growth kept
        "coupling": c.coupling_class,             # strong | weak | none
    }


def _hypothesis_row(c, intervention_label) -> dict:
    """A verified regulatory/mixed design -- a lead, not a buildable construct."""
    genes = next((o.genes for o in c.gene_options if o.genes), ())
    return {
        "k": c.cardinality,
        "label": intervention_label(c.reaction_actions),
        # `action` per intervention carries the type (upregulation / downregulation
        # / knockout); engine provenance (source) is not reported.
        "interventions": [{"reaction_id": a.reaction_id, "action": a.action,
                           "forced_flux": round(a.forced_flux, 4)}
                          for a in c.reaction_actions],
        "putative_genes": list(genes),
        "guaranteed_floor": c.verified_floor,
        "max_ceiling": c.verified_ceiling,
        "retained_mu": c.retained_mu,
        "coupling": c.coupling_class,
        "translation_status": c.translation_status,
    }


def rank_design_interventions(source: dict, target: str, biomass: str, *,
                              substrate: str | None = None,
                              methods: Sequence[str] = _DESIGN_METHODS,
                              require_buildable: bool = True,
                              parallel: bool = False,
                              params: dict | None = None,
                              discovery_source: dict | None = None,
                              db=None) -> dict:
    """Run the unified design engine over `source` and return a design-focused payload.

    The engine pools the selected discovery `methods` (OptForce and/or OptKnock)
    into one verification/ranking core and, with `require_buildable=True` (the
    default), ranks only knockout designs that have a concrete, verified
    gene-deletion construct. Regulatory/mixed designs are returned separately as
    hypotheses (leads), never ranked.

    `source`           model-resolution dict (repo / custom / session).
    `target` / `biomass`  product and growth reaction ids (both required).
    `substrate`        optional substrate-uptake id; auto-detected when omitted.
    `methods`          discovery engines to pool (default both).
    `require_buildable` gate the ranked set to buildable gene-KO designs (default True).
    `parallel`         run the discovery engines concurrently (one model copy each).
    `params`           PARTIAL tuning dict (`max_knockouts`, `enable_k2`, `gamma`,
                       `top_k`, `backend`, ...); absent fields take engine defaults,
                       unknown keys are rejected.
    `discovery_source` optional reduced discovery model (same resolution dict).

    Returns a dict whose headline is ``designs`` -- the ranked, buildable designs
    reduced to the fields that define a design (see :func:`_design_row`): the
    reaction knockouts, the gene deletion that realizes them, and the guaranteed
    production floor / ceiling / retained growth / coupling. Also returns the
    self-describing ``table`` + download ``csv`` (design columns). The regulatory
    leads are returned both as the ``hypotheses`` rows and as their own
    ``hypotheses_table`` + ``hypotheses_csv`` grid (parallel to the designs' grid),
    alongside the wild-type ``baseline``, the ``methods`` actually run, and a
    ``rejected`` summary (count + reasons -- e.g. how many verified KO designs were
    dropped as not-buildable). Model resolution/overrides, the RAM-only
    job slot, and typed base-model errors behave exactly as in the single-engine
    entries.
    """
    import src.design_engine as _de

    methods = tuple(methods)
    bad = [m for m in methods if m not in _DESIGN_METHODS]
    if bad:
        raise ToolError(f"unknown discovery method(s): {bad}",
                        hint=f"methods must be a subset of {list(_DESIGN_METHODS)}")
    params = dict(params or {})
    unknown = set(params) - _DESIGN_PARAM_KEYS
    if unknown:
        raise ToolError(f"unknown design param(s): {sorted(unknown)}",
                        hint=f"allowed params: {sorted(_DESIGN_PARAM_KEYS)}")

    with _CACHE.job_slot():
        model = resolve_overridden(source, db=db)
        disc = None
        if discovery_source is not None:
            disc = resolve_overridden(discovery_source, db=db)
        result = _de.rank_intervention_candidates(
            model, target, biomass, substrate, discovery_model=disc,
            methods=methods, require_buildable=require_buildable,
            parallel=parallel, **params)

    label = _de.intervention_label
    rej_reasons: dict[str, int] = {}
    for c in result.rejected:
        key = c.rejection_reason or c.import_loop_status
        rej_reasons[key] = rej_reasons.get(key, 0) + 1

    return {
        "designs": [_design_row(c, label) for c in result.accepted],
        "hypotheses": [_hypothesis_row(c, label) for c in result.hypotheses],
        "table": _de.ranking_table(result),
        "csv": _de.ranking_csv(result),
        # Regulatory/mixed leads rendered as their own grid (+ download CSV),
        # parallel to `table`/`csv` for the ranked designs.
        "hypotheses_table": _de.hypotheses_table(result),
        "hypotheses_csv": _de.hypotheses_csv(result),
        "baseline": {
            "mu_wt": result.mu_wt,
            "floor_wt": result.floor_wt,
            "ceiling_wt": result.ceiling_wt,
            "substrate_rxn": result.substrate_rxn,
        },
        "methods": list(result.methods),
        "require_buildable": require_buildable,
        "rejected": {"n": len(result.rejected), "reasons": rej_reasons},
        "pareto_objectives": list(result.pareto_objectives),
        "extraction_fallback": result.extraction_fallback,
        "extraction_failure_reason": result.extraction_failure_reason,
        "timings": result.timings,
        "notes": list(result.notes),
        "params": params,
    }
