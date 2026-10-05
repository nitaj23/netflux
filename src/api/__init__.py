"""NetFlux API - a thin FastAPI routing layer over the backend.

Owns request/response schemas, validation, error mapping and DB injection;
implements no business logic, delegating to ``src.tool`` and ``src.queries``.

    from src.api import app             # the ASGI application
    from src.api.app import create_app  # or build a fresh one (tests)

Endpoints (under /api/v1): GET /health; GET /models[/{id}]; GET /conditions[/{id}];
GET /reactions/{bigg}/gems; GET /metabolites/{bigg}/gems; POST /gems; POST /mfg;
POST /pathways; POST /designs.
"""
from src.api.app import app, create_app

__all__ = ["app", "create_app"]
