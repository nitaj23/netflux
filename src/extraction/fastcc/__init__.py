"""
src.extraction.fastcc

FASTCC flux-consistency boundary (algorithm) + the FASTCC-derived reduced GEMs.

This package replaces the former single-module ``extraction/fastcc.py`` -- the
public API is re-exported unchanged, so existing imports keep working:

    from src.extraction.fastcc import fastcc, fastcc_from_cobra, FastccResult

Derived models
--------------
``models/<key>.json`` are FASTCC-reduced copies of the genome-scale models: the
full GEM with all flux-BLOCKED reactions removed, leaving only the consistent
set. They load exactly like a GEM and are usable as drop-in reduced models:

    from src.extraction.fastcc import load_fastcc_model
    model = load_fastcc_model("iML1515")      # cobra.Model, consistent reactions only

``models/metadata.json`` catalogs them (organism, source GEM, epsilon, sizes).
Regenerate with:  ``.venv/Scripts/python.exe -m scripts.build_fastcc_models``
"""
from __future__ import annotations

import json
from pathlib import Path

from .fastcc import FastccResult, fastcc, fastcc_from_cobra, _lp7

__all__ = ["FastccResult", "fastcc", "fastcc_from_cobra", "_lp7",
           "MODELS_DIR", "load_fastcc_model", "fastcc_metadata", "available_models"]

MODELS_DIR = Path(__file__).resolve().parent / "models"


def fastcc_metadata() -> dict:
    """Catalog of derived models (keyed by model id)."""
    meta = MODELS_DIR / "metadata.json"
    return json.loads(meta.read_text()) if meta.exists() else {}


def available_models() -> list[str]:
    """Model keys with a generated FASTCC-reduced JSON."""
    return sorted(p.stem for p in MODELS_DIR.glob("*.json") if p.stem != "metadata")


def load_fastcc_model(key: str):
    """Load the FASTCC-reduced cobra model for `key` (like a genome-scale model)."""
    import cobra
    path = MODELS_DIR / f"{key}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"no FASTCC-reduced model for {key!r} in {MODELS_DIR}. "
            f"available: {available_models()}  "
            f"(build with: python -m scripts.build_fastcc_models)")
    return cobra.io.load_json_model(str(path))
