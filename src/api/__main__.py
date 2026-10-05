"""Run the API with ``python -m src.api``.

A thin convenience wrapper around uvicorn; honours HOST / PORT / API_RELOAD
environment variables. For production use uvicorn/gunicorn directly against
``src.api.app:app``.
"""
from __future__ import annotations

import os


def main() -> None:
    import uvicorn

    uvicorn.run(
        "src.api.app:app",
        host=os.getenv("HOST", "127.0.0.1"),
        port=int(os.getenv("PORT", "8000")),
        reload=os.getenv("API_RELOAD", "").lower() in ("1", "true"),
    )


if __name__ == "__main__":
    main()
