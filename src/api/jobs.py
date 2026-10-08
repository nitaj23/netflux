"""Bounded in-process job runner for async /designs.

A single ThreadPoolExecutor caps concurrency (default 4 workers == at most 4
design jobs running at once; extras stay ``queued`` until a worker frees up).
Job state is RAM-only and ephemeral -- nothing is persisted, matching the
tool-layer invariant. No Celery/Redis/broker: this is deliberately the smallest
thing that gives queued/running/completed/failed with typed-error retention.
"""
from __future__ import annotations

import os
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any, Callable

# Env-configurable so deployments (compose, Cloud Run) can tune concurrency to
# the instance's RAM without a code change. Defaults to 4 to preserve prior behaviour.
MAX_CONCURRENT_DESIGN_JOBS = int(os.getenv("NETFLUX_MAX_DESIGN_JOBS", "4"))


def default_error_payload(exc: BaseException) -> dict:
    """Serialise a failed job's exception into the {error, detail, hint} contract."""
    return {
        "error": type(exc).__name__,
        "detail": str(exc) or type(exc).__name__,
        "hint": getattr(exc, "hint", None),
    }


@dataclass
class Job:
    id: str
    status: str = "queued"           # queued | running | completed | failed
    result: Any = None
    error: dict | None = None


class JobManager:
    """Thread-safe registry + bounded executor. One instance per process."""

    def __init__(self, *, max_workers: int = MAX_CONCURRENT_DESIGN_JOBS,
                 error_serializer: Callable[[BaseException], dict] = default_error_payload):
        self._executor = ThreadPoolExecutor(
            max_workers=max_workers, thread_name_prefix="design-job")
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()
        self._serialize_error = error_serializer

    def submit(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> str:
        job = Job(id=uuid.uuid4().hex)
        with self._lock:
            self._jobs[job.id] = job

        def _run() -> None:
            with self._lock:
                job.status = "running"
            try:
                result = fn(*args, **kwargs)
            except BaseException as exc:  # noqa: BLE001 - retain any typed error
                with self._lock:
                    job.error = self._serialize_error(exc)
                    job.status = "failed"
            else:
                with self._lock:
                    job.result = result
                    job.status = "completed"

        self._executor.submit(_run)
        return job.id

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)

    def shutdown(self, *, wait: bool = False) -> None:
        self._executor.shutdown(wait=wait)
