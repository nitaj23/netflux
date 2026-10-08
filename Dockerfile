# NetFlux backend image. Used by both docker-compose (local) and Cloud Run (prod).
# The frontend is NOT built here: locally it runs from a node image via compose,
# in production it is built by Vercel.
FROM python:3.12-slim

# build-essential covers the few scientific wheels that may compile from source.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install deps first so the layer caches across source-only changes.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# App source + committed data (gem_models/, graphs/).
COPY . .

# Build the read-only SQLite DB from the committed GEMs/graphs at build time.
# mfghub.db is gitignored, so it is always regenerated here, never shipped in.
RUN python -m src.build_db --reset

# Cloud Run injects $PORT (8080); compose overrides it to 8000. The CMD honors
# whichever is set, so the same image is portable. Single worker: each worker
# loads the full solver stack, so more than one risks OOM on small instances.
ENV PORT=8080
CMD ["sh", "-c", "uvicorn src.api.app:app --host 0.0.0.0 --port ${PORT} --workers 1"]
