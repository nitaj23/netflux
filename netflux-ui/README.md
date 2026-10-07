# NetFlux UI

React + TypeScript + Vite frontend for the NetFlux API (`/api/v1`).

## Run

```bash
# backend (from netflux/): build the DB once, then serve
python -m src.build_db
uvicorn src.api.app:app --port 8000

# frontend (from netflux/netflux-ui/)
npm install
npm run dev          # http://localhost:5173, proxies /api -> 127.0.0.1:8000
```

Env: `NETFLUX_API` (dev-proxy target, default `http://127.0.0.1:8000`);
`VITE_API_BASE` (API base for builds served from another origin; that origin
must be listed in the server's `NETFLUX_CORS_ORIGINS`).

## Layout

- `src/api/` typed client, `{error, detail, hint}` envelope -> `ApiError`, one function per endpoint
- `src/components/ParamPanel.tsx` Compound / Organism / Substrate / Flux / Objective form, repo and uploaded-GEM routes
- `src/components/ResultsPane.tsx` Subgraph tabs, Interventions tab, download menu
- `src/components/SubgraphView.tsx` SVG pathway map (pan/zoom, details sidebar); layout in `src/lib/pathwayLayout.ts`, a port of the Escher-style map in `scripts/pathway_map_template.html` plus subsystem lanes
- `src/components/InterventionsView.tsx` job polling and result tables

Conventions: never sends `condition_id`; Objective is sent as `overrides.objective_rxn` and as `biomass`;
Substrate goes to `/designs` only; flux bounds apply to the Substrate exchange (or the model's default
carbon source when Substrate is empty).
