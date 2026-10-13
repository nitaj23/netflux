# NetFlux

NetFlux is a web tool Built by iGEM MIT-MAHE 2026 for metabolic engineering.

## Description

A genome-scale metabolic model (GEM) contains thousands of reactions. Finding the ones that turn a carbon source into a target product, and deciding which genes to change, normally takes several separate tools and a lot of compute. NetFlux does both in one place.

- **Pathway extraction (SCOPE).** NetFlux represents metabolism as a Mass Flow Graph, in which each edge records how much metabolite mass flows from one reaction to the next. SCOPE removes reactions that can never carry flux (FASTCC), finds the flux state that maximises your product, and isolates the source-to-product corridor (GUDIE: personalised PageRank plus flux coupling). A network of thousands of reactions shrinks to a pathway small enough to read.
- **Strain design.** The design engine pools proposals from OptForce and OptKnock and checks each one on the full model. It translates every reaction knockout into a buildable gene deletion and ranks the designs by guaranteed production, retained growth, robustness and maximum production. Regulatory designs (up- or down-regulation) are returned separately as leads.
- **Models.** Eleven curated models are included (table below). You can also upload your own model in COBRA-JSON or SBML format. Uploads are kept in memory for 30 minutes and are never saved.

| Model | Organism |
| --- | --- |
| iML1515 | *Escherichia coli* K-12 MG1655 |
| iJN1026 | *Bacillus subtilis* |
| iJN1463 | *Pseudomonas putida* KT2440 |
| iCN1361 | *Cupriavidus necator* H16 |
| iAF987 | *Geobacter metallireducens* GS-15 |
| iCTH669 | *Clostridium thermocellum* |
| iNF517 | *Lactococcus lactis* subsp. *cremoris* MG1363 |
| iJN678 | *Synechocystis* sp. PCC 6803 |
| iAF692 | *Methanosarcina barkeri* |
| iMM904 | *Saccharomyces cerevisiae* S288C |
| iLB1027 | *Phaeodactylum tricornutum* |

The backend is a FastAPI service written in Python using COBRApy. The frontend is a React and TypeScript app built with Vite.

## Installation

### Option 1: Docker (recommended)

Requires [Docker](https://docs.docker.com/get-docker/) with Docker Compose.

```bash
git clone https://gitlab.igem.org/2026/software/mit-mahe/netflux.git
cd netflux
docker compose up --build
```

Then open http://localhost:5173. The first build takes several minutes because it installs the solvers and builds the model database.

### Option 2: Run locally

Requires Python 3.12 and Node.js 20 or newer.

**Backend**

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m src.build_db --reset   # builds mfghub.db from the models in gem_models/
python -m src.api                # serves the API on http://127.0.0.1:8000
```

On Debian or Ubuntu, install `libexpat1` first if `import cobra` fails with a `libexpat.so.1` error.

**Frontend** (in a second terminal)

```bash
cd netflux-ui
npm install
npm run dev                      # http://localhost:5173, proxies /api to the backend
```

## Usage

### In the browser

1. **Compound:** type the product you want, e.g. `succinate`.
2. **Organism:** pick one of the organisms that contain it. Those that already make it are listed first. To use your own model instead, click **Upload** and choose a `.json` or `.xml` file.
3. Substrate, flux bounds and objective fill in automatically. Change them if you need to.
4. Click **Generate**. The pathway appears on the right. Reactions are arrows and metabolites are circles, and the highlighted line is the main route from substrate to product. Hover over a reaction to see its flux, genes and coupling. Each Generate opens a new tab, so you can compare settings side by side.
5. Open **Interventions** to run strain design. It takes up to a few minutes, then shows the ranked knockouts with the genes to delete.
6. **Download** saves the pathway as `.graphml` or `.json`, the source model, or the interventions as `.csv`.

### Through the API

Interactive API documentation is at http://localhost:8000/docs. All endpoints are under `/api/v1`.

```bash
# Find the organisms that contain a compound
curl "http://localhost:8000/api/v1/compounds?q=succ"
curl "http://localhost:8000/api/v1/compounds/succ/organisms"

# Look up the E. coli model's id (replace <id> below with it)
curl "http://localhost:8000/api/v1/models?organism=coli"

# Extract the production pathway
curl -X POST http://localhost:8000/api/v1/pathways \
  -H "Content-Type: application/json" \
  -d '{"source": {"strategy": "repo", "model_id": <id>}, "target": "EX_succ_e", "include_graphml": true}'

# Start a strain-design job, then poll it until "status" is "completed"
curl -X POST http://localhost:8000/api/v1/designs \
  -H "Content-Type: application/json" \
  -d '{"source": {"strategy": "repo", "model_id": <id>}, "target": "EX_succ_e", "biomass": "BIOMASS_Ec_iML1515_core_75p37M"}'
curl http://localhost:8000/api/v1/jobs/<job_id>
```

To use an uploaded model, send it to `POST /api/v1/gems/upload` and use `{"strategy": "custom", "gem_token": "<token>"}` as the `source`.

## Files

```
netflux/
├── src/                    # Python backend
│   ├── api/                # FastAPI app, routes, schemas, async job runner
│   ├── extraction/         # SCOPE: MFG → FASTCC → GUDIE → FVA pathway extraction
│   ├── design_engine/      # OptForce + OptKnock discovery, verification, ranking
│   ├── graph_ops.py        # Mass Flow Graph construction, FBA/pFBA
│   ├── ingest.py           # model ingest pipeline
│   ├── build_db.py         # builds the model database from gem_models/
│   ├── schema.py / db.py / queries.py   # database layer
│   └── tool.py             # on-demand compute entry points
├── netflux-ui/             # React + Vite browser interface
├── gem_models/             # curated genome-scale models (COBRA-JSON)
├── graphs/                 # generated Mass Flow Graph artifacts
├── scripts/                # model-building, validation, and evaluation scripts
├── Dockerfile              # backend container image
├── docker-compose.yaml     # full-stack local deployment
├── requirements.txt        # Python dependencies (pinned)
└── mfghub.db               # prebuilt model database (SQLite)

```

## Contributing

The most valuable contribution you can make is **a new genome-scale metabolic model** for our curated database. NetFlux is only as broad as the organisms it knows about, and we actively welcome well-curated GEMs — especially for organisms relevant to biomanufacturing that are not yet covered.

To submit a model:

1. Make sure the model is a valid **COBRA-JSON** or **SBML/XML** file that loads in COBRApy, with a defined biomass objective and exchange reactions.
2. Confirm it grows and carries flux under at least one medium condition — you can check this quickly by uploading it through the UI or `POST /api/v1/gems` and running an FBA.
3. Open an issue with the model file (or a link to it), its source/reference, the organism, and any medium conditions it was curated for.

We review each submission for consistency before adding it to the repository database. Please include a citation or provenance for the model so it can be properly attributed.

You can also send us data at mitmahe2026software@gmail.com


## Authors and acknowledgment

NetFlux was built by **Team MIT-MAHE 2026** for the iGEM competition. We thank our instructors, advisors, and the wider iGEM community for their support.

NetFlux stands on a body of open methods and tools from the metabolic-modelling and network-science communities. We gratefully acknowledge the following key works, which the pipeline directly builds on:

- Beguerisse-Díaz, M., Bosque, G., Oyarzún, D., Picó, J., & Barahona, M. (2018). Flux-dependent graphs for metabolic networks. *npj Systems Biology and Applications, 4*(1), 32. https://doi.org/10.1038/s41540-018-0067-y
- Burgard, A. P., Nikolaev, E. V., Schilling, C. H., & Maranas, C. D. (2004). Flux coupling analysis of genome-scale metabolic network reconstructions. *Genome Research, 14*(2), 301–312. https://doi.org/10.1101/gr.1926504
- Burgard, A. P., Pharkya, P., & Maranas, C. D. (2003). OptKnock: A bilevel programming framework for identifying gene knockout strategies for microbial strain optimization. *Biotechnology and Bioengineering, 84*(6), 647–657. https://doi.org/10.1002/bit.10803
- Ebrahim, A., Lerman, J. A., Palsson, B. O., & Hyduke, D. R. (2013). COBRApy: Constraints-based reconstruction and analysis for Python. *BMC Systems Biology, 7*(1), 74. https://doi.org/10.1186/1752-0509-7-74
- Hasibi, R., Michoel, T., & Oyarzún, D. A. (2024). Integration of graph neural networks and genome-scale metabolic models for predicting gene essentiality. *npj Systems Biology and Applications, 10*(1), 24. https://doi.org/10.1038/s41540-024-00348-2
- Orth, J. D., Thiele, I., & Palsson, B. Ø. (2010). What is flux balance analysis? *Nature Biotechnology, 28*(3), 245–248. https://doi.org/10.1038/nbt.1614
- Ranganathan, S., Suthers, P. F., & Maranas, C. D. (2010). OptForce: An optimization procedure for identifying all genetic manipulations leading to targeted overproductions. *PLoS Computational Biology, 6*(4), e1000744. https://doi.org/10.1371/journal.pcbi.1000744
- Silva, M. I., Aparício, D., Malveiro, B., Ascensão, J. T., & Bizarro, P. (2021). GUDIE: A flexible, user-defined method to extract subgraphs of interest from large graphs. *arXiv*. https://doi.org/10.48550/arXiv.2108.09200
- Tong, H., Faloutsos, C., & Pan, J. (2006). Fast random walk with restart and its applications. *Proceedings of the Sixth International Conference on Data Mining (ICDM)*, 613–622. https://doi.org/10.1109/ICDM.2006.70
- Vlassis, N., Pacheco, M. P., & Sauter, T. (2014). Fast reconstruction of compact context-specific metabolic network models. *PLoS Computational Biology, 10*(1), e1003424. https://doi.org/10.1371/journal.pcbi.1003424


## License
This repository is licensed under the [Apache License 2.0](LICENSE)
