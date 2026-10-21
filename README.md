---
title: Steroid Atlas
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
---

# Steroid Atlas — Web

A **Next.js + deck.gl** frontend for the [Steroid Atlas](https://github.com/jinichlab/steroid_atlas):
an interactive UMAP of **14,089 steroid-interacting proteins**, **676 steroid
molecules**, and **2,976 natural + synthetic steroid entries**, curated from
UniProt, Rhea, ChEBI, RefSeq, and hand-audited 2024–2026 literature
recruitments.

Companion to the marimo notebook viewer in `jinichlab/steroid_atlas` — same
underlying data, redesigned as a clone-and-deploy web app with faster WebGL
rendering, richer hover cards, structure thumbnails on every steroid tile,
and a multi-column search bar.

## Features

- **WebGL UMAP** (deck.gl `ScatterplotLayer` + `IconLayer`) — 14 k points
  in a single frame, smooth pan + zoom, no jank
- **Three views**: *Protein centric* (82 clusters), *Steroid centric*
  (8 family-labeled clusters), *Natural + synthetic* (2 groups)
- **Zoom-conditional cluster labels** at each centroid (family stem +
  cluster id), just like the marimo Vega version
- **Marimo-style hover tooltip**: name · accession · gene · organism ·
  cluster · EC, with a `★ NEW` badge for literature-recruited proteins
- **Free-text search** across name, gene, accession, GO, EC, ChEBI, Rhea,
  keyword, and sequence — matches stay lit on the map, everything else
  greys out, and the pool below the map narrows to the matches
- **Pan / Select toggle** (top-left of the plot) — drag = pan (default) or
  drag = lasso rectangle; the last drawn rectangle **persists** on the map
  for reference until you clear it
- **Steroid catalogue tiles** with pre-rendered 2D structure PNGs (RDKit
  → 2,441 unique files) — click a tile → per-compound protein list → tick
  any rows to open their full detail cards side-by-side
- **Multi-select detail cards** with EC badges, Rhea links, GO / Keyword
  chips, PubMed refs, UniProt + AlphaFold links, and a 3-column grid of
  interacting-steroid thumbnails per protein
- **Spatially-aware distinct color palette** — clusters that are neighbors
  on the map get maximally-contrasting hues (greedy graph coloring)
- **Newly-recruited proteins render as stars**, tinted by their cluster
  color so you still see which cluster they belong to

## Quickstart

```bash
git clone https://github.com/jinichlab/steroid_atlas_web.git
cd steroid_atlas_web/web
npm install
npm run dev -- --port 3000
```

## Chatbot / RAG

"Ask the Atlas" (floating widget, bottom-right) answers questions grounded in a
prebuilt corpus of ChEBI compound records and UniProt protein literature. On an
Explore page it also sees your current selection (lasso / search / cluster), so
"what are these?" works — the selected items are shared with the chat via
`web/src/lib/selection-context.tsx` and shown as a "N selected" chip in the widget.

Two processes:

1. **Retrieval sidecar** — a local Python service that owns the 2 GB FAISS index.
   ```bash
   npm run dev:rag          # from web/; loads the index (~20 s), serves 127.0.0.1:8000
   ```
2. **Next.js** — `npm run dev` as usual. `/api/chat` embeds each question, asks
   the sidecar for the top passages, and prompts OpenAI with them + a
   cite-your-sources instruction. If the sidecar is down the chat still answers,
   flagged as ungrounded.

Setup:

- `cp web/.env.example web/.env.local` and set `OPENAI_API_KEY`.
- The vector store lives in `rag_data/rag_store/` (`index.faiss` + `catalog.jsonl`,
  ~2.2 GB, gitignored — copy it in out of band). Details and a systemd unit for
  deployment: [`rag/README.md`](rag/README.md).
- Needs a persistent host (~4 GB RAM) — **not** Vercel serverless.

## Routes

| URL                  | Page                                                 |
|----------------------|------------------------------------------------------|
| `/`                  | Landing — stats + three view cards                   |
| `/explore/protein`   | UMAP + catalogue for the 14,089-protein view         |
| `/explore/molecule`  | UMAP + catalogue for the 676-steroid view            |
| `/explore/natsyn`    | UMAP + catalogue for the 2,976 nat + syn entries     |

## Repo layout

```
web/                     Next.js web app — run npm commands from here
  src/app/               pages + API routes (/api/chat, /api/config)
  src/components/        UmapPlot (deck.gl), ChatWidget, API-key UI
  src/lib/               palette, types, key mode, RAG client
  public/atlas/          generated web data: JSON + 2,441 structure PNGs
rag/                     "Ask the Atlas" retrieval service + corpus/index builders
data/                    canonical CSVs (proteins, molecules, natural + synthetic)
  provenance/            paper-audit record for the literature-recruited entries
analysis/                every paper analysis — see analysis/README.md
  scripts/  results/  figures/       reproducible from data/
  historical/                        how proteins.csv was built (record only)
scripts/                 reproduce_analysis.sh, build_atlas_data.py,
                         build_structures.py, deploy_space.sh, start_space.sh
docs/                    methods + changelog
rag_data/                gitignored: 2.1 GB index + raw corpus
Dockerfile               Hugging Face Space image (web/ + rag/)
requirements.txt         Python deps for analysis/, scripts/, rag/
```

## Rebuilding the data

Everything downstream of `data/*.csv` — analyses, figures and the web app's
`web/public/atlas/` — rebuilds with one command:

```bash
pip install -r requirements.txt
sh scripts/reproduce_analysis.sh          # ~2 min; add --no-web to skip web/public/atlas/
```

To rebuild only the web data: `python3 scripts/build_atlas_data.py` and
`python3 scripts/build_structures.py` (needs RDKit).

## UMAP interactions

- **drag** — pan the map (default). Toggle to **Select** in the top-left
  to switch drag → lasso rectangle
- **scroll** — zoom in / out. Cluster labels appear at moderate zoom
- **click a dot** — open the summary card with UniProt / AlphaFold links
- **hover** — rich tooltip with name / accession / gene / organism /
  cluster / EC
- **type in the search bar** — instantly filter the map + the catalogue
- **cluster sidebar** — click any chip to highlight that cluster's dots;
  the family name is next to the id, so you always know what you picked
- **catalogue tiles** — click a steroid tile → panel of proteins that act
  on it → tick any protein to append its detail card to the stack below

## Stack

- Next.js 14 (App Router) · TypeScript · Tailwind CSS
- deck.gl (ScatterplotLayer + IconLayer + TextLayer) for WebGL scatter
- lucide-react (Pan / Select icons)
- Node.js 20 LTS

## Related

- [`jinichlab/steroid_atlas`](https://github.com/jinichlab/steroid_atlas) —
  canonical data, marimo visualizer, analysis scripts. All data changes
  happen there first; this repo re-imports via `scripts/build_atlas_data.py`.

## License

MIT
