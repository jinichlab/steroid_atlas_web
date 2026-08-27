# Steroid Atlas — Web

A **Next.js + deck.gl** frontend for the [Steroid Atlas](https://github.com/jinichlab/steroid_atlas):
an interactive UMAP of **14,089 steroid-interacting proteins**, **677 steroid
molecules**, and **2,889 natural + synthetic steroid entries**, curated from
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
git clone git@github.com:jinichlab/steroid_atlas_web.git
cd steroid_atlas_web
npm install
npm run dev -- --port 3000
```

Open http://localhost:3000. On a remote machine, tunnel the port from your
laptop:

```bash
ssh -N -L 3000:localhost:3000 <user>@<server>
```

## Chatbot / RAG

"Ask the Atlas" (floating widget, bottom-right) answers questions grounded in a
prebuilt corpus of ChEBI compound records and UniProt protein literature. On an
Explore page it also sees your current selection (lasso / search / cluster), so
"what are these?" works — the selected items are shared with the chat via
`src/lib/selection-context.tsx` and shown as a "N selected" chip in the widget.

Two processes:

1. **Retrieval sidecar** — a local Python service that owns the 2 GB FAISS index.
   ```bash
   npm run dev:rag          # loads the index (~20 s), serves 127.0.0.1:8000
   ```
2. **Next.js** — `npm run dev` as usual. `/api/chat` embeds each question, asks
   the sidecar for the top passages, and prompts OpenAI with them + a
   cite-your-sources instruction. If the sidecar is down the chat still answers,
   flagged as ungrounded.

Setup:

- `cp .env.example .env.local` and set `OPENAI_API_KEY`.
- The vector store lives in `rag_data/rag_store/` (`index.faiss` + `catalog.jsonl`,
  ~2.2 GB, gitignored — copy it in out of band). Details and a systemd unit for
  deployment: [`scripts/rag/README.md`](scripts/rag/README.md).
- Needs a persistent host (~4 GB RAM) — **not** Vercel serverless.

## Routes

| URL                  | Page                                                 |
|----------------------|------------------------------------------------------|
| `/`                  | Landing — stats + three view cards                   |
| `/explore/protein`   | UMAP + catalogue for the 14,089-protein view         |
| `/explore/molecule`  | UMAP + catalogue for the 677-steroid view            |
| `/explore/natsyn`    | UMAP + catalogue for the 2,889 nat + syn entries     |

## Repo layout

```
src/
  app/
    page.tsx                        Landing page
    explore/[view]/                 Three views share one dynamic route
      page.tsx                      Loads per-view cluster metadata
      explore-client.tsx            All interaction logic (client component)
  components/
    UmapPlot.tsx                    deck.gl scatter + icon + text layers
  lib/
    palette.ts                      Golden-angle palette + spatial re-mapping
    types.ts                        Protein / Molecule / ClusterMeta types
public/
  atlas/
    proteins.json                   14,089 rows  (~23 MB)
    molecules.json                    677 rows   (~3 MB)
    natsyn.json                     2,889 rows   (~4 MB)
    protein_clusters.json              82 clusters (fingerprint stems + top GO)
    molecule_clusters.json              8 clusters (Bile acids, Estrogens, …)
    natsyn_clusters.json                2 clusters (natural / synthetic)
    structures/                    ~2,441 pre-rendered 260×260 PNGs
    structures_index.json          name / CHEBI → PNG filename lookup
    summary.json                   Landing-page counts
scripts/
  build_atlas_data.py              CSV → JSON pre-processing
  build_structures.py              RDKit renders 2D structure PNGs
  rag/rag_server.py                RAG retrieval sidecar (FAISS + OpenAI embeddings)
rag_data/rag_store/                index.faiss + catalog.jsonl (gitignored, ~2.2 GB)
```

## Rebuilding the data

The bundled JSON + PNGs are regenerated from the marimo repo's CSVs whenever
they change. Clone both repos side-by-side and run:

```bash
# from steroid_atlas_web/
LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 scripts/build_atlas_data.py
LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 scripts/build_structures.py
```

`build_atlas_data.py` reads from `../steroid_atlas/data/` and writes to
`public/atlas/`; `build_structures.py` needs RDKit (available in the
project's miniconda env).

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
