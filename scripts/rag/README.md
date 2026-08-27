# RAG retrieval sidecar

A small FastAPI service that loads the prebuilt FAISS index + chunk catalog once
and answers similarity searches. The Next.js `/api/chat` route calls it before
prompting the chat model, so answers are grounded in the Steroid Atlas corpus.

The web app never touches the 2 GB index directly — it only speaks HTTP to this
service (`src/lib/rag.ts`). Swapping in a hosted vector DB later means
reimplementing `POST /search`, nothing else.

The HTTP layer is Python stdlib (`http.server`); the only third-party imports
are `faiss`, `numpy` and `openai`, all already in the project conda env.

## Data

Lives in `rag_data/rag_store/` (gitignored — copy it to the server out of band):

| file | size | notes |
|---|---|---|
| `index.faiss` | ~2.0 GB | `IndexFlatIP`, dim **3072**, **168,739** vectors |
| `catalog.jsonl` | ~222 MB | 168,739 lines, line _N_ ↔ vector _N_ |

Dim 3072 ⇒ the index was built with OpenAI `text-embedding-3-large`. Query
embeddings **must** use the same model (`EMBEDDING_MODEL`) — the server asserts
`len(embedding) == index.d` and returns 500 on mismatch rather than silently
returning garbage.

## Run (dev)

```bash
cd /path/to/steroid_atlas_web
export $(grep OPENAI_API_KEY .env.local)      # sidecar does NOT read .env.local
LD_LIBRARY_PATH=$HOME/miniconda3/lib \
RAG_STORE_DIR=$PWD/rag_data/rag_store \
$HOME/miniconda3/bin/python3 scripts/rag/rag_server.py
```

or just `npm run dev:rag` (same command, wired in `package.json`).

Startup reads ~2.2 GB off disk and holds ~2.3 GB RSS. Expect **10–40 s** before:

```
[rag] loaded index: ntotal=168739 dim=3072 model=text-embedding-3-large threads=4 (23.4s)
[rag] serving on http://127.0.0.1:8000
```

The process only starts listening after the index finishes loading, so a
connection refused during startup is expected; `/health` returns 503 only in the
brief window the socket is up but state is not.

## Environment

| var | default | purpose |
|---|---|---|
| `OPENAI_API_KEY` | — | **required**; used for query embeddings |
| `RAG_STORE_DIR` | `<repo>/rag_data/rag_store` | index + catalog location |
| `EMBEDDING_MODEL` | `text-embedding-3-large` | must yield dim 3072 |
| `RAG_SERVER_HOST` | `127.0.0.1` | keep on loopback unless firewalled |
| `RAG_SERVER_PORT` | `8000` | uvicorn port |
| `RAG_FAISS_THREADS` | `4` | `faiss.omp_set_num_threads` |

## API

### `GET /health`
`200 {"status":"ok","ntotal":168739,"dim":3072,"model":"text-embedding-3-large"}`
· `503 {"status":"loading"}` while the index loads.

### `POST /search`
Request:
```json
{ "query": "What is 2-hydroxyestrone?", "k": 6, "kind": null }
```
`kind` is optional — `"molecule"` or `"protein"` to restrict results (the server
over-fetches 5×k then filters).

Response:
```json
{
  "results": [
    { "rank": 1, "score": 0.83, "kind": "molecule",
      "paper": "2-hydroxyestrone", "section": "compound record",
      "chebi": "CHEBI:1156", "accession": null, "text": "Compound: 2-hydroxyestrone ..." }
  ],
  "model": "text-embedding-3-large",
  "took_ms": 42
}
```

Errors: `422` bad body · `502` embedding call failed · `503` still loading ·
`500` `EMBEDDING_MODEL` dim mismatch.

> If `index.faiss` was not L2-normalised at build time, `score` values come back
> `> 1`. Ranking is still correct (inner product is monotonic in cosine for
> fixed-norm query). The server logs the observed range on the first search.

## Smoke test

```bash
curl -s localhost:8000/health | jq
curl -s localhost:8000/search -H 'content-type: application/json' \
  -d '{"query":"What is 2-hydroxyestrone?","k":4}' | jq '.results[] | {rank,score,paper}'
curl -s localhost:8000/search -H 'content-type: application/json' \
  -d '{"query":"enzymes that hydroxylate estrone","k":4,"kind":"protein"}' | jq '.results[].kind'
```

## Deploy (lab server / persistent host — NOT Vercel)

Needs ~4 GB RAM and ~2.3 GB disk for the store. Run the Next app on the same host
so `RAG_SERVER_URL` stays loopback.

`/etc/systemd/system/atlas-rag.service`:

```ini
[Unit]
Description=Steroid Atlas RAG sidecar
After=network-online.target

[Service]
Type=simple
User=atlas
WorkingDirectory=/srv/steroid_atlas_web
Environment=LD_LIBRARY_PATH=/opt/miniconda3/lib
Environment=RAG_STORE_DIR=/srv/steroid_atlas_web/rag_data/rag_store
Environment=EMBEDDING_MODEL=text-embedding-3-large
EnvironmentFile=/srv/steroid_atlas_web/.env.rag   # OPENAI_API_KEY=...
ExecStart=/opt/miniconda3/bin/python3 /srv/steroid_atlas_web/scripts/rag/rag_server.py
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now atlas-rag
curl -s localhost:8000/health
```

If the sidecar ever binds a non-loopback address, add a shared-secret header
check to `/search` and firewall the port to the web host only.
