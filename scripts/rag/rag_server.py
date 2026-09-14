"""Steroid Atlas RAG retrieval sidecar.

Loads the prebuilt FAISS index + chunk catalog once, then answers similarity
searches over HTTP. The Next.js `/api/chat` route calls `POST /search` to pull
grounding passages before it prompts the chat model.

Stdlib-only (http.server) on purpose — the only third-party imports are faiss,
numpy and openai, which the project's conda env already has. No FastAPI/uvicorn
to install.

Run:
    LD_LIBRARY_PATH=~/miniconda3/lib \
    RAG_STORE_DIR=$PWD/rag_data/rag_store \
    OPENAI_API_KEY=sk-... \
    ~/miniconda3/bin/python3 scripts/rag/rag_server.py

See scripts/rag/README.md for the full contract and a systemd unit.
"""
from __future__ import annotations

import json
import os
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import faiss
import numpy as np
import openai

# ── Config ───────────────────────────────────────────────────────────────
_REPO_ROOT = Path(__file__).resolve().parents[2]
RAG_STORE_DIR = Path(os.getenv("RAG_STORE_DIR", _REPO_ROOT / "rag_data" / "rag_store"))
INDEX_PATH = RAG_STORE_DIR / "index.faiss"
CATALOG_PATH = RAG_STORE_DIR / "catalog.jsonl"

EMBEDDING_MODEL = os.getenv("EMBEDDING_MODEL", "text-embedding-3-large")
FAISS_THREADS = int(os.getenv("RAG_FAISS_THREADS", "4"))
HOST = os.getenv("RAG_SERVER_HOST", "127.0.0.1")
PORT = int(os.getenv("RAG_SERVER_PORT", "8000"))

_VALID_KINDS = {"molecule", "protein"}

# ── State (populated by load()) ──────────────────────────────────────────
_index: faiss.Index | None = None
_catalog: list[str] = []          # raw JSONL lines, index i ↔ vector i
_default_client: openai.OpenAI | None = None  # from server env, if any
_scores_logged = False


def load() -> None:
    global _index, _catalog, _default_client

    if not INDEX_PATH.exists():
        raise FileNotFoundError(f"FAISS index not found: {INDEX_PATH}")
    if not CATALOG_PATH.exists():
        raise FileNotFoundError(f"Catalog not found: {CATALOG_PATH}")

    faiss.omp_set_num_threads(FAISS_THREADS)

    t0 = time.time()
    print(f"[rag] reading index {INDEX_PATH} ...", flush=True)
    _index = faiss.read_index(str(INDEX_PATH))
    print(f"[rag] reading catalog {CATALOG_PATH} ...", flush=True)
    with CATALOG_PATH.open("r", encoding="utf-8") as fh:
        _catalog = fh.readlines()

    if _index.ntotal != len(_catalog):
        raise RuntimeError(
            f"index/catalog drift: index.ntotal={_index.ntotal} "
            f"catalog lines={len(_catalog)} — rebuild the store"
        )

    # A server-wide key is optional: each /search request may instead carry
    # its own `api_key` (forwarded from the visitor's browser). The index
    # itself needs no OpenAI access to load.
    if os.getenv("OPENAI_API_KEY"):
        _default_client = openai.OpenAI()
        key_note = "server key configured"
    else:
        _default_client = None
        key_note = "no server key — each request must supply api_key"

    print(
        f"[rag] loaded index: ntotal={_index.ntotal} dim={_index.d} "
        f"model={EMBEDDING_MODEL} threads={FAISS_THREADS} ({key_note}) "
        f"({time.time() - t0:.1f}s)",
        flush=True,
    )


class MissingApiKey(Exception):
    """Neither the request nor the sidecar's own env has an OpenAI key."""


def _client_for(api_key: str | None) -> openai.OpenAI:
    if api_key:
        return openai.OpenAI(api_key=api_key)
    if _default_client is not None:
        return _default_client
    raise MissingApiKey(
        "no OpenAI API key available (pass api_key in the request, or set "
        "OPENAI_API_KEY on the sidecar)"
    )


def _embed(query: str, api_key: str | None) -> np.ndarray:
    client = _client_for(api_key)
    emb = (
        client.embeddings.create(model=EMBEDDING_MODEL, input=[query])
        .data[0]
        .embedding
    )
    if _index is not None and len(emb) != _index.d:
        raise ValueError(
            f"embedding dim {len(emb)} != index dim {_index.d}; "
            f"EMBEDDING_MODEL={EMBEDDING_MODEL} is wrong for this index"
        )
    v = np.asarray(emb, dtype="float32")[None, :]
    faiss.normalize_L2(v)
    return v


def do_search(query: str, k: int, kind: str | None, api_key: str | None) -> dict:
    global _scores_logged
    assert _index is not None

    t0 = time.time()
    v = _embed(query, api_key)

    over = k * 5 if kind else k
    over = min(over, _index.ntotal)
    distances, ids = _index.search(v, over)

    if not _scores_logged:
        row = distances[0]
        print(
            f"[rag] first search score range: min={row.min():.3f} max={row.max():.3f} "
            f"(>1 means the index was not L2-normalised at build time; ranking is "
            f"unaffected)",
            flush=True,
        )
        _scores_logged = True

    results: list[dict] = []
    for score, i in zip(distances[0], ids[0]):
        if i < 0:
            continue
        row = json.loads(_catalog[i])
        if kind and row.get("kind") != kind:
            continue
        results.append(
            {
                "rank": len(results) + 1,
                "score": float(score),
                "kind": row.get("kind", ""),
                "paper": row.get("paper") or "",
                "section": row.get("section") or "",
                "chebi": row.get("chebi"),
                "accession": row.get("accession"),
                "text": row.get("text") or "",
            }
        )
        if len(results) >= k:
            break

    return {
        "results": results,
        "model": EMBEDDING_MODEL,
        "took_ms": int((time.time() - t0) * 1000),
    }


# ── HTTP ─────────────────────────────────────────────────────────────────
class Handler(BaseHTTPRequestHandler):
    # HTTP/1.0 + explicit close: one request per connection. Avoids keep-alive
    # races between the stdlib server and undici's connection pool (the Next.js
    # route's fetch), which showed up as spurious "operation aborted" timeouts.
    protocol_version = "HTTP/1.0"

    def _send(self, status: int, payload: dict) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.close_connection = True
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Connection", "close")
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt: str, *args) -> None:  # quieter default logging
        print(f"[rag] {self.address_string()} {fmt % args}", flush=True)

    def do_GET(self) -> None:
        if self.path.rstrip("/") == "/health":
            if _index is None:
                self._send(503, {"status": "loading"})
            else:
                self._send(
                    200,
                    {
                        "status": "ok",
                        "ntotal": _index.ntotal,
                        "dim": _index.d,
                        "model": EMBEDDING_MODEL,
                    },
                )
            return
        self._send(404, {"error": "not found"})

    def do_POST(self) -> None:
        if self.path.rstrip("/") != "/search":
            self._send(404, {"error": "not found"})
            return
        if _index is None:
            self._send(503, {"error": "index still loading"})
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
            body = json.loads(self.rfile.read(length) or b"{}")
        except (ValueError, json.JSONDecodeError):
            self._send(422, {"error": "invalid JSON body"})
            return

        query = (body.get("query") or "").strip()
        if not query:
            self._send(422, {"error": "missing 'query'"})
            return

        try:
            k = int(body.get("k", 6))
        except (TypeError, ValueError):
            k = 6
        k = max(1, min(k, 50))

        kind = body.get("kind")
        if kind is not None and kind not in _VALID_KINDS:
            self._send(422, {"error": "kind must be 'molecule', 'protein' or null"})
            return

        # Per-request key (forwarded from the visitor's browser via Next.js).
        # Falls back to the sidecar's own OPENAI_API_KEY when omitted.
        api_key = (body.get("api_key") or "").strip() or None

        try:
            self._send(200, do_search(query, k, kind, api_key))
        except MissingApiKey as exc:
            self._send(401, {"error": str(exc)})
        except ValueError as exc:  # dim mismatch
            self._send(500, {"error": str(exc)})
        except Exception as exc:  # noqa: BLE001 — embedding / search failure
            print(f"[rag] search failed: {exc}", flush=True)
            self._send(502, {"error": f"search failed: {exc}"})


def main() -> None:
    load()
    server = ThreadingHTTPServer((HOST, PORT), Handler)
    print(f"[rag] serving on http://{HOST}:{PORT}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[rag] shutting down", flush=True)
        server.shutdown()


if __name__ == "__main__":
    main()
