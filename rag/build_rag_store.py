"""
build_rag_store.py — turn the compound + protein corpus into a searchable FAISS index.

Reads:
  · rag_data/rag_train/small_molecule/data_txt/*.txt   — small_molecule_scraper.py
  · rag_data/rag_train/small_molecule/_status.csv
  · rag_data/rag_train/proteins/data_txt/*.txt         — protein_scraper.py (.txt files)

Writes:
  · rag_data/rag_store/index.faiss     — the vectors
  · rag_data/rag_store/catalog.jsonl   — the text each vector came from

Vector i in the index corresponds to line i in the catalog. That holds only
because IndexFlatIP preserves insertion order — if you ever switch to IVF or
HNSW, add explicit IDs with IndexIDMap or your citations will point at the
wrong document.

    python rag/build_rag_store.py
    python rag/build_rag_store.py --rebuild     # ignore existing index
"""

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import faiss
import numpy as np
from dotenv import load_dotenv
from openai import OpenAI

ROOT = Path(__file__).resolve().parents[1]  # repo root (rag/ → ..)
MOLECULE_DIR = ROOT / "rag_data" / "rag_train" / "small_molecule" / "data_txt"
PROTEIN_DIR = ROOT / "rag_data" / "rag_train" / "proteins" / "data_txt"
STATUS_CSV = ROOT / "rag_data" / "rag_train" / "small_molecule" / "_status.csv"
STORE = ROOT / "rag_data" / "rag_store"

MODEL = "text-embedding-3-large"
DIM = 3072                      # must match MODEL
CHUNK = 1200                    # characters per chunk
OVERLAP = 200                   # characters shared between neighbours
MIN_CHUNK = 200                 # drop headers, footers, stub lines
BATCH = 64                      # texts per embedding request

load_dotenv(ROOT / ".env.local")


# ------------------------- chunking -------------------------

def chunk(text):
    text = " ".join(text.split())
    step = CHUNK - OVERLAP
    return [text[i:i + CHUNK] for i in range(0, len(text), step)]


# ------------------------- collection -----------------------

def load_names():
    """ChEBI slug -> readable compound name, so citations aren't bare IDs."""
    names = {}
    status = STATUS_CSV
    if not status.exists():
        return names
    with open(status, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            slug = (row.get("chebi_id") or "").replace(":", "_")
            if slug and row.get("name"):
                names[slug] = row["name"]
    return names


def collect_compounds():
    if not MOLECULE_DIR.exists():
        sys.exit(f"No {MOLECULE_DIR}. Run scripts/build_corpus.py first.")
    if not PROTEIN_DIR.exists():
        sys.exit(f"No {PROTEIN_DIR}. Run scripts/protein_scraper.py first.")

    names = load_names()
    records = []
    tiny = 0

    sources = [
        (MOLECULE_DIR, "molecule"),
        (PROTEIN_DIR, "protein"),
    ]

    total_files = 0
    for directory, kind in sources:
        files = sorted(directory.glob("*.txt"))
        total_files += len(files)
        for txt in files:
            body = txt.read_text(encoding="utf-8")
            if len(body) < MIN_CHUNK:
                tiny += 1
                continue
            label = names.get(txt.stem, txt.stem) if kind == "molecule" else txt.stem
            for piece in chunk(body):
                if len(piece) < MIN_CHUNK:
                    continue
                records.append({
                    "paper": label,
                    "section": "compound record" if kind == "molecule" else "protein record",
                    "chebi": txt.stem.replace("_", ":") if kind == "molecule" else None,
                    "accession": txt.stem if kind == "protein" else None,
                    "kind": kind,
                    "page": None,
                    "text": piece,
                })

    print(f"  {total_files} files -> {len(records)} chunks"
          + (f" ({tiny} files too short to index)" if tiny else ""))
    return records


# ------------------------- embedding ------------------------

def embed(texts, client):
    vectors = []
    for i in range(0, len(texts), BATCH):
        batch = texts[i:i + BATCH]
        for attempt in range(4):
            try:
                resp = client.embeddings.create(model=MODEL, input=batch)
                vectors += [d.embedding for d in resp.data]
                break
            except Exception as e:
                if attempt == 3:
                    raise
                wait = 2 ** attempt
                print(f"\n    retry in {wait}s ({e})")
                time.sleep(wait)
        print(f"  embedded {min(i + BATCH, len(texts))}/{len(texts)}",
              end="\r", flush=True)
    print()
    arr = np.array(vectors, dtype="float32")
    # Normalize so inner product == cosine similarity.
    arr /= np.linalg.norm(arr, axis=1, keepdims=True) + 1e-12
    return arr


# ------------------------- main -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true",
                    help="overwrite an existing index")
    args = ap.parse_args()

    key = os.getenv("OPENAI_API_KEY")
    if not key:
        sys.exit(f"No OPENAI_API_KEY. Expected it in {ROOT / '.env'}")

    if (STORE / "index.faiss").exists() and not args.rebuild:
        sys.exit(f"{STORE / 'index.faiss'} already exists. "
                 f"Use --rebuild to replace it.")

    print(f"Reading {MOLECULE_DIR} and {PROTEIN_DIR} ...")
    records = collect_compounds()
    if not records:
        sys.exit("Nothing to index — every file was empty or too short. "
                 "Check the scraper output before spending on embeddings.")

    total_chars = sum(len(r["text"]) for r in records)
    est_cost = (total_chars / 4) / 1_000_000 * 0.13    # ~4 chars/token
    median = sorted(len(r["text"]) for r in records)[len(records) // 2]
    print(f"\n{len(records):,} chunks · {total_chars:,} chars · "
          f"median chunk {median:,} · est. cost ${est_cost:.2f}")
    if input("proceed? [y/N] ").strip().lower() != "y":
        sys.exit("aborted")

    client = OpenAI(api_key=key)
    print("\nEmbedding...")
    vectors = embed([r["text"] for r in records], client)

    if vectors.shape[1] != DIM:
        sys.exit(f"Model returned {vectors.shape[1]} dims, expected {DIM}. "
                 f"Update DIM at the top of this file.")

    STORE.mkdir(parents=True, exist_ok=True)
    index = faiss.IndexFlatIP(DIM)
    index.add(vectors)
    faiss.write_index(index, str(STORE / "index.faiss"))

    with open(STORE / "catalog.jsonl", "w", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    print(f"\nwrote {index.ntotal:,} vectors -> {STORE}")
    print("Restart the RAG service (npm run dev:rag) to pick it up.")


if __name__ == "__main__":
    main()