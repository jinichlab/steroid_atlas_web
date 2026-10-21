#!/bin/sh
# Hugging Face Space entrypoint: Next.js on :7860, RAG sidecar on 127.0.0.1:8000.
# The web app starts immediately; the 2.1 GB index downloads and loads in the
# background (chat answers ungrounded until it's ready).
if [ -d /data ] && [ -w /data ]; then STORE=/data/rag_store; else STORE=/tmp/rag_store; fi
export RAG_STORE_DIR="$STORE"
(
  python3 -c "from huggingface_hub import snapshot_download; snapshot_download('${RAG_DATASET:-akuppal/steroid-atlas-rag}', repo_type='dataset', local_dir='$STORE')" \
  && exec python3 rag/rag_server.py
) &
exec node server.js
