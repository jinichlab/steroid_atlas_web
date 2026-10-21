#!/bin/sh
# Upload the web app to the Hugging Face Space, which rebuilds it from the
# Dockerfile. Needs `hf auth login` with write access to the Space.
#
#   sh scripts/deploy_space.sh                      # jinichlab/steroid-atlas
#   SPACE=you/your-space sh scripts/deploy_space.sh
#
# Only what the Docker build needs is uploaded — never .env*, node_modules or
# the RAG index (the Space downloads that from its private dataset at startup;
# see scripts/start_space.sh). The analysis side (data/, analysis/, …) stays in git only.
#
# The Space is made an exact mirror of what's uploaded: `--delete "*"` removes
# any remote file not in this upload (Hugging Face never deletes .gitattributes),
# so files from old layouts or a mistaken upload can't linger.
set -eu
cd "$(dirname "$0")/.."
SPACE="${SPACE:-jinichlab/steroid-atlas}"

# Refuse to run from anything but the steroid_atlas_web repo root.
for f in Dockerfile web/package.json rag/rag_server.py scripts/start_space.sh; do
  [ -e "$f" ] || { echo "deploy_space.sh: $f not found — not the steroid_atlas_web repo root" >&2; exit 1; }
done

hf upload "$SPACE" . . --repo-type space \
  --exclude "web/node_modules/*" "web/.next/*" "web/.env*" "web/*.tsbuildinfo" ".env*" \
            "rag_data/*" ".git/*" "data/*" "analysis/*" "docs/*" ".venv/*" \
  --delete "*" \
  --commit-message "${MSG:-Deploy from $(git rev-parse --short HEAD 2>/dev/null || echo local)}"
