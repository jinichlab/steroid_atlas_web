FROM node:20-bookworm-slim AS build
WORKDIR /app
COPY package.json package-lock.json ./
RUN npm ci
COPY . .
RUN npm run build

FROM node:20-bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH
COPY scripts/rag/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt huggingface_hub
WORKDIR /app
COPY --from=build --chown=node /app/.next/standalone ./
COPY --from=build --chown=node /app/.next/static ./.next/static
COPY --from=build --chown=node /app/public ./public
COPY --chown=node scripts/rag ./scripts/rag
COPY --chown=node start.sh ./start.sh
RUN chown node /app && chmod +x start.sh
USER node
ENV NODE_ENV=production PORT=7860 HOSTNAME=0.0.0.0 HF_HOME=/tmp/hf
EXPOSE 7860
CMD ["./start.sh"]
