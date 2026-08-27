/**
 * HTTP boundary to the RAG retrieval sidecar (scripts/rag/rag_server.py).
 *
 * `retrieve()` never throws — on any failure (sidecar down, timeout, bad
 * response) it returns `{ hits: [], ok: false }` so the chat route can fall
 * back to an ungrounded answer. Swapping in a hosted vector DB later means
 * reimplementing this file and nothing else.
 */

export type RagKind = "molecule" | "protein";

export type RagHit = {
  rank: number;
  score: number;
  kind: RagKind;
  paper: string;
  section: string;
  chebi: string | null;
  accession: string | null;
  text: string;
};

export type RagResult = { hits: RagHit[]; ok: boolean };

const SERVER_URL = process.env.RAG_SERVER_URL || "http://127.0.0.1:8000";
const TIMEOUT_MS = Number(process.env.RAG_TIMEOUT_MS) || 8000;
const MAX_CHUNK_CHARS = Number(process.env.RAG_MAX_CHUNK_CHARS) || 1600;
const MAX_CONTEXT_CHARS = Number(process.env.RAG_MAX_CONTEXT_CHARS) || 12000;

type SearchResponse = {
  results?: Array<Partial<RagHit>>;
};

export async function retrieve(
  query: string,
  opts?: { k?: number; kind?: RagKind },
): Promise<RagResult> {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), TIMEOUT_MS);

  try {
    const res = await fetch(`${SERVER_URL}/search`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        query,
        k: opts?.k ?? 6,
        kind: opts?.kind ?? null,
      }),
      signal: controller.signal,
      cache: "no-store",
    });

    if (!res.ok) {
      console.warn(`[rag] sidecar responded ${res.status}`);
      return { hits: [], ok: false };
    }

    const data = (await res.json()) as SearchResponse;
    const hits: RagHit[] = (data.results ?? []).map((r, i) => ({
      rank: r.rank ?? i + 1,
      score: r.score ?? 0,
      kind: (r.kind as RagKind) ?? "molecule",
      paper: r.paper ?? "",
      section: r.section ?? "",
      chebi: r.chebi ?? null,
      accession: r.accession ?? null,
      text: r.text ?? "",
    }));

    return { hits, ok: true };
  } catch (err) {
    const reason = err instanceof Error ? err.message : String(err);
    console.warn(`[rag] retrieve failed: ${reason}`);
    return { hits: [], ok: false };
  } finally {
    clearTimeout(timer);
  }
}

/** "2-hydroxyestrone (CHEBI:1156)" for molecules, "UniProt W8DLN2" for proteins. */
export function sourceLabel(hit: RagHit): string {
  if (hit.kind === "protein" && hit.accession) {
    return `UniProt ${hit.accession}`;
  }
  const name = hit.paper || hit.chebi || "unknown";
  return hit.chebi ? `${name} (${hit.chebi})` : name;
}

/** Numbered context block for the system prompt. Truncates per-chunk and overall. */
export function formatContext(hits: RagHit[]): string {
  const blocks: string[] = [];
  let total = 0;

  for (const hit of hits) {
    const header =
      hit.kind === "protein"
        ? `[${hit.rank}] protein · ${sourceLabel(hit)} · ${hit.section}`
        : `[${hit.rank}] molecule · ${sourceLabel(hit)} · ${hit.section}`;

    let body = hit.text.replace(/\s+/g, " ").trim();
    if (body.length > MAX_CHUNK_CHARS) {
      body = `${body.slice(0, MAX_CHUNK_CHARS)}…`;
    }

    const block = `${header}\n${body}`;
    if (total + block.length > MAX_CONTEXT_CHARS) break;
    blocks.push(block);
    total += block.length + 2;
  }

  return blocks.join("\n\n");
}
