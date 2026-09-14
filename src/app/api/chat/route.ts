import { NextResponse } from "next/server";
import OpenAI from "openai";
import { retrieve, formatContext, sourceLabel } from "@/lib/rag";

export const runtime = "nodejs";

type Message = { role: "user" | "assistant"; content: string };

type SelectionItem = {
  kind?: string;
  name?: string;
  chebi?: string | null;
  smiles?: string | null;
  accession?: string | null;
  gene?: string | null;
  organism?: string | null;
  cluster?: number | string | null;
};

type Selection = {
  view?: string | null;
  origin?: string | null;
  items?: SelectionItem[];
  total?: number;
};

const SYSTEM_PROMPT =
  "You are 'Ask the Atlas', an assistant for the Steroid Atlas web app — an interactive UMAP atlas of steroid-metabolizing enzymes and their small-molecule substrates. " +
  "Answer questions about steroid compounds, metabolizing proteins/enzymes, and how they relate. Be concise and precise. If you are unsure, say so.";

const GROUNDING_INSTRUCTIONS = [
  "Answer using ONLY the numbered CONTEXT passages below. Each is a chunk from the",
  "Steroid Atlas knowledge base (compound records with ChEBI definitions and chemical",
  "properties, or protein records from UniProt literature).",
  '- Cite every claim with the bracket number(s) it came from, e.g. "Estrone is a 17-oxo steroid [1][3]."',
  "- If the CONTEXT does not contain the answer, say so plainly and do not guess.",
  '- Prefer ChEBI definitions / chemical properties for "what is" questions; prefer protein records for enzyme/substrate questions.',
  "- Keep answers concise. Do not dump raw passage text.",
].join("\n");

const DEGRADED_NOTE =
  "The knowledge base is temporarily unavailable, so answer from general knowledge and state that this response is not grounded in the Steroid Atlas corpus.";

const NO_CONTEXT_NOTE =
  "No relevant passages were found in the knowledge base for this question. Say the Atlas corpus has nothing specific on it, then optionally give a brief general answer flagged as such.";

const MAX_SELECTION_ITEMS = 80;

/** Human-readable block describing what the user has selected on the map. */
function formatSelection(sel: Selection): string {
  const items = (sel.items ?? []).slice(0, MAX_SELECTION_ITEMS);
  if (!items.length) return "";

  const originText =
    sel.origin === "lasso"
      ? "lasso-selected a region of"
      : sel.origin === "search"
        ? "filtered"
        : sel.origin === "cluster"
          ? "highlighted a cluster on"
          : "selected items on";
  const view = sel.view ? `the ${sel.view} map` : "the map";
  const total = sel.total ?? items.length;
  const shown =
    total > items.length ? ` (showing the first ${items.length})` : "";

  const lines = items.map((it) => {
    if (it.kind === "protein") {
      const bits = [it.accession, it.name, it.gene && `gene ${it.gene}`, it.organism]
        .filter(Boolean)
        .join(" · ");
      return `- ${bits}`;
    }
    const bits = [it.name, it.chebi, it.smiles && `SMILES ${it.smiles}`]
      .filter(Boolean)
      .join(" · ");
    return `- ${bits}`;
  });

  return [
    `The user has ${originText} ${view}: ${total} item(s)${shown}.`,
    "Treat this as the subject of questions like \"these\", \"what I selected\", \"this set\".",
    "SELECTION:",
    ...lines,
  ].join("\n");
}

export async function POST(req: Request) {
  // The key lives only on the server (.env.local, set via /api/config or by
  // hand). It is never accepted from the request body.
  const apiKey = process.env.OPENAI_API_KEY?.trim();
  if (!apiKey) {
    return NextResponse.json(
      {
        error:
          "Chat is disabled: no OpenAI API key is configured. Add one on the home page.",
      },
      { status: 503 },
    );
  }

  let body: { query?: string; history?: Message[]; selection?: Selection };
  try {
    body = await req.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON body." }, { status: 400 });
  }

  const query = body.query?.trim();
  if (!query) {
    return NextResponse.json({ error: "Missing 'query'." }, { status: 400 });
  }

  const history = (body.history ?? [])
    .filter((m) => (m.role === "user" || m.role === "assistant") && typeof m.content === "string")
    .slice(-10);

  const selectionBlock = body.selection ? formatSelection(body.selection) : "";
  const selectionNames = (body.selection?.items ?? [])
    .slice(0, 20)
    .map((it) => it.name)
    .filter(Boolean)
    .join(", ");

  // ── Retrieval ──────────────────────────────────────────────────────────
  // Fold the selected item names into the retrieval query so passages about
  // them surface even when the user just says "these".
  const retrievalQuery = selectionNames
    ? `${query}\nRelevant items: ${selectionNames}`
    : query;

  const ragEnabled = process.env.RAG_ENABLED !== "false";
  const { hits, ok } = ragEnabled
    ? await retrieve(retrievalQuery, {
        k: Number(process.env.RAG_TOP_K) || 6,
        // Forward the key so the sidecar picks up a newly-saved one without
        // needing its own env or a restart.
        apiKey,
      })
    : { hits: [], ok: true };

  let systemContent: string;
  if (hits.length > 0) {
    systemContent = `${SYSTEM_PROMPT}\n\n${GROUNDING_INSTRUCTIONS}\n\nCONTEXT:\n${formatContext(hits)}`;
  } else if (!ok) {
    console.warn("[rag] sidecar unavailable, answering without retrieval");
    systemContent = `${SYSTEM_PROMPT}\n\n${DEGRADED_NOTE}`;
  } else {
    systemContent = `${SYSTEM_PROMPT}\n\n${NO_CONTEXT_NOTE}`;
  }

  if (selectionBlock) {
    systemContent = `${systemContent}\n\n${selectionBlock}`;
  }

  const sources = hits.map((h) => ({
    n: h.rank,
    label: sourceLabel(h),
    kind: h.kind,
    chebi: h.chebi,
    accession: h.accession,
    score: Number(h.score.toFixed(3)),
  }));

  const openai = new OpenAI({ apiKey });

  try {
    const completion = await openai.chat.completions.create({
      model: process.env.OPENAI_MODEL || "gpt-4o-mini",
      temperature: 0.2,
      messages: [
        { role: "system", content: systemContent },
        ...history,
        { role: "user", content: query },
      ],
    });

    const answer = completion.choices[0]?.message?.content?.trim() ?? "";
    return NextResponse.json({ answer, sources });
  } catch (err) {
    const message = err instanceof Error ? err.message : "Unknown error calling OpenAI.";
    return NextResponse.json({ error: message }, { status: 502 });
  }
}
