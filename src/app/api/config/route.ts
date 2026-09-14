import { NextResponse } from "next/server";
import { writeEnvVar, isSafeEnvValue } from "@/lib/env-file";

export const runtime = "nodejs";
export const dynamic = "force-dynamic";

const KEY = "OPENAI_API_KEY";

/** Whether saving a key from the UI is allowed. Set ALLOW_KEY_SETUP=false on a
 *  publicly-reachable deployment — otherwise any visitor could replace the key. */
const setupAllowed = process.env.ALLOW_KEY_SETUP !== "false";

function status() {
  const key = process.env[KEY]?.trim();
  return {
    hasKey: !!key,
    // Enough to recognise which key is configured, never enough to use it.
    hint: key ? `····${key.slice(-4)}` : null,
    canEdit: setupAllowed,
  };
}

/** Is a key configured? The key itself is never sent to the browser. */
export async function GET() {
  return NextResponse.json(status());
}

/** Save a key to .env.local (and to the running process). */
export async function POST(req: Request) {
  if (!setupAllowed) {
    return NextResponse.json(
      { error: "Key setup is disabled on this deployment (ALLOW_KEY_SETUP=false)." },
      { status: 403 },
    );
  }

  let body: { apiKey?: string };
  try {
    body = await req.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON body." }, { status: 400 });
  }

  const apiKey = body.apiKey?.trim() ?? "";
  if (!isSafeEnvValue(apiKey)) {
    return NextResponse.json(
      { error: "Provide a single-line, non-empty API key." },
      { status: 400 },
    );
  }
  if (!apiKey.startsWith("sk-")) {
    return NextResponse.json(
      { error: "That doesn't look like an OpenAI key — it should start with 'sk-'." },
      { status: 400 },
    );
  }

  try {
    await writeEnvVar(KEY, apiKey);
  } catch (err) {
    const message = err instanceof Error ? err.message : "Unknown error.";
    return NextResponse.json(
      { error: `Could not write .env.local: ${message}` },
      { status: 500 },
    );
  }

  return NextResponse.json(status());
}

/** Remove the saved key. */
export async function DELETE() {
  if (!setupAllowed) {
    return NextResponse.json(
      { error: "Key setup is disabled on this deployment (ALLOW_KEY_SETUP=false)." },
      { status: 403 },
    );
  }
  try {
    await writeEnvVar(KEY, null);
  } catch (err) {
    const message = err instanceof Error ? err.message : "Unknown error.";
    return NextResponse.json(
      { error: `Could not write .env.local: ${message}` },
      { status: 500 },
    );
  }
  return NextResponse.json(status());
}
