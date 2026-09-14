"use client";

import { useState } from "react";
import { useApiKey } from "@/lib/api-key-context";

/** Save / change / remove the server's OpenAI key. Used by ApiKeyGate on the
 *  home page and, compactly, inside ChatWidget when no key is configured. */
export default function ApiKeyForm({ compact = false }: { compact?: boolean }) {
  const { hasKey, hint, canEdit, saveKey, clearKey } = useApiKey();
  const [input, setInput] = useState("");
  const [editing, setEditing] = useState(false);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function save() {
    if (!input.trim() || busy) return;
    setBusy(true);
    setError(null);
    const err = await saveKey(input.trim());
    setBusy(false);
    if (err) {
      setError(err);
      return;
    }
    setInput("");
    setEditing(false);
  }

  async function remove() {
    setBusy(true);
    setError(null);
    const err = await clearKey();
    setBusy(false);
    if (err) setError(err);
  }

  const text = compact ? "text-xs" : "text-sm";

  if (!canEdit) {
    return (
      <div className={`${text} text-slate-600`}>
        {hasKey
          ? `🔑 Key configured by the site admin (${hint}).`
          : "🔑 No key configured. Ask the site admin to set OPENAI_API_KEY."}
      </div>
    );
  }

  if (hasKey && !editing) {
    return (
      <div className={`flex flex-wrap items-center gap-2 ${text}`}>
        <span className="text-slate-600">🔑 Key saved ({hint})</span>
        <button
          className="rounded border border-slate-300 bg-white px-2 py-0.5 text-slate-700 hover:bg-slate-100 disabled:opacity-50"
          onClick={() => setEditing(true)}
          disabled={busy}
        >
          change
        </button>
        <button
          className="rounded border border-slate-300 bg-white px-2 py-0.5 text-slate-700 hover:bg-slate-100 disabled:opacity-50"
          onClick={remove}
          disabled={busy}
        >
          remove
        </button>
        {error && <span className="text-red-600">{error}</span>}
      </div>
    );
  }

  return (
    <div>
      <div className="flex flex-wrap items-center gap-2">
        <input
          type="password"
          autoComplete="off"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && save()}
          placeholder="sk-..."
          disabled={busy}
          className={`min-w-[200px] flex-1 rounded-md border border-slate-300 bg-white px-3 py-1.5 font-mono shadow-sm focus:border-sky-500 focus:outline-none focus:ring-2 focus:ring-sky-100 ${text}`}
        />
        <button
          onClick={save}
          disabled={!input.trim() || busy}
          className={`rounded-md bg-black text-white disabled:opacity-40 ${
            compact ? "px-2 py-1.5 text-xs" : "px-3 py-1.5 text-sm"
          }`}
        >
          {busy ? "Saving…" : "Save"}
        </button>
        {hasKey && (
          <button
            className="text-xs text-slate-500 hover:text-slate-800"
            onClick={() => {
              setEditing(false);
              setInput("");
              setError(null);
            }}
          >
            cancel
          </button>
        )}
      </div>
      {error && <div className="mt-1 text-xs text-red-600">{error}</div>}
    </div>
  );
}
