"use client";

import { useState } from "react";
import { MessageCircle, Minus } from "lucide-react";
import { useSelection } from "@/lib/selection-context";
import { useApiKey } from "@/lib/api-key-context";
import ApiKeyForm from "@/components/ApiKeyForm";

type Source = { n: number; label: string };
type Message = {
  role: "user" | "assistant";
  content: string;
  sources?: Source[];
};

export default function ChatWidget() {
  const { selection } = useSelection();
  const { hasKey } = useApiKey();
  const [isOpen, setIsOpen] = useState(false);
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState("");
  const [isLoading, setIsLoading] = useState(false);

  async function sendMessage() {
    if (!input.trim() || isLoading || !hasKey) return;
    const query = input.trim();
    const history = messages;
    const userMsg: Message = { role: "user", content: query };
    setMessages((prev) => [...prev, userMsg]);
    setInput("");
    setIsLoading(true);

    try {
      const res = await fetch("/api/chat", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          query,
          history,
          // Only the chat-relevant fields — not the (potentially large)
          // cross-view bridge arrays.
          selection: selection.items.length
            ? {
                view: selection.view,
                origin: selection.origin,
                items: selection.items,
                total: selection.total,
              }
            : undefined,
        }),
      });
      const data = await res.json();
      const content = res.ok
        ? data.answer || "(empty response)"
        : `Error: ${data.error || res.statusText}`;
      const sources: Source[] | undefined = res.ok ? data.sources : undefined;
      setMessages((prev) => [...prev, { role: "assistant", content, sources }]);
    } catch (err) {
      const message = err instanceof Error ? err.message : "Request failed.";
      setMessages((prev) => [
        ...prev,
        { role: "assistant", content: `Error: ${message}` },
      ]);
    } finally {
      setIsLoading(false);
    }
  }

  if (!isOpen) {
    return (
      <button
        onClick={() => setIsOpen(true)}
        className="fixed bottom-6 right-6 z-50 rounded-full bg-black p-4 text-white shadow-lg hover:bg-neutral-800"
        aria-label="Open chat"
      >
        <MessageCircle size={24} />
      </button>
    );
  }

  return (
    <div className="fixed bottom-6 right-6 z-50 flex h-[500px] w-[350px] flex-col rounded-xl border bg-white shadow-2xl">
      <div className="flex items-center justify-between border-b p-3">
        <span className="font-medium">Ask the Atlas</span>
        <div className="flex items-center gap-2">
          {selection.items.length > 0 && (
            <span
              className="rounded-full bg-sky-100 px-2 py-0.5 text-xs font-medium text-sky-700"
              title={`${selection.total} ${selection.view ?? ""} item(s) selected on the map are shared with the chat`}
            >
              {selection.total} selected
            </span>
          )}
          <button onClick={() => setIsOpen(false)} aria-label="Minimize chat">
            <Minus size={18} />
          </button>
        </div>
      </div>

      {!hasKey && (
        <div className="border-b bg-amber-50 p-2 text-xs text-amber-800">
          <div className="mb-1">
            🔑 Chat is disabled — no OpenAI API key is configured.
          </div>
          <ApiKeyForm compact />
        </div>
      )}

      <div className="flex-1 space-y-2 overflow-y-auto p-3">
        {messages.map((m, i) => (
          <div
            key={i}
            className={m.role === "user" ? "text-right" : "text-left"}
          >
            <span
              className={`inline-block whitespace-pre-wrap rounded-lg px-3 py-2 text-sm ${
                m.role === "user" ? "bg-black text-white" : "bg-neutral-100"
              }`}
            >
              {m.content}
            </span>
            {m.sources && m.sources.length > 0 && (
              <ul className="mt-1 space-y-0.5 text-xs text-neutral-500">
                {m.sources.map((s) => (
                  <li key={s.n}>
                    [{s.n}] {s.label}
                  </li>
                ))}
              </ul>
            )}
          </div>
        ))}
        {isLoading && (
          <div className="text-left">
            <span className="inline-block rounded-lg bg-neutral-100 px-3 py-2 text-sm text-neutral-500">
              Thinking…
            </span>
          </div>
        )}
      </div>

      <div className="flex items-center gap-2 border-t p-3">
        <input
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && sendMessage()}
          placeholder={
            hasKey
              ? "Ask about a compound or protein..."
              : "Add an API key to start chatting"
          }
          disabled={!hasKey}
          className="flex-1 rounded-md border px-3 py-2 text-sm disabled:bg-neutral-100"
        />
        <button
          onClick={sendMessage}
          disabled={isLoading || !hasKey}
          className="rounded-md bg-black px-3 py-2 text-sm text-white disabled:opacity-50"
        >
          Send
        </button>
      </div>
    </div>
  );
}
