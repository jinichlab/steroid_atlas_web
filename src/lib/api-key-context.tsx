"use client";

import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useState,
  type ReactNode,
} from "react";

/**
 * Tracks the OpenAI key that powers "Ask the Atlas", in one of two modes
 * (decided by the server, see src/lib/key-mode.ts):
 *
 *  - "server": the key lives on the server (.env.local). The browser only
 *    ever sees a boolean and a last-4-characters hint; saving a key POSTs it
 *    to /api/config, which persists it to .env.local.
 *  - "user": each visitor uses their own key. It is kept in this tab's
 *    sessionStorage (gone when the tab closes), never sent to /api/config, and
 *    attached to each chat request via `chatHeaders()`.
 */

export type KeyMode = "server" | "user";

export type KeyStatus = {
  mode: KeyMode;
  hasKey: boolean;
  hint: string | null;
  canEdit: boolean;
};

type ApiKeyContextValue = KeyStatus & {
  /** False until the first /api/config response lands. */
  ready: boolean;
  /** Returns null on success, or an error message. */
  saveKey: (key: string) => Promise<string | null>;
  clearKey: () => Promise<string | null>;
  /** Extra headers for /api/chat — carries the visitor's key in "user" mode. */
  chatHeaders: () => Record<string, string>;
};

const ApiKeyContext = createContext<ApiKeyContextValue | null>(null);

const EMPTY: KeyStatus = { mode: "server", hasKey: false, hint: null, canEdit: true };
const STORAGE_KEY = "steroid-atlas:openai-key";
const USER_KEY_HEADER = "x-openai-key"; // mirrors src/lib/key-mode.ts

// sessionStorage can throw (privacy modes, blocked storage) — degrade to
// "no remembered key" rather than breaking the page.
function readStored(): string {
  try {
    return sessionStorage.getItem(STORAGE_KEY) ?? "";
  } catch {
    return "";
  }
}
function writeStored(key: string | null) {
  try {
    if (key) sessionStorage.setItem(STORAGE_KEY, key);
    else sessionStorage.removeItem(STORAGE_KEY);
  } catch {
    // key still works for this page view via state
  }
}

const hintFor = (key: string) => (key ? `····${key.slice(-4)}` : null);

export function ApiKeyProvider({ children }: { children: ReactNode }) {
  const [status, setStatus] = useState<KeyStatus>(EMPTY);
  const [ready, setReady] = useState(false);
  const [userKey, setUserKey] = useState("");

  useEffect(() => {
    let alive = true;
    fetch("/api/config", { cache: "no-store" })
      .then((r) => r.json())
      .then((d: Partial<KeyStatus>) => {
        if (!alive) return;
        const mode: KeyMode = d.mode === "user" ? "user" : "server";
        setStatus({ ...EMPTY, ...d, mode });
        if (mode === "user") setUserKey(readStored());
      })
      .catch(() => {})
      .finally(() => alive && setReady(true));
    return () => {
      alive = false;
    };
  }, []);

  const send = useCallback(
    async (method: "POST" | "DELETE", key?: string): Promise<string | null> => {
      try {
        const res = await fetch("/api/config", {
          method,
          headers: { "Content-Type": "application/json" },
          body: key === undefined ? undefined : JSON.stringify({ apiKey: key }),
        });
        const data = await res.json();
        if (!res.ok) return data.error || res.statusText;
        setStatus((s) => ({ ...s, ...(data as KeyStatus) }));
        return null;
      } catch (err) {
        return err instanceof Error ? err.message : "Request failed.";
      }
    },
    [],
  );

  const isUser = status.mode === "user";

  const saveKey = useCallback(
    async (key: string) => {
      if (!isUser) return send("POST", key);
      if (!key.startsWith("sk-")) {
        return "That doesn't look like an OpenAI key — it should start with 'sk-'.";
      }
      setUserKey(key);
      writeStored(key);
      return null;
    },
    [isUser, send],
  );

  const clearKey = useCallback(async () => {
    if (!isUser) return send("DELETE");
    setUserKey("");
    writeStored(null);
    return null;
  }, [isUser, send]);

  const chatHeaders = useCallback(
    (): Record<string, string> =>
      isUser && userKey ? { [USER_KEY_HEADER]: userKey } : {},
    [isUser, userKey],
  );

  const value: ApiKeyContextValue = isUser
    ? {
        mode: "user",
        hasKey: !!userKey,
        hint: hintFor(userKey),
        canEdit: true,
        ready,
        saveKey,
        clearKey,
        chatHeaders,
      }
    : { ...status, ready, saveKey, clearKey, chatHeaders };

  return <ApiKeyContext.Provider value={value}>{children}</ApiKeyContext.Provider>;
}

/** Safe outside a provider — reports "no key, still loading nothing". */
export function useApiKey(): ApiKeyContextValue {
  return (
    useContext(ApiKeyContext) ?? {
      ...EMPTY,
      ready: true,
      saveKey: async () => "No provider.",
      clearKey: async () => "No provider.",
      chatHeaders: () => ({}),
    }
  );
}
