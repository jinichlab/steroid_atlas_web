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
 * Tracks whether the *server* has an OpenAI API key configured (in .env.local).
 *
 * The key value itself never reaches the browser — only a boolean and a
 * last-4-characters hint. Saving a key POSTs it to /api/config, which persists
 * it to .env.local so it survives a restart.
 */

export type KeyStatus = {
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
};

const ApiKeyContext = createContext<ApiKeyContextValue | null>(null);

const EMPTY: KeyStatus = { hasKey: false, hint: null, canEdit: true };

export function ApiKeyProvider({ children }: { children: ReactNode }) {
  const [status, setStatus] = useState<KeyStatus>(EMPTY);
  const [ready, setReady] = useState(false);

  useEffect(() => {
    let alive = true;
    fetch("/api/config", { cache: "no-store" })
      .then((r) => r.json())
      .then((d: KeyStatus) => alive && setStatus(d))
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
        setStatus(data as KeyStatus);
        return null;
      } catch (err) {
        return err instanceof Error ? err.message : "Request failed.";
      }
    },
    [],
  );

  const saveKey = useCallback((key: string) => send("POST", key), [send]);
  const clearKey = useCallback(() => send("DELETE"), [send]);

  return (
    <ApiKeyContext.Provider value={{ ...status, ready, saveKey, clearKey }}>
      {children}
    </ApiKeyContext.Provider>
  );
}

/** Safe outside a provider — reports "no key, still loading nothing". */
export function useApiKey(): ApiKeyContextValue {
  return (
    useContext(ApiKeyContext) ?? {
      ...EMPTY,
      ready: true,
      saveKey: async () => "No provider.",
      clearKey: async () => "No provider.",
    }
  );
}
