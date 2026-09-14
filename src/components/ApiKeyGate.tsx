"use client";

import ApiKeyForm from "@/components/ApiKeyForm";
import { useApiKey } from "@/lib/api-key-context";

/** Home-page box for the OpenAI key that powers "Ask the Atlas". The atlas
 *  itself needs no key, so this never blocks the rest of the page. */
export default function ApiKeyGate() {
  const { hasKey, ready } = useApiKey();
  if (!ready) return null; // avoid a flash before /api/config responds

  return (
    <div
      className={`mb-10 rounded-xl border p-4 ${
        hasKey ? "border-slate-200 bg-white" : "border-amber-200 bg-amber-50"
      }`}
    >
      <div className="mb-2 text-sm font-semibold text-slate-900">
        {hasKey
          ? "🔑 Ask the Atlas is enabled"
          : "🔑 Add an OpenAI API key to enable Ask the Atlas"}
      </div>
      <p className="mb-3 text-sm text-slate-600">
        {hasKey ? (
          <>
            The chat widget (bottom-right) is ready to use. The key is stored in
            this server&rsquo;s <code className="text-xs">.env.local</code>, so
            it persists across restarts.
          </>
        ) : (
          <>
            The chat widget is disabled until a key is set. Paste one below and
            it&rsquo;s saved to this server&rsquo;s{" "}
            <code className="text-xs">.env.local</code> — you only do this once.
            Get a key at{" "}
            <a
              href="https://platform.openai.com/api-keys"
              target="_blank"
              rel="noreferrer"
              className="underline hover:text-sky-700"
            >
              platform.openai.com/api-keys
            </a>
            .
          </>
        )}
      </p>
      <ApiKeyForm />
    </div>
  );
}
