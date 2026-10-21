/**
 * Where the OpenAI key for "Ask the Atlas" comes from.
 *
 *  - "server" (default): one key on the server (.env.local or a deployment
 *    secret), shared by every visitor.
 *  - "user" (API_KEY_MODE=user): each visitor pastes their own key. It lives
 *    only in their browser tab and rides along with each /api/chat request in
 *    USER_KEY_HEADER; the server uses it for that request and never stores it.
 *    Any server key is ignored, so a public deployment can't spend it.
 */
export type KeyMode = "server" | "user";

export const USER_KEY_HEADER = "x-openai-key";

export function keyMode(): KeyMode {
  return process.env.API_KEY_MODE === "user" ? "user" : "server";
}
