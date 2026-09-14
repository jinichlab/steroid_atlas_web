/**
 * Read/write a single variable in the project's `.env.local`, so a key entered
 * in the UI survives a server restart.
 *
 * Server-only (uses fs). `.env.local` is gitignored, so the key is never
 * committed.
 */
import { promises as fs } from "fs";
import path from "path";

const ENV_PATH = path.join(process.cwd(), ".env.local");

/** Reject anything that could break out of `NAME=value` into another line. */
export function isSafeEnvValue(value: string): boolean {
  return value.length > 0 && !/[\r\n]/.test(value);
}

/**
 * Set (or, with `null`, remove) `name` in .env.local, preserving every other
 * line — comments and blank lines included. Also updates `process.env` so the
 * change takes effect without a restart.
 */
export async function writeEnvVar(
  name: string,
  value: string | null,
): Promise<void> {
  let existing = "";
  try {
    existing = await fs.readFile(ENV_PATH, "utf8");
  } catch {
    // no .env.local yet — we'll create it
  }

  const lines = existing ? existing.split("\n") : [];
  const matches = (line: string) =>
    line.trim().replace(/^export\s+/, "").startsWith(`${name}=`);

  const kept = lines.filter((line) => !matches(line));
  if (value !== null) kept.push(`${name}=${value}`);

  // Collapse trailing blank lines, then end with exactly one newline.
  while (kept.length && kept[kept.length - 1].trim() === "") kept.pop();
  const next = kept.join("\n") + "\n";

  await fs.writeFile(ENV_PATH, next, { encoding: "utf8", mode: 0o600 });
  try {
    await fs.chmod(ENV_PATH, 0o600); // pre-existing file keeps its old mode otherwise
  } catch {
    // best-effort (e.g. unusual filesystems)
  }

  if (value === null) delete process.env[name];
  else process.env[name] = value;
}
