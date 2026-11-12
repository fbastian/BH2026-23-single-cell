/**
 * Page state in the URL query string, so that a link reopens the same view (filters, search,
 * pivot axes and drill-down, ...). Each component reads its parameters at start and writes
 * them back when they change, with history.replaceState (no new history entry per click).
 *
 * Lists are comma-separated, each value URI-encoded, so that a comma inside a value cannot
 * break a list; ':', '/', '=' and '|' are left readable (e.g. tissue=UBERON:0002048).
 */

/** Raw (still encoded) value of each parameter, in URL order. */
function rawParams(): Map<string, string> {
  const out = new Map<string, string>();
  for (const part of window.location.search.replace(/^\?/, "").split("&")) {
    if (!part) continue;
    const i = part.indexOf("=");
    out.set(decode(i < 0 ? part : part.slice(0, i)), i < 0 ? "" : part.slice(i + 1));
  }
  return out;
}

function decode(raw: string): string {
  try {
    return decodeURIComponent(raw.replace(/\+/g, " "));
  } catch {
    return raw;
  }
}

function encode(value: string): string {
  return encodeURIComponent(value)
    .replace(/%3A/gi, ":")
    .replace(/%2F/gi, "/")
    .replace(/%3D/gi, "=")
    .replace(/%7C/gi, "|");
}

/** A parameter of the URL, or null if absent. */
export function urlParam(name: string): string | null {
  const raw = rawParams().get(name);
  return raw == null ? null : decode(raw);
}

/** A comma-separated list parameter of the URL ([] if absent). */
export function urlList(name: string): string[] {
  const raw = rawParams().get(name);
  return raw ? raw.split(",").map(decode).filter(Boolean) : [];
}

/** Set parameters of the URL (a list is written comma-separated); null or empty removes one. */
export function setUrlParams(updates: Record<string, string | string[] | null | undefined>): void {
  const params = rawParams();
  for (const [key, value] of Object.entries(updates)) {
    if (value == null || value === "" || (Array.isArray(value) && !value.length)) params.delete(key);
    else params.set(key, Array.isArray(value) ? value.map(encode).join(",") : encode(value));
  }
  const search = [...params].map(([k, v]) => `${encode(k)}=${v}`).join("&");
  const next = `${window.location.pathname}${search ? `?${search}` : ""}${window.location.hash}`;
  if (next !== `${window.location.pathname}${window.location.search}${window.location.hash}`) {
    window.history.replaceState(window.history.state, "", next);
  }
}
