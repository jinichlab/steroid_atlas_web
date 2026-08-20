/**
 * Cluster palette — matches the Python visualizer exactly so users switching
 * between the marimo and Next.js apps see the same colors.
 *
 * Algorithm: golden-angle hue rotation combined with alternating lightness
 * (0.62 / 0.36) and a 3-way saturation cycle (0.95 / 0.72 / 0.55). Yields
 * 82 pairwise-distinct colors where no two adjacent IDs share both hue and
 * lightness.
 */

function hslToRgb(h: number, s: number, l: number): [number, number, number] {
  const k = (n: number) => (n + h * 12) % 12;
  const a = s * Math.min(l, 1 - l);
  const f = (n: number) =>
    l - a * Math.max(-1, Math.min(k(n) - 3, 9 - k(n), 1));
  return [f(0), f(8), f(4)];
}

function rgbToHex([r, g, b]: [number, number, number]): string {
  const to = (x: number) =>
    Math.round(x * 255)
      .toString(16)
      .padStart(2, "0");
  return `#${to(r)}${to(g)}${to(b)}`;
}

export function clusterPalette(n: number): string[] {
  const out: string[] = [];
  for (let i = 0; i < n; i++) {
    const h = ((i * 137.508) / 360) % 1;
    const l = i % 2 === 0 ? 0.62 : 0.36;
    const s = [0.95, 0.72, 0.55][i % 3];
    out.push(rgbToHex(hslToRgb(h, s, l)));
  }
  return out;
}

export function clusterRgb(n: number): [number, number, number][] {
  const out: [number, number, number][] = [];
  for (let i = 0; i < n; i++) {
    const h = ((i * 137.508) / 360) % 1;
    const l = i % 2 === 0 ? 0.62 : 0.36;
    const s = [0.95, 0.72, 0.55][i % 3];
    const [r, g, b] = hslToRgb(h, s, l);
    out.push([Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)]);
  }
  return out;
}

/**
 * Perceptual distance between two RGB colors — cheap-and-good enough
 * (sum of squared component differences).
 */
function colorDist(a: [number, number, number], b: [number, number, number]) {
  const dr = a[0] - b[0], dg = a[1] - b[1], db = a[2] - b[2];
  return dr * dr + dg * dg + db * db;
}

/**
 * Assign a distinct-color mapping keyed by cluster id, given each cluster's
 * (cx, cy) centroid on the plot. For every cluster we pick the palette
 * color that maximizes the minimum perceptual distance to its spatial
 * neighbors' already-assigned colors — a greedy graph-colouring pass.
 *
 * This fixes the “two adjacent clusters end up as similar shades of green”
 * problem: even if their palette ids happen to be close on the golden
 * angle, they get shuffled apart in the final assignment.
 *
 * @param centroids   array of {id, cx, cy}
 * @param nColorsHint how many raw colors to draw from the base palette
 * @returns           Map<clusterId, [r,g,b]>
 */
export function distinctPaletteForCentroids(
  centroids: { id: number; cx: number; cy: number }[],
  nColorsHint = 200,
): Map<number, [number, number, number]> {
  const n = centroids.length;
  const pool = clusterRgb(Math.max(nColorsHint, n * 2));
  const assigned = new Map<number, [number, number, number]>();

  // Pre-compute each cluster's k-nearest neighbors (k=6).
  const K = Math.min(6, n - 1);
  const neighbors = new Map<number, number[]>();
  for (const c of centroids) {
    const others = centroids
      .filter((o) => o.id !== c.id)
      .map((o) => ({
        id: o.id,
        d: (o.cx - c.cx) ** 2 + (o.cy - c.cy) ** 2,
      }))
      .sort((a, b) => a.d - b.d)
      .slice(0, K)
      .map((o) => o.id);
    neighbors.set(c.id, others);
  }

  // Order clusters by size of neighborhood (most-constrained first). Since
  // we made a symmetric K, this is just by id — but we sort anyway to
  // stay deterministic.
  const order = [...centroids].sort((a, b) => a.id - b.id);

  const used = new Set<number>();
  for (const c of order) {
    const nbrs = neighbors.get(c.id) ?? [];
    const nbrColors: [number, number, number][] = nbrs
      .map((nid) => assigned.get(nid))
      .filter((x): x is [number, number, number] => !!x);

    let bestIdx = -1;
    let bestScore = -Infinity;
    for (let i = 0; i < pool.length; i++) {
      if (used.has(i)) continue;
      // Score = min distance to any already-picked neighbor color.
      let minDist = Infinity;
      for (const nb of nbrColors) {
        const d = colorDist(pool[i], nb);
        if (d < minDist) minDist = d;
      }
      if (minDist > bestScore) {
        bestScore = minDist;
        bestIdx = i;
      }
    }
    if (bestIdx === -1) bestIdx = pool.findIndex((_, i) => !used.has(i));
    if (bestIdx === -1) bestIdx = 0;
    used.add(bestIdx);
    assigned.set(c.id, pool[bestIdx]);
  }
  return assigned;
}
