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
 * Categorical palette for views with ≤ 8 clusters. Checked across ALL pairs
 * (not just neighbours) on a white surface: normal-vision OKLab ΔE ≥ 18 and
 * protan / deutan / tritan ΔE ≥ 9, so every cluster is distinguishable,
 * including for colour-blind readers. #ccaa44 and #11aaaa sit below 3:1
 * contrast on white — fine here because every cluster is also named by its
 * on-plot label and the sidebar legend.
 */
const SMALL_CATEGORICAL = [
  "#4400ee", // indigo
  "#ff0066", // rose
  "#11aaaa", // teal
  "#ccaa44", // ochre
  "#008300", // green
  "#ff00ff", // magenta
  "#773388", // plum
  "#8877cc", // lavender
];

function hexToRgb(h: string): [number, number, number] {
  return [
    parseInt(h.slice(1, 3), 16),
    parseInt(h.slice(3, 5), 16),
    parseInt(h.slice(5, 7), 16),
  ];
}

type Lab = [number, number, number];

/** sRGB (0–255) → OKLab, a space where Euclidean distance ≈ perceived difference. */
function rgbToOklab([r, g, b]: [number, number, number]): Lab {
  const lin = (c: number) => {
    const x = c / 255;
    return x <= 0.04045 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4;
  };
  const [R, G, B] = [lin(r), lin(g), lin(b)];
  const l = Math.cbrt(0.4122214708 * R + 0.5363325363 * G + 0.0514459929 * B);
  const m = Math.cbrt(0.2119034982 * R + 0.6806995451 * G + 0.1073969566 * B);
  const s = Math.cbrt(0.0883024619 * R + 0.2817188376 * G + 0.6299787005 * B);
  return [
    0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s,
  ];
}

/** Squared OKLab distance. */
function colorDist(a: Lab, b: Lab) {
  const dl = a[0] - b[0], da = a[1] - b[1], db = a[2] - b[2];
  return dl * dl + da * da + db * db;
}

/**
 * n mutually-distinct colors (Glasbey-style): farthest-point sampling in
 * OKLab over a dense sRGB grid. Each pick maximises its minimum distance to
 * every color already chosen, so no two clusters anywhere on the map share
 * a near-identical shade. Candidates are limited to mid lightness with some
 * chroma, so points stay visible on both the light and dark themes and never
 * read as the grey used for unlabeled proteins.
 */
function glasbeyRgb(n: number): [number, number, number][] {
  const cands: [number, number, number][] = [];
  const labs: Lab[] = [];
  for (let r = 0; r <= 255; r += 15)
    for (let g = 0; g <= 255; g += 15)
      for (let b = 0; b <= 255; b += 15) {
        const lab = rgbToOklab([r, g, b]);
        const chroma = Math.hypot(lab[1], lab[2]);
        if (lab[0] < 0.42 || lab[0] > 0.85 || chroma < 0.07) continue;
        cands.push([r, g, b]);
        labs.push(lab);
      }

  // Seed the distances with colors to stay away from: the unlabeled grey
  // plus the white / near-black plot backgrounds.
  const avoid = ([[190, 190, 190], [255, 255, 255], [20, 20, 20]] as const).map(
    (c) => rgbToOklab([c[0], c[1], c[2]]),
  );
  const minD = labs.map((lab) => Math.min(...avoid.map((a) => colorDist(lab, a))));

  const out: [number, number, number][] = [];
  for (let k = 0; k < n && k < cands.length; k++) {
    let best = 0;
    for (let i = 1; i < minD.length; i++) if (minD[i] > minD[best]) best = i;
    out.push(cands[best]);
    const picked = labs[best];
    for (let i = 0; i < minD.length; i++) {
      const d = colorDist(labs[i], picked);
      if (d < minD[i]) minD[i] = d;
    }
  }
  return out;
}

/**
 * Assign a distinct-color mapping keyed by cluster id, given each cluster's
 * (cx, cy) centroid on the plot. The pool is exactly n Glasbey colors, so
 * every pair of clusters is globally distinct; on top of that, a greedy
 * graph-colouring pass hands each cluster the remaining color farthest
 * (in OKLab) from its spatial neighbours' colors, so the closest-looking
 * pairs never sit next to each other on the map.
 *
 * @param centroids   array of {id, cx, cy}
 * @returns           Map<clusterId, [r,g,b]>
 */
export function distinctPaletteForCentroids(
  centroids: { id: number; cx: number; cy: number }[],
): Map<number, [number, number, number]> {
  const n = centroids.length;

  // Few clusters (steroid view = 8, natsyn = 2): every cluster gets its own
  // fixed categorical color, assigned by cluster id. The greedy pass below
  // only keeps *map neighbours* apart, so with a small view it could hand
  // two non-adjacent clusters near-identical shades.
  if (n <= SMALL_CATEGORICAL.length) {
    const ids = centroids.map((c) => c.id).sort((a, b) => a - b);
    return new Map(ids.map((id, i) => [id, hexToRgb(SMALL_CATEGORICAL[i])]));
  }

  const pool = glasbeyRgb(n);
  const poolLab = pool.map(rgbToOklab);
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
    const nbrColors: Lab[] = nbrs
      .map((nid) => assigned.get(nid))
      .filter((x): x is [number, number, number] => !!x)
      .map(rgbToOklab);

    let bestIdx = -1;
    let bestScore = -Infinity;
    for (let i = 0; i < pool.length; i++) {
      if (used.has(i)) continue;
      // Score = min distance to any already-picked neighbor color.
      let minDist = Infinity;
      for (const nb of nbrColors) {
        const d = colorDist(poolLab[i], nb);
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
