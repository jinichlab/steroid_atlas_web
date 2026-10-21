"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import Link from "next/link";
import dynamic from "next/dynamic";
import type {
  ViewKind,
  Protein,
  Molecule,
  NatsynEntry,
  ClusterMeta,
} from "@/lib/types";
import type { UmapPoint, ClusterCentroidLabel } from "@/components/UmapPlot";
import { clusterPalette, distinctPaletteForCentroids } from "@/lib/palette";
import { useSelection, type SelectionItem } from "@/lib/selection-context";

// Cap on how many selected items we hand to the chat assistant — keeps the
// /api/chat payload and prompt size bounded on huge lasso selections.
const MAX_SELECTION_ITEMS = 80;

const UmapPlot = dynamic(() => import("@/components/UmapPlot"), { ssr: false });

type Row = Protein | Molecule | NatsynEntry;

interface Props {
  kind: ViewKind;
  clusters: ClusterMeta[];
  /** High-dimensional k-means clusters (analysis/scripts/19_highdim_kmeans.py);
   *  empty when not built — the clustering toggle then stays hidden. */
  clustersHd: ClusterMeta[];
}

/** "umap": k-means on the 2D map (published); "hd": k-means in the full
 *  embedding (ProtT5 for proteins, Morgan fingerprints for steroids). */
type ClusterMode = "umap" | "hd";

/** A row's cluster id under the active clustering; -1 = no high-dim label. */
function clusterIdOf(r: Row, mode: ClusterMode): number {
  const raw = mode === "hd" ? (r as any).cluster_hd : (r as any).cluster;
  if (typeof raw === "number") return raw;
  if (typeof raw === "string") return raw === "natural" ? 0 : 1; // natsyn groups
  return mode === "hd" ? -1 : 0;
}

function clusterText(r: Row, mode: ClusterMode): string {
  const c = clusterIdOf(r, mode);
  const label = mode === "hd" ? "Embedding cluster" : "Cluster";
  return c < 0 ? `${label} —` : `${label} ${c + 1}`;
}

const DATA_URL: Record<ViewKind, string> = {
  protein: "/atlas/proteins.json",
  molecule: "/atlas/molecules.json",
  natsyn: "/atlas/natsyn.json",
};
const VIEW_LABELS: Record<ViewKind, string> = {
  protein: "Protein centric",
  molecule: "Steroid centric",
  natsyn: "Natural + synthetic",
};

function rowKey(kind: ViewKind, r: Row): string {
  if (kind === "protein") return (r as Protein).accession;
  if (kind === "molecule") return (r as Molecule).compound_name;
  return `${(r as NatsynEntry).compound_name}::${(r as NatsynEntry).chebi_id}`;
}
function rowLabel(kind: ViewKind, r: Row): string {
  if (kind === "protein") {
    const p = r as Protein;
    const nm = (p.protein_names ?? "").slice(0, 80);
    const org = (p.organism ?? "").slice(0, 40);
    const acc = p.accession ? ` [${p.accession}]` : "";
    return org ? `${nm}${acc} · ${org}` : `${nm}${acc}`;
  }
  return ((r as Molecule).compound_name ?? "").slice(0, 80);
}

function esc(s: string): string {
  return String(s ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

/** Marimo-style multi-line HTML tooltip: name / accession · gene · organism /
 *  EC + ChEBI counts + cluster + newly-recruited badge. */
function tooltipHtmlFor(kind: ViewKind, r: Row, mode: ClusterMode): string {
  if (kind === "protein") {
    const p = r as Protein;
    const name = esc((p.protein_names ?? "").slice(0, 90));
    const gene = esc(p.gene_names ?? "");
    const org = esc(p.organism ?? "");
    const acc = esc(p.accession ?? "");
    const ecs = esc((p.ec_numbers ?? "").split(";").slice(0, 2).join(" · "));
    const cluster = clusterText(p, mode);
    const star = p.is_literature_recruited
      ? `<span style="color:#f59e0b;font-weight:700;">★ NEW</span> `
      : "";
    return `
      <div style="font-size:12px;line-height:1.35;max-width:340px;">
        <div style="font-weight:600;color:currentColor;">${name}</div>
        <div style="margin-top:3px;font-size:11px;opacity:.85;">
          ${star}<span style="font-family:ui-monospace,monospace;">${acc}</span>
          ${gene ? ` · ${gene}` : ""}
        </div>
        <div style="margin-top:2px;font-size:11px;opacity:.75;">
          <em>${org}</em>
        </div>
        <div style="margin-top:2px;font-size:11px;opacity:.75;">
          ${cluster}${ecs ? ` · EC ${ecs}` : ""}
        </div>
      </div>`;
  }
  const m = r as Molecule;
  const name = esc((m.compound_name ?? "").slice(0, 90));
  const chebi = esc(m.chebi_id ?? "");
  return `
    <div style="font-size:12px;line-height:1.35;max-width:320px;">
      <div style="font-weight:600;">${name}</div>
      ${chebi ? `<div style="margin-top:3px;font-family:ui-monospace,monospace;font-size:11px;opacity:.75;">CHEBI:${chebi}</div>` : ""}
    </div>`;
}
function splitList(s: string): string[] {
  if (!s) return [];
  return Array.from(
    new Set(s.replace(/;/g, "\n").split("\n").map((x) => x.trim()).filter(Boolean)),
  );
}

/** `molecules.json` stores interacting proteins as a Python-list string:
 *  "['A0A016SKM3', 'A0A016VXS0', …]". */
function parseAccessionList(s: string): string[] {
  if (!s) return [];
  return s
    .replace(/[[\]'"]/g, "")
    .split(",")
    .map((x) => x.trim())
    .filter(Boolean);
}

/** "CHEBI:17263", "17263.0", 17263 → "17263"; anything else → "". */
function chebiNum(c: string | number | null | undefined): string {
  const s = String(c ?? "")
    .trim()
    .toLowerCase()
    .replace(/^chebi:/, "");
  return /^\d+(\.\d+)?$/.test(s) ? String(parseInt(s, 10)) : "";
}

type StructIndex = Record<string, string>;

export default function ExploreClient({ kind, clusters: clustersUmap, clustersHd }: Props) {
  // Clustering toggle (protein + steroid views, when high-dim labels exist).
  const hasHd = kind !== "natsyn" && clustersHd.length > 0;
  const [clusterMode, setClusterMode] = useState<ClusterMode>("umap");
  const mode: ClusterMode = hasHd ? clusterMode : "umap";
  const clusters = mode === "hd" ? clustersHd : clustersUmap;
  const { selection: crossSel, setSelection, clearSelection } = useSelection();
  const [rows, setRows] = useState<Row[]>([]);
  const [loading, setLoading] = useState(true);
  const [structIdx, setStructIdx] = useState<StructIndex>({});
  const [highlight, setHighlight] = useState<number | null>(null);
  const [selectionSet, setSelectionSet] = useState<Set<string> | null>(null);
  const [focusedPoint, setFocusedPoint] = useState<Row | null>(null);
  // Multi-select for the detail cards — clicking a row in the protein table
  // (or "add to details" on the focused-point card) toggles that protein here.
  const [focusedProteins, setFocusedProteins] = useState<Protein[]>([]);
  const toggleProtein = (p: Protein) =>
    setFocusedProteins((cur) => {
      const i = cur.findIndex((x) => x.accession === p.accession);
      if (i >= 0) return [...cur.slice(0, i), ...cur.slice(i + 1)];
      return [...cur, p];
    });
  const isProteinFocused = (acc: string) =>
    focusedProteins.some((x) => x.accession === acc);
  const clearFocusedProteins = () => setFocusedProteins([]);
  // Protein view: the interacting-steroid tile the user clicked in a detail card.
  const [focusedSteroid, setFocusedSteroid] = useState<{
    name: string;
    chebi: string;
  } | null>(null);
  const [legendOpen, setLegendOpen] = useState(true);
  // Free-text search across names / accessions / gene / GO / EC / keywords /
  // ChEBI / compound / sequence.
  const [search, setSearch] = useState("");
  // A bumping counter that tells UmapPlot to erase its sticky rectangle
  // when the user hits the top-of-plot "clear" button.
  const [lassoClearSignal, setLassoClearSignal] = useState(0);
  // The "🧬 Steroids in this…" results panel sits far below the plot; after a
  // search (attention is at the top search bar) the user has no cue it appeared.
  const resultsRef = useRef<HTMLDivElement | null>(null);
  const scrollToResults = () =>
    resultsRef.current?.scrollIntoView({ behavior: "smooth", block: "start" });
  // The steroid detail card renders below the structure grid, out of view of
  // the tile that was clicked — scroll to it after a tile click (not after a
  // map-dot click, which would yank the user away from the plot).
  const steroidCardRef = useRef<HTMLDivElement | null>(null);
  const scrollToSteroidCard = useRef(false);

  useEffect(() => {
    let alive = true;
    setLoading(true);
    fetch(DATA_URL[kind])
      .then((r) => r.json())
      .then((d: Row[]) => {
        if (!alive) return;
        setRows(d);
        setLoading(false);
      })
      .catch(() => setLoading(false));
    return () => {
      alive = false;
    };
  }, [kind]);

  useEffect(() => {
    let alive = true;
    fetch("/atlas/structures_index.json")
      .then((r) => r.json())
      .then((d: StructIndex) => alive && setStructIdx(d))
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, []);

  const structUrl = (
    name: string,
    chebi?: string | number | null,
  ): string | null => {
    // molecules.json compound_name has a leading space; its chebi_id is a JSON
    // float (17263.0) that JSON.parse turns into 17263, while the index keys it
    // as "chebi:17263.0". Normalise both, and ignore junk chebi ("[]", "nan").
    const nm = (name || "").trim().toLowerCase();
    if (nm && structIdx[nm]) return `/atlas/structures/${structIdx[nm]}`;

    const c = String(chebi ?? "")
      .trim()
      .toLowerCase()
      .replace(/^chebi:/, "");
    if (/^\d+(\.\d+)?$/.test(c)) {
      const n = String(parseInt(c, 10));
      for (const key of [`chebi:${c}`, `chebi:${n}`, `chebi:${n}.0`]) {
        if (structIdx[key]) return `/atlas/structures/${structIdx[key]}`;
      }
    }
    return null;
  };

  const points: UmapPoint[] = useMemo(
    () =>
      rows.map((r) => {
        const clusterNum = clusterIdOf(r, mode);
        return {
          x: (r as any).umap_1,
          y: (r as any).umap_2,
          cluster: clusterNum,
          is_star: (r as any).is_literature_recruited ?? 0,
          label: rowLabel(kind, r),
          key: rowKey(kind, r),
          tooltipHtml: tooltipHtmlFor(kind, r, mode),
        };
      }),
    [rows, kind, mode],
  );

  // ── Search — filters into a Set<key> that feeds the same "selection" slot
  //  the lasso uses, so search matches highlight on the map + populate the
  //  catalogue below the same way.
  const searchMatchKeys = useMemo<Set<string> | null>(() => {
    const q = search.trim().toLowerCase();
    if (!q) return null;
    // Multi-column contains-any-token match
    const tokens = q.split(/\s+/).filter(Boolean);
    const match = (r: Row): boolean => {
      const parts: string[] = [];
      if (kind === "protein") {
        const p = r as Protein;
        parts.push(
          p.accession, p.entry_name, p.protein_names, p.gene_names,
          p.organism, p.ec_numbers, p.rhea_reactions,
          p.interacting_chebi_ids, p.interacting_compounds,
          p.go_ids, p.go_labels, p.keyword_ids, p.keyword_labels,
        );
      } else if (kind === "molecule") {
        const m = r as Molecule;
        parts.push(m.compound_name, m.chebi_id, m.smiles,
                   m.interacting_protein_accessions);
      } else {
        const n = r as NatsynEntry;
        parts.push(n.compound_name, n.chebi_id, n.smiles, n.protein_entries);
      }
      const blob = parts.join("\n").toLowerCase();
      return tokens.every((t) => blob.includes(t));
    };
    const s = new Set<string>();
    for (const r of rows) if (match(r)) s.add(rowKey(kind, r));
    return s;
  }, [rows, search, kind]);

  const focusedPointKey = focusedPoint ? rowKey(kind, focusedPoint) : null;

  useEffect(() => {
    if (!focusedPoint || !scrollToSteroidCard.current) return;
    scrollToSteroidCard.current = false;
    steroidCardRef.current?.scrollIntoView({ behavior: "smooth", block: "start" });
  }, [focusedPoint]);

  // Cluster centroids drive on-plot labels AND the distinct-color palette
  // (so neighbors on the map don't end up as similar shades of green).
  const centroids = useMemo(() => {
    const byCid = new Map<number, { sx: number; sy: number; n: number }>();
    for (const p of points) {
      if (p.cluster < 0) continue; // no high-dim label
      const e = byCid.get(p.cluster);
      if (e) {
        e.sx += p.x;
        e.sy += p.y;
        e.n += 1;
      } else {
        byCid.set(p.cluster, { sx: p.x, sy: p.y, n: 1 });
      }
    }
    return Array.from(byCid.entries()).map(([id, v]) => ({
      id,
      cx: v.sx / v.n,
      cy: v.sy / v.n,
    }));
  }, [points]);

  const paletteMap = useMemo(() => {
    if (!centroids.length) return new Map<number, [number, number, number]>();
    const m = distinctPaletteForCentroids(centroids);
    m.set(-1, [190, 190, 190]); // proteins without a high-dim label: grey
    return m;
  }, [centroids]);

  const paletteFallback = useMemo(
    () => clusterPalette(Math.max(clusters.length, 2)),
    [clusters],
  );
  const swatchFor = (cid: number): string => {
    const m = paletteMap.get(cid);
    if (m) return `rgb(${m[0]},${m[1]},${m[2]})`;
    return paletteFallback[cid] ?? "#999";
  };

  const clusterLabels: ClusterCentroidLabel[] = useMemo(
    () =>
      clusters.map((c) => ({
        cluster: c.id,
        display_id: c.display_id,
        name: c.name || "",
        n: c.n,
      })),
    [clusters],
  );

  const hasLocalSelection =
    (searchMatchKeys?.size ?? 0) > 0 ||
    (selectionSet?.size ?? 0) > 0 ||
    highlight != null;

  // rowKey → interacting protein accessions, for the molecule view (built once).
  const molAccByKey = useMemo(() => {
    if (kind === "protein") return null;
    const map = new Map<string, string[]>();
    for (const r of rows) {
      // natsyn.json stores the same Python-list string under a different key.
      const raw =
        kind === "natsyn"
          ? (r as NatsynEntry).protein_entries
          : (r as Molecule).interacting_protein_accessions;
      map.set(rowKey(kind, r), parseAccessionList(raw ?? ""));
    }
    return map;
  }, [rows, kind]);

  // Protein view: compound → proteins that interact with it, inverted from
  // each protein's interacting_compounds / interacting_chebi_ids.
  const proteinsByCompound = useMemo(() => {
    if (kind !== "protein") return null;
    const map = new Map<string, Protein[]>();
    const add = (k: string, p: Protein) => {
      const cur = map.get(k);
      if (cur) cur.push(p);
      else map.set(k, [p]);
    };
    for (const r of rows) {
      const p = r as Protein;
      for (const nm of splitList(p.interacting_compounds ?? ""))
        add(`n:${nm.toLowerCase()}`, p);
      for (const c of splitList(p.interacting_chebi_ids ?? "")) {
        const n = chebiNum(c);
        if (n) add(`c:${n}`, p);
      }
    }
    return map;
  }, [rows, kind]);

  // Molecule / natsyn views: accession → Protein, fetched lazily the first
  // time a steroid card needs protein names.
  const [proteinByAcc, setProteinByAcc] = useState<Map<string, Protein> | null>(
    null,
  );
  const proteinFetchStarted = useRef(false);
  useEffect(() => {
    if (kind === "protein" || !focusedPoint || proteinFetchStarted.current)
      return;
    proteinFetchStarted.current = true;
    fetch(DATA_URL.protein)
      .then((r) => r.json())
      .then((d: Protein[]) =>
        setProteinByAcc(new Map(d.map((p) => [p.accession, p]))),
      )
      .catch(() => {
        proteinFetchStarted.current = false;
      });
  }, [kind, focusedPoint]);

  // ── Cross-view highlight ──────────────────────────────────────────────
  // The keys this view should light up given a selection made in another
  // view (protein ↔ molecule via interacting-protein accessions). Purely
  // derived — never written back to the context.
  const incomingKeys = useMemo<Set<string> | null>(() => {
    // A live local selection overrides anything carried in from another view;
    // returning null here also keeps `pool` referentially stable (no re-render
    // loop with the context-publish effect below).
    if (hasLocalSelection) return null;
    if (!rows.length || !crossSel.view) return null;

    const isCrossView = crossSel.view !== kind;
    const accSet = new Set(crossSel.proteinAccessions);
    const s = new Set<string>();

    if (kind === "protein") {
      for (const r of rows) {
        if (accSet.has((r as Protein).accession)) s.add(rowKey(kind, r));
      }
    } else if (crossSel.view === kind) {
      // returning to the view that owns the selection — restore it exactly
      const mk = new Set(crossSel.moleculeKeys);
      for (const r of rows) {
        const k = rowKey(kind, r);
        if (mk.has(k)) s.add(k);
      }
    } else if (kind === "molecule" && crossSel.view === "protein" && molAccByKey) {
      for (const r of rows) {
        const k = rowKey(kind, r);
        if ((molAccByKey.get(k) ?? []).some((a) => accSet.has(a))) s.add(k);
      }
    } else {
      return null; // e.g. natsyn ← protein: no interaction data
    }

    if (s.size) return s;
    // Nothing here matches the other view's selection. For a real cross-view
    // carry-over of a non-empty selection, return the EMPTY set so the plot
    // greys everything ("interacts with nothing in this view"); for a failed
    // same-view restore, fall back to null (leave the plot untouched).
    return isCrossView && crossSel.total > 0 ? s : null;
  }, [rows, kind, crossSel, molAccByKey, hasLocalSelection]);

  // Cross-view highlight is showing (no local edit in this view yet).
  const showingIncoming = !!incomingKeys && crossSel.view !== kind;
  const incomingEmpty = showingIncoming && incomingKeys?.size === 0;

  const pool = useMemo(() => {
    if (searchMatchKeys && searchMatchKeys.size) {
      return rows.filter((r) => searchMatchKeys.has(rowKey(kind, r)));
    }
    if (selectionSet && selectionSet.size) {
      return rows.filter((r) => selectionSet.has(rowKey(kind, r)));
    }
    if (highlight != null) {
      return rows.filter((r) => clusterIdOf(r, mode) === highlight);
    }
    if (incomingKeys && incomingKeys.size) {
      return rows.filter((r) => incomingKeys.has(rowKey(kind, r)));
    }
    return [];
  }, [rows, highlight, selectionSet, searchMatchKeys, incomingKeys, kind, mode]);

  // ── Publish a *local* selection to the shared context (chat + other views).
  //  Incoming cross-view highlights are display-only and never written back.
  useEffect(() => {
    if (!hasLocalSelection) return;
    const origin =
      searchMatchKeys && searchMatchKeys.size
        ? "search"
        : selectionSet && selectionSet.size
          ? "lasso"
          : "cluster";
    const items: SelectionItem[] = pool
      .slice(0, MAX_SELECTION_ITEMS)
      .map((r) => {
        if (kind === "protein") {
          const p = r as Protein;
          return {
            kind: "protein",
            name: p.protein_names || p.accession,
            accession: p.accession,
            gene: p.gene_names || null,
            organism: p.organism || null,
            cluster: clusterIdOf(p, mode),
          };
        }
        const m = r as Molecule & NatsynEntry;
        return {
          kind: "molecule",
          name: m.compound_name,
          chebi: m.chebi_id || null,
          smiles: m.smiles || null,
          cluster: clusterIdOf(m, mode),
        };
      });
    const proteinAccessions =
      kind === "protein"
        ? (pool as Protein[]).map((p) => p.accession)
        : Array.from(
            new Set(
              pool.flatMap((r) =>
                molAccByKey?.get(rowKey(kind, r)) ?? [],
              ),
            ),
          );
    const moleculeKeys =
      kind === "protein" ? [] : pool.map((r) => rowKey(kind, r));
    setSelection({
      view: kind,
      origin,
      items,
      total: pool.length,
      proteinAccessions,
      moleculeKeys,
    });
  }, [
    pool,
    kind,
    hasLocalSelection,
    searchMatchKeys,
    selectionSet,
    molAccByKey,
    setSelection,
    mode,
  ]);

  // Clear the open detail cards when the user jumps to a different cluster.
  useEffect(() => {
    clearFocusedProteins();
    setFocusedSteroid(null);
  }, [highlight]);

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <div className="mx-auto max-w-[1500px] px-4 py-6">
        {/* Header */}
        <div className="mb-3 flex flex-wrap items-center gap-3">
          <Link href="/" className="text-sm text-slate-500 hover:text-slate-900">
            ← Atlas home
          </Link>
          <div className="ml-2 text-xl font-semibold tracking-tight text-slate-900">
            {VIEW_LABELS[kind]}
          </div>
          <div className="text-sm text-slate-500">
            {loading ? "loading…" : `${rows.length.toLocaleString()} entries`}
          </div>
          <div className="ml-auto flex gap-2 text-xs">
            {(["protein", "molecule", "natsyn"] as ViewKind[]).map((v) => (
              <Link
                key={v}
                href={`/explore/${v}`}
                className={`rounded-md border px-3 py-1 transition ${
                  v === kind
                    ? "border-sky-600 bg-sky-600 text-white"
                    : "border-slate-300 bg-white text-slate-700 hover:bg-slate-100"
                }`}
              >
                {VIEW_LABELS[v]}
              </Link>
            ))}
          </div>
        </div>

        {/* Search bar */}
        <div className="mb-3 flex flex-wrap items-center gap-2">
          <input
            type="text"
            placeholder={
              kind === "protein"
                ? "Search — name, gene, accession, GO, EC, ChEBI, Rhea, keyword, sequence…"
                : kind === "molecule"
                  ? "Search — compound name, ChEBI, SMILES…"
                  : "Search — compound name, ChEBI, SMILES…"
            }
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter") scrollToResults();
            }}
            className="w-full max-w-2xl rounded-md border border-slate-300 bg-white px-3 py-1.5 text-sm shadow-sm focus:border-sky-500 focus:outline-none focus:ring-2 focus:ring-sky-100"
          />
          {search && (
            <div className="flex items-center gap-2 text-xs text-slate-500">
              {searchMatchKeys && searchMatchKeys.size > 0 ? (
                <button
                  className="rounded border border-sky-300 bg-sky-50 px-2 py-0.5 font-medium text-sky-700 hover:bg-sky-100"
                  onClick={scrollToResults}
                >
                  ⌄ {searchMatchKeys.size.toLocaleString()}{" "}
                  {kind === "protein" ? "proteins" : "structures"}
                </button>
              ) : (
                <span>no matches</span>
              )}
              <button
                className="rounded border border-slate-300 bg-white px-2 py-0.5 text-slate-700 hover:bg-slate-100"
                onClick={() => setSearch("")}
              >
                ✕ clear
              </button>
            </div>
          )}
        </div>

        {/* Above-plot legend */}
        <div className="mb-2 flex flex-wrap items-center gap-4 text-xs text-slate-600">
          <span className="inline-flex items-center gap-1">
            <span className="inline-block h-3 w-3 rounded-full bg-cyan-500" />{" "}
            existing
          </span>
          <span className="inline-flex items-center gap-1">
            <svg width="14" height="14" viewBox="-1 -1 2 2" className="inline-block">
              <path
                d="M0,-0.95 L0.22,-0.31 L0.9,-0.31 L0.36,0.12 L0.56,0.78 L0,0.38 L-0.56,0.78 L-0.36,0.12 L-0.9,-0.31 L-0.22,-0.31 Z"
                fill="currentColor"
                className="text-cyan-500"
              />
            </svg>
            newly recruited (colored by cluster)
          </span>
          <span className="inline-flex items-center gap-1">
            <span className="inline-block h-3 w-3 rounded-full bg-cyan-500 ring-2 ring-rose-500" />{" "}
            selected dot
          </span>
          {(highlight != null ||
            focusedPoint ||
            selectionSet ||
            search ||
            // any carried-over selection, incl. one restored when returning
            // to the view that made it (not just cross-view highlights)
            !!incomingKeys) && (
            <button
              className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-slate-700 hover:bg-slate-100"
              onClick={() => {
                setHighlight(null);
                setSelectionSet(null);
                setFocusedPoint(null);
                clearFocusedProteins();
                setSearch("");
                setLassoClearSignal((n) => n + 1);
                clearSelection();
              }}
            >
              ✕ clear
            </button>
          )}
        </div>

        {showingIncoming && (
          <div
            className={`mb-3 flex flex-wrap items-center gap-2 rounded-md border px-3 py-2 text-xs ${
              incomingEmpty
                ? "border-amber-200 bg-amber-50 text-amber-800"
                : "border-sky-200 bg-sky-50 text-sky-800"
            }`}
          >
            <span>
              {incomingEmpty ? (
                <>
                  None of the {kind === "protein" ? "proteins" : "steroids"} in
                  this view interact with the{" "}
                  <b>
                    {crossSel.total.toLocaleString()}{" "}
                    {crossSel.view === "protein" ? "proteins" : "steroids"}
                  </b>{" "}
                  selected in {VIEW_LABELS[crossSel.view as ViewKind]} view —
                  everything is greyed out.
                </>
              ) : (
                <>
                  Highlighting{" "}
                  <b>
                    {(incomingKeys?.size ?? 0).toLocaleString()}{" "}
                    {kind === "protein" ? "proteins" : "steroids"}
                  </b>{" "}
                  that interact with the{" "}
                  <b>
                    {crossSel.total.toLocaleString()}{" "}
                    {crossSel.view === "protein" ? "proteins" : "steroids"}
                  </b>{" "}
                  selected in {VIEW_LABELS[crossSel.view as ViewKind]} view.
                  Lasso or search to refine.
                </>
              )}
            </span>
          </div>
        )}

        <div className="grid gap-4 md:grid-cols-[minmax(0,1fr)_260px]">
          <div>
            <UmapPlot
              points={points}
              clusterLabels={clusterLabels}
              highlightCluster={highlight}
              matchKeys={
                searchMatchKeys && searchMatchKeys.size
                  ? searchMatchKeys
                  : selectionSet && selectionSet.size
                    ? selectionSet
                    : // incomingKeys may be an EMPTY set on purpose (cross-view
                      // selection that matches nothing here) → grey everything.
                      (incomingKeys ?? null)
              }
              focusedKey={focusedPointKey}
              darkMode={false}
              paletteOverride={paletteMap}
              clearLassoSignal={lassoClearSignal}
              onClickPoint={(p) => {
                const r = rows.find((r) => rowKey(kind, r) === p.key);
                if (r) setFocusedPoint(r);
              }}
              onLassoSelect={(pts) => {
                const s = new Set(pts.map((p) => p.key!).filter(Boolean));
                setSelectionSet(s.size ? s : null);
                setHighlight(null);
                setFocusedPoint(null);
                clearFocusedProteins();
              }}
              height={640}
            />
            {focusedPoint && kind === "protein" && (
              <div className="mt-3 rounded-xl border-2 border-sky-500 bg-sky-50 p-4">
                <div className="text-xs uppercase tracking-wide text-sky-700">
                  {kind === "protein" ? "Selected protein" : "Selected"}
                </div>
                <div className="mt-1 text-lg font-semibold text-slate-900">
                  {kind === "protein"
                    ? (focusedPoint as Protein).protein_names
                    : (focusedPoint as Molecule).compound_name}
                </div>
                {kind === "protein" && (
                  <div className="mt-1 text-sm text-slate-600">
                    {(focusedPoint as Protein).gene_names} ·{" "}
                    <em>{(focusedPoint as Protein).organism}</em> ·{" "}
                    {(focusedPoint as Protein).length_aa} aa ·{" "}
                    {clusterText(focusedPoint, mode)}
                    {(focusedPoint as Protein).is_literature_recruited
                      ? " · ★ newly recruited"
                      : ""}
                  </div>
                )}
                <div className="mt-2 flex flex-wrap gap-2 text-sm">
                  {kind === "protein" && (
                    <>
                      <a
                        className="rounded border border-sky-500 px-2 py-1 text-sky-700 hover:bg-sky-100"
                        href={`https://www.uniprot.org/uniprotkb/${(focusedPoint as Protein).accession}/entry`}
                        target="_blank"
                      >
                        UniProt
                      </a>
                      <a
                        className="rounded border border-emerald-500 px-2 py-1 text-emerald-700 hover:bg-emerald-100"
                        href={`https://alphafold.ebi.ac.uk/entry/${(focusedPoint as Protein).accession}`}
                        target="_blank"
                      >
                        AlphaFold
                      </a>
                      <button
                        className="rounded border border-slate-400 px-2 py-1 text-slate-700 hover:bg-slate-100"
                        onClick={() =>
                          toggleProtein(focusedPoint as Protein)
                        }
                      >
                        {isProteinFocused((focusedPoint as Protein).accession)
                          ? "✓ In details"
                          : "＋ Add to details ↓"}
                      </button>
                    </>
                  )}
                  <button
                    className="ml-auto text-xs text-slate-500 hover:text-slate-900"
                    onClick={() => setFocusedPoint(null)}
                  >
                    clear
                  </button>
                </div>
              </div>
            )}
          </div>

          {/* Cluster sidebar */}
          <aside className="rounded-xl border border-slate-200 bg-white shadow-sm">
            {hasHd && (
              <div className="border-b border-slate-100 p-3">
                <div className="mb-1.5 text-xs font-medium text-slate-500">
                  Clustering
                </div>
                <div className="flex rounded-md border border-slate-200 p-0.5 text-xs">
                  {(
                    [
                      ["umap", kind === "molecule" ? "2D clustering" : "2D map"],
                      ["hd", kind === "molecule" ? "Fingerprint clustering" : "Full embedding"],
                    ] as const
                  ).map(([m, label]) => (
                    <button
                      key={m}
                      onClick={() => {
                        setClusterMode(m);
                        setHighlight(null); // cluster ids differ between modes
                      }}
                      className={`flex-1 rounded px-2 py-1 transition ${
                        mode === m
                          ? "bg-slate-900 text-white"
                          : "text-slate-600 hover:bg-slate-100"
                      }`}
                    >
                      {label}
                    </button>
                  ))}
                </div>
                <div className="mt-1.5 text-[11px] leading-snug text-slate-500">
                  {mode === "hd"
                    ? `k-means in the full ${kind === "protein" ? "ProtT5 embedding" : "fingerprint space"}, same k. Dots stay at their map positions.`
                    : "k-means on the 2D map coordinates (published clusters)."}
                </div>
              </div>
            )}
            <button
              onClick={() => setLegendOpen((v) => !v)}
              className="flex w-full items-center justify-between border-b border-slate-100 p-3 text-left"
            >
              <div>
                <div className="text-sm font-semibold text-slate-900">
                  Clusters ({clusters.length})
                </div>
                <div className="text-xs text-slate-500">
                  {legendOpen ? "click any row to highlight" : "click to expand"}
                </div>
              </div>
              <span className="text-slate-500">{legendOpen ? "▾" : "▸"}</span>
            </button>
            {legendOpen && (
              <div className="max-h-[540px] overflow-y-auto p-2">
                {clusters.map((c) => {
                  const isPick = highlight === c.id;
                  const clr = swatchFor(c.id);
                  return (
                    <button
                      key={c.id}
                      onClick={() => {
                        setHighlight(isPick ? null : c.id);
                        clearFocusedProteins();
                      }}
                      className={`mb-1 flex w-full items-center gap-2 rounded-md border px-2 py-1.5 text-left text-[12px] transition ${
                        isPick
                          ? "border-sky-500 bg-sky-50 text-slate-900"
                          : "border-slate-200 bg-white text-slate-700 hover:bg-slate-50"
                      }`}
                    >
                      <span
                        className="inline-block h-4 w-4 flex-shrink-0 rounded-sm ring-1 ring-black/10"
                        style={{ background: clr }}
                      />
                      <span className="w-6 flex-shrink-0 font-mono text-[11px] text-slate-500">
                        {c.display_id}
                      </span>
                      <span className="flex-1 truncate">
                        {c.name || "(unnamed)"}
                      </span>
                      <span className="flex-shrink-0 text-[10px] text-slate-500">
                        n={c.n.toLocaleString()}
                      </span>
                    </button>
                  );
                })}
              </div>
            )}
          </aside>
        </div>

        {/* Protein view — table of the selected / matched proteins. Click a
            row to open (or close) that protein's full detail card + its
            interacting-steroid structures, in the "Detail cards" section
            below — no second click needed. */}
        {kind === "protein" && pool.length > 0 && (
          <div className="mt-6" ref={resultsRef}>
            <div className="mb-3 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                {searchMatchKeys && searchMatchKeys.size
                  ? `🔬 ${pool.length.toLocaleString()} proteins matching “${search.trim()}”`
                  : "🔬 Proteins in this selection"}
              </div>
              <div className="text-sm text-slate-500">
                click a row to open its detail + structures below · click again
                to close
              </div>
              {focusedProteins.length > 0 && (
                <button
                  className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                  onClick={clearFocusedProteins}
                >
                  ✕ close all {focusedProteins.length} open
                </button>
              )}
            </div>
            <div className="max-h-[420px] overflow-y-auto rounded-lg border border-slate-200 bg-white">
              {(pool as Protein[]).slice(0, 300).map((p) => {
                const isPick = isProteinFocused(p.accession);
                const nSteroids = splitList(
                  p.interacting_compounds ?? "",
                ).length;
                return (
                  <button
                    key={p.accession}
                    onClick={() => toggleProtein(p)}
                    className={`flex w-full flex-wrap items-baseline gap-x-2 gap-y-0.5 border-b border-slate-100 px-3 py-2 text-left text-[12px] transition ${
                      isPick ? "bg-sky-100" : "hover:bg-slate-50"
                    }`}
                  >
                    <span
                      className={`inline-block h-3.5 w-3.5 flex-shrink-0 self-center rounded-sm border transition ${
                        isPick
                          ? "border-sky-600 bg-sky-600"
                          : "border-slate-400 bg-white"
                      }`}
                    >
                      {isPick && (
                        <span className="block text-center text-[9px] leading-[13px] text-white">
                          ✓
                        </span>
                      )}
                    </span>
                    <span className="font-mono font-semibold text-sky-700">
                      {p.accession}
                    </span>
                    {p.is_literature_recruited ? (
                      <span className="text-amber-600" title="newly recruited">
                        ★
                      </span>
                    ) : null}
                    <span className="text-slate-900">
                      {(p.protein_names || "").slice(0, 70)}
                    </span>
                    <span className="text-slate-500">
                      {p.gene_names ? `· ${p.gene_names.slice(0, 18)} ` : ""}·{" "}
                      {(p.organism || "").slice(0, 34)} ·{" "}
                      {clusterText(p, mode)}
                      {nSteroids
                        ? ` · ${nSteroids} steroid${nSteroids !== 1 ? "s" : ""}`
                        : ""}
                      {p.pubmed_count ? ` · ${p.pubmed_count} refs` : ""}
                    </span>
                  </button>
                );
              })}
            </div>
            {pool.length > 300 && (
              <div className="mt-2 text-xs text-slate-500">
                …showing first 300 of {pool.length.toLocaleString()}
              </div>
            )}
          </div>
        )}

        {/* Molecule / natsyn view — steroid tile grid (each is its OWN
            structure). Shows up when a cluster is picked, the user drags to
            select a region, or a search matches. */}
        {(kind === "molecule" || kind === "natsyn") && pool.length > 0 && (
          <div className="mt-6" ref={resultsRef}>
            <div className="mb-3 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                {searchMatchKeys && searchMatchKeys.size
                  ? `🧬 ${pool.length.toLocaleString()} structures matching “${search.trim()}”`
                  : "🧬 Steroids in this selection"}
              </div>
              <div className="text-sm text-slate-500">
                {pool.length.toLocaleString()} compounds
              </div>
            </div>
            <div className="grid gap-3 sm:grid-cols-2 md:grid-cols-3 xl:grid-cols-4">
              {(pool as (Molecule | NatsynEntry)[])
                .slice(0, 60)
                .map((m) => {
                  const nm = m.compound_name;
                  const cb = m.chebi_id;
                  const src = structUrl(nm, cb);
                  const isPick = focusedPoint
                    ? rowKey(kind, focusedPoint) ===
                      rowKey(kind, m as any)
                    : false;
                  return (
                    <button
                      key={rowKey(kind, m as any)}
                      onClick={() => {
                        setFocusedPoint(m as any);
                        scrollToSteroidCard.current = true;
                      }}
                      className={`overflow-hidden rounded-xl border bg-white p-3 text-left text-xs transition ${
                        isPick
                          ? "border-sky-500 ring-2 ring-sky-400"
                          : "border-slate-200 hover:border-slate-300 hover:shadow-md"
                      }`}
                    >
                      <div
                        className="flex items-center justify-center rounded-lg bg-white"
                        style={{ height: 200, padding: 6 }}
                      >
                        {src ? (
                          // eslint-disable-next-line @next/next/no-img-element
                          <img
                            src={src}
                            alt={nm}
                            className="h-full w-full object-contain"
                          />
                        ) : (
                          <span className="text-[11px] text-slate-400">
                            no structure
                          </span>
                        )}
                      </div>
                      <div className="mt-2 line-clamp-2 min-h-[2.4em] text-center text-sm font-medium text-slate-900">
                        {nm}
                      </div>
                      {cb && (
                        <div className="mt-0.5 text-center font-mono text-[10px] text-slate-500">
                          CHEBI:{cb}
                        </div>
                      )}
                    </button>
                  );
                })}
            </div>
            {pool.length > 60 && (
              <div className="mt-2 text-xs text-slate-500">
                …and {pool.length - 60} more compounds
              </div>
            )}
          </div>
        )}

        {/* Molecule / natsyn view — the clicked steroid's detail card sits
            below the whole structure grid. */}
        {focusedPoint && kind !== "protein" && (
          <div ref={steroidCardRef} className="mt-6 scroll-mt-4">
            <SteroidDetailCard
              name={(focusedPoint as Molecule).compound_name}
              chebi={(focusedPoint as Molecule).chebi_id}
              smiles={(focusedPoint as Molecule).smiles}
              structUrl={structUrl}
              accessions={molAccByKey?.get(rowKey(kind, focusedPoint)) ?? []}
              proteinByAcc={proteinByAcc}
              onClose={() => setFocusedPoint(null)}
            />
          </div>
        )}

        {focusedProteins.length > 0 && (
          <div className="mt-6 space-y-4">
            <div className="flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                📄 Detail cards ({focusedProteins.length})
              </div>
              <button
                className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                onClick={clearFocusedProteins}
              >
                ✕ close all
              </button>
            </div>
            {focusedProteins.map((p) => (
              <ProteinDetailCard
                key={p.accession}
                protein={p}
                structUrl={structUrl}
                onClose={() => toggleProtein(p)}
                onPickSteroid={(name, chebi) => setFocusedSteroid({ name, chebi })}
                clusterLabel={clusterText(p, mode)}
              />
            ))}
          </div>
        )}

        {kind === "protein" && focusedSteroid && (
          <div className="mt-6">
            <SteroidDetailCard
              name={focusedSteroid.name}
              chebi={focusedSteroid.chebi}
              structUrl={structUrl}
              proteins={
                proteinsByCompound?.get(
                  `n:${focusedSteroid.name.trim().toLowerCase()}`,
                ) ??
                proteinsByCompound?.get(`c:${chebiNum(focusedSteroid.chebi)}`) ??
                []
              }
              isProteinOpen={isProteinFocused}
              onPickProtein={toggleProtein}
              onClose={() => setFocusedSteroid(null)}
            />
          </div>
        )}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
function ProteinDetailCard({
  protein,
  structUrl,
  onClose,
  onPickSteroid,
  clusterLabel,
}: {
  protein: Protein;
  structUrl: (n: string, c?: string) => string | null;
  onClose: () => void;
  onPickSteroid: (name: string, chebi: string) => void;
  clusterLabel: string;
}) {
  const compounds = splitList(protein.interacting_compounds ?? "");
  const chebis = splitList(protein.interacting_chebi_ids ?? "");
  const rheas = splitList(protein.rhea_reactions ?? "");
  const gos = splitList(protein.go_labels ?? "").slice(0, 8);
  const kws = splitList(protein.keyword_labels ?? "").slice(0, 8);
  const pubmeds = splitList(protein.pubmed_ids ?? "").slice(0, 5);

  return (
    <div className="mt-6 rounded-xl border border-slate-200 bg-white p-5 shadow-md">
      <div className="mb-3 flex flex-wrap items-baseline gap-3">
        <div className="text-lg font-bold text-slate-900">
          {protein.protein_names}
        </div>
        {protein.is_literature_recruited ? (
          <span className="rounded-full bg-amber-100 px-2 py-0.5 text-[10px] font-semibold text-amber-800">
            ★ NEW
          </span>
        ) : null}
        <button
          className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
          onClick={onClose}
        >
          ✕ close
        </button>
      </div>
      <div className="mb-3 text-sm text-slate-600">
        <span className="font-mono text-sky-700">{protein.accession}</span> ·{" "}
        {protein.gene_names || "—"} · <em>{protein.organism}</em> ·{" "}
        {protein.length_aa} aa · {clusterLabel}
      </div>

      <div className="mb-4 flex flex-wrap gap-2 text-sm">
        <a
          href={`https://www.uniprot.org/uniprotkb/${protein.accession}/entry`}
          target="_blank"
          className="rounded border border-sky-500 px-2 py-1 text-sky-700 hover:bg-sky-50"
        >
          UniProt
        </a>
        <a
          href={`https://alphafold.ebi.ac.uk/entry/${protein.accession}`}
          target="_blank"
          className="rounded border border-emerald-500 px-2 py-1 text-emerald-700 hover:bg-emerald-50"
        >
          AlphaFold
        </a>
        {protein.paper_url && protein.paper_url.startsWith("http") && (
          <a
            href={protein.paper_url}
            target="_blank"
            className="rounded border border-amber-500 px-2 py-1 text-amber-700 hover:bg-amber-50"
          >
            Source paper
          </a>
        )}
      </div>

      <div className="grid gap-4 md:grid-cols-2">
        <div>
          {protein.ec_numbers && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">EC:</span>{" "}
              {splitList(protein.ec_numbers.replace(/,/g, ";"))
                .slice(0, 8)
                .map((e) => (
                  <span
                    key={e}
                    className="mr-1 inline-block rounded bg-slate-100 px-1.5 py-0.5 font-mono text-[11px] text-sky-700"
                  >
                    {e}
                  </span>
                ))}
            </div>
          )}
          {rheas.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">Rhea:</span>{" "}
              {rheas.slice(0, 6).map((r) => (
                <a
                  key={r}
                  href={`https://www.rhea-db.org/rhea/${r}`}
                  target="_blank"
                  className="mr-1 inline-block rounded bg-slate-100 px-1.5 py-0.5 font-mono text-[11px] text-sky-700 hover:bg-slate-200"
                >
                  Rhea:{r}
                </a>
              ))}
            </div>
          )}
          {gos.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">GO:</span>{" "}
              {gos.map((g) => (
                <span
                  key={g}
                  className="mr-1 inline-block rounded bg-indigo-50 px-1.5 py-0.5 text-[11px] text-indigo-800"
                >
                  {g.slice(0, 40)}
                </span>
              ))}
            </div>
          )}
          {kws.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">Keywords:</span>{" "}
              {kws.map((k) => (
                <span
                  key={k}
                  className="mr-1 inline-block rounded bg-amber-50 px-1.5 py-0.5 text-[11px] text-amber-800"
                >
                  {k.slice(0, 30)}
                </span>
              ))}
            </div>
          )}
          {pubmeds.length > 0 && (
            <div className="mb-2 text-sm text-slate-700">
              <span className="font-semibold text-slate-900">
                PubMed ({protein.pubmed_count}):
              </span>{" "}
              {pubmeds.map((pm) => (
                <a
                  key={pm}
                  href={`https://pubmed.ncbi.nlm.nih.gov/${pm}`}
                  target="_blank"
                  className="mr-1 inline-block rounded bg-emerald-50 px-1.5 py-0.5 text-[11px] text-emerald-800 hover:bg-emerald-100"
                >
                  PMID {pm}
                </a>
              ))}
            </div>
          )}
        </div>
        <div>
          <div className="mb-2 text-sm font-semibold text-slate-900">
            Interacting steroids ({compounds.length})
          </div>
          <div className="grid grid-cols-3 gap-2">
            {compounds.slice(0, 12).map((nm, i) => {
              const cb = chebis[i] ?? "";
              const src = structUrl(nm, cb);
              return (
                <button
                  key={nm + i}
                  onClick={() => onPickSteroid(nm, cb)}
                  title="Show proteins that interact with this steroid"
                  className="overflow-hidden rounded-md border border-slate-200 bg-white p-2 text-center text-[10px] transition hover:border-sky-400 hover:shadow-md"
                >
                  <div
                    className="flex items-center justify-center"
                    style={{ height: 110 }}
                  >
                    {src ? (
                      // eslint-disable-next-line @next/next/no-img-element
                      <img
                        src={src}
                        alt={nm}
                        className="max-h-full max-w-full object-contain"
                      />
                    ) : (
                      <span className="text-slate-400">no img</span>
                    )}
                  </div>
                  <div className="mt-1 line-clamp-2 min-h-[2.2em] text-slate-700">
                    {nm.slice(0, 40)}
                  </div>
                </button>
              );
            })}
            {compounds.length > 12 && (
              <div className="col-span-3 text-center text-[11px] text-slate-500">
                …and {compounds.length - 12} more
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
/** Steroid detail: structure, ChEBI / SMILES, and the proteins it interacts
 *  with. Pass `proteins` when they're already resolved (protein view), or
 *  `accessions` + `proteinByAcc` (molecule / natsyn views, where the protein
 *  table loads lazily — bare accessions show until it arrives). */
function SteroidDetailCard({
  name,
  chebi,
  smiles,
  structUrl,
  proteins,
  accessions,
  proteinByAcc,
  isProteinOpen,
  onPickProtein,
  onClose,
}: {
  name: string;
  chebi?: string | number | null;
  smiles?: string | null;
  structUrl: (n: string, c?: string) => string | null;
  proteins?: Protein[];
  accessions?: string[];
  proteinByAcc?: Map<string, Protein> | null;
  isProteinOpen?: (acc: string) => boolean;
  onPickProtein?: (p: Protein) => void;
  onClose: () => void;
}) {
  const LIMIT = 50;
  const cb = chebiNum(chebi);
  const src = structUrl(name, cb);
  const list: { acc: string; p?: Protein }[] = proteins
    ? proteins.map((p) => ({ acc: p.accession, p }))
    : (accessions ?? []).map((acc) => ({ acc, p: proteinByAcc?.get(acc) }));
  const loadingNames = !proteins && list.length > 0 && !proteinByAcc;

  return (
    <div className="mt-3 rounded-xl border-2 border-sky-500 bg-white p-5 shadow-md">
      <div className="mb-3 flex flex-wrap items-baseline gap-3">
        <div className="text-xs uppercase tracking-wide text-sky-700">
          Selected steroid
        </div>
        <div className="text-lg font-bold text-slate-900">{name.trim()}</div>
        <button
          className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
          onClick={onClose}
        >
          ✕ close
        </button>
      </div>

      <div className="grid gap-4 md:grid-cols-[220px_1fr]">
        <div>
          <div
            className="flex items-center justify-center rounded-lg border border-slate-200 bg-white"
            style={{ height: 200, padding: 6 }}
          >
            {src ? (
              // eslint-disable-next-line @next/next/no-img-element
              <img
                src={src}
                alt={name}
                className="max-h-full max-w-full object-contain"
              />
            ) : (
              <span className="text-[11px] text-slate-400">no structure</span>
            )}
          </div>
          {cb && (
            <a
              href={`https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:${cb}`}
              target="_blank"
              className="mt-2 inline-block rounded border border-sky-500 px-2 py-1 text-xs text-sky-700 hover:bg-sky-50"
            >
              ChEBI:{cb}
            </a>
          )}
          {smiles && (
            <div
              className="mt-2 break-all font-mono text-[10px] text-slate-500"
              title={smiles}
            >
              {smiles.length > 160 ? `${smiles.slice(0, 160)}…` : smiles}
            </div>
          )}
        </div>

        <div>
          <div className="mb-2 text-sm font-semibold text-slate-900">
            Interacting proteins ({list.length.toLocaleString()})
            {loadingNames && (
              <span className="ml-2 text-xs font-normal text-slate-400">
                loading names…
              </span>
            )}
            {onPickProtein && list.length > 0 && (
              <span className="ml-2 text-xs font-normal text-slate-500">
                click one to open its detail card
              </span>
            )}
          </div>
          {list.length === 0 ? (
            <div className="text-sm text-slate-500">
              No known interacting proteins.
            </div>
          ) : (
            <div className="max-h-[320px] overflow-y-auto rounded-lg border border-slate-200">
              {list.slice(0, LIMIT).map(({ acc, p }) => {
                const open = isProteinOpen?.(acc) ?? false;
                const body = (
                  <>
                    <span className="font-mono font-semibold text-sky-700">
                      {acc}
                    </span>
                    {p && (
                      <>
                        <span className="text-slate-900">
                          {(p.protein_names || "").slice(0, 70)}
                        </span>
                        <span className="text-slate-500">
                          · <em>{(p.organism || "").slice(0, 34)}</em> · cl{" "}
                          {(p.cluster ?? 0) + 1}
                        </span>
                      </>
                    )}
                  </>
                );
                return (
                  <div
                    key={acc}
                    className={`flex items-baseline gap-2 border-b border-slate-100 px-3 py-1.5 text-[12px] ${
                      open ? "bg-sky-100" : ""
                    }`}
                  >
                    {onPickProtein && p ? (
                      <button
                        onClick={() => onPickProtein(p)}
                        className="flex flex-1 flex-wrap items-baseline gap-x-2 text-left hover:underline"
                      >
                        {body}
                      </button>
                    ) : (
                      <div className="flex flex-1 flex-wrap items-baseline gap-x-2">
                        {body}
                      </div>
                    )}
                    <a
                      href={`https://www.uniprot.org/uniprotkb/${acc}/entry`}
                      target="_blank"
                      className="flex-shrink-0 text-[11px] text-sky-700 hover:underline"
                    >
                      UniProt ↗
                    </a>
                  </div>
                );
              })}
            </div>
          )}
          {list.length > LIMIT && (
            <div className="mt-1 text-xs text-slate-500">
              …and {(list.length - LIMIT).toLocaleString()} more
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
