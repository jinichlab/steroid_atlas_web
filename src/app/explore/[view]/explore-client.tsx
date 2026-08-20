"use client";

import { useEffect, useMemo, useState } from "react";
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

const UmapPlot = dynamic(() => import("@/components/UmapPlot"), { ssr: false });

type Row = Protein | Molecule | NatsynEntry;

interface Props {
  kind: ViewKind;
  clusters: ClusterMeta[];
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
function tooltipHtmlFor(kind: ViewKind, r: Row): string {
  if (kind === "protein") {
    const p = r as Protein;
    const name = esc((p.protein_names ?? "").slice(0, 90));
    const gene = esc(p.gene_names ?? "");
    const org = esc(p.organism ?? "");
    const acc = esc(p.accession ?? "");
    const ecs = esc((p.ec_numbers ?? "").split(";").slice(0, 2).join(" · "));
    const cluster = ((p.cluster ?? 0) + 1).toString();
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
          Cluster ${cluster}${ecs ? ` · EC ${ecs}` : ""}
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

type StructIndex = Record<string, string>;

export default function ExploreClient({ kind, clusters }: Props) {
  const [rows, setRows] = useState<Row[]>([]);
  const [loading, setLoading] = useState(true);
  const [structIdx, setStructIdx] = useState<StructIndex>({});
  const [highlight, setHighlight] = useState<number | null>(null);
  const [selectionSet, setSelectionSet] = useState<Set<string> | null>(null);
  const [focusedPoint, setFocusedPoint] = useState<Row | null>(null);
  const [focusedCompound, setFocusedCompound] = useState<string | null>(null);
  // Multi-select for the bottom detail cards — clicking a row in the
  // per-compound list toggles that protein in the set; the "See full
  // detail" button on the focused-point card adds that one protein too.
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
  const [legendOpen, setLegendOpen] = useState(true);
  // Free-text search across names / accessions / gene / GO / EC / keywords /
  // ChEBI / compound / sequence.
  const [search, setSearch] = useState("");
  // A bumping counter that tells UmapPlot to erase its sticky rectangle
  // when the user hits the top-of-plot "clear" button.
  const [lassoClearSignal, setLassoClearSignal] = useState(0);

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

  const structUrl = (name: string, chebi?: string): string | null => {
    const nm = (name || "").toLowerCase();
    if (nm && structIdx[nm]) return `/atlas/structures/${structIdx[nm]}`;
    if (chebi && structIdx[`chebi:${chebi}`])
      return `/atlas/structures/${structIdx[`chebi:${chebi}`]}`;
    return null;
  };

  const points: UmapPoint[] = useMemo(
    () =>
      rows.map((r) => {
        const raw = (r as any).cluster;
        const clusterNum =
          typeof raw === "number"
            ? raw
            : typeof raw === "string"
              ? raw === "natural"
                ? 0
                : 1
              : 0;
        return {
          x: (r as any).umap_1,
          y: (r as any).umap_2,
          cluster: clusterNum,
          is_star: (r as any).is_literature_recruited ?? 0,
          label: rowLabel(kind, r),
          key: rowKey(kind, r),
          tooltipHtml: tooltipHtmlFor(kind, r),
        };
      }),
    [rows, kind],
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

  // Cluster centroids drive on-plot labels AND the distinct-color palette
  // (so neighbors on the map don't end up as similar shades of green).
  const centroids = useMemo(() => {
    const byCid = new Map<number, { sx: number; sy: number; n: number }>();
    for (const p of points) {
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
    return distinctPaletteForCentroids(centroids);
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

  const pool = useMemo(() => {
    if (searchMatchKeys && searchMatchKeys.size) {
      return rows.filter((r) => searchMatchKeys.has(rowKey(kind, r)));
    }
    if (selectionSet && selectionSet.size) {
      return rows.filter((r) => selectionSet.has(rowKey(kind, r)));
    }
    if (highlight != null) {
      return rows.filter((r) => Number((r as any).cluster) === highlight);
    }
    return [];
  }, [rows, highlight, selectionSet, searchMatchKeys, kind]);

  const compoundTiles = useMemo(() => {
    if (kind !== "protein" || !pool.length) return [];
    const byCompound = new Map<string, { chebis: Set<string>; prots: Protein[] }>();
    for (const r of pool as Protein[]) {
      const names = splitList(r.interacting_compounds ?? "");
      const chebis = splitList(r.interacting_chebi_ids ?? "");
      names.forEach((nm, i) => {
        const cb = chebis[i] ?? "";
        const entry = byCompound.get(nm);
        if (entry) {
          entry.prots.push(r);
          if (cb) entry.chebis.add(cb);
        } else {
          byCompound.set(nm, {
            prots: [r],
            chebis: new Set(cb ? [cb] : []),
          });
        }
      });
    }
    return Array.from(byCompound.entries())
      .map(([name, v]) => ({
        name,
        chebi: v.chebis.values().next().value ?? "",
        prots: v.prots,
      }))
      .sort((a, b) => b.prots.length - a.prots.length);
  }, [pool, kind]);

  const focusedCompoundProts = useMemo(() => {
    if (!focusedCompound) return [];
    const t = compoundTiles.find((t) => t.name === focusedCompound);
    return t?.prots ?? [];
  }, [compoundTiles, focusedCompound]);

  useEffect(() => {
    clearFocusedProteins();
  }, [focusedCompound, highlight]);

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
            className="w-full max-w-2xl rounded-md border border-slate-300 bg-white px-3 py-1.5 text-sm shadow-sm focus:border-sky-500 focus:outline-none focus:ring-2 focus:ring-sky-100"
          />
          {search && (
            <div className="text-xs text-slate-500">
              {searchMatchKeys
                ? `${searchMatchKeys.size.toLocaleString()} matches`
                : ""}
              <button
                className="ml-2 rounded border border-slate-300 bg-white px-2 py-0.5 text-slate-700 hover:bg-slate-100"
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
          {(highlight != null || focusedPoint || selectionSet || search) && (
            <button
              className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-slate-700 hover:bg-slate-100"
              onClick={() => {
                setHighlight(null);
                setSelectionSet(null);
                setFocusedPoint(null);
                setFocusedCompound(null);
                clearFocusedProteins();
                setSearch("");
                setLassoClearSignal((n) => n + 1);
              }}
            >
              ✕ clear
            </button>
          )}
        </div>

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
                    : null
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
                setFocusedCompound(null);
                clearFocusedProteins();
              }}
              height={640}
            />
            {focusedPoint && (
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
                    {(focusedPoint as Protein).length_aa} aa · Cluster{" "}
                    {((focusedPoint as Protein).cluster ?? 0) + 1}
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
                        setFocusedCompound(null);
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

        {/* Steroid catalogue */}
        {kind === "protein" && compoundTiles.length > 0 && (
          <div className="mt-6">
            <div className="mb-3 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                🧬 Steroids in this pool
              </div>
              <div className="text-sm text-slate-500">
                {compoundTiles.length.toLocaleString()} compounds across{" "}
                {pool.length.toLocaleString()} proteins
              </div>
              {focusedCompound && (
                <button
                  className="ml-auto rounded-md border border-slate-300 bg-white px-2 py-1 text-xs text-slate-700 hover:bg-slate-100"
                  onClick={() => setFocusedCompound(null)}
                >
                  ✕ clear steroid focus
                </button>
              )}
            </div>
            <div className="grid gap-3 sm:grid-cols-2 md:grid-cols-3 xl:grid-cols-4">
              {compoundTiles.slice(0, 60).map((t) => {
                const isPick = focusedCompound === t.name;
                const src = structUrl(t.name, t.chebi);
                return (
                  <button
                    key={t.name}
                    onClick={() =>
                      setFocusedCompound(isPick ? null : t.name)
                    }
                    className={`overflow-hidden rounded-xl border bg-white p-3 text-left text-xs transition ${
                      isPick
                        ? "border-sky-500 ring-2 ring-sky-400"
                        : "border-slate-200 hover:border-slate-300 hover:shadow-md"
                    }`}
                  >
                    <div className="flex items-center justify-center rounded-lg bg-white"
                         style={{ height: 210, padding: 6 }}>
                      {src ? (
                        // eslint-disable-next-line @next/next/no-img-element
                        <img
                          src={src}
                          alt={t.name}
                          className="h-full w-full object-contain"
                          style={{ imageRendering: "auto" }}
                        />
                      ) : (
                        <span className="text-[11px] text-slate-400">
                          no structure
                        </span>
                      )}
                    </div>
                    <div className="mt-2 line-clamp-2 min-h-[2.4em] text-sm font-medium text-slate-900 text-center">
                      {t.name}
                    </div>
                    {t.chebi && (
                      <div className="mt-0.5 text-center font-mono text-[10px] text-slate-500">
                        CHEBI:{t.chebi}
                      </div>
                    )}
                    <div className="mt-1 flex justify-center">
                      <span className="rounded-full bg-indigo-50 px-2 py-0.5 text-[11px] font-medium text-indigo-700">
                        {t.prots.length.toLocaleString()} protein
                        {t.prots.length !== 1 ? "s" : ""}
                      </span>
                    </div>
                    {isPick && (
                      <div className="mt-1 text-center text-[11px] font-semibold text-sky-700">
                        ◉ FOCUSED
                      </div>
                    )}
                  </button>
                );
              })}
            </div>
            {compoundTiles.length > 60 && (
              <div className="mt-2 text-xs text-slate-500">
                …and {compoundTiles.length - 60} more compounds
              </div>
            )}

            {focusedCompound && focusedCompoundProts.length > 0 && (
              <div className="mt-4 rounded-xl border-2 border-sky-400 bg-sky-50 p-4">
                <div className="mb-2 flex flex-wrap items-baseline gap-2 text-sm text-sky-800">
                  <b className="text-slate-900">{focusedCompound}</b>
                  <span className="text-slate-500">
                    ({focusedCompoundProts.length.toLocaleString()})
                  </span>
                  <span className="text-slate-600">
                    — tick any protein below to add its full detail card
                  </span>
                  {focusedProteins.length > 0 && (
                    <button
                      className="ml-auto rounded border border-slate-300 bg-white px-2 py-0.5 text-xs text-slate-700 hover:bg-slate-100"
                      onClick={clearFocusedProteins}
                    >
                      ✕ clear all {focusedProteins.length} details
                    </button>
                  )}
                </div>
                <div className="max-h-[380px] overflow-y-auto rounded-md border border-slate-200 bg-white">
                  {focusedCompoundProts.slice(0, 200).map((p) => {
                    const isPick = isProteinFocused(p.accession);
                    return (
                      <button
                        key={p.accession}
                        onClick={() => toggleProtein(p)}
                        className={`flex w-full flex-wrap items-baseline gap-2 border-b border-slate-100 px-3 py-2 text-left text-[12px] transition ${
                          isPick ? "bg-sky-100" : "hover:bg-slate-50"
                        }`}
                      >
                        <span
                          className={`inline-block h-3.5 w-3.5 flex-shrink-0 rounded-sm border transition ${
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
                        <span className="text-slate-900">
                          {(p.protein_names || "").slice(0, 60)}
                        </span>
                        <span className="text-slate-500">
                          {p.gene_names ? `· ${p.gene_names.slice(0, 20)}` : ""}{" "}
                          · {(p.organism || "").slice(0, 40)}
                          {p.pubmed_count ? ` · ${p.pubmed_count} refs` : ""}
                        </span>
                      </button>
                    );
                  })}
                </div>
                {focusedCompoundProts.length > 200 && (
                  <div className="mt-2 text-xs text-slate-500">
                    …showing first 200 of {focusedCompoundProts.length}
                  </div>
                )}
              </div>
            )}
          </div>
        )}

        {/* Molecule / natsyn view — steroid tile grid (each is its OWN
            structure). Shows up when a cluster is picked or the user
            drags to select a region. */}
        {(kind === "molecule" || kind === "natsyn") && pool.length > 0 && (
          <div className="mt-6">
            <div className="mb-3 flex flex-wrap items-baseline gap-3">
              <div className="text-lg font-semibold text-slate-900">
                🧬 Steroids in this selection
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
                      onClick={() => setFocusedPoint(m as any)}
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
              />
            ))}
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
}: {
  protein: Protein;
  structUrl: (n: string, c?: string) => string | null;
  onClose: () => void;
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
        {protein.length_aa} aa · Cluster {(protein.cluster ?? 0) + 1}
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
                <div
                  key={nm + i}
                  className="overflow-hidden rounded-md border border-slate-200 bg-white p-2 text-center text-[10px]"
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
                </div>
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
