"use client";

import { useMemo, useRef, useState, useCallback, useEffect } from "react";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer, IconLayer, TextLayer } from "@deck.gl/layers";
import { OrthographicView, PickingInfo } from "@deck.gl/core";
import { clusterRgb } from "@/lib/palette";
import { MousePointer2, Hand } from "lucide-react";

export interface UmapPoint {
  x: number;
  y: number;
  cluster: number;
  is_star?: number;
  label?: string;
  key?: string;
  /** Optional per-point HTML block shown on hover. Overrides `label`. */
  tooltipHtml?: string;
}

export interface ClusterCentroidLabel {
  cluster: number;
  display_id: number;
  name: string;
  n: number;
}

interface UmapPlotProps {
  points: UmapPoint[];
  clusterLabels?: ClusterCentroidLabel[];
  highlightCluster?: number | null;
  /** When set (non-null), only points whose key is in the set stay colored;
   *  everyone else greys out. Used for free-text search matches + lasso
   *  selection. */
  matchKeys?: Set<string> | null;
  focusedKey?: string | null;
  /** When true, the last drawn lasso rectangle stays on the plot until the
   *  next drag or an external clear. Defaults to true. */
  lassoPersistent?: boolean;
  /** Set to true to erase any drawn lasso — bump this counter on external
   *  "clear" events so the effect fires. */
  clearLassoSignal?: number;
  onClickPoint?: (p: UmapPoint) => void;
  onLassoSelect?: (points: UmapPoint[]) => void;
  height?: number;
  darkMode?: boolean;
  paletteOverride?: Map<number, [number, number, number]>;
}

// White star SVG, tinted per-instance via getColor (mask: true).
const STAR_SVG = `data:image/svg+xml;utf8,${encodeURIComponent(
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="-1 -1 2 2" width="32" height="32">` +
    `<path d="M0,-0.95 L0.22,-0.31 L0.9,-0.31 L0.36,0.12 L0.56,0.78 L0,0.38 L-0.56,0.78 L-0.36,0.12 L-0.9,-0.31 L-0.22,-0.31 Z" ` +
    `fill="#ffffff"/></svg>`,
)}`;
const ICON_MAPPING = {
  star: { x: 0, y: 0, width: 32, height: 32, anchorY: 16, anchorX: 16, mask: true },
} as const;

export default function UmapPlot({
  points,
  clusterLabels = [],
  highlightCluster = null,
  matchKeys = null,
  focusedKey = null,
  lassoPersistent = true,
  clearLassoSignal = 0,
  onClickPoint,
  onLassoSelect,
  height = 620,
  darkMode = true,
  paletteOverride,
}: UmapPlotProps) {
  // ── Bounds → initial view ─────────────────────────────────────────────
  const bounds = useMemo(() => {
    if (!points.length) return { cx: 0, cy: 0, zoom: 6 };
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const p of points) {
      if (p.x < minX) minX = p.x;
      if (p.x > maxX) maxX = p.x;
      if (p.y < minY) minY = p.y;
      if (p.y > maxY) maxY = p.y;
    }
    const cx = (minX + maxX) / 2;
    const cy = (minY + maxY) / 2;
    const span = Math.max(maxX - minX, maxY - minY) || 1;
    const zoom = Math.log2(900 / (span * 1.15));
    return { cx, cy, zoom };
  }, [points]);

  const [viewState, setViewState] = useState({
    target: [0, 0, 0] as [number, number, number],
    zoom: 6,
  });
  useEffect(() => {
    setViewState({ target: [bounds.cx, bounds.cy, 0], zoom: bounds.zoom });
  }, [bounds]);

  // ── Mode toggle: pan vs select ─────────────────────────────────────────
  const [selectMode, setSelectMode] = useState(false);

  // ── Palette ─────────────────────────────────────────────────────────────
  const nClusters = useMemo(
    () => Math.max(...points.map((p) => p.cluster), 0) + 1,
    [points],
  );
  const rgbFallback = useMemo(() => clusterRgb(Math.max(nClusters, 2)), [nClusters]);
  const colorFor = useCallback(
    (cid: number): [number, number, number] => {
      if (paletteOverride) {
        const c = paletteOverride.get(cid);
        if (c) return c;
      }
      return rgbFallback[cid] ?? [30, 100, 180];
    },
    [paletteOverride, rgbFallback],
  );

  // ── Lasso (only active when select-mode is on) ────────────────────────
  // Two rectangles: the live one being drawn (`liveLasso`) and the last
  // committed one (`stickyLasso`) that persists on-plot until the next
  // drag or an external clear signal.
  const [liveLasso, setLiveLasso] = useState<
    { x0: number; y0: number; x1: number; y1: number } | null
  >(null);
  const [stickyLasso, setStickyLasso] = useState<
    { x0: number; y0: number; x1: number; y1: number } | null
  >(null);
  const dragging = useRef(false);
  const deckRef = useRef<DeckGL | null>(null);

  // External clear signal — parent bumps this counter to erase the sticky.
  useEffect(() => {
    setStickyLasso(null);
    setLiveLasso(null);
  }, [clearLassoSignal]);

  const startLasso = useCallback(
    (e: React.PointerEvent<HTMLDivElement>) => {
      if (!selectMode) return;
      if (e.button !== 0) return;
      const rect = (e.currentTarget as HTMLDivElement).getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;
      setLiveLasso({ x0: x, y0: y, x1: x, y1: y });
      setStickyLasso(null); // starting a new drag clears the old sticky
      dragging.current = true;
    },
    [selectMode],
  );
  const moveLasso = useCallback((e: React.PointerEvent<HTMLDivElement>) => {
    if (!dragging.current) return;
    const rect = (e.currentTarget as HTMLDivElement).getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    setLiveLasso((l) => (l ? { ...l, x1: x, y1: y } : null));
  }, []);
  const endLasso = useCallback(
    (_e: React.PointerEvent<HTMLDivElement>) => {
      if (!dragging.current || !liveLasso) {
        dragging.current = false;
        setLiveLasso(null);
        return;
      }
      dragging.current = false;
      const x0 = Math.min(liveLasso.x0, liveLasso.x1);
      const y0 = Math.min(liveLasso.y0, liveLasso.y1);
      const x1 = Math.max(liveLasso.x0, liveLasso.x1);
      const y1 = Math.max(liveLasso.y0, liveLasso.y1);
      const w = x1 - x0;
      const h = y1 - y0;
      if (w * h < 25) {
        setLiveLasso(null);
        return;
      }
      // Persist the rectangle on the plot for reference until the next drag.
      if (lassoPersistent) setStickyLasso({ x0, y0, x1, y1 });
      setLiveLasso(null);
      if (!deckRef.current || !onLassoSelect) return;
      const picks = (deckRef.current as any).pickObjects({
        x: x0,
        y: y0,
        width: w,
        height: h,
        layerIds: ["scatter", "stars"],
      });
      const seen = new Set<string>();
      const out: UmapPoint[] = [];
      for (const info of picks) {
        const obj = (info as any).object as UmapPoint | undefined;
        if (!obj) continue;
        const k = obj.key ?? `${obj.x},${obj.y},${obj.cluster}`;
        if (!seen.has(k)) {
          seen.add(k);
          out.push(obj);
        }
      }
      onLassoSelect(out);
    },
    [liveLasso, onLassoSelect, lassoPersistent],
  );

  // ── Cluster centroids for label layer ─────────────────────────────────
  const centroids = useMemo(() => {
    if (!points.length || !clusterLabels.length) return [];
    const byCid = new Map<number, number[]>();
    points.forEach((p, i) => {
      const arr = byCid.get(p.cluster);
      if (arr) arr.push(i);
      else byCid.set(p.cluster, [i]);
    });
    const nameById = new Map<number, ClusterCentroidLabel>();
    clusterLabels.forEach((c) => nameById.set(c.cluster, c));
    const out: { x: number; y: number; text: string; cluster: number }[] = [];
    const MIN_N = 30;
    for (const [cid, idxs] of byCid) {
      if (idxs.length < MIN_N) continue;
      const xs = idxs.map((i) => points[i].x).sort((a, b) => a - b);
      const ys = idxs.map((i) => points[i].y).sort((a, b) => a - b);
      const mid = Math.floor(idxs.length / 2);
      const meta = nameById.get(cid);
      if (!meta) continue;
      // Show the FULL family name (with a small safety cap so runaway
      // labels don't overflow the plot). Marimo used ≤60 chars.
      const short = (meta.name || "").slice(0, 60);
      const text = short
        ? `${meta.display_id} · ${short}`
        : String(meta.display_id);
      out.push({ x: xs[mid], y: ys[mid], text, cluster: cid });
    }
    return out;
  }, [points, clusterLabels]);
  // Show cluster labels only once the user has zoomed in noticeably.
  const showLabels = viewState.zoom > bounds.zoom + 1.4;

  // ── Layers ──────────────────────────────────────────────────────────────
  const scatterPoints = useMemo(() => points.filter((p) => !p.is_star), [points]);
  const starPoints = useMemo(() => points.filter((p) => p.is_star), [points]);

  const layers = useMemo(() => {
    const highlight = highlightCluster;
    const grey: [number, number, number, number] = darkMode
      ? [80, 100, 130, 90]
      : [210, 220, 230, 90];
    const scatter = new ScatterplotLayer<UmapPoint>({
      id: "scatter",
      data: scatterPoints,
      getPosition: (d: UmapPoint) => [d.x, d.y, 0],
      getFillColor: (d: UmapPoint) => {
        if (matchKeys && d.key && !matchKeys.has(d.key)) return grey;
        if (highlight != null && d.cluster !== highlight) return grey;
        const c = colorFor(d.cluster);
        return [c[0], c[1], c[2], 220];
      },
      getRadius: () => 2.6,
      getLineColor: (d: UmapPoint) =>
        focusedKey && d.key === focusedKey
          ? [224, 20, 76]
          : darkMode
            ? [15, 23, 42]
            : [255, 255, 255],
      getLineWidth: (d: UmapPoint) =>
        focusedKey && d.key === focusedKey ? 2 : 0.4,
      lineWidthUnits: "pixels",
      radiusUnits: "pixels",
      stroked: true,
      pickable: true,
      onClick: (info: PickingInfo) => {
        if (selectMode) return true; // ignore clicks in select mode
        const p = info.object as UmapPoint | undefined;
        if (p && onClickPoint) onClickPoint(p);
        return true;
      },
      updateTriggers: {
        getFillColor: [highlightCluster, matchKeys, colorFor, darkMode],
        getLineColor: [focusedKey, darkMode],
        getLineWidth: [focusedKey],
      },
    });
    const stars = new IconLayer<UmapPoint>({
      id: "stars",
      data: starPoints,
      getPosition: (d: UmapPoint) => [d.x, d.y, 0],
      iconAtlas: STAR_SVG,
      iconMapping: ICON_MAPPING,
      getIcon: () => "star",
      sizeUnits: "pixels",
      getSize: () => 22,
      getColor: (d: UmapPoint) => {
        if (matchKeys && d.key && !matchKeys.has(d.key)) return grey;
        if (highlight != null && d.cluster !== highlight) return grey;
        const c = colorFor(d.cluster);
        return [c[0], c[1], c[2], 255];
      },
      pickable: true,
      onClick: (info: PickingInfo) => {
        if (selectMode) return true;
        const p = info.object as UmapPoint | undefined;
        if (p && onClickPoint) onClickPoint(p);
        return true;
      },
      updateTriggers: {
        getColor: [highlightCluster, matchKeys, colorFor, darkMode],
      },
    });
    const labelData = showLabels ? centroids : [];
    const labels = new TextLayer<{ x: number; y: number; text: string; cluster: number }>({
      id: "cluster-labels",
      data: labelData,
      getPosition: (d) => [d.x, d.y, 0],
      getText: (d) => d.text,
      getSize: 13,
      getColor: darkMode ? [255, 255, 255, 240] : [15, 23, 42, 240],
      background: true,
      getBackgroundColor: darkMode ? [15, 23, 42, 200] : [255, 255, 255, 230],
      backgroundPadding: [4, 2, 4, 2],
      getBorderColor: [148, 163, 184, 180],
      getBorderWidth: 0.5,
      fontFamily: "-apple-system,BlinkMacSystemFont,sans-serif",
      fontWeight: 600,
      sizeUnits: "pixels",
      getTextAnchor: "middle",
      getAlignmentBaseline: "center",
      pickable: false,
    });
    return [scatter, stars, labels];
  }, [
    scatterPoints,
    starPoints,
    colorFor,
    highlightCluster,
    matchKeys,
    focusedKey,
    onClickPoint,
    darkMode,
    centroids,
    showLabels,
    selectMode,
  ]);

  const bg = darkMode ? "bg-slate-900" : "bg-white";
  const border = darkMode ? "border-slate-700" : "border-slate-200";

  return (
    <div
      className={`relative overflow-hidden rounded-xl border ${border} ${bg} shadow-sm`}
      style={{ height }}
      onPointerDown={startLasso}
      onPointerMove={moveLasso}
      onPointerUp={endLasso}
      onPointerCancel={endLasso}
    >
      <DeckGL
        ref={deckRef}
        views={new OrthographicView({ flipY: false })}
        viewState={viewState}
        controller={{
          dragPan: !selectMode, // in select mode, we handle the drag as lasso
          dragRotate: false,
          scrollZoom: true,
          doubleClickZoom: true,
        }}
        onViewStateChange={(e: any) => setViewState(e.viewState)}
        layers={layers}
        style={{ position: "absolute", inset: 0 }}
        getCursor={({ isDragging }: any) =>
          selectMode ? "crosshair" : isDragging ? "grabbing" : "grab"
        }
        getTooltip={({ object }: any) => {
          const p = object as UmapPoint | undefined;
          if (!p) return null;
          const bg = darkMode ? "#0f172a" : "#ffffff";
          const fg = darkMode ? "#f1f5f9" : "#111827";
          const bd = darkMode ? "#334155" : "#cbd5e1";
          const inner =
            p.tooltipHtml ??
            `<div style="font-size:11.5px;line-height:1.35;max-width:320px;
              white-space:normal;">${(p.label ?? "").replace(/</g, "&lt;")}</div>`;
          return {
            html: inner,
            style: {
              background: bg,
              color: fg,
              border: `1px solid ${bd}`,
              borderRadius: "6px",
              padding: "8px 10px",
              boxShadow: "0 2px 8px rgba(0,0,0,0.18)",
              pointerEvents: "none",
              maxWidth: "360px",
            },
          };
        }}
      />
      {stickyLasso && !liveLasso && (
        <div
          className="pointer-events-none absolute border-2 border-dashed border-sky-500 bg-sky-200/15"
          style={{
            left: stickyLasso.x0,
            top: stickyLasso.y0,
            width: stickyLasso.x1 - stickyLasso.x0,
            height: stickyLasso.y1 - stickyLasso.y0,
          }}
        />
      )}
      {liveLasso && (
        <div
          className="pointer-events-none absolute border-2 border-sky-500 bg-sky-200/30"
          style={{
            left: Math.min(liveLasso.x0, liveLasso.x1),
            top: Math.min(liveLasso.y0, liveLasso.y1),
            width: Math.abs(liveLasso.x1 - liveLasso.x0),
            height: Math.abs(liveLasso.y1 - liveLasso.y0),
          }}
        />
      )}

      {/* Pan / select toggle */}
      <div className="absolute top-2 left-2 flex gap-1 rounded-md border border-slate-300 bg-white/95 p-0.5 shadow-sm">
        <button
          title="Pan mode — drag pans the map"
          onClick={() => setSelectMode(false)}
          className={`flex items-center gap-1 rounded px-2 py-1 text-[11px] transition ${
            !selectMode
              ? "bg-sky-600 text-white"
              : "text-slate-700 hover:bg-slate-100"
          }`}
        >
          <Hand size={12} /> Pan
        </button>
        <button
          title="Select mode — drag draws a rectangle to select"
          onClick={() => setSelectMode(true)}
          className={`flex items-center gap-1 rounded px-2 py-1 text-[11px] transition ${
            selectMode
              ? "bg-sky-600 text-white"
              : "text-slate-700 hover:bg-slate-100"
          }`}
        >
          <MousePointer2 size={12} /> Select
        </button>
      </div>

      <div
        className={`pointer-events-none absolute bottom-2 right-2 rounded-md px-3 py-1 text-[11px] shadow-sm ${
          darkMode
            ? "bg-slate-800/85 text-slate-200"
            : "bg-white/90 text-slate-600"
        }`}
      >
        {selectMode
          ? "drag to draw a select rectangle · scroll to zoom"
          : "drag to pan · scroll to zoom in for cluster labels · click a dot"}
      </div>
    </div>
  );
}
