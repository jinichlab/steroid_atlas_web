"""Extended k-means k-selection sweep on the current 35,349-protein atlas.

The original silhouette sweep in refit_full_umap.py was hard-capped at k=100
and returned argmax at k=95 — right at the boundary. This script re-runs the
sweep over k in [10, 200] on the deployed atlas coordinates to check whether
silhouette peaks beyond 100 or plateaus.

Reads:  ../data/proteins.csv   (uses the umap_1, umap_2 columns)
Writes: ../analysis/kselection_extended.tsv           per-k silhouette + inertia + Davies-Bouldin
        ../analysis/kselection_extended_report.txt    human-readable summary
"""
from __future__ import annotations

import time
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
OUT_TSV = output("kselection_extended.tsv")
OUT_REPORT = output("kselection_extended_report.txt")

K_MIN, K_MAX = 10, 200
SIL_SAMPLE = 5000            # subsample for silhouette computation (expensive on all pairs)
RANDOM_STATE = 42


def main() -> int:
    print(f"Loading {IN}...")
    df = pd.read_csv(IN, low_memory=False)
    coords = df[["umap_1", "umap_2"]].dropna().to_numpy()
    n = len(coords)
    print(f"  {n:,} points  (dropped {len(df)-n} rows with missing UMAP)")

    sample_size = min(SIL_SAMPLE, n)
    print(f"\nSweeping k in [{K_MIN}, {K_MAX}]  (silhouette sample = {sample_size:,})")
    print(f"{'k':>5}{'silhouette':>14}{'inertia':>16}{'davies-bouldin':>18}{'t(s)':>7}")

    records = []
    for k in range(K_MIN, K_MAX + 1):
        t0 = time.time()
        km = KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=10).fit(coords)
        sil = silhouette_score(coords, km.labels_,
                               sample_size=sample_size, random_state=RANDOM_STATE)
        inertia = km.inertia_
        # Davies-Bouldin uses all points (fast enough here since 2D)
        db = davies_bouldin_score(coords, km.labels_)
        dt = time.time() - t0
        records.append({"k": k, "silhouette": sil, "inertia": inertia, "davies_bouldin": db})
        print(f"{k:>5}{sil:>14.4f}{inertia:>16.0f}{db:>18.4f}{dt:>7.1f}")

    out = pd.DataFrame(records)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    best_sil_k = int(out.loc[out["silhouette"].idxmax(), "k"])
    best_sil = float(out["silhouette"].max())
    best_db_k = int(out.loc[out["davies_bouldin"].idxmin(), "k"])
    best_db = float(out["davies_bouldin"].min())

    lines = []
    lines.append("=== Extended k-selection sweep on 35,349-protein atlas UMAP ===\n")
    lines.append(f"Search range: k ∈ [{K_MIN}, {K_MAX}]")
    lines.append(f"Silhouette computed on {sample_size:,} sampled points per k")
    lines.append(f"Davies-Bouldin computed on all points")
    lines.append("")
    lines.append(f"Silhouette argmax:      k = {best_sil_k}   (silhouette = {best_sil:.4f})")
    lines.append(f"Davies-Bouldin argmin:  k = {best_db_k}   (D-B = {best_db:.4f})")
    lines.append("")
    lines.append("Silhouette values near the argmax (compare with the currently deployed k=95):")
    subset = out[(out["k"] >= max(K_MIN, best_sil_k - 15)) & (out["k"] <= min(K_MAX, best_sil_k + 15))]
    for _, r in subset.iterrows():
        marker = "  ← current k=95" if int(r["k"]) == 95 else ("  ← argmax" if int(r["k"]) == best_sil_k else "")
        lines.append(f"  k={int(r['k']):>4}  sil={r['silhouette']:.4f}  DB={r['davies_bouldin']:.4f}{marker}")
    OUT_REPORT.write_text("\n".join(lines) + "\n")

    print()
    print(f"Silhouette argmax:      k = {best_sil_k}  (silhouette = {best_sil:.4f})")
    print(f"Davies-Bouldin argmin:  k = {best_db_k}   (D-B = {best_db:.4f})")
    print(f"Currently deployed:     k = 95           (silhouette = {out.loc[out['k']==95, 'silhouette'].iloc[0]:.4f})")
    print()
    print(f"Wrote {OUT_TSV.name}, {OUT_REPORT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
