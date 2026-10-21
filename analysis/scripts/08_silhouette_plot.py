"""Publication figure — silhouette + Davies-Bouldin vs k.

Reads the extended k-selection sweep and produces a two-panel figure that
justifies the deployed k choice.  The paper narrative: silhouette climbs
modestly beyond k≈95 (∼10% relative gain from k=95 → argmax at k=185) but
Davies-Bouldin keeps improving through the whole tested range — we chose
k=95 for interpretability and disclose the sweep transparently.

Reads:  ../analysis/kselection_extended.tsv    (produced by 03_kmeans_kselection_extended.py)
Writes: ../analysis/silhouette_vs_k.png        300 dpi PNG for the paper
        ../analysis/silhouette_vs_k.pdf        vector PDF for typesetting
"""
from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = output("kselection_extended.tsv")
OUT_PNG = output("silhouette_vs_k.png")
OUT_PDF = output("silhouette_vs_k.pdf")

DEPLOYED_K = 95   # chosen for interpretability
LINE = "#1F4B99"  # deep blue — magnitude plot uses one hue
DEPLOYED = "#0E7490"
ARGMAX = "#D97706"


def main() -> int:
    df = pd.read_csv(IN, sep="\t")
    df = df.sort_values("k").reset_index(drop=True)

    argmax_k = int(df.loc[df["silhouette"].idxmax(), "k"])
    argmax_sil = float(df["silhouette"].max())
    deployed_sil = float(df.loc[df["k"] == DEPLOYED_K, "silhouette"].iloc[0])
    argmin_db_k = int(df.loc[df["davies_bouldin"].idxmin(), "k"])
    argmin_db = float(df["davies_bouldin"].min())
    deployed_db = float(df.loc[df["k"] == DEPLOYED_K, "davies_bouldin"].iloc[0])

    # ── Figure layout ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    fig.suptitle(
        f"k-means selection on the 2D UMAP  (k ∈ [{int(df['k'].min())}, {int(df['k'].max())}], "
        f"n = 35,117 proteins)",
        fontsize=12, y=1.03,
    )

    for ax in axes:
        ax.set_xlabel("Number of clusters k", fontsize=10)
        ax.grid(True, alpha=0.25, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(25))
        ax.tick_params(labelsize=9)

    # ── Panel A: silhouette (higher = better) ────────────────────────────────
    ax = axes[0]
    ax.plot(df["k"], df["silhouette"], color=LINE, linewidth=1.6, zorder=2)
    ax.set_ylabel("Silhouette score (higher = better)", fontsize=10)
    ax.set_title("A. Silhouette", fontsize=11, loc="left", pad=8)

    # Give ourselves headroom so annotations don't collide with the curve
    _y_lo, _y_hi = ax.get_ylim()
    ax.set_ylim(_y_lo - 0.03, _y_hi + 0.04)

    # Mark deployed k=95 — annotation goes BELOW-RIGHT into whitespace under the curve
    ax.scatter([DEPLOYED_K], [deployed_sil], s=55, color=DEPLOYED,
               edgecolor="white", linewidth=1.4, zorder=3)
    ax.annotate(
        f"k = {DEPLOYED_K}  (deployed)\nsilhouette = {deployed_sil:.3f}",
        xy=(DEPLOYED_K, deployed_sil), xytext=(DEPLOYED_K + 8, deployed_sil - 0.09),
        fontsize=8.5, color=DEPLOYED, ha="left",
        arrowprops=dict(arrowstyle="-", color=DEPLOYED, lw=0.7),
    )

    # Mark silhouette argmax — annotation goes UP-LEFT into whitespace above the plateau
    ax.scatter([argmax_k], [argmax_sil], s=55, color=ARGMAX,
               edgecolor="white", linewidth=1.4, zorder=3)
    ax.annotate(
        f"k = {argmax_k}  (argmax)\nsilhouette = {argmax_sil:.3f}",
        xy=(argmax_k, argmax_sil), xytext=(argmax_k - 8, argmax_sil + 0.035),
        fontsize=8.5, color=ARGMAX, ha="right",
        arrowprops=dict(arrowstyle="-", color=ARGMAX, lw=0.7),
    )

    # ── Panel B: Davies-Bouldin (lower = better) ────────────────────────────
    ax = axes[1]
    ax.plot(df["k"], df["davies_bouldin"], color=LINE, linewidth=1.6, zorder=2)
    ax.set_ylabel("Davies-Bouldin index (lower = better)", fontsize=10)
    ax.set_title("B. Davies-Bouldin", fontsize=11, loc="left", pad=8)

    # Headroom on both ends so annotations don't clip
    _y_lo, _y_hi = ax.get_ylim()
    ax.set_ylim(_y_lo - 0.05, _y_hi + 0.03)

    # Panel B curve is monotonically decreasing, so the upper-right region is
    # empty — place both annotations well above the curve to avoid overlap.
    ax.scatter([DEPLOYED_K], [deployed_db], s=55, color=DEPLOYED,
               edgecolor="white", linewidth=1.4, zorder=3)
    ax.annotate(
        f"k = {DEPLOYED_K}  (deployed)\nD-B = {deployed_db:.3f}",
        xy=(DEPLOYED_K, deployed_db),
        xytext=(DEPLOYED_K + 25, deployed_db + 0.13),
        fontsize=8.5, color=DEPLOYED, ha="left",
        arrowprops=dict(arrowstyle="-", color=DEPLOYED, lw=0.7),
    )

    ax.scatter([argmin_db_k], [argmin_db], s=55, color=ARGMAX,
               edgecolor="white", linewidth=1.4, zorder=3)
    ax.annotate(
        f"k = {argmin_db_k}  (argmin)\nD-B = {argmin_db:.3f}",
        xy=(argmin_db_k, argmin_db),
        xytext=(argmin_db_k - 40, argmin_db + 0.15),
        fontsize=8.5, color=ARGMAX, ha="right",
        arrowprops=dict(arrowstyle="-", color=ARGMAX, lw=0.7),
    )

    for f in (OUT_PNG, OUT_PDF):
        fig.savefig(f, dpi=300, bbox_inches="tight")
        print(f"Wrote {f.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
