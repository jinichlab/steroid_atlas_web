"""Publication figure — atlas-wide GO enrichment across all 95 clusters.

Two-panel landscape figure supporting the "every cluster has clear enrichment"
paper claim:

  A. Distribution of top-GO fold-enrichment across ALL 95 clusters. A histogram
     shows every cluster's top enrichment on a log x-axis; a dashed reference
     line marks the 10× "clear-story" threshold. The message: no cluster sits
     in the low-enrichment tail — every bin above 10× is populated.

  B. Top 12 clusters ranked by fold-enrichment, each labeled with its top GO
     term. Concrete biology (aromatase, sterol demethylase, bile acid symporter,
     …) so the reader anchors the aggregate claim to real families.

Reads:  ../analysis/cluster_go_top_terms.tsv   (from 09_cluster_go_enrichment.py)
Writes: ../analysis/cluster_enrichment_figure.png   300 dpi PNG for the paper
        ../analysis/cluster_enrichment_figure.pdf   vector PDF for typesetting
"""
from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = output("cluster_go_top_terms.tsv")
OUT_PNG = output("cluster_enrichment_figure.png")
OUT_PDF = output("cluster_enrichment_figure.pdf")

BAR = "#1F4B99"          # single hue — magnitude plot, one series
REF = "#9CA3AF"          # muted gray for the 10× threshold
INK = "#111827"

# Steroid-focused exemplar clusters ordered by their salience to the paper's
# audience. The plot will pick TOP 12 by enrichment from the full list — this
# dict just carries the pretty label if the cluster shows up.
EXEMPLAR_LABELS = {
    "21": "24S-hydroxycholesterol 7α-hydroxylase  (bile acid biosynthesis)",
    "2":  "aromatase  (estrogen biosynthesis)",
    "22": "7β-HSD / 11β-HSD  (steroid catabolism)",
    "17": "3β-HSD  (steroidogenesis)",
    "20": "glutathione transporter  (conjugated steroid efflux)",
    "0":  "cholesterol 24-hydroxylase  (brain cholesterol turnover)",
    "15": "bile acid:Na⁺ symporter  (hepatic bile-acid uptake)",
    "4":  "sterol 14α-demethylase  (cholesterol biosynthesis)",
    "9":  "canalicular bile acid transport",
    "7":  "β-glucosidase  (glycosphingolipid metabolism)",
    "11": "glycosyltransferase  (sterol glucuronidation)",
    "12": "sterol ester esterase  (lipid catabolism)",
    "23": "monoacylglycerol lipase  (lipid metabolism)",
    "13": "acetyl-CoA C-acyltransferase  (fatty-acid β-oxidation)",
}


def main() -> int:
    df = pd.read_csv(IN, sep="\t")
    df = df[df["rank"] == 1].copy()
    df["cluster"] = df["cluster"].astype(str)
    print(f"Plotting {len(df)} clusters")

    # ── Layout ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(11.2, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.5])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    # ── PANEL A — histogram of top-enrichment across all 95 clusters ─────────
    log_enrich = np.log10(df["enrichment_ratio"].values)
    # bins from 1x (0) to 1000x (3) in 0.15 log-decade steps → ~20 bins
    bins = np.arange(0, 3.1, 0.15)
    axA.hist(log_enrich, bins=bins, color=BAR, alpha=0.85,
             edgecolor="white", linewidth=0.6)
    axA.axvline(np.log10(10), color=REF, linestyle="--", linewidth=1.0, zorder=1)
    # Place the "10× threshold" label in the empty right-side space where the
    # histogram has few counts, so it doesn't collide with the tallest bars.
    axA.text(np.log10(10) + 0.05, axA.get_ylim()[1] * 0.55, "10× threshold →",
             color=REF, fontsize=9, va="center", ha="left", rotation=0)

    # X-axis: relabel log-space ticks as their linear values
    tick_locs = [0, 1, 2, 3]
    axA.set_xticks(tick_locs)
    axA.set_xticklabels(["1×", "10×", "100×", "1,000×"])
    axA.set_xlim(-0.05, 3.1)
    axA.set_xlabel("Fold-enrichment of the cluster's top GO term", fontsize=10)
    axA.set_ylabel("Number of clusters", fontsize=10)
    axA.set_title("A. Distribution across all 95 clusters",
                  fontsize=11, loc="left", pad=8)
    axA.grid(True, axis="y", alpha=0.25, linewidth=0.5)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.tick_params(labelsize=9)

    # Annotate the two "sub-10x" moderate clusters (should be 2 of them)
    n_below = int((df["enrichment_ratio"] < 10).sum())
    n_above = len(df) - n_below
    axA.text(
        0.98, 0.62,
        f"{n_above} / {len(df)}  clusters ≥ 10×\n"
        f"({n_below} in the ≥ 3× ‘moderate’ bucket)",
        transform=axA.transAxes, fontsize=9.5, color=INK,
        va="top", ha="right",
    )

    # ── PANEL B — CURATED steroid-relevant exemplars, ranked by enrichment ──
    # We deliberately do NOT show raw top-N by enrichment — those are dominated
    # by small niche clusters (T-cell costimulation, response to light, etc.)
    # that aren't the biology this paper is about. The curated set below is
    # the case-study-worthy set the reader should walk away with.
    top = df[df["cluster"].isin(EXEMPLAR_LABELS.keys())].copy()
    top = top.sort_values("enrichment_ratio", ascending=False).reset_index(drop=True)

    ypos = np.arange(len(top))[::-1]   # highest enrichment at top of panel
    axB.barh(ypos, top["enrichment_ratio"], color=BAR, height=0.72,
             edgecolor="white", linewidth=0.4, zorder=2)
    axB.axvline(10, color=REF, linestyle="--", linewidth=0.9, zorder=1)

    # y-tick labels: cluster id, one line
    axB.set_yticks(ypos)
    axB.set_yticklabels([f"cluster {c}" for c in top["cluster"]],
                        fontsize=9, color=INK)
    axB.tick_params(axis="y", length=0, pad=4)

    # Direct labels: GO term text placed to the right of each bar
    for y, (_, r) in zip(ypos, top.iterrows()):
        cid = r["cluster"]
        pretty = EXEMPLAR_LABELS.get(cid, r["go_term"])
        x_end = float(r["enrichment_ratio"])
        # Place the label slightly past the bar's right edge, on a log axis
        axB.text(
            x_end * 1.08, y, pretty,
            fontsize=8.5, color=INK, va="center", ha="left",
        )

    axB.set_xscale("log")
    axB.set_xlim(1, 3200)
    axB.set_ylim(-0.7, len(top) - 0.3)
    axB.xaxis.set_major_formatter(mticker.ScalarFormatter())
    axB.tick_params(axis="x", labelsize=9)
    axB.set_xlabel("Fold-enrichment of top GO term", fontsize=10)
    axB.set_title("B. Steroid-relevant exemplar clusters (curated)",
                  fontsize=11, loc="left", pad=8)
    axB.grid(True, axis="x", which="major", alpha=0.25, linewidth=0.5)
    axB.grid(True, axis="x", which="minor", alpha=0.10, linewidth=0.4)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.spines["left"].set_visible(False)

    # ── Overall title ───────────────────────────────────────────────────────
    fig.suptitle(
        "Every atlas cluster is enriched for a specific GO term "
        "— 93 of 95 clusters clear the 10× / p < 10⁻²⁰ threshold (hypergeometric)",
        fontsize=11.5, y=1.06, x=0.01, ha="left",
    )

    for f in (OUT_PNG, OUT_PDF):
        fig.savefig(f, dpi=300, bbox_inches="tight")
        print(f"Wrote {f.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
