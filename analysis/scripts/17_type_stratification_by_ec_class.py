"""Refined Panel C — protein-class stratification, with the Enzyme category
subdivided by EC first-digit (EC1 oxidoreductases, EC2 transferases, EC3
hydrolases, EC4 lyases, EC5 isomerases, EC6 ligases, EC7 translocases).

For each cluster:
  - if dominant protein type is 'Enzyme', re-label by dominant EC class
    (first digit of the most-common EC number in the cluster)
  - else keep the original type (Transporter / Receptor / Other)

Then re-render Panel C with the expanded category axis.

Outputs (under analysis/):
    figure_X_C_type_stratification_ec.{pdf,png}
    cluster_type_composition_ec.tsv           per-cluster refined labels
"""
from __future__ import annotations

import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402

EC_LABELS = {
    "1": "EC 1\nOxidoreductases",
    "2": "EC 2\nTransferases",
    "3": "EC 3\nHydrolases",
    "4": "EC 4\nLyases",
    "5": "EC 5\nIsomerases",
    "6": "EC 6\nLigases",
    "7": "EC 7\nTranslocases",
}

EC_SHORT = {
    "1": "EC 1", "2": "EC 2", "3": "EC 3", "4": "EC 4",
    "5": "EC 5", "6": "EC 6", "7": "EC 7",
}

# Color palette by category (harmonize with Fig X)
CAT_COLORS = {
    "EC 1": "#1e40af",   # deep blue — oxidoreductases (dominant class)
    "EC 2": "#0891b2",   # cyan — transferases
    "EC 3": "#059669",   # green — hydrolases
    "EC 4": "#ca8a04",   # gold — lyases
    "EC 5": "#e11d48",   # rose — isomerases
    "EC 6": "#a21caf",   # purple — ligases
    "EC 7": "#6b7280",   # slate — translocases
    "Transporter": "#16a34a",
    "Receptor":    "#ea580c",
    "Other":       "#a855f7",
}
NULL_COLOR = "#7f1d1d"
RANDOM_MEDIAN = 0.20

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Load ────────────────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["cluster"] = prot["cluster"].astype(int)
prot["ec_numbers"] = prot["ec_numbers"].fillna("").astype(str)
prot["keyword_labels"] = prot["keyword_labels"].fillna("").astype(str)

sig = pd.read_csv(output("cluster_tanimoto_significance.tsv"), sep="\t")
comp = pd.read_csv(output("cluster_type_composition.tsv"), sep="\t")
tbl = sig.merge(comp[["cluster", "dominant_type", "n_proteins"]],
                on="cluster", how="left")

# ── For each protein, extract its dominant EC class (first digit) ────────
_EC_RE = re.compile(r"(\d)\.\d+\.\d+\.[\d\-]+")


def ec_first_digit(s: str) -> str:
    """Return the most-common first-digit EC class in the string, or ''."""
    hits = _EC_RE.findall(s)
    if not hits:
        return ""
    # Most common (proteins can list multiple ECs; take dominant)
    from collections import Counter
    return Counter(hits).most_common(1)[0][0]


prot["ec_class"] = prot["ec_numbers"].map(ec_first_digit)

# ── Per-cluster: dominant EC class (only for enzyme-dominant clusters) ──
rows = []
for _, r in tbl.iterrows():
    cid = int(r["cluster"])
    sub = prot[prot["cluster"] == cid]
    if r["dominant_type"] == "Enzyme":
        # Dominant EC class among proteins that have an EC
        classes = sub[sub["ec_class"] != ""]["ec_class"]
        if not classes.empty:
            top_class, top_n = classes.value_counts().index[0], \
                               classes.value_counts().iloc[0]
            frac = top_n / len(classes)
            refined = EC_SHORT.get(top_class, "")
        else:
            refined, frac = "Enzyme (no EC)", 0.0
    else:
        refined = r["dominant_type"]
        frac = float("nan")
    rows.append(dict(
        cluster=cid,
        family=r["family"],
        original_type=r["dominant_type"],
        refined_category=refined,
        refined_frac=round(frac, 3) if frac == frac else np.nan,
        observed=r["observed"],
        n_proteins=r["n_proteins"],
    ))
refined = pd.DataFrame(rows)
refined.to_csv(output("cluster_type_composition_ec.tsv"),
               sep="\t", index=False)

# Category breakdown
cat_counts = refined["refined_category"].value_counts()
print("Refined category breakdown across the 82 clusters:")
for k in ["EC 1", "EC 2", "EC 3", "EC 4", "EC 5", "EC 6", "EC 7",
         "Transporter", "Receptor", "Other", "Mixed"]:
    n = int(cat_counts.get(k, 0))
    if n:
        print(f"  {k:14s}  {n}")

# ── Figure ─────────────────────────────────────────────────────────────
# Order: EC classes first (1→7), then non-enzymes
order = [k for k in ["EC 1", "EC 2", "EC 3", "EC 4", "EC 5", "EC 6", "EC 7"]
         if cat_counts.get(k, 0) > 0]
order += [k for k in ["Transporter", "Receptor", "Other"]
          if cat_counts.get(k, 0) > 0]

data = []
labels = []
for cat in order:
    vals = refined[refined["refined_category"] == cat]["observed"].dropna().values
    if len(vals) == 0:
        continue
    data.append(vals)
    n = len(vals)
    if cat.startswith("EC"):
        pretty = EC_LABELS.get(cat.split()[-1], cat)
    else:
        pretty = cat
    labels.append(f"{pretty}\n(n={n})")

fig, ax = plt.subplots(figsize=(13, 6.5))
pos = np.arange(len(data))
parts = ax.violinplot(data, positions=pos, widths=0.75,
                      showmedians=True, showextrema=False)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(CAT_COLORS.get(order[i], "#94a3b8"))
    body.set_alpha(0.78)
    body.set_edgecolor("black")
    body.set_linewidth(0.5)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.4)

# Jittered dots for individual clusters
for i, d in enumerate(data):
    x = np.random.default_rng(i).normal(i, 0.055, size=len(d))
    ax.scatter(x, d, s=22, color="black", alpha=0.6, zorder=3,
               edgecolor="white", linewidth=0.4)

# Divider between enzyme classes and non-enzymes
n_ec = sum(1 for k in order if k.startswith("EC"))
if n_ec < len(order):
    ax.axvline(n_ec - 0.5, color="black", linewidth=1, linestyle=":")
    ax.text((n_ec - 1) / 2, 0.83, "ENZYMES (split by EC class)",
            ha="center", fontsize=10, fontweight="bold", color="#334155")
    ax.text(n_ec + (len(order) - n_ec - 1) / 2, 0.83, "NON-ENZYMES",
            ha="center", fontsize=10, fontweight="bold", color="#334155")

ax.axhline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.2,
           label=f"Random-pair baseline = {RANDOM_MEDIAN:.2f}")
ax.set_xticks(pos)
ax.set_xticklabels(labels)
ax.set_ylabel("Cluster mean substrate Tanimoto (ECFP4)")
ax.set_title("Substrate chemistry stratifies by EC class and protein type",
             pad=10)
ax.set_ylim(0.15, 0.88)
ax.grid(axis="y", alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig(output("figure_X_C_type_stratification_ec.pdf"))
fig.savefig(output("figure_X_C_type_stratification_ec.png"),
            dpi=300, bbox_inches="tight")
print("\nWrote figure_X_C_type_stratification_ec.{pdf,png}")

# ── Also list which specific enzyme clusters belong to each EC class ──
print("\nEnzyme clusters by EC class:")
for cat in ["EC 1", "EC 2", "EC 3", "EC 4", "EC 5", "EC 6", "EC 7"]:
    rows_ = refined[refined["refined_category"] == cat]
    if rows_.empty:
        continue
    print(f"\n  {cat}: {len(rows_)} clusters")
    for _, r in rows_.sort_values("observed", ascending=False).iterrows():
        print(f"    C{int(r['cluster'])+1:2d}  T={r['observed']:.3f}  "
              f"n={int(r['n_proteins']):4d}  {str(r['family'])[:60]}")

print("\nDone.")
