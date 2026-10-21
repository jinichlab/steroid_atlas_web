"""Publication Figure X — cluster coherence + protein-class stratification.

4 panels on one page:

  A) UMAP overview: all 14,089 proteins colored by cluster (scene-setter).
  B) Cluster coherence: 82 clusters as horizontal bars, sorted by observed
     substrate Tanimoto. Bars colored by FDR significance (dark = q<0.05,
     grey = ns). Vertical dashed line at atlas random-pair baseline.
  C) Protein-class stratification: violins of cluster median Tanimoto for
     each of Enzyme / Transporter / Receptor / Other.
  D) Enzyme specialists vs promiscuous: horizontal violins of pairwise
     Tanimoto for top-10 tightest and bottom-10 loosest enzyme clusters,
     with family names on the y-axis.

Outputs:
    analysis/figure_X_coherence_and_type.pdf
    analysis/figure_X_coherence_and_type.png    (dpi=300)
"""
from __future__ import annotations

import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except ImportError as e:
    raise SystemExit("rdkit not installed. Run under the conda python.") from e

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402

# ── Palette ──────────────────────────────────────────────────────────────
TYPE_COLORS = {
    "Enzyme":      "#2563eb",  # blue
    "Receptor":    "#ea580c",  # orange
    "Transporter": "#16a34a",  # green
    "Other":       "#a855f7",  # purple
    "Mixed":       "#64748b",  # slate
}
SIG_COLOR   = "#0369a1"   # deep blue for significant
NSIG_COLOR  = "#cbd5e1"   # light slate for non-significant
NULL_COLOR  = "#7f1d1d"   # dark red for null baseline

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "figure.dpi": 120,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Load atlas + prior analysis outputs ──────────────────────────────────
print("Loading data ...")
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["cluster"] = prot["cluster"].astype(int)

sig = pd.read_csv(output("cluster_tanimoto_significance.tsv"), sep="\t")
comp = pd.read_csv(output("cluster_type_composition.tsv"), sep="\t")

# Merge significance + protein-type composition + family name
tbl = sig.merge(comp[["cluster", "dominant_type", "n_proteins"]],
                on="cluster", how="left")

print(f"  atlas: {len(prot):,} proteins, {prot['cluster'].nunique()} clusters")
print(f"  sig table: {len(tbl):,} rows")

# Random pair Tanimoto baseline (from earlier permutation runs it's ~0.20)
RANDOM_MEDIAN = 0.20

# Substrate pool sizes for the enzyme extremes panel (need raw Tanimoto arrays)
print("Precomputing molecule fingerprints ...")
mol = pd.read_csv(ROOT / "data" / "molecules.csv", low_memory=False)
mol = mol[["compound_name", "chebi_id", "smiles"]].copy()
mol = mol[mol["smiles"].astype(str).str.strip().ne("")]
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol.dropna(subset=["chebi_id"])
fp_by_chebi = {}
for _, r in mol.iterrows():
    m = Chem.MolFromSmiles(str(r["smiles"]))
    if m is None:
        continue
    fp_by_chebi[int(r["chebi_id"])] = AllChem.GetMorganFingerprintAsBitVect(
        m, 2, nBits=2048
    )

_NUM_RE = re.compile(r"\d+")
prot["chebi_set"] = prot["interacting_chebi_ids"].fillna("").astype(str).map(
    lambda s: {int(x) for x in _NUM_RE.findall(s)}
)


def cluster_sims(cid: int) -> np.ndarray:
    sub = prot[prot["cluster"] == cid]
    all_c: set[int] = set()
    for s in sub["chebi_set"]:
        all_c.update(s)
    idxs = sorted(c for c in all_c if c in fp_by_chebi)
    if len(idxs) < 2:
        return np.array([], dtype=np.float32)
    fps = [fp_by_chebi[c] for c in idxs]
    sims = []
    for i in range(len(idxs) - 1):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    return np.asarray(sims, dtype=np.float32)


# ── Set up figure ────────────────────────────────────────────────────────
print("Composing figure ...")
fig = plt.figure(figsize=(13, 11))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.2], height_ratios=[1.0, 1.2],
              hspace=0.32, wspace=0.28)


def panel_label(ax, letter, y=1.03, x=-0.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom", ha="left")


# ── Panel A: UMAP overview ───────────────────────────────────────────────
print("  Panel A: UMAP overview ...")
axA = fig.add_subplot(gs[0, 0])
# Use a large discrete-ish palette
n_clusters = int(prot["cluster"].max()) + 1
cmap = plt.get_cmap("gist_ncar", n_clusters)
axA.scatter(prot["umap_1"], prot["umap_2"],
            c=prot["cluster"], cmap=cmap, s=1.4,
            alpha=0.75, linewidth=0)
axA.set_xlabel("UMAP-1")
axA.set_ylabel("UMAP-2")
axA.set_title(f"{len(prot):,} proteins across {n_clusters} k-means clusters",
              pad=6)
# Small annotation with atlas totals
axA.text(0.02, 0.98,
         f"14,089 steroid-interacting proteins\n"
         f"82 clusters (k selected by silhouette)",
         transform=axA.transAxes, fontsize=8, va="top",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor="black", linewidth=0.4, alpha=0.9))
axA.grid(alpha=0.25)
panel_label(axA, "A")

# ── Panel B: Cluster-coherence bar plot ─────────────────────────────────
print("  Panel B: coherence bars ...")
axB = fig.add_subplot(gs[0, 1])
tblB = tbl.dropna(subset=["observed", "z_score"]).copy()
tblB = tblB.sort_values("observed", ascending=True)  # smallest at bottom → biggest at top
tblB["is_sig"] = tblB["q_specialist"] < 0.05
colors = [SIG_COLOR if s else NSIG_COLOR for s in tblB["is_sig"]]
axB.barh(np.arange(len(tblB)), tblB["observed"], color=colors,
         edgecolor="black", linewidth=0.2, height=0.85)
axB.axvline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.2,
            label=f"random-pair baseline = {RANDOM_MEDIAN:.2f}")
axB.set_yticks([])
axB.set_ylabel(f"77 clusters (with ≥2 substrates), sorted by Tanimoto")
axB.set_xlabel("Mean pairwise substrate Tanimoto (ECFP4)")
n_sig = int(tblB["is_sig"].sum())
n_tot = len(tblB)
axB.set_title(f"Cluster substrate coherence  ·  "
              f"{n_sig}/{n_tot} clusters significant (BH-FDR q<0.05)",
              pad=6)

# Legend
from matplotlib.patches import Patch
legend_h = [
    Patch(facecolor=SIG_COLOR, edgecolor="black", label="q < 0.05 (coherent)"),
    Patch(facecolor=NSIG_COLOR, edgecolor="black", label="n.s."),
    Patch(facecolor="none", edgecolor=NULL_COLOR, label=f"Random baseline"),
]
axB.legend(handles=legend_h, loc="lower right", framealpha=0.95, fontsize=8)
axB.set_xlim(0, 0.85)
axB.grid(axis="x", alpha=0.3)
panel_label(axB, "B")

# ── Panel C: violins of cluster medians by protein type ─────────────────
print("  Panel C: type stratification ...")
axC = fig.add_subplot(gs[1, 0])
tblC = tbl.dropna(subset=["observed", "dominant_type"]).copy()
type_order = ["Enzyme", "Transporter", "Receptor", "Other", "Mixed"]
data_by_type = [tblC.loc[tblC["dominant_type"] == t, "observed"].values
                for t in type_order]
counts_by_type = [len(d) for d in data_by_type]

parts = axC.violinplot(
    data_by_type,
    positions=np.arange(len(type_order)),
    widths=0.75,
    showmedians=True, showextrema=False,
)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(TYPE_COLORS[type_order[i]])
    body.set_alpha(0.75)
    body.set_edgecolor("black")
    body.set_linewidth(0.5)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.2)

# Overlay dots for individual clusters
for i, d in enumerate(data_by_type):
    x = np.random.default_rng(i).normal(i, 0.06, size=len(d))
    axC.scatter(x, d, s=14, color="black", alpha=0.55, zorder=3,
                edgecolor="white", linewidth=0.3)

axC.axhline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1)
axC.set_xticks(np.arange(len(type_order)))
axC.set_xticklabels([f"{t}\n(n={c})" for t, c in zip(type_order, counts_by_type)])
axC.set_ylabel("Cluster mean substrate Tanimoto (ECFP4)")
axC.set_title("Substrate chemistry stratifies by protein class", pad=6)
axC.set_ylim(0, 0.85)
axC.grid(axis="y", alpha=0.3)
panel_label(axC, "C")

# ── Panel D: enzyme extremes ─────────────────────────────────────────────
print("  Panel D: enzyme specialists vs promiscuous ...")
axD = fig.add_subplot(gs[1, 1])
enz = tbl[(tbl["dominant_type"] == "Enzyme") &
          tbl["observed"].notna()].copy()
enz = enz.merge(prot.groupby("cluster")
                    .apply(lambda s: len({c for x in s["chebi_set"] for c in x
                                          if c in fp_by_chebi}))
                    .rename("n_substrates_here").reset_index(),
                on="cluster", how="left")
enz = enz[enz["n_substrates_here"] >= 5].sort_values("observed", ascending=False)

TOP_K = 8
top_enz = enz.head(TOP_K).sort_values("observed", ascending=True)
bot_enz = enz.tail(TOP_K).sort_values("observed", ascending=True)
show = pd.concat([bot_enz, top_enz])  # bottom-to-top: loose → tight

sims_arrays = []
labels = []
for _, r in show.iterrows():
    sims = cluster_sims(int(r["cluster"]))
    if sims.size == 0:
        continue
    sims_arrays.append(sims)
    fam = str(r["family"])[:36]
    labels.append(f"{fam}\n(n_prot={int(r['n_proteins'])}, "
                  f"n_sub={int(r['n_substrates_here'])})")

pos = np.arange(len(sims_arrays))
parts = axD.violinplot(sims_arrays, positions=pos, vert=False,
                        widths=0.85, showmedians=True, showextrema=False)
# Half get one color, half get another — split at midpoint
midpoint = TOP_K   # bottom TOP_K are loose; top TOP_K are tight
for i, body in enumerate(parts["bodies"]):
    if i < TOP_K:
        body.set_facecolor("#f97316")   # orange = promiscuous
    else:
        body.set_facecolor("#0369a1")   # blue = specialist
    body.set_alpha(0.75)
    body.set_edgecolor("black")
    body.set_linewidth(0.5)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.2)

axD.axvline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1,
            label=f"random-pair baseline")

# Divider line between promiscuous and specialist blocks
axD.axhline(TOP_K - 0.5, color="black", linewidth=0.6, linestyle=":")
axD.text(0.82, TOP_K - 1, "◄ BROAD-SPECIFICITY", fontsize=8, color="#f97316",
         fontweight="bold", va="center")
axD.text(0.82, TOP_K + 0.2, "◄ SPECIALIST", fontsize=8, color="#0369a1",
         fontweight="bold", va="center")

axD.set_yticks(pos)
axD.set_yticklabels(labels, fontsize=7)
axD.set_xlabel("Pairwise substrate Tanimoto (ECFP4)")
axD.set_title(f"Enzyme clusters — top {TOP_K} specialists vs top {TOP_K} broad-specificity",
              pad=6)
axD.set_xlim(0, 1)
axD.grid(axis="x", alpha=0.3)
axD.legend(loc="lower right", fontsize=8)
panel_label(axD, "D")

# ── Save ────────────────────────────────────────────────────────────────
fig.suptitle(
    "Figure X · Cluster substrate coherence and protein-class stratification",
    fontsize=13, fontweight="bold", y=0.995,
)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(output("figure_X_coherence_and_type.pdf"))
fig.savefig(output("figure_X_coherence_and_type.png"),
            dpi=300, bbox_inches="tight")
print(f"\nWrote figure_X_coherence_and_type.{{pdf,png}}")
print("Done.")
