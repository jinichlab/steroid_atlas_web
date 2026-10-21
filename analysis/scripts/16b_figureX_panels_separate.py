"""Render the 4 panels of Figure X as SEPARATE standalone figures, so each
can be polished independently before being composed together at the end.

Outputs:
    figure_X_A_umap.{pdf,png}
    figure_X_B_coherence_bars.{pdf,png}
    figure_X_C_type_stratification.{pdf,png}
    figure_X_D_enzyme_extremes.{pdf,png}
"""
from __future__ import annotations

import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except ImportError as e:
    raise SystemExit("rdkit not installed. Run under the conda python.") from e

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402

# ── Global palette ──────────────────────────────────────────────────────
TYPE_COLORS = {
    "Enzyme":      "#2563eb",
    "Receptor":    "#ea580c",
    "Transporter": "#16a34a",
    "Other":       "#a855f7",
    "Mixed":       "#64748b",
}
SIG_COLOR = "#0369a1"
NSIG_COLOR = "#cbd5e1"
NULL_COLOR = "#7f1d1d"
SPECIALIST_COLOR = "#0369a1"
BROAD_COLOR = "#f97316"
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

# ── Load data ────────────────────────────────────────────────────────────
print("Loading atlas + analysis tables ...")
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot["cluster"] = prot["cluster"].astype(int)

sig = pd.read_csv(output("cluster_tanimoto_significance.tsv"), sep="\t")
comp = pd.read_csv(output("cluster_type_composition.tsv"), sep="\t")
tbl = sig.merge(comp[["cluster", "dominant_type", "n_proteins"]],
                on="cluster", how="left")

# Fingerprints for D
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


# ═════════════════════════════════════════════════════════════════════════
# Panel A — UMAP overview
# ═════════════════════════════════════════════════════════════════════════
print("A: UMAP overview ...")
fig, ax = plt.subplots(figsize=(8, 7))
n_clusters = int(prot["cluster"].max()) + 1
cmap = plt.get_cmap("gist_ncar", n_clusters)
ax.scatter(prot["umap_1"], prot["umap_2"],
           c=prot["cluster"], cmap=cmap, s=2.0,
           alpha=0.8, linewidth=0)
ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.set_title(f"14,089 steroid-interacting proteins  ·  "
             f"{n_clusters} k-means clusters", pad=10)
ax.text(0.02, 0.98,
        "k selected by silhouette + composite score in [50, 95]",
        transform=ax.transAxes, fontsize=9, va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                  edgecolor="black", linewidth=0.5, alpha=0.9))
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(output("figure_X_A_umap.pdf"))
fig.savefig(output("figure_X_A_umap.png"), dpi=300, bbox_inches="tight")
plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════
# Panel B — Coherence bars
# ═════════════════════════════════════════════════════════════════════════
print("B: coherence bars ...")
fig, ax = plt.subplots(figsize=(9, 10))
tblB = tbl.dropna(subset=["observed", "z_score"]).copy()
tblB = tblB.sort_values("observed", ascending=True)
tblB["is_sig"] = tblB["q_specialist"] < 0.05
colors = [SIG_COLOR if s else NSIG_COLOR for s in tblB["is_sig"]]
ypos = np.arange(len(tblB))
ax.barh(ypos, tblB["observed"], color=colors,
        edgecolor="black", linewidth=0.2, height=0.85)
ax.axvline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.4,
           label=f"Random-pair baseline = {RANDOM_MEDIAN:.2f}")

# Sparse label every ~10 rows so reader can see WHICH clusters
label_every = max(1, len(tblB) // 12)
ytick_pos, ytick_labels = [], []
for i, (_, row) in enumerate(tblB.iterrows()):
    if i % label_every == 0 or i == len(tblB) - 1:
        fam = str(row["family"])[:34]
        ytick_pos.append(i)
        ytick_labels.append(f"C{int(row['cluster'])+1} · {fam}")
ax.set_yticks(ytick_pos)
ax.set_yticklabels(ytick_labels, fontsize=8)

ax.set_xlabel("Mean pairwise substrate Tanimoto (ECFP4)")
n_sig = int(tblB["is_sig"].sum())
n_tot = len(tblB)
ax.set_title(f"Cluster substrate coherence  ·  "
             f"{n_sig}/{n_tot} clusters significantly coherent "
             f"(BH-FDR q<0.05)", pad=10)

# Annotate the block
mid_y = len(tblB) * 0.62
ax.annotate(f"n = {n_sig}", xy=(0.45, mid_y),
            fontsize=16, fontweight="bold", color="white",
            ha="center", va="center")

from matplotlib.patches import Patch
legend_h = [
    Patch(facecolor=SIG_COLOR, edgecolor="black",
          label=f"q < 0.05 (coherent) — {n_sig} clusters"),
    Patch(facecolor=NSIG_COLOR, edgecolor="black",
          label=f"n.s. — {n_tot - n_sig} clusters"),
    Patch(facecolor="none", edgecolor=NULL_COLOR,
          label=f"Random-pair baseline = {RANDOM_MEDIAN:.2f}"),
]
ax.legend(handles=legend_h, loc="lower right", framealpha=0.95)
ax.set_xlim(0, 0.85)
ax.set_ylim(-0.5, len(tblB) - 0.5)
ax.grid(axis="x", alpha=0.3)
fig.tight_layout()
fig.savefig(output("figure_X_B_coherence_bars.pdf"))
fig.savefig(output("figure_X_B_coherence_bars.png"), dpi=300, bbox_inches="tight")
plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════
# Panel C — Type stratification
# ═════════════════════════════════════════════════════════════════════════
print("C: type stratification ...")
fig, ax = plt.subplots(figsize=(9, 6.5))
tblC = tbl.dropna(subset=["observed", "dominant_type"]).copy()
# Drop 'Mixed' — it's a curiosity, not a claim.
type_order = ["Enzyme", "Transporter", "Receptor", "Other"]
tblC = tblC[tblC["dominant_type"].isin(type_order)]
data_by_type = [tblC.loc[tblC["dominant_type"] == t, "observed"].values
                for t in type_order]
counts = [len(d) for d in data_by_type]

parts = ax.violinplot(
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
    parts["cmedians"].set_linewidth(1.4)

# Jittered dots
for i, d in enumerate(data_by_type):
    x = np.random.default_rng(i).normal(i, 0.055, size=len(d))
    ax.scatter(x, d, s=22, color="black", alpha=0.6, zorder=3,
               edgecolor="white", linewidth=0.4)

ax.axhline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.2,
           label=f"Random-pair baseline = {RANDOM_MEDIAN:.2f}")
ax.set_xticks(np.arange(len(type_order)))
ax.set_xticklabels([f"{t}\n(n={c} clusters)" for t, c in zip(type_order, counts)])
ax.set_ylabel("Cluster mean substrate Tanimoto (ECFP4)")
ax.set_title("Substrate chemistry stratifies by protein class", pad=10)
ax.set_ylim(0.15, 0.85)
ax.grid(axis="y", alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig(output("figure_X_C_type_stratification.pdf"))
fig.savefig(output("figure_X_C_type_stratification.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════
# Panel D — Enzyme specialists vs broad-specificity
# ═════════════════════════════════════════════════════════════════════════
print("D: enzyme extremes ...")
enz = tbl[(tbl["dominant_type"] == "Enzyme") &
          tbl["observed"].notna()].copy()
n_sub_map = (prot.groupby("cluster")
                 .apply(lambda s: len({c for x in s["chebi_set"] for c in x
                                       if c in fp_by_chebi}))
                 .rename("n_substrates_here").reset_index())
enz = enz.merge(n_sub_map, on="cluster", how="left")
enz = enz[enz["n_substrates_here"] >= 5].sort_values("observed", ascending=False)

TOP_K = 8
top_enz = enz.head(TOP_K).sort_values("observed", ascending=True)
bot_enz = enz.tail(TOP_K).sort_values("observed", ascending=True)
show = pd.concat([bot_enz, top_enz])

sims_arrays, labels = [], []
for _, r in show.iterrows():
    sims = cluster_sims(int(r["cluster"]))
    if sims.size == 0:
        continue
    sims_arrays.append(sims)
    fam = str(r["family"])[:40]
    labels.append(
        f"C{int(r['cluster'])+1} · {fam}  "
        f"(n_prot={int(r['n_proteins'])}, n_sub={int(r['n_substrates_here'])})"
    )

fig, ax = plt.subplots(figsize=(11, 8))
pos = np.arange(len(sims_arrays))
parts = ax.violinplot(sims_arrays, positions=pos, vert=False,
                      widths=0.85, showmedians=True, showextrema=False)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(BROAD_COLOR if i < TOP_K else SPECIALIST_COLOR)
    body.set_alpha(0.78)
    body.set_edgecolor("black")
    body.set_linewidth(0.5)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.4)

ax.axvline(RANDOM_MEDIAN, color=NULL_COLOR, linestyle="--", linewidth=1.2,
           label=f"Random-pair baseline = {RANDOM_MEDIAN:.2f}")

# Divider between the two blocks
ax.axhline(TOP_K - 0.5, color="black", linewidth=1.0, linestyle=":")
ax.text(1.02, TOP_K - TOP_K / 2 - 0.5,
        "BROAD-SPECIFICITY", fontsize=10, color=BROAD_COLOR,
        fontweight="bold", va="center", rotation=270,
        transform=ax.get_yaxis_transform())
ax.text(1.02, TOP_K + TOP_K / 2 - 0.5,
        "SPECIALIST", fontsize=10, color=SPECIALIST_COLOR,
        fontweight="bold", va="center", rotation=270,
        transform=ax.get_yaxis_transform())

ax.set_yticks(pos)
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("Pairwise substrate Tanimoto (ECFP4)")
ax.set_title(f"Enzyme clusters — top {TOP_K} specialists vs top {TOP_K} broad-specificity",
             pad=10)
ax.set_xlim(0, 1)
ax.grid(axis="x", alpha=0.3)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(output("figure_X_D_enzyme_extremes.pdf"))
fig.savefig(output("figure_X_D_enzyme_extremes.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

print("\nWrote 4 standalone figures:")
for name in ("A_umap", "B_coherence_bars", "C_type_stratification",
             "D_enzyme_extremes"):
    print(f"  analysis/figure_X_{name}.pdf/.png")
print("Done.")
