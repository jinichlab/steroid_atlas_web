"""Publication figure: within-cluster pairwise Tanimoto for every protein cluster.

Two panels:
  (a) All 82 clusters, horizontal violin plot, sorted top-to-bottom from
      most-uniform-chemistry (highest median Tanimoto) to most-diverse
      (lowest median), with family names on the y-axis. Each violin shows
      the full shape of the pairwise Tanimoto distribution WITHIN that
      cluster's substrate pool.
  (b) Zoom: top-15 tightest and top-15 most-diverse clusters on either
      side, with larger violins and per-cluster N + median annotation.

The atlas-wide random-pair Tanimoto median is drawn as a dashed reference
line so the reader can see where "random chemistry" would sit.

Run under the project's miniconda:
    LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 \\
        analysis/12_within_cluster_tanimoto_figure.py
"""
from __future__ import annotations

import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except ImportError as e:
    raise SystemExit("rdkit not installed. Run under the conda python.") from e

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
PROT_CSV = ROOT / "data" / "proteins.csv"
MOL_CSV = ROOT / "data" / "molecules.csv"
FP_TSV = output("cluster_fingerprints.tsv")

# ── Load fingerprints ────────────────────────────────────────────────────
print(f"Loading {MOL_CSV.name} ...")
mol = pd.read_csv(MOL_CSV, low_memory=False)
mol = mol[["compound_name", "chebi_id", "smiles"]].copy()
mol = mol[mol["smiles"].astype(str).str.strip().ne("")]
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol.dropna(subset=["chebi_id"])

fps_list = []
chebis_list = []
for _, r in mol.iterrows():
    m = Chem.MolFromSmiles(str(r["smiles"]))
    if m is None:
        continue
    fps_list.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
    chebis_list.append(int(r["chebi_id"]))
n_pool = len(fps_list)
idx_of_chebi = {c: i for i, c in enumerate(chebis_list)}
print(f"  {n_pool:,} unique steroids fingerprinted")

# Precompute full pairwise matrix (fast lookup for every cluster)
print("Precomputing full pairwise Tanimoto matrix ...")
T = np.zeros((n_pool, n_pool), dtype=np.float32)
for i in range(n_pool):
    T[i, :] = np.asarray(
        DataStructs.BulkTanimotoSimilarity(fps_list[i], fps_list),
        dtype=np.float32,
    )
T = np.maximum(T, T.T)
np.fill_diagonal(T, 1.0)

# Random-pair baseline
rng = np.random.default_rng(20260824)
n_sample = 20000
ii = rng.integers(0, n_pool, size=n_sample)
jj = rng.integers(0, n_pool, size=n_sample)
mask = ii != jj
random_pair_sims = T[ii[mask], jj[mask]]
random_median = float(np.median(random_pair_sims))
print(f"  random-pair median = {random_median:.3f}")

# ── Per-cluster substrate lists → pairwise Tanimoto arrays ───────────────
print(f"Loading {PROT_CSV.name} ...")
prot = pd.read_csv(PROT_CSV, low_memory=False)
prot = prot[["accession", "cluster", "interacting_chebi_ids"]].copy()
prot["cluster"] = prot["cluster"].astype(int)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)
_NUM_RE = re.compile(r"\d+")
prot["chebi_set"] = prot["interacting_chebi_ids"].map(
    lambda s: {int(x) for x in _NUM_RE.findall(s)}
)

name_by_cluster: dict[int, str] = {}
if FP_TSV.exists():
    fp_df = pd.read_csv(FP_TSV, sep="\t")
    for _, r in fp_df.iterrows():
        name_by_cluster[int(r["cluster"])] = str(r["dominant_stem"])[:55]


def pair_sims_for(cluster_id: int, sub: pd.DataFrame) -> np.ndarray:
    """Upper-triangle pairwise Tanimoto values for this cluster's substrates."""
    all_c: set[int] = set()
    for s in sub["chebi_set"]:
        all_c.update(s)
    idxs = [idx_of_chebi[c] for c in all_c if c in idx_of_chebi]
    idxs = sorted(idxs)
    if len(idxs) < 2:
        return np.array([], dtype=np.float32)
    sub_T = T[np.ix_(idxs, idxs)]
    iu = np.triu_indices(len(idxs), k=1)
    return sub_T[iu]


records = []
for cid, sub in prot.groupby("cluster"):
    sims = pair_sims_for(int(cid), sub)
    if sims.size == 0:
        continue
    records.append({
        "cluster": int(cid),
        "family": name_by_cluster.get(int(cid), f"C{int(cid)+1}"),
        "n_substrates": len({c for s in sub["chebi_set"] for c in s
                             if c in idx_of_chebi}),
        "sims": sims,
        "median": float(np.median(sims)),
        "mean": float(sims.mean()),
    })
print(f"  {len(records)} clusters with ≥2 fingerprinted substrates")

# ── PANEL A: horizontal violin, ALL clusters, sorted by median ───────────
records.sort(key=lambda r: r["median"], reverse=True)

# Colour by median Tanimoto: cool for tight, warm for diverse
cmap = plt.get_cmap("coolwarm_r")
med_vals = np.array([r["median"] for r in records])
norm = plt.Normalize(vmin=0.15, vmax=0.75)

fig_h = max(9, 0.24 * len(records))
fig, ax = plt.subplots(figsize=(11, fig_h))

positions = np.arange(len(records))
parts = ax.violinplot(
    [r["sims"] for r in records],
    positions=positions,
    vert=False,
    widths=0.85,
    showmedians=True,
    showextrema=False,
)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(cmap(norm(records[i]["median"])))
    body.set_edgecolor("black")
    body.set_alpha(0.85)
    body.set_linewidth(0.4)
if "cmedians" in parts:
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.4)

# y-tick labels: family (n=#substrates)
labels = [f"{r['family']}  (n={r['n_substrates']})" for r in records]
ax.set_yticks(positions)
ax.set_yticklabels(labels, fontsize=7)
ax.invert_yaxis()  # tightest at top

ax.axvline(random_median, color="#7f1d1d", linestyle="--", linewidth=1,
           label=f"Atlas random-pair median = {random_median:.2f}")
ax.set_xlim(0, 1)
ax.set_xlabel("Pairwise Tanimoto (ECFP4) — within cluster")
ax.set_title(
    f"Within-cluster substrate chemical similarity  "
    f"({len(records)} of 82 clusters with ≥2 fingerprinted substrates, "
    f"sorted by median)"
)
ax.grid(axis="x", alpha=0.3)
ax.legend(loc="lower right", fontsize=9)

# Small colorbar for the color→median mapping
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.02)
cb.set_label("Cluster median Tanimoto", fontsize=8)
cb.ax.tick_params(labelsize=8)

fig.tight_layout()
fig.savefig(output("within_cluster_tanimoto_all.pdf"))
fig.savefig(output("within_cluster_tanimoto_all.png"), dpi=180,
            bbox_inches="tight")
print("  wrote within_cluster_tanimoto_all.{pdf,png}")

# ── PANEL B: zoom on the extremes ────────────────────────────────────────
TOP_K = 15
tight = records[:TOP_K]
loose = records[-TOP_K:][::-1]

fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharex=True)

for ax_, group, title, base_color in (
    (axes[0], tight, f"Top {TOP_K} most CHEMICALLY UNIFORM clusters",
     "#0369a1"),
    (axes[1], loose, f"Top {TOP_K} most CHEMICALLY DIVERSE clusters",
     "#c2410c"),
):
    pos = np.arange(len(group))
    parts = ax_.violinplot(
        [r["sims"] for r in group],
        positions=pos,
        vert=False,
        widths=0.85,
        showmedians=True,
        showextrema=True,
    )
    for body in parts["bodies"]:
        body.set_facecolor(base_color)
        body.set_alpha(0.55)
        body.set_edgecolor("black")
        body.set_linewidth(0.5)
    for k in ("cmedians", "cmaxes", "cmins", "cbars"):
        if k in parts:
            parts[k].set_color("black")
            parts[k].set_linewidth(1.0)
    labels = [f"{r['family']}\n(n={r['n_substrates']}, "
              f"median={r['median']:.2f})" for r in group]
    ax_.set_yticks(pos)
    ax_.set_yticklabels(labels, fontsize=8)
    ax_.invert_yaxis()
    ax_.axvline(random_median, color="#7f1d1d", linestyle="--", linewidth=1)
    ax_.set_xlim(0, 1)
    ax_.set_xlabel("Pairwise Tanimoto (ECFP4)")
    ax_.set_title(title, fontsize=11)
    ax_.grid(axis="x", alpha=0.3)

axes[0].text(random_median + 0.01, len(tight) - 0.5,
             f"random-pair\nmedian = {random_median:.2f}",
             fontsize=8, color="#7f1d1d", va="top")

fig.suptitle(
    "Substrate chemistry distributions — cluster extremes",
    fontsize=13, y=1.005,
)
fig.tight_layout()
fig.savefig(output("within_cluster_tanimoto_extremes.pdf"))
fig.savefig(output("within_cluster_tanimoto_extremes.png"), dpi=180,
            bbox_inches="tight")
print("  wrote within_cluster_tanimoto_extremes.{pdf,png}")

print("\nDone.")
