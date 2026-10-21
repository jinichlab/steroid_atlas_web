"""Intra-cluster Tanimoto analysis for the steroid-molecule UMAP clusters.

For each of the 8 molecule-view clusters:
  * pairwise Tanimoto similarity of Morgan (ECFP4) fingerprints
  * per-cluster stats: n, mean, median, std, min, max
  * per-molecule "centroid similarity" = mean Tanimoto to all other members
    of the same cluster (an entry-level measure of how chemically central
    the molecule is inside its cluster)
  * background: random-pair Tanimoto distribution across ALL molecules,
    as a null to compare against

Outputs (all under analysis/):
    intra_cluster_tanimoto.tsv                per-cluster summary
    per_molecule_tanimoto.tsv                 per-molecule centrality
    intra_cluster_tanimoto_distribution.pdf   box plot per cluster
    intra_vs_random_tanimoto.pdf              KDE overlay
    intra_cluster_tanimoto_report.md          plain-text summary

Run under the project's miniconda so RDKit is on the path:
    LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 \\
        analysis/09_intra_cluster_tanimoto.py
"""
from __future__ import annotations

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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
DATA = ROOT / "data" / "molecules.csv"

RNG_SEED = 20260824

# ── Load molecules ────────────────────────────────────────────────────────
print(f"Loading {DATA.name} ...")
df = pd.read_csv(DATA, low_memory=False)
print(f"  {len(df):,} rows")

need = ["compound_name", "smiles", "cluster"]
missing = [c for c in need if c not in df.columns]
if missing:
    raise SystemExit(f"molecules.csv missing required columns: {missing}")

df = df[need + (["chebi_id"] if "chebi_id" in df.columns else [])].copy()
df = df[df["smiles"].astype(str).str.strip().ne("")]
df["cluster"] = df["cluster"].astype(int)
print(f"  {len(df):,} rows with a SMILES string")

# ── Morgan (ECFP4) fingerprints ───────────────────────────────────────────
print("Computing ECFP4 fingerprints (radius=2, nBits=2048) ...")
mols, fps, keep_idx = [], [], []
for i, smi in enumerate(df["smiles"].astype(str)):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    mols.append(m)
    fps.append(fp)
    keep_idx.append(i)

df = df.iloc[keep_idx].reset_index(drop=True)
print(f"  {len(fps):,} fingerprints; dropped {len(keep_idx) - len(fps):,} unparseable rows")


def pairwise_tanimoto(fp_list):
    """Return dense n×n numpy array of Tanimoto similarities (diagonal = 1)."""
    n = len(fp_list)
    M = np.ones((n, n), dtype=np.float32)
    for i in range(n):
        # RDKit's BulkTanimotoSimilarity returns [sim(fp_i, fp_j) for j in fps]
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list)
        M[i, :] = np.asarray(sims, dtype=np.float32)
    return M


# ── Per-cluster stats ─────────────────────────────────────────────────────
print("Computing per-cluster pairwise Tanimoto ...")
rows = []
per_mol_rows = []
for cid, sub in df.groupby("cluster"):
    idx = sub.index.to_list()
    fp_sub = [fps[i] for i in idx]
    n = len(fp_sub)
    if n < 2:
        rows.append({
            "cluster": int(cid),
            "n": n,
            "mean_tanimoto": np.nan,
            "median_tanimoto": np.nan,
            "std_tanimoto": np.nan,
            "min_tanimoto": np.nan,
            "max_tanimoto": np.nan,
        })
        # Single-member cluster: centrality is undefined.
        per_mol_rows.append({
            "compound_name": sub["compound_name"].iloc[0],
            "chebi_id": sub.get("chebi_id", pd.Series([""] * n)).iloc[0]
                       if "chebi_id" in sub.columns else "",
            "cluster": int(cid),
            "n_in_cluster": n,
            "mean_tanimoto_to_cluster": np.nan,
            "nearest_neighbor_tanimoto_in_cluster": np.nan,
        })
        continue

    M = pairwise_tanimoto(fp_sub)
    # Upper-triangle pairs only (exclude diagonal + duplicates)
    iu = np.triu_indices(n, k=1)
    pair_sims = M[iu]

    rows.append({
        "cluster": int(cid),
        "n": n,
        "mean_tanimoto": float(np.mean(pair_sims)),
        "median_tanimoto": float(np.median(pair_sims)),
        "std_tanimoto": float(np.std(pair_sims)),
        "min_tanimoto": float(np.min(pair_sims)),
        "max_tanimoto": float(np.max(pair_sims)),
    })

    # Per-molecule: mean sim to all OTHER members (exclude self)
    np.fill_diagonal(M, np.nan)
    mean_to_cluster = np.nanmean(M, axis=1)
    # Nearest-neighbor sim within cluster (excluding self)
    nn_in_cluster = np.nanmax(M, axis=1)

    for j, i in enumerate(idx):
        row_ = df.iloc[i]
        per_mol_rows.append({
            "compound_name": row_["compound_name"],
            "chebi_id": row_.get("chebi_id", ""),
            "cluster": int(cid),
            "n_in_cluster": n,
            "mean_tanimoto_to_cluster": float(mean_to_cluster[j]),
            "nearest_neighbor_tanimoto_in_cluster": float(nn_in_cluster[j]),
        })

stats = pd.DataFrame(rows).sort_values("cluster")
per_mol = pd.DataFrame(per_mol_rows).sort_values(
    ["cluster", "mean_tanimoto_to_cluster"], ascending=[True, False]
)

# Round for readability in the TSVs
for c in stats.columns:
    if c not in ("cluster", "n"):
        stats[c] = stats[c].astype(float).round(4)
for c in ("mean_tanimoto_to_cluster", "nearest_neighbor_tanimoto_in_cluster"):
    if c in per_mol.columns:
        per_mol[c] = per_mol[c].astype(float).round(4)

stats.to_csv(output("intra_cluster_tanimoto.tsv"), sep="\t", index=False)
per_mol.to_csv(output("per_molecule_tanimoto.tsv"), sep="\t", index=False)
print(f"  wrote intra_cluster_tanimoto.tsv ({len(stats)} clusters)")
print(f"  wrote per_molecule_tanimoto.tsv ({len(per_mol)} rows)")

# ── Null / background: random-pair Tanimoto across ALL molecules ─────────
print("Sampling random-pair Tanimoto for background distribution ...")
rng = np.random.default_rng(RNG_SEED)
n_all = len(fps)
n_samples = min(20_000, n_all * (n_all - 1) // 2)
i_idx = rng.integers(0, n_all, size=n_samples)
j_idx = rng.integers(0, n_all, size=n_samples)
mask = i_idx != j_idx
i_idx, j_idx = i_idx[mask], j_idx[mask]
random_sims = np.array([
    DataStructs.TanimotoSimilarity(fps[i], fps[j])
    for i, j in zip(i_idx, j_idx)
], dtype=np.float32)
print(f"  {len(random_sims):,} random pairs sampled")

# ── Box plot per cluster ─────────────────────────────────────────────────
print("Rendering box plot per cluster ...")
cluster_pair_lists = {}
for cid, sub in df.groupby("cluster"):
    idx = sub.index.to_list()
    fp_sub = [fps[i] for i in idx]
    n = len(fp_sub)
    if n < 2:
        cluster_pair_lists[int(cid)] = np.array([])
        continue
    M = pairwise_tanimoto(fp_sub)
    iu = np.triu_indices(n, k=1)
    cluster_pair_lists[int(cid)] = M[iu]

fig, ax = plt.subplots(figsize=(9, 4.5))
data = [cluster_pair_lists[c] for c in sorted(cluster_pair_lists)]
labels = [f"C{c+1}\n(n={len(df[df['cluster']==c])})" for c in sorted(cluster_pair_lists)]
bp = ax.boxplot(
    data,
    labels=labels,
    showfliers=False,
    patch_artist=True,
    boxprops=dict(facecolor="#93c5fd", edgecolor="#1e40af"),
    medianprops=dict(color="#1e3a8a", linewidth=2),
    whiskerprops=dict(color="#1e40af"),
    capprops=dict(color="#1e40af"),
)
ax.axhline(np.median(random_sims), color="#dc2626", linestyle="--",
           linewidth=1, label=f"random-pair median = {np.median(random_sims):.2f}")
ax.set_ylabel("Pairwise Tanimoto (ECFP4)")
ax.set_title("Intra-cluster chemical similarity — steroid molecule clusters")
ax.set_ylim(0, 1)
ax.grid(axis="y", alpha=0.3)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(output("intra_cluster_tanimoto_distribution.pdf"))
fig.savefig(output("intra_cluster_tanimoto_distribution.png"), dpi=180)
print(f"  wrote intra_cluster_tanimoto_distribution.{{pdf,png}}")

# ── KDE-style comparison: intra vs random ────────────────────────────────
print("Rendering intra-vs-random KDE ...")
all_intra = np.concatenate([v for v in cluster_pair_lists.values() if v.size])

fig, ax = plt.subplots(figsize=(7, 4))
bins = np.linspace(0, 1, 51)
ax.hist(random_sims, bins=bins, density=True, alpha=0.5,
        color="#94a3b8", label=f"Random pairs (n={len(random_sims):,})")
ax.hist(all_intra, bins=bins, density=True, alpha=0.6,
        color="#0369a1", label=f"Intra-cluster pairs (n={len(all_intra):,})")
ax.set_xlabel("Tanimoto similarity (ECFP4)")
ax.set_ylabel("Density")
ax.set_title("Intra-cluster chemical similarity is above random background")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(output("intra_vs_random_tanimoto.pdf"))
fig.savefig(output("intra_vs_random_tanimoto.png"), dpi=180)
print(f"  wrote intra_vs_random_tanimoto.{{pdf,png}}")

# ── Markdown summary ─────────────────────────────────────────────────────
print("Writing plain-text report ...")
lines = [
    "# Intra-cluster Tanimoto — steroid molecule clusters",
    "",
    f"- Total molecules with valid SMILES: **{len(df):,}**",
    f"- Random-pair Tanimoto: mean **{random_sims.mean():.3f}**, "
    f"median **{np.median(random_sims):.3f}** "
    f"(n = {len(random_sims):,} sampled pairs)",
    f"- All intra-cluster pairs: mean **{all_intra.mean():.3f}**, "
    f"median **{np.median(all_intra):.3f}** "
    f"(n = {len(all_intra):,} pairs)",
    "",
    "## Per-cluster",
    "",
    "| Cluster | n | mean | median | std | min | max |",
    "|---|---|---|---|---|---|---|",
]
for _, row in stats.iterrows():
    lines.append(
        f"| C{int(row['cluster'])+1} | {int(row['n'])} | "
        f"{row['mean_tanimoto']:.3f} | {row['median_tanimoto']:.3f} | "
        f"{row['std_tanimoto']:.3f} | {row['min_tanimoto']:.3f} | "
        f"{row['max_tanimoto']:.3f} |"
    )
lines += [
    "",
    "**Reading:** higher mean = the cluster is chemically tight (its "
    "members share substructure). Big gap between mean and min = the "
    "cluster contains some outliers RDKit sees as chemically distant "
    "despite sharing the UMAP neighborhood.",
    "",
    "## Files",
    "- `intra_cluster_tanimoto.tsv` — per-cluster summary",
    "- `per_molecule_tanimoto.tsv` — per-molecule centrality inside its cluster",
    "- `intra_cluster_tanimoto_distribution.{pdf,png}` — box plot",
    "- `intra_vs_random_tanimoto.{pdf,png}` — KDE overlay vs null",
]
(output("intra_cluster_tanimoto_report.md")).write_text(
    "\n".join(lines), encoding="utf-8"
)
print(f"  wrote intra_cluster_tanimoto_report.md")

print("\nDone.")
