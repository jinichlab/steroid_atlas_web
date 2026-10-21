"""Statistical significance of intra-cluster substrate Tanimoto — per cluster.

For each of the 82 protein clusters, ask: **is the observed intra-cluster
mean pairwise Tanimoto significantly higher (specialist) or lower
(promiscuous) than what you'd get by drawing the same number of substrates
uniformly at random from the full molecule pool?**

Method
------
1. Precompute the full 589×589 Tanimoto matrix over every steroid molecule
   we have a fingerprint for (fast lookup for every downstream permutation).
2. For each cluster c with N_c substrates:
     - observed = mean pairwise Tanimoto over the N_c substrates it actually
       has (upper-triangle only)
     - null distribution = draw N_c random substrates from the full pool,
       compute the same mean, repeat N_PERM times
     - two-sided empirical p-value + Z-score
3. Benjamini-Hochberg FDR across the 82 clusters, separately in each
   direction (specialist vs promiscuous).

Outputs (under analysis/):
    cluster_tanimoto_significance.tsv         one row per cluster (82)
    cluster_tanimoto_significance_volcano.pdf/png
    cluster_tanimoto_significance_report.md   plain-text summary

Run under the project's miniconda:
    LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 \\
        analysis/11_cluster_tanimoto_significance.py
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
PROT_CSV = ROOT / "data" / "proteins.csv"
MOL_CSV = ROOT / "data" / "molecules.csv"
FP_TSV = output("cluster_fingerprints.tsv")
N_PERM = 2000
RNG_SEED = 20260824

# ── Fingerprints for every steroid we can ─────────────────────────────────
print(f"Loading {MOL_CSV.name} ...")
mol = pd.read_csv(MOL_CSV, low_memory=False)
mol = mol[["compound_name", "chebi_id", "smiles"]].copy()
mol = mol[mol["smiles"].astype(str).str.strip().ne("")]
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol.dropna(subset=["chebi_id"])

chebi_ids: list[int] = []
fps = []
for _, row in mol.iterrows():
    m = Chem.MolFromSmiles(str(row["smiles"]))
    if m is None:
        continue
    fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
    chebi_ids.append(int(row["chebi_id"]))
n_pool = len(fps)
idx_of_chebi = {c: i for i, c in enumerate(chebi_ids)}
print(f"  built fingerprints for {n_pool:,} unique ChEBI IDs")

# ── Precompute the full n_pool × n_pool Tanimoto matrix ─────────────────
print("Precomputing full pairwise Tanimoto matrix ...")
T = np.zeros((n_pool, n_pool), dtype=np.float32)
for i in range(n_pool):
    row = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
    T[i, :] = np.asarray(row, dtype=np.float32)
# Symmetrize (guard against tiny float error)
T = np.maximum(T, T.T)
np.fill_diagonal(T, 1.0)
print(f"  T shape = {T.shape}  ({T.nbytes / 1e6:.2f} MB)")


def mean_pairwise(indices: np.ndarray) -> float:
    """Mean of the upper-triangular Tanimoto values for a set of indices."""
    n = len(indices)
    if n < 2:
        return float("nan")
    sub = T[np.ix_(indices, indices)]
    iu = np.triu_indices(n, k=1)
    return float(sub[iu].mean())


# ── Cluster → substrate index list ────────────────────────────────────────
print(f"Loading {PROT_CSV.name} ...")
prot = pd.read_csv(PROT_CSV, low_memory=False)
prot = prot[["accession", "cluster", "interacting_chebi_ids"]].copy()
prot["cluster"] = prot["cluster"].astype(int)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)

_NUM_RE = re.compile(r"\d+")
def parse_chebi(s: str) -> set[int]:
    return {int(x) for x in _NUM_RE.findall(s)}
prot["chebi_set"] = prot["interacting_chebi_ids"].map(parse_chebi)


def cluster_indices(sub: pd.DataFrame) -> np.ndarray:
    """Unique substrate indices for one protein cluster."""
    all_c: set[int] = set()
    for s in sub["chebi_set"]:
        all_c.update(s)
    idxs = [idx_of_chebi[c] for c in all_c if c in idx_of_chebi]
    return np.asarray(sorted(idxs), dtype=np.int32)


# ── Permutation test per cluster ─────────────────────────────────────────
print(f"Running permutation test (N_PERM={N_PERM}) ...")
rng = np.random.default_rng(RNG_SEED)
rows = []
for cid, sub in prot.groupby("cluster"):
    idxs = cluster_indices(sub)
    n = len(idxs)
    obs = mean_pairwise(idxs) if n >= 2 else float("nan")

    if n < 2:
        rows.append(dict(
            cluster=int(cid), n_substrates=n,
            observed=obs, null_mean=np.nan, null_std=np.nan,
            z_score=np.nan,
            p_specialist=np.nan, p_promiscuous=np.nan, p_two_sided=np.nan,
        ))
        continue

    # Draw N_PERM random samples of size n from the pool
    null = np.empty(N_PERM, dtype=np.float32)
    for k in range(N_PERM):
        pick = rng.choice(n_pool, size=n, replace=False)
        null[k] = mean_pairwise(pick)

    null_mean = float(null.mean())
    null_std = float(null.std(ddof=1))
    z = (obs - null_mean) / null_std if null_std > 0 else np.nan
    # +1 to numerator + denominator gives conservative empirical p (avoids p=0)
    p_specialist = float((np.sum(null >= obs) + 1) / (N_PERM + 1))
    p_promiscuous = float((np.sum(null <= obs) + 1) / (N_PERM + 1))
    p_two_sided = 2 * min(p_specialist, p_promiscuous)

    rows.append(dict(
        cluster=int(cid), n_substrates=int(n),
        observed=round(obs, 4),
        null_mean=round(null_mean, 4),
        null_std=round(null_std, 4),
        z_score=round(z, 3),
        p_specialist=round(p_specialist, 5),
        p_promiscuous=round(p_promiscuous, 5),
        p_two_sided=round(p_two_sided, 5),
    ))

out = pd.DataFrame(rows).sort_values("cluster")


# ── Benjamini-Hochberg FDR across the 82 clusters (per direction) ────────
def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Return q-values from Benjamini-Hochberg."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    # Monotone non-decreasing
    q = np.minimum.accumulate(q[::-1])[::-1]
    q_out = np.empty_like(q)
    q_out[order] = np.clip(q, 0, 1)
    return q_out


mask = out["p_specialist"].notna()
out.loc[mask, "q_specialist"] = bh_fdr(out.loc[mask, "p_specialist"].to_numpy())
out.loc[mask, "q_promiscuous"] = bh_fdr(out.loc[mask, "p_promiscuous"].to_numpy())
out["q_specialist"] = out["q_specialist"].round(5)
out["q_promiscuous"] = out["q_promiscuous"].round(5)

# ── Attach family names from cluster_fingerprints.tsv ────────────────────
name_by_cluster: dict[int, str] = {}
if FP_TSV.exists():
    fp_df = pd.read_csv(FP_TSV, sep="\t")
    if "cluster" in fp_df.columns and "dominant_stem" in fp_df.columns:
        for _, r in fp_df.iterrows():
            name_by_cluster[int(r["cluster"])] = str(r["dominant_stem"])[:70]
out["family"] = out["cluster"].map(lambda c: name_by_cluster.get(int(c), "(unnamed)"))

# Reorder columns for readability
out = out[[
    "cluster", "family", "n_substrates",
    "observed", "null_mean", "null_std", "z_score",
    "p_specialist", "q_specialist",
    "p_promiscuous", "q_promiscuous",
    "p_two_sided",
]]

tsv_path = output("cluster_tanimoto_significance.tsv")
out.to_csv(tsv_path, sep="\t", index=False)
print(f"  wrote {tsv_path.name}")

# ── Volcano plot: z-score vs -log10 q ───────────────────────────────────
print("Rendering volcano plot ...")
valid = out.dropna(subset=["z_score"])
# Choose the smaller q per row (dominant direction)
q_dir = np.where(
    valid["z_score"] > 0,
    valid["q_specialist"].fillna(1.0),
    valid["q_promiscuous"].fillna(1.0),
)
neglogq = -np.log10(np.clip(q_dir, 1e-4, 1.0))

fig, ax = plt.subplots(figsize=(9, 6))
colors = np.where(
    (valid["z_score"] > 0) & (q_dir < 0.05), "#0369a1",
    np.where((valid["z_score"] < 0) & (q_dir < 0.05), "#c2410c", "#94a3b8"),
)
ax.scatter(valid["z_score"], neglogq, c=colors, s=60, alpha=0.85,
           edgecolor="black", linewidth=0.4)
ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8,
           label="FDR q = 0.05")
ax.axvline(0, color="black", linewidth=0.4)
ax.set_xlabel("Z-score (observed − null) / null-std")
ax.set_ylabel("−log₁₀(BH-FDR q)")
ax.set_title("Substrate coherence per cluster — permutation test vs random pool")

# Annotate top hits in each direction
tag_top = valid.assign(_q=q_dir).nlargest(6, "z_score").head(6)
tag_bot = valid.assign(_q=q_dir).nsmallest(6, "z_score").head(6)
for _, r in pd.concat([tag_top, tag_bot]).iterrows():
    nlq = -np.log10(np.clip(r["_q"], 1e-4, 1.0))
    ax.annotate(f"C{int(r['cluster'])+1}",
                (r["z_score"], nlq),
                fontsize=8, xytext=(4, 3), textcoords="offset points")

# Legend for colors
from matplotlib.patches import Patch
handles = [
    Patch(facecolor="#0369a1", edgecolor="black", label="Specialist (q<0.05)"),
    Patch(facecolor="#c2410c", edgecolor="black", label="Promiscuous (q<0.05)"),
    Patch(facecolor="#94a3b8", edgecolor="black", label="Not significant"),
]
ax.legend(handles=handles, loc="upper left", fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(output("cluster_tanimoto_significance_volcano.pdf"))
fig.savefig(output("cluster_tanimoto_significance_volcano.png"), dpi=180)
print("  wrote cluster_tanimoto_significance_volcano.{pdf,png}")

# ── Report ───────────────────────────────────────────────────────────────
print("Writing report ...")
n_valid = int(mask.sum())
n_sig_spec = int(((out["z_score"] > 0) & (out["q_specialist"] < 0.05)).sum())
n_sig_prom = int(((out["z_score"] < 0) & (out["q_promiscuous"] < 0.05)).sum())

top_spec = out[(out["z_score"] > 0) & (out["q_specialist"] < 0.05)] \
    .sort_values("z_score", ascending=False).head(10)
top_prom = out[(out["z_score"] < 0) & (out["q_promiscuous"] < 0.05)] \
    .sort_values("z_score", ascending=True).head(10)

lines = [
    "# Cluster substrate-Tanimoto significance",
    "",
    "## Method",
    "",
    f"For each cluster we take the observed mean pairwise Tanimoto over its "
    f"unique substrate set, then draw **{N_PERM:,} random substrate samples "
    f"of the same size** from the full molecule pool (n = {n_pool}) and "
    f"compute the same mean each time to build a null distribution. The "
    f"empirical p-value is `(#null ≥ observed + 1) / (N_perm + 1)` "
    f"(specialist direction) or `(#null ≤ observed + 1) / (N_perm + 1)` "
    f"(promiscuous direction). We correct each direction separately by "
    f"Benjamini-Hochberg FDR.",
    "",
    "**Z-score sign convention:** positive = observed > null (specialist), "
    "negative = observed < null (promiscuous). Larger magnitudes mean "
    "the cluster is farther from what random sampling would produce.",
    "",
    "## Summary",
    "",
    f"- Testable clusters (≥2 substrates): **{n_valid} of 82**",
    f"- Significantly SPECIALIST (q_specialist < 0.05, z > 0): **{n_sig_spec}**",
    f"- Significantly PROMISCUOUS (q_promiscuous < 0.05, z < 0): **{n_sig_prom}**",
    "",
    "---",
    "",
    "## Top 10 significantly SPECIALIST protein families",
    "*(substrates are chemically MORE similar than random pools of the same size)*",
    "",
    "| Protein family | # substrates | Observed | Null mean | Z | q |",
    "|---|---|---|---|---|---|",
]
for _, r in top_spec.iterrows():
    lines.append(
        f"| **{r['family']}** | {int(r['n_substrates'])} | "
        f"{r['observed']:.3f} | {r['null_mean']:.3f} | {r['z_score']:+.2f} | "
        f"{r['q_specialist']:.1e} |"
    )

lines += [
    "",
    "## Top 10 significantly PROMISCUOUS protein families",
    "*(substrates are chemically MORE diverse than random pools of the same size)*",
    "",
    "| Protein family | # substrates | Observed | Null mean | Z | q |",
    "|---|---|---|---|---|---|",
]
for _, r in top_prom.iterrows():
    lines.append(
        f"| **{r['family']}** | {int(r['n_substrates'])} | "
        f"{r['observed']:.3f} | {r['null_mean']:.3f} | {r['z_score']:+.2f} | "
        f"{r['q_promiscuous']:.1e} |"
    )

lines += [
    "",
    "---",
    "",
    "## Files",
    "- `cluster_tanimoto_significance.tsv` — full 82-row table with z-scores + "
    "q-values in both directions.",
    "- `cluster_tanimoto_significance_volcano.{pdf,png}` — Z vs −log₁₀(q).",
]
(output("cluster_tanimoto_significance_report.md")).write_text(
    "\n".join(lines), encoding="utf-8"
)
print("  wrote cluster_tanimoto_significance_report.md")

print("\nDone.")
