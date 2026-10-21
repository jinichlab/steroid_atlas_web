"""Substrate chemical-diversity for every protein cluster (all 82).

For each of the 82 protein-UMAP clusters:
  * gather the union of interacting steroid substrates (via ChEBI IDs) across
    every protein in the cluster
  * compute pairwise Tanimoto of ECFP4 fingerprints over that substrate set
  * report the intra-cluster substrate diversity (mean/median/std/min/max)
  * also compute the analogous per-PROTEIN measure ("substrate range of this
    one enzyme") — mean pairwise Tanimoto over the compounds the enzyme
    itself acts on

Reads:
    ../data/proteins.csv         (14,089 rows w/ interacting_chebi_ids + cluster)
    ../data/molecules.csv        (677 rows w/ chebi_id + smiles)

Writes (under analysis/):
    protein_cluster_substrate_tanimoto.tsv     one row per protein cluster (82)
    per_protein_substrate_tanimoto.tsv         one row per protein w/ ≥2 substrates
    protein_cluster_substrate_boxplot.pdf/png  box plot per cluster (top-40 by size)
    protein_cluster_substrate_scatter.pdf/png  substrate-count vs mean-Tanimoto
    protein_cluster_substrate_report.md        plain-text summary

Run under the project's miniconda:
    LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 \\
        analysis/10_protein_cluster_substrate_tanimoto.py
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

# ── Load molecule SMILES + build a ChEBI → fingerprint dictionary ─────────
print(f"Loading {MOL_CSV.name} ...")
mol = pd.read_csv(MOL_CSV, low_memory=False)
mol = mol[["compound_name", "chebi_id", "smiles"]].copy()
mol = mol[mol["smiles"].astype(str).str.strip().ne("")]
mol["chebi_id"] = pd.to_numeric(mol["chebi_id"], errors="coerce").astype("Int64")
mol = mol.dropna(subset=["chebi_id"])
print(f"  {len(mol):,} molecules with a ChEBI ID + SMILES")

fp_by_chebi: dict[int, object] = {}
name_by_chebi: dict[int, str] = {}
for _, row in mol.iterrows():
    m = Chem.MolFromSmiles(str(row["smiles"]))
    if m is None:
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    cid = int(row["chebi_id"])
    fp_by_chebi[cid] = fp
    name_by_chebi[cid] = str(row["compound_name"])
print(f"  built fingerprints for {len(fp_by_chebi):,} ChEBI IDs")

# ── Load proteins + parse their interacting_chebi_ids column ─────────────
print(f"Loading {PROT_CSV.name} ...")
prot = pd.read_csv(PROT_CSV, low_memory=False)
need = ["accession", "protein_names", "gene_names", "organism", "cluster",
        "interacting_chebi_ids"]
missing = [c for c in need if c not in prot.columns]
if missing:
    raise SystemExit(f"proteins.csv missing required columns: {missing}")
prot = prot[need].copy()
prot["cluster"] = prot["cluster"].astype(int)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)
print(f"  {len(prot):,} proteins across {prot['cluster'].nunique()} clusters")

_NUM_RE = re.compile(r"\d+")


def parse_chebi_list(s: str) -> set[int]:
    """Extract every integer ChEBI id from a semicolon/newline-delimited string."""
    if not s:
        return set()
    return {int(x) for x in _NUM_RE.findall(s)}


prot["chebi_set"] = prot["interacting_chebi_ids"].map(parse_chebi_list)

# Discard ChEBI ids we have no fingerprint for (either not in molecules.csv
# or unparseable SMILES).
def filter_known(s: set[int]) -> set[int]:
    return {c for c in s if c in fp_by_chebi}


prot["chebi_set_known"] = prot["chebi_set"].map(filter_known)
prot["n_known_substrates"] = prot["chebi_set_known"].map(len)
print(f"  {(prot['n_known_substrates'] >= 1).sum():,} proteins map to ≥1 known substrate")
print(f"  {(prot['n_known_substrates'] >= 2).sum():,} proteins map to ≥2 known substrates")

# ── Tanimoto helpers ─────────────────────────────────────────────────────

def pairwise_tanimoto(chebis: list[int]) -> np.ndarray:
    """Upper-triangle pairwise Tanimoto vector over the given ChEBI list."""
    n = len(chebis)
    if n < 2:
        return np.array([], dtype=np.float32)
    fps = [fp_by_chebi[c] for c in chebis]
    sims = []
    for i in range(n - 1):
        row = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        sims.extend(row)
    return np.asarray(sims, dtype=np.float32)


def stats_of(sims: np.ndarray) -> dict:
    if sims.size == 0:
        return dict(mean=np.nan, median=np.nan, std=np.nan, min=np.nan, max=np.nan)
    return dict(
        mean=float(sims.mean()),
        median=float(np.median(sims)),
        std=float(sims.std()),
        min=float(sims.min()),
        max=float(sims.max()),
    )

# ── Per-protein substrate range ──────────────────────────────────────────
print("Computing per-protein substrate Tanimoto (proteins with ≥2 substrates)...")
per_protein_rows = []
for _, row in prot.iterrows():
    chebis = sorted(row["chebi_set_known"])
    if len(chebis) < 2:
        continue
    sims = pairwise_tanimoto(chebis)
    s = stats_of(sims)
    per_protein_rows.append({
        "accession": row["accession"],
        "protein_names": row["protein_names"],
        "gene_names": row["gene_names"],
        "organism": row["organism"],
        "cluster": int(row["cluster"]),
        "n_substrates": len(chebis),
        "mean_tanimoto": round(s["mean"], 4),
        "median_tanimoto": round(s["median"], 4),
        "std_tanimoto": round(s["std"], 4),
        "min_tanimoto": round(s["min"], 4),
        "max_tanimoto": round(s["max"], 4),
    })
per_protein = pd.DataFrame(per_protein_rows).sort_values(
    ["cluster", "mean_tanimoto"], ascending=[True, False]
)
per_protein.to_csv(output("per_protein_substrate_tanimoto.tsv"),
                   sep="\t", index=False)
print(f"  wrote per_protein_substrate_tanimoto.tsv ({len(per_protein):,} proteins)")

# ── Per-cluster substrate diversity ──────────────────────────────────────
print("Computing per-cluster substrate Tanimoto (union of substrates)...")
cluster_rows = []
cluster_sim_arrays: dict[int, np.ndarray] = {}
for cid, sub in prot.groupby("cluster"):
    # Union of ChEBI IDs across all proteins in this cluster (only ones we
    # have fingerprints for).
    all_chebis: set[int] = set()
    for s in sub["chebi_set_known"]:
        all_chebis.update(s)
    chebis = sorted(all_chebis)
    n_prot_with_sub = int((sub["n_known_substrates"] >= 1).sum())

    sims = pairwise_tanimoto(chebis)
    cluster_sim_arrays[int(cid)] = sims
    st = stats_of(sims)
    cluster_rows.append({
        "cluster": int(cid),
        "display_id": int(cid) + 1,
        "n_proteins_in_cluster": int(len(sub)),
        "n_proteins_with_known_substrate": n_prot_with_sub,
        "n_unique_substrates": len(chebis),
        "n_pairs": int(sims.size),
        "mean_tanimoto": round(st["mean"], 4),
        "median_tanimoto": round(st["median"], 4),
        "std_tanimoto": round(st["std"], 4),
        "min_tanimoto": round(st["min"], 4),
        "max_tanimoto": round(st["max"], 4),
    })
cluster_stats = pd.DataFrame(cluster_rows).sort_values("cluster")
cluster_stats.to_csv(output("protein_cluster_substrate_tanimoto.tsv"),
                     sep="\t", index=False)
print(f"  wrote protein_cluster_substrate_tanimoto.tsv "
      f"({len(cluster_stats)} clusters)")

# ── Figure 1: box plot per cluster (top-40 by protein count) ─────────────
print("Rendering box plot for the 40 largest clusters ...")
top40 = cluster_stats.nlargest(40, "n_proteins_in_cluster").sort_values("cluster")
data = []
labels = []
for _, row in top40.iterrows():
    cid = int(row["cluster"])
    sims = cluster_sim_arrays.get(cid, np.array([]))
    if sims.size == 0:
        continue
    data.append(sims)
    labels.append(f"C{cid+1}\n(np={row['n_proteins_in_cluster']}, "
                  f"ns={int(row['n_unique_substrates'])})")

fig, ax = plt.subplots(figsize=(14, 5))
ax.boxplot(
    data,
    tick_labels=labels,
    showfliers=False,
    patch_artist=True,
    boxprops=dict(facecolor="#93c5fd", edgecolor="#1e40af"),
    medianprops=dict(color="#1e3a8a", linewidth=1.6),
    whiskerprops=dict(color="#1e40af"),
    capprops=dict(color="#1e40af"),
)
ax.set_ylabel("Pairwise Tanimoto (ECFP4)")
ax.set_title(
    "Substrate chemical diversity per protein cluster — 40 largest clusters"
)
ax.set_ylim(0, 1)
ax.grid(axis="y", alpha=0.3)
for lbl in ax.get_xticklabels():
    lbl.set_rotation(0)
    lbl.set_fontsize(7)
fig.tight_layout()
fig.savefig(output("protein_cluster_substrate_boxplot.pdf"))
fig.savefig(output("protein_cluster_substrate_boxplot.png"), dpi=180)
print("  wrote protein_cluster_substrate_boxplot.{pdf,png}")

# ── Figure 2: substrate count vs mean Tanimoto scatter ───────────────────
print("Rendering substrate-count vs mean-Tanimoto scatter ...")
valid = cluster_stats.dropna(subset=["mean_tanimoto"])
fig, ax = plt.subplots(figsize=(7.5, 5))
sizes = np.clip(valid["n_proteins_in_cluster"] * 0.6, 30, 800)
sc = ax.scatter(
    valid["n_unique_substrates"], valid["mean_tanimoto"],
    s=sizes, c=valid["n_proteins_in_cluster"],
    cmap="viridis", edgecolor="black", linewidth=0.4, alpha=0.85,
)
plt.colorbar(sc, ax=ax, label="# proteins in cluster")
# Annotate a handful of extreme points
n = len(valid)
for _, r in valid.nlargest(6, "mean_tanimoto").iterrows():
    ax.annotate(f"C{int(r['cluster'])+1}",
                (r["n_unique_substrates"], r["mean_tanimoto"]),
                fontsize=8, xytext=(4, 3), textcoords="offset points")
for _, r in valid.nsmallest(6, "mean_tanimoto").iterrows():
    ax.annotate(f"C{int(r['cluster'])+1}",
                (r["n_unique_substrates"], r["mean_tanimoto"]),
                fontsize=8, xytext=(4, -8), textcoords="offset points")
ax.set_xlabel("# unique substrate steroids (with known fingerprint)")
ax.set_ylabel("Mean pairwise Tanimoto of substrates")
ax.set_title(
    "Protein clusters by substrate breadth vs chemical tightness"
)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(output("protein_cluster_substrate_scatter.pdf"))
fig.savefig(output("protein_cluster_substrate_scatter.png"), dpi=180)
print("  wrote protein_cluster_substrate_scatter.{pdf,png}")

# ── Size-adjusted Tanimoto (subsample every cluster to the SAME N) ────────
# The concern: a cluster with only 2 substrates that happen to be similar
# looks artificially "tight" versus a cluster with 100 substrates. To make
# clusters comparable, we resample each cluster's substrate pool down to
# SUBSAMPLE_N substrates, compute mean pairwise Tanimoto, and average over
# many bootstrap draws. Clusters with fewer than SUBSAMPLE_N substrates get
# reported as NaN in the adjusted column (they're not comparable).
SUBSAMPLE_N = 10
N_BOOTSTRAP = 200
print(f"Bootstrapped size-adjusted Tanimoto: subsample={SUBSAMPLE_N}, "
      f"bootstrap={N_BOOTSTRAP} ...")
rng = np.random.default_rng(20260824)

adj_mean: dict[int, float] = {}
for cid, sub in prot.groupby("cluster"):
    all_chebis = set()
    for s in sub["chebi_set_known"]:
        all_chebis.update(s)
    chebis = sorted(all_chebis)
    if len(chebis) < SUBSAMPLE_N:
        adj_mean[int(cid)] = float("nan")
        continue
    means = []
    for _ in range(N_BOOTSTRAP):
        pick = rng.choice(len(chebis), size=SUBSAMPLE_N, replace=False)
        sims = pairwise_tanimoto([chebis[i] for i in pick])
        means.append(sims.mean())
    adj_mean[int(cid)] = float(np.mean(means))

cluster_stats["mean_tanimoto_adj"] = cluster_stats["cluster"].map(
    lambda c: round(adj_mean[int(c)], 4) if not np.isnan(adj_mean[int(c)]) else np.nan
)
cluster_stats.to_csv(output("protein_cluster_substrate_tanimoto.tsv"),
                     sep="\t", index=False)
print(f"  wrote protein_cluster_substrate_tanimoto.tsv "
      f"(now includes size-adjusted mean_tanimoto_adj)")

# ── Report ───────────────────────────────────────────────────────────────
print("Writing plain-text report ...")

# Bring in cluster_fingerprints so we can name each protein cluster.
# We want the PROTEIN FAMILY NAME to be the primary column, not the raw ID.
fp_file = output("cluster_fingerprints.tsv")
name_by_cluster: dict[int, str] = {}
if fp_file.exists():
    fp_df = pd.read_csv(fp_file, sep="\t")
    if "cluster" in fp_df.columns and "dominant_stem" in fp_df.columns:
        for _, r in fp_df.iterrows():
            name_by_cluster[int(r["cluster"])] = str(r["dominant_stem"])[:70]


def fam(cid: int) -> str:
    return name_by_cluster.get(cid, "(unnamed)")


N_TOP = 10

lines = [
    "# Substrate chemical diversity — 82 protein clusters",
    "",
    "## What this measures",
    "",
    "For each protein cluster we look at every steroid its member enzymes are "
    "known to bind or transform, and ask: **how chemically similar are those "
    "steroids to each other?**",
    "",
    "The similarity score is **Tanimoto over ECFP4 fingerprints** — a "
    "widely-used cheminformatics measure that compares which chemical "
    "substructures two molecules share.",
    "",
    "**Scale: 0 (completely different chemistry) → 1 (identical molecules).** "
    "**Higher = the enzymes in the cluster act on a narrower / more uniform "
    "chemical class (specialists).** Lower = they act on a chemically diverse "
    "range of steroids (broad-specificity / promiscuous enzymes).",
    "",
    "## Two flavors reported",
    "",
    f"- **Raw mean Tanimoto** — mean over ALL pairs in the cluster's substrate "
    f"pool. Size-invariant in principle, but with N=2–5 substrates a lucky "
    f"pair can inflate the score.",
    f"- **Size-adjusted mean Tanimoto** — for every cluster with ≥{SUBSAMPLE_N} "
    f"substrates, we draw {N_BOOTSTRAP} random subsamples of exactly "
    f"{SUBSAMPLE_N} substrates each, compute the mean pairwise Tanimoto on "
    f"each subsample, and average. **This is the fair cross-cluster ranking.**",
    "",
    "Coverage: **{n_have:,} of {n_tot:,}** proteins have at least one substrate "
    "we could fingerprint. **{n_adj}** of 82 clusters have ≥{ss} substrates and "
    "are comparable in the size-adjusted view.".format(
        n_have=(prot['n_known_substrates'] >= 1).sum(),
        n_tot=len(prot),
        n_adj=int(cluster_stats['mean_tanimoto_adj'].notna().sum()),
        ss=SUBSAMPLE_N,
    ),
    "",
    "---",
    "",
    f"## Top {N_TOP} SPECIALIST protein families",
    "*(size-adjusted mean Tanimoto — narrowest substrate chemistry)*",
    "",
    "| Protein family | # proteins | # substrates | Tanimoto (adj) |",
    "|---|---|---|---|",
]
adj_valid = cluster_stats.dropna(subset=["mean_tanimoto_adj"])
top = adj_valid.nlargest(N_TOP, "mean_tanimoto_adj")
for _, r in top.iterrows():
    lines.append(
        f"| **{fam(int(r['cluster']))}** | "
        f"{int(r['n_proteins_in_cluster'])} | "
        f"{int(r['n_unique_substrates'])} | "
        f"{r['mean_tanimoto_adj']:.3f} |"
    )

lines += [
    "",
    f"## Top {N_TOP} BROAD-SPECIFICITY protein families",
    "*(size-adjusted mean Tanimoto — widest substrate chemistry)*",
    "",
    "| Protein family | # proteins | # substrates | Tanimoto (adj) |",
    "|---|---|---|---|",
]
bot = adj_valid.nsmallest(N_TOP, "mean_tanimoto_adj")
for _, r in bot.iterrows():
    lines.append(
        f"| **{fam(int(r['cluster']))}** | "
        f"{int(r['n_proteins_in_cluster'])} | "
        f"{int(r['n_unique_substrates'])} | "
        f"{r['mean_tanimoto_adj']:.3f} |"
    )

lines += [
    "",
    "---",
    "",
    "## Files",
    "- `protein_cluster_substrate_tanimoto.tsv` — one row per cluster (82). "
    "Columns include both `mean_tanimoto` (raw) and `mean_tanimoto_adj` "
    "(subsampled).",
    "- `per_protein_substrate_tanimoto.tsv` — one row per protein with ≥2 "
    "known substrates.",
    "- `protein_cluster_substrate_boxplot.{pdf,png}` — box plot of raw "
    "pairwise Tanimoto for the 40 largest clusters.",
    "- `protein_cluster_substrate_scatter.{pdf,png}` — substrate breadth vs "
    "chemical tightness.",
]
# Hand-written section appended after the generated report (kept in the
# script so re-runs reproduce it rather than dropping it).
SUBSTRATE_PRODUCT_CAVEAT = '''

## Caveat: substrate / product mixing for enzyme clusters

For catalytic enzymes, the substrate pool used here contains both the
substrates and the products of each Rhea reaction the enzyme is
annotated with. Because most steroid-modifying reactions add or remove a
single functional group, substrate and product typically share ≥90% of
the ECFP4 bits set — their pairwise Tanimoto is often 0.9–1.0, so each
enzyme-catalyzed reaction contributes one near-duplicate pair to the
pool. The estimated impact is a **+0.02–0.05** inflation of the mean
Tanimoto for enzyme clusters, roughly uniform across clusters. Impact on
receptor / transporter clusters is zero (they only bind, they do not
turn substrate over into product).

**What is preserved:** cross-cluster ranking (which cluster is more
chemically coherent than another).

**What is subtly biased:** absolute Tanimoto magnitudes for enzyme
clusters, and direct enzyme-vs-non-enzyme comparisons.

A more rigorous revision would re-fetch each Rhea reaction and separate
input from output ChEBIs. See the corresponding note in
`within_cluster_tanimoto_by_type_report.md` for the full rationale.'''

(output("protein_cluster_substrate_report.md")).write_text(
    "\n".join(lines) + SUBSTRATE_PRODUCT_CAVEAT, encoding="utf-8"
)
print("  wrote protein_cluster_substrate_report.md")

print("\nDone.")
