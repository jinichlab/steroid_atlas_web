"""Cytochrome P450 superfamily landscape across the atlas.

For every cluster, count how many proteins are P450s (either the string
"cytochrome P450" appears in `protein_names`, or `gene_names` matches
the standard `CYP<digits>[A-Z]?\\d*` nomenclature). Clusters where P450s
are the majority are the "P450 clusters" of the atlas.

For every P450 cluster:
  * dominant CYP family (CYP3 vs CYP7 vs CYP46, etc.)
  * substrate diversity (mean pairwise Tanimoto ECFP4)
  * # unique substrates
  * top organisms
  * UMAP centroid

Outputs (under analysis/):
    cyp_cluster_landscape.tsv               per-P450-cluster summary
    cyp_family_summary.tsv                  per-CYP-family aggregate
    cyp_landscape_umap.pdf/png              UMAP colored by CYP family
    cyp_landscape_scatter.pdf/png           tanimoto vs #substrates
                                             colored by family
    cyp_landscape_report.md                 plain-text narrative summary
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

CYP_MAJORITY_THRESHOLD = 0.50  # cluster is a "CYP cluster" if ≥50% are P450s
DOMINANT_FAM_THRESHOLD = 0.30  # cluster gets a family label if ≥30% share it

# ── Load fingerprints ────────────────────────────────────────────────────
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
print(f"Loaded {len(fp_by_chebi):,} substrate fingerprints")

# ── Load proteins ────────────────────────────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot = prot[[
    "accession", "protein_names", "gene_names", "organism",
    "length_aa", "ec_numbers",
    "interacting_chebi_ids", "umap_1", "umap_2", "cluster",
]].copy()
prot["cluster"] = prot["cluster"].astype(int)
prot["protein_names"] = prot["protein_names"].fillna("").astype(str)
prot["gene_names"] = prot["gene_names"].fillna("").astype(str)
prot["ec_numbers"] = prot["ec_numbers"].fillna("").astype(str)
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)
print(f"Loaded {len(prot):,} proteins across {prot['cluster'].nunique()} clusters")

_NUM_RE = re.compile(r"\d+")
prot["chebi_set"] = prot["interacting_chebi_ids"].map(
    lambda s: {int(x) for x in _NUM_RE.findall(s)}
)

# ── Detect P450 proteins + extract CYP family ────────────────────────────
_CYP_RE = re.compile(r"CYP(\d+)([A-Z])?(\d*)")


def is_p450(row) -> bool:
    """Protein is a P450 if any of: name contains 'cytochrome P450', or gene
    matches the CYP nomenclature."""
    n = row["protein_names"].lower()
    if "cytochrome p450" in n or "cytochrome p-450" in n:
        return True
    if _CYP_RE.search(row["gene_names"].upper()):
        return True
    return False


def cyp_family(row) -> str:
    """Return CYP family label ('CYP3', 'CYP7', 'CYP46', ...) or ''."""
    m = _CYP_RE.search(row["gene_names"].upper())
    if m:
        return f"CYP{m.group(1)}"
    # Try to pull a family digit from protein_names (e.g. 'Cytochrome P450 46A1')
    m2 = re.search(r"cytochrome\s+p[-\s]?450\s+(\d+)", row["protein_names"].lower())
    if m2:
        return f"CYP{m2.group(1)}"
    return ""


prot["is_cyp"] = prot.apply(is_p450, axis=1)
prot["cyp_family"] = prot.apply(cyp_family, axis=1)
n_cyps = int(prot["is_cyp"].sum())
n_named = int((prot["cyp_family"] != "").sum())
print(f"  {n_cyps:,} P450s in the atlas ({100*n_cyps/len(prot):.1f}%)")
print(f"  {n_named:,} with a resolvable CYP family designation")

top_fams = prot[prot["cyp_family"] != ""]["cyp_family"].value_counts().head(15)
print("  Top 15 CYP families:")
for f, c in top_fams.items():
    print(f"    {f:8s}  {c:,d}")

# ── Per-cluster CYP composition ──────────────────────────────────────────
print("\nComputing per-cluster CYP composition ...")
cyp_cluster_rows = []
for cid, sub in prot.groupby("cluster"):
    n = len(sub)
    n_cyp = int(sub["is_cyp"].sum())
    frac_cyp = n_cyp / n if n else 0.0
    fam_counts = sub[sub["cyp_family"] != ""]["cyp_family"].value_counts()
    dom_fam, dom_frac = "", 0.0
    if not fam_counts.empty:
        dom_fam = fam_counts.index[0]
        dom_frac = fam_counts.iloc[0] / n

    # substrate diversity: mean pairwise Tanimoto within this cluster
    all_chebis: set[int] = set()
    for s in sub["chebi_set"]:
        all_chebis.update(s)
    known = sorted(c for c in all_chebis if c in fp_by_chebi)
    if len(known) >= 2:
        fps = [fp_by_chebi[c] for c in known]
        sims = []
        for i in range(len(known) - 1):
            sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
        mean_tanimoto = float(np.mean(sims))
    else:
        mean_tanimoto = float("nan")

    cyp_cluster_rows.append(dict(
        cluster=int(cid),
        display_id=int(cid) + 1,
        n_proteins=n,
        n_cyps=n_cyp,
        pct_cyps=round(100 * frac_cyp, 1),
        dominant_family=dom_fam,
        dominant_family_pct=round(100 * dom_frac, 1),
        n_substrates=len(known),
        mean_tanimoto=round(mean_tanimoto, 4)
        if not np.isnan(mean_tanimoto) else np.nan,
        umap_cx=round(sub["umap_1"].mean(), 3),
        umap_cy=round(sub["umap_2"].mean(), 3),
    ))

cyp_df = pd.DataFrame(cyp_cluster_rows)
cyp_df["is_cyp_cluster"] = cyp_df["pct_cyps"] >= 100 * CYP_MAJORITY_THRESHOLD

n_cyp_clusters = int(cyp_df["is_cyp_cluster"].sum())
print(f"  {n_cyp_clusters} of 82 clusters are majority-P450 "
      f"(≥{int(100*CYP_MAJORITY_THRESHOLD)}% CYPs)")

# Save
cyp_df.sort_values(["is_cyp_cluster", "pct_cyps"],
                   ascending=[False, False]).to_csv(
    output("cyp_cluster_landscape.tsv"), sep="\t", index=False)
print(f"  wrote cyp_cluster_landscape.tsv")

# ── Per-CYP-family aggregate ─────────────────────────────────────────────
p450_only = prot[prot["cyp_family"] != ""]
fam_rows = []
for fam, sub in p450_only.groupby("cyp_family"):
    if len(sub) < 5:
        continue
    # which clusters?
    cl_counts = sub["cluster"].value_counts()
    fam_rows.append(dict(
        cyp_family=fam,
        n_proteins=len(sub),
        n_clusters_present=int((cl_counts >= 3).sum()),
        top_cluster=f"C{int(cl_counts.index[0]) + 1}",
        top_cluster_pct=round(100 * cl_counts.iloc[0] / len(sub), 1),
    ))
fam_df = pd.DataFrame(fam_rows).sort_values("n_proteins", ascending=False)
fam_df.to_csv(output("cyp_family_summary.tsv"), sep="\t", index=False)
print(f"  wrote cyp_family_summary.tsv ({len(fam_df)} families with ≥5 proteins)")

# ── Figure 1: UMAP with CYP clusters highlighted + coloured by family ───
print("\nRendering UMAP CYP landscape ...")
fig, ax = plt.subplots(figsize=(10, 8))

# Background: all proteins in light grey
ax.scatter(prot["umap_1"], prot["umap_2"], s=1.5, c="#e2e8f0",
           alpha=0.5, linewidth=0)

# Foreground: P450 clusters, coloured by dominant family
cyp_clusters = cyp_df[cyp_df["is_cyp_cluster"]].copy()
# Top 8 families get named colours; others become 'Other CYP'
top8 = cyp_clusters["dominant_family"].value_counts().head(8).index.to_list()
palette = plt.get_cmap("tab10")
color_for = {fam: palette(i) for i, fam in enumerate(top8)}

for _, row in cyp_clusters.iterrows():
    cid = int(row["cluster"])
    sub = prot[prot["cluster"] == cid]
    fam = row["dominant_family"]
    col = color_for.get(fam, "#94a3b8")
    ax.scatter(sub["umap_1"], sub["umap_2"], s=6, c=[col],
               alpha=0.85, edgecolor="black", linewidth=0.15,
               label=fam if fam in top8 else None)
    ax.text(row["umap_cx"], row["umap_cy"], f"C{int(row['display_id'])}",
            fontsize=7, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.15",
                      facecolor="white", edgecolor="black",
                      alpha=0.85, linewidth=0.4))

# One legend entry per family
handles, labels = ax.get_legend_handles_labels()
seen = set(); dedup_h, dedup_l = [], []
for h, l in zip(handles, labels):
    if l in seen: continue
    seen.add(l); dedup_h.append(h); dedup_l.append(l)
ax.legend(dedup_h, dedup_l, loc="upper left", fontsize=8,
          title="Dominant CYP family", title_fontsize=9)

ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.set_title(
    f"Cytochrome P450 clusters across the Steroid Atlas UMAP\n"
    f"({n_cyps:,} P450 proteins forming {n_cyp_clusters} majority-P450 "
    f"clusters of the {prot['cluster'].nunique()})"
)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(output("cyp_landscape_umap.pdf"))
fig.savefig(output("cyp_landscape_umap.png"), dpi=180, bbox_inches="tight")
print("  wrote cyp_landscape_umap.{pdf,png}")

# ── Figure 2: substrate breadth vs chemistry tightness ──────────────────
print("Rendering substrate breadth vs chemistry tightness ...")
fig, ax = plt.subplots(figsize=(9, 6))
plot_df = cyp_clusters.dropna(subset=["mean_tanimoto"])
for fam in top8:
    fam_rows = plot_df[plot_df["dominant_family"] == fam]
    ax.scatter(fam_rows["n_substrates"], fam_rows["mean_tanimoto"],
               s=fam_rows["n_proteins"] * 0.6 + 30,
               c=[color_for[fam]], alpha=0.85,
               edgecolor="black", linewidth=0.4, label=fam)
other = plot_df[~plot_df["dominant_family"].isin(top8)]
if not other.empty:
    ax.scatter(other["n_substrates"], other["mean_tanimoto"],
               s=other["n_proteins"] * 0.6 + 30,
               c="#94a3b8", alpha=0.6,
               edgecolor="black", linewidth=0.4, label="Other CYP")

# Label each cluster with its ID
for _, r in plot_df.iterrows():
    ax.annotate(f"C{int(r['display_id'])}",
                (r["n_substrates"], r["mean_tanimoto"]),
                fontsize=7, xytext=(4, 3), textcoords="offset points")

ax.set_xlabel("# unique substrates (fingerprinted)")
ax.set_ylabel("Mean pairwise substrate Tanimoto")
ax.set_title("P450 cluster landscape — substrate breadth × chemistry tightness")
ax.legend(loc="upper right", fontsize=8, title="Dominant CYP family",
          title_fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(output("cyp_landscape_scatter.pdf"))
fig.savefig(output("cyp_landscape_scatter.png"), dpi=180, bbox_inches="tight")
print("  wrote cyp_landscape_scatter.{pdf,png}")

# ── Report ───────────────────────────────────────────────────────────────
lines = [
    "# Cytochrome P450 landscape across the Steroid Atlas",
    "",
    "## Motivation",
    "",
    "The cytochrome P450 (CYP) superfamily is the backbone of steroid "
    "biosynthesis, catabolism, and drug metabolism. Different CYP families "
    "have wildly different substrate specificities — CYP19A1 aromatizes a "
    "narrow set of androgens, CYP7A1 hydroxylates cholesterol at C7, "
    "CYP3A4 metabolises hundreds of xenobiotic and endogenous steroids. "
    "**If the atlas is capturing real functional biology, P450 subfamilies "
    "should split across many clusters, not collapse into one.**",
    "",
    "## Coverage",
    "",
    f"- **P450 proteins in atlas:** {n_cyps:,} of {len(prot):,} "
    f"({100*n_cyps/len(prot):.1f}%)",
    f"- **Proteins with a resolvable CYP<family> designation:** {n_named:,}",
    f"- **Majority-P450 clusters (≥{int(100*CYP_MAJORITY_THRESHOLD)}% CYPs):** "
    f"{n_cyp_clusters} of {prot['cluster'].nunique()}",
    "",
    "## Top CYP families in the atlas",
    "",
    "| Family | # proteins | # clusters (with ≥3 members) | Top cluster |",
    "|---|---|---|---|",
]
for _, r in fam_df.head(20).iterrows():
    lines.append(
        f"| **{r['cyp_family']}** | {int(r['n_proteins']):,} | "
        f"{int(r['n_clusters_present'])} | {r['top_cluster']} "
        f"({r['top_cluster_pct']:.0f}%) |"
    )

lines += [
    "",
    "## P450 clusters — the atlas landscape",
    "",
    f"({n_cyp_clusters} clusters where ≥{int(100*CYP_MAJORITY_THRESHOLD)}% "
    "of proteins are P450s, sorted by substrate diversity)",
    "",
    "| Cluster | Dominant family | # proteins | # substrates | mean Tanimoto |",
    "|---|---|---|---|---|",
]
for _, r in cyp_clusters.sort_values("mean_tanimoto",
                                     ascending=False).iterrows():
    lines.append(
        f"| C{int(r['display_id'])} | "
        f"{r['dominant_family']} ({r['dominant_family_pct']:.0f}%) | "
        f"{int(r['n_proteins'])} | {int(r['n_substrates'])} | "
        f"{r['mean_tanimoto']:.3f} |"
    )

lines += [
    "",
    "## Interpretation",
    "",
    "The P450 superfamily is not a single blob on the atlas — it fans out "
    f"across {n_cyp_clusters} clusters. Some observations:",
    "",
    "- Multiple clusters can share the same CYP-family designation "
    "(e.g. multiple CYP3 or CYP46 clusters). These are typically "
    "**taxonomic variants** of the same enzyme — vertebrate paralogs, or "
    "orthologs in fish vs mammal vs reptile. The ProtT5 embedding "
    "distinguishes them at the sequence level even though their substrate "
    "chemistry may be near-identical (see the Cholesterol 24-hydroxylase "
    "case study, `case_cholesterol24hydroxylase.md`).",
    "",
    "- Some P450 clusters are **narrow specialists** (high mean Tanimoto, "
    "few substrates) — e.g. CYP85 (brassinosteroid biosynthesis) or the "
    "steroidogenic CYP19 / CYP17 / CYP21 / CYP11 clusters.",
    "",
    "- Others are **broad-specificity metabolizers** (low mean Tanimoto, "
    "many substrates) — CYP3A being the canonical example.",
    "",
    "The two figures (`cyp_landscape_umap` and `cyp_landscape_scatter`) "
    "show this fan-out visually and cluster-by-cluster, respectively.",
    "",
    "## Files",
    "- `cyp_cluster_landscape.tsv` — every cluster with its CYP composition + "
    "substrate stats",
    "- `cyp_family_summary.tsv` — one row per CYP<family> with proteins + "
    "clusters",
    "- `cyp_landscape_umap.{pdf,png}` — UMAP scatter with CYP clusters "
    "coloured by dominant family",
    "- `cyp_landscape_scatter.{pdf,png}` — substrate breadth vs Tanimoto for "
    "the CYP clusters",
]
(output("cyp_landscape_report.md")).write_text(
    "\n".join(lines), encoding="utf-8"
)
print("  wrote cyp_landscape_report.md")

print("\nDone.")
