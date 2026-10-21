"""Case study: WHY do 5 different UMAP clusters all get named
'Cholesterol 24-hydroxylase' by the fingerprint script?

Compare C6, C18, C26, C48, C76 across:
  - taxonomic composition (top organisms, kingdom hint)
  - protein length distribution (soluble vs membrane vs multi-domain hint)
  - UMAP centroid distance (ProtT5 embedding divergence)
  - substrate ChEBI overlap between cluster pairs (Jaccard)
  - within- vs across-cluster mean substrate Tanimoto

Outputs (under analysis/):
    case_cholesterol24hydroxylase.md              plain-text report
    case_cholesterol24hydroxylase.pdf/png         4-panel comparison figure
    case_cholesterol24hydroxylase_composition.tsv per-cluster breakdown
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

# The 5 Cholesterol 24-hydroxylase clusters (0-indexed → display_id = c+1)
CASE_CLUSTERS = [5, 17, 25, 47, 75]   # display: C6, C18, C26, C48, C76
CASE_LABEL = "Cholesterol 24-hydroxylase"

# ── Load fingerprints for substrate Tanimoto ──────────────────────────────
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

# ── Load proteins for the 5 case clusters ────────────────────────────────
prot = pd.read_csv(ROOT / "data" / "proteins.csv", low_memory=False)
prot = prot[prot["cluster"].isin(CASE_CLUSTERS)].copy()
prot["cluster"] = prot["cluster"].astype(int)
prot["length_aa"] = prot["length_aa"].astype(int)
prot["organism"] = prot["organism"].astype(str).fillna("")
prot["interacting_chebi_ids"] = prot["interacting_chebi_ids"].fillna("").astype(str)
print(f"Loaded {len(prot):,} proteins across {len(CASE_CLUSTERS)} clusters")

_NUM_RE = re.compile(r"\d+")


def parse_chebi(s: str) -> set[int]:
    return {int(x) for x in _NUM_RE.findall(s)}


prot["chebi_set"] = prot["interacting_chebi_ids"].map(parse_chebi)


# Simple domain-of-life hint from the organism string
def kingdom_hint(org: str) -> str:
    """Loose kingdom guess from UniProt organism string."""
    o = org.lower()
    # Very common animal / mammal patterns
    for tok in ("homo ", "mus ", "rattus ", "bos ", "sus ", "canis ",
                "felis ", "gallus ", "danio ", "xenopus ", "macaca "):
        if tok in o:
            return "Metazoa (mammal/vertebrate)"
    if "sapiens" in o:
        return "Metazoa (mammal/vertebrate)"
    # Bacterial hints
    for tok in ("escherichia", "staphyloco", "streptoco", "bacillus",
                "mycobact", "clostrid", "lactoba", "salmonella",
                "pseudomonas", "vibrio", "corynebac", "listeria",
                "helicobact", "enterococc"):
        if tok in o:
            return "Bacteria"
    # Fungal hints
    for tok in ("saccharomy", "candida", "aspergill", "penicill",
                "neurospora", "schizosac", "kluyveromy", "yarrowia"):
        if tok in o:
            return "Fungi"
    # Plant hints
    for tok in ("arabidopsis", "oryza ", "zea ", "solanum",
                "nicotiana", "vitis", "medicago"):
        if tok in o:
            return "Viridiplantae"
    # Broad taxa often present in UniProt organism strings
    if any(t in o for t in ("bacterium", "phage", "virus")):
        return "Bacteria" if "bacterium" in o else "Virus"
    return "Other/unknown"


prot["kingdom"] = prot["organism"].map(kingdom_hint)

# ── Per-cluster composition ──────────────────────────────────────────────
rows = []
for cid in CASE_CLUSTERS:
    sub = prot[prot["cluster"] == cid]
    # top organisms
    top_orgs = sub["organism"].value_counts().head(5)
    kingdom_frac = sub["kingdom"].value_counts(normalize=True).round(3)
    rows.append(dict(
        cluster=cid + 1,
        n_proteins=len(sub),
        median_length=int(sub["length_aa"].median()),
        length_p10=int(np.percentile(sub["length_aa"], 10)),
        length_p90=int(np.percentile(sub["length_aa"], 90)),
        umap_cx=round(sub["umap_1"].mean(), 3),
        umap_cy=round(sub["umap_2"].mean(), 3),
        top_kingdom=kingdom_frac.index[0] if not kingdom_frac.empty else "",
        top_kingdom_pct=(int(100 * kingdom_frac.iloc[0])
                         if not kingdom_frac.empty else 0),
        top_organism=top_orgs.index[0] if not top_orgs.empty else "",
        top_organism_n=int(top_orgs.iloc[0]) if not top_orgs.empty else 0,
        substrate_pool_size=len({c for s in sub["chebi_set"] for c in s
                                 if c in fp_by_chebi}),
    ))
comp = pd.DataFrame(rows)
comp.to_csv(output("case_cholesterol24hydroxylase_composition.tsv"),
            sep="\t", index=False)
print(comp.to_string(index=False))

# ── UMAP centroid pairwise distances ─────────────────────────────────────
print("\nUMAP centroid pairwise distances (proxy for ProtT5 divergence):")
cx = comp["umap_cx"].to_numpy()
cy = comp["umap_cy"].to_numpy()
ids = comp["cluster"].to_numpy()
n_c = len(ids)
umap_dist = np.zeros((n_c, n_c))
for i in range(n_c):
    for j in range(n_c):
        umap_dist[i, j] = np.hypot(cx[i] - cx[j], cy[i] - cy[j])
umap_df = pd.DataFrame(umap_dist,
                       index=[f"C{i}" for i in ids],
                       columns=[f"C{i}" for i in ids]).round(2)
print(umap_df)

# ── Substrate overlap (Jaccard) between cluster pairs ───────────────────
cluster_chebis = {}
for cid in CASE_CLUSTERS:
    sub = prot[prot["cluster"] == cid]
    all_c = set().union(*sub["chebi_set"].to_list())
    cluster_chebis[cid + 1] = {c for c in all_c if c in fp_by_chebi}

print("\nSubstrate ChEBI-set Jaccard overlap between clusters:")
jac = np.zeros((n_c, n_c))
for i, ci in enumerate(ids):
    for j, cj in enumerate(ids):
        A = cluster_chebis[ci]
        B = cluster_chebis[cj]
        u = len(A | B)
        jac[i, j] = len(A & B) / u if u else 0.0
jac_df = pd.DataFrame(jac,
                      index=[f"C{i}" for i in ids],
                      columns=[f"C{i}" for i in ids]).round(2)
print(jac_df)

# ── Mean pairwise Tanimoto: within-cluster vs cross-cluster ─────────────
print("\nSubstrate Tanimoto: within-cluster diagonal, cross-cluster off-diagonal")
tanimat = np.zeros((n_c, n_c))
for i, ci in enumerate(ids):
    A = list(cluster_chebis[ci])
    for j, cj in enumerate(ids):
        B = list(cluster_chebis[cj])
        if not A or not B:
            continue
        sims = []
        for a in A:
            row = DataStructs.BulkTanimotoSimilarity(
                fp_by_chebi[a], [fp_by_chebi[b] for b in B]
            )
            if i == j:
                # within-cluster: take upper-triangle only
                for k, b in enumerate(B):
                    if b != a:
                        sims.append(row[k])
            else:
                sims.extend(row)
        tanimat[i, j] = float(np.mean(sims)) if sims else 0.0
tani_df = pd.DataFrame(tanimat,
                       index=[f"C{i}" for i in ids],
                       columns=[f"C{i}" for i in ids]).round(3)
print(tani_df)

# ── Figure: 4 panels ─────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Panel 1: length distribution per cluster
ax = axes[0, 0]
for cid in CASE_CLUSTERS:
    sub = prot[prot["cluster"] == cid]
    ax.hist(sub["length_aa"], bins=40, alpha=0.55, label=f"C{cid+1} (n={len(sub)})")
ax.set_xlabel("Protein length (aa)")
ax.set_ylabel("count")
ax.set_title(f"Length distribution\n"
             f"{CASE_LABEL} — 5 clusters")
ax.legend(fontsize=9)
ax.set_xlim(0, min(2000, prot["length_aa"].quantile(0.99)))
ax.grid(alpha=0.3)

# Panel 2: kingdom composition (stacked bar)
ax = axes[0, 1]
counts = prot.groupby(["cluster", "kingdom"]).size().unstack(fill_value=0)
counts = counts.reindex(CASE_CLUSTERS)
pcts = counts.div(counts.sum(axis=1), axis=0) * 100
king_order = pcts.sum().sort_values(ascending=False).index.to_list()
pcts = pcts[king_order]
bottoms = np.zeros(len(CASE_CLUSTERS))
kingdom_color = {
    "Metazoa (mammal/vertebrate)": "#0369a1",
    "Bacteria": "#16a34a",
    "Fungi": "#a855f7",
    "Viridiplantae": "#eab308",
    "Virus": "#ec4899",
    "Other/unknown": "#94a3b8",
}
xpos = np.arange(len(CASE_CLUSTERS))
for king in king_order:
    ax.bar(xpos, pcts[king].values, bottom=bottoms,
           color=kingdom_color.get(king, "#94a3b8"),
           edgecolor="black", linewidth=0.4, label=king)
    bottoms += pcts[king].values
ax.set_xticks(xpos)
ax.set_xticklabels([f"C{c+1}" for c in CASE_CLUSTERS])
ax.set_ylabel("% of proteins")
ax.set_title("Kingdom composition per cluster")
ax.set_ylim(0, 100)
ax.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.35),
          ncol=2)
ax.grid(axis="y", alpha=0.3)

# Panel 3: substrate Jaccard heatmap
ax = axes[1, 0]
im = ax.imshow(jac, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(n_c))
ax.set_yticks(range(n_c))
ax.set_xticklabels([f"C{i}" for i in ids])
ax.set_yticklabels([f"C{i}" for i in ids])
for i in range(n_c):
    for j in range(n_c):
        v = jac[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color="white" if v > 0.5 else "black", fontsize=9)
ax.set_title("Substrate overlap (Jaccard)")
plt.colorbar(im, ax=ax, shrink=0.75).set_label("Jaccard")

# Panel 4: cross-cluster mean substrate Tanimoto heatmap
ax = axes[1, 1]
im = ax.imshow(tanimat, cmap="coolwarm_r", vmin=0.15, vmax=0.6, aspect="auto")
ax.set_xticks(range(n_c))
ax.set_yticks(range(n_c))
ax.set_xticklabels([f"C{i}" for i in ids])
ax.set_yticklabels([f"C{i}" for i in ids])
for i in range(n_c):
    for j in range(n_c):
        v = tanimat[i, j]
        color = "white" if abs(v - 0.375) > 0.15 else "black"
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color=color, fontsize=9)
ax.set_title("Mean pairwise substrate Tanimoto\n"
             "(diagonal = within, off-diagonal = across)")
plt.colorbar(im, ax=ax, shrink=0.75).set_label("Mean Tanimoto")

fig.suptitle(
    f"Case study: 5 UMAP clusters all called '{CASE_LABEL}'",
    fontsize=13, y=0.995,
)
fig.tight_layout()
fig.savefig(output("case_cholesterol24hydroxylase.pdf"))
fig.savefig(output("case_cholesterol24hydroxylase.png"),
            dpi=180, bbox_inches="tight")
print("\nWrote case_cholesterol24hydroxylase.{pdf,png}")

# ── Report ───────────────────────────────────────────────────────────────
lines = [
    f"# Case study — {CASE_LABEL} (5 clusters)",
    "",
    "## Motivation",
    "",
    "Five separate UMAP clusters (C6, C18, C26, C48, C76 — total 1,353 "
    "proteins) all get labelled `Cholesterol 24-hydroxylase` by the "
    "dominant-name fingerprint. Why does k-means put them in five distinct "
    "clusters rather than one?",
    "",
    "**Short answer:** ProtT5 embeddings encode the amino-acid sequence, "
    "so proteins that share a *recommended enzyme name* but belong to "
    "different sequence families / organisms / lengths land in different "
    "embedding-space clusters. The dominant-name label collapses that "
    "diversity into a single string, but the underlying biology is 5 "
    "distinct enzyme populations.",
    "",
    "## Per-cluster composition",
    "",
    "| Cluster | n proteins | median length (aa) | top kingdom | top organism |",
    "|---|---|---|---|---|",
]
for _, r in comp.iterrows():
    lines.append(
        f"| C{int(r['cluster'])} | {int(r['n_proteins'])} | "
        f"{int(r['median_length'])} ({int(r['length_p10'])}–{int(r['length_p90'])}) | "
        f"{r['top_kingdom']} ({int(r['top_kingdom_pct'])}%) | "
        f"{r['top_organism']} ({int(r['top_organism_n'])}) |"
    )

lines += [
    "",
    "## UMAP centroid distances (ProtT5 embedding proxy)",
    "",
    "```",
    umap_df.to_string(),
    "```",
    "",
    "Higher distance = more divergent embedding = more sequence-level "
    "difference. Clusters that share a name but sit far apart on the UMAP "
    "are the strongest evidence that the name is convergent, not "
    "reflecting a single sequence family.",
    "",
    "## Substrate ChEBI overlap between clusters (Jaccard)",
    "",
    "```",
    jac_df.to_string(),
    "```",
    "",
    "## Mean pairwise substrate Tanimoto (within & cross)",
    "",
    "```",
    tani_df.to_string(),
    "```",
    "",
    "Diagonal = within-cluster mean pairwise Tanimoto. Off-diagonal = "
    "mean pairwise Tanimoto between substrates of cluster i and substrates "
    "of cluster j. If the diagonal >> off-diagonal for a pair, that pair "
    "of clusters really does act on chemically distinct substrate sets.",
    "",
    "## Files",
    "- `case_cholesterol24hydroxylase_composition.tsv` — the 5-row summary "
    "table with lengths, organisms, kingdoms, UMAP centroids.",
    "- `case_cholesterol24hydroxylase.{pdf,png}` — 4-panel comparison "
    "figure (length distribution, kingdom composition, substrate Jaccard, "
    "substrate cross-Tanimoto).",
]
(output("case_cholesterol24hydroxylase.md")).write_text("\n".join(lines),
                                                           encoding="utf-8")
print("Wrote case_cholesterol24hydroxylase.md")
print("\nDone.")
