# Cytochrome P450 landscape across the Steroid Atlas

## Motivation

The cytochrome P450 (CYP) superfamily is the backbone of steroid biosynthesis, catabolism, and drug metabolism. Different CYP families have wildly different substrate specificities — CYP19A1 aromatizes a narrow set of androgens, CYP7A1 hydroxylates cholesterol at C7, CYP3A4 metabolises hundreds of xenobiotic and endogenous steroids. **If the atlas is capturing real functional biology, P450 subfamilies should split across many clusters, not collapse into one.**

## Coverage

- **P450 proteins in atlas:** 2,601 of 14,089 (18.5%)
- **Proteins with a resolvable CYP<family> designation:** 2,592
- **Majority-P450 clusters (≥50% CYPs):** 12 of 82

## Top CYP families in the atlas

| Family | # proteins | # clusters (with ≥3 members) | Top cluster |
|---|---|---|---|
| **CYP46** | 1,202 | 9 | C6 (59%) |
| **CYP1** | 504 | 5 | C20 (88%) |
| **CYP39** | 303 | 3 | C9 (82%) |
| **CYP51** | 211 | 5 | C29 (78%) |
| **CYP17** | 180 | 5 | C46 (72%) |
| **CYP3** | 78 | 3 | C57 (65%) |
| **CYP11** | 23 | 1 | C48 (100%) |
| **CYP2** | 17 | 2 | C24 (65%) |
| **CYP19** | 13 | 1 | C64 (100%) |
| **CYP85** | 12 | 2 | C65 (42%) |
| **CYP7** | 8 | 1 | C69 (100%) |
| **CYP90** | 7 | 1 | C47 (100%) |
| **CYP4** | 6 | 1 | C39 (100%) |

## P450 clusters — the atlas landscape

(12 clusters where ≥50% of proteins are P450s, sorted by substrate diversity)

| Cluster | Dominant family | # proteins | # substrates | mean Tanimoto |
|---|---|---|---|---|
| C47 | CYP90 (64%) | 11 | 29 | 0.453 |
| C26 | CYP46 (100%) | 200 | 20 | 0.434 |
| C6 | CYP46 (100%) | 713 | 26 | 0.378 |
| C29 | CYP51 (100%) | 164 | 12 | 0.368 |
| C20 | CYP1 (99%) | 445 | 26 | 0.319 |
| C48 | CYP46 (60%) | 199 | 66 | 0.301 |
| C46 | CYP17 (90%) | 144 | 25 | 0.295 |
| C34 | CYP39 (82%) | 56 | 38 | 0.291 |
| C71 | CYP39 (50%) | 16 | 17 | 0.282 |
| C9 | CYP39 (93%) | 268 | 39 | 0.277 |
| C18 | CYP46 (69%) | 182 | 44 | 0.272 |
| C57 | CYP3 (70%) | 73 | 42 | 0.263 |

## Interpretation

The P450 superfamily is not a single blob on the atlas — it fans out across 12 clusters. Some observations:

- Multiple clusters can share the same CYP-family designation (e.g. multiple CYP3 or CYP46 clusters). These are typically **taxonomic variants** of the same enzyme — vertebrate paralogs, or orthologs in fish vs mammal vs reptile. The ProtT5 embedding distinguishes them at the sequence level even though their substrate chemistry may be near-identical (see the Cholesterol 24-hydroxylase case study, `case_cholesterol24hydroxylase.md`).

- Some P450 clusters are **narrow specialists** (high mean Tanimoto, few substrates) — e.g. CYP85 (brassinosteroid biosynthesis) or the steroidogenic CYP19 / CYP17 / CYP21 / CYP11 clusters.

- Others are **broad-specificity metabolizers** (low mean Tanimoto, many substrates) — CYP3A being the canonical example.

The two figures (`cyp_landscape_umap` and `cyp_landscape_scatter`) show this fan-out visually and cluster-by-cluster, respectively.

## Files
- `cyp_cluster_landscape.tsv` — every cluster with its CYP composition + substrate stats
- `cyp_family_summary.tsv` — one row per CYP<family> with proteins + clusters
- `cyp_landscape_umap.{pdf,png}` — UMAP scatter with CYP clusters coloured by dominant family
- `cyp_landscape_scatter.{pdf,png}` — substrate breadth vs Tanimoto for the CYP clusters