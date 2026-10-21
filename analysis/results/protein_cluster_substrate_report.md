# Substrate chemical diversity — 82 protein clusters

## What this measures

For each protein cluster we look at every steroid its member enzymes are known to bind or transform, and ask: **how chemically similar are those steroids to each other?**

The similarity score is **Tanimoto over ECFP4 fingerprints** — a widely-used cheminformatics measure that compares which chemical substructures two molecules share.

**Scale: 0 (completely different chemistry) → 1 (identical molecules).** **Higher = the enzymes in the cluster act on a narrower / more uniform chemical class (specialists).** Lower = they act on a chemically diverse range of steroids (broad-specificity / promiscuous enzymes).

## Two flavors reported

- **Raw mean Tanimoto** — mean over ALL pairs in the cluster's substrate pool. Size-invariant in principle, but with N=2–5 substrates a lucky pair can inflate the score.
- **Size-adjusted mean Tanimoto** — for every cluster with ≥10 substrates, we draw 200 random subsamples of exactly 10 substrates each, compute the mean pairwise Tanimoto on each subsample, and average. **This is the fair cross-cluster ranking.**

Coverage: **12,825 of 14,089** proteins have at least one substrate we could fingerprint. **58** of 82 clusters have ≥10 substrates and are comparable in the size-adjusted view.

---

## Top 10 SPECIALIST protein families
*(size-adjusted mean Tanimoto — narrowest substrate chemistry)*

| Protein family | # proteins | # substrates | Tanimoto (adj) |
|---|---|---|---|
| **Bile salt export pump** | 83 | 10 | 0.581 |
| **Acyl-coenzyme A diphosphatase NUDT19** | 252 | 17 | 0.572 |
| **G-protein coupled estrogen receptor 1** | 159 | 10 | 0.492 |
| **Cytochrome P450 85A1** | 11 | 29 | 0.453 |
| **Cholesterol 24-hydroxylase** | 200 | 20 | 0.437 |
| **Non-lysosomal glucosylceramidase** | 50 | 10 | 0.421 |
| **Glucosylceramidase** | 46 | 11 | 0.403 |
| **Caveolin** | 77 | 12 | 0.400 |
| **Sterol O-acyltransferase 1** | 46 | 61 | 0.400 |
| **Cholesterol 24-hydroxylase** | 713 | 26 | 0.375 |

## Top 10 BROAD-SPECIFICITY protein families
*(size-adjusted mean Tanimoto — widest substrate chemistry)*

| Protein family | # proteins | # substrates | Tanimoto (adj) |
|---|---|---|---|
| **Oxysterol-binding protein-related protein 8** | 110 | 65 | 0.252 |
| **11-beta-hydroxysteroid dehydrogenase 1** | 204 | 87 | 0.252 |
| **G-protein coupled estrogen receptor 1** | 66 | 27 | 0.258 |
| **3-hydroxyacyl-CoA dehydrogenase type-2** | 63 | 108 | 0.259 |
| **Cytochrome P450 3A** | 73 | 42 | 0.262 |
| **Cholesterol 24-hydroxylase** | 182 | 44 | 0.268 |
| **Hormone-sensitive lipase** | 299 | 39 | 0.271 |
| **Bile acid receptor** | 82 | 29 | 0.273 |
| **Carboxylic ester hydrolase** | 198 | 99 | 0.277 |
| **Carboxylic ester hydrolase** | 368 | 11 | 0.278 |

---

## Files
- `protein_cluster_substrate_tanimoto.tsv` — one row per cluster (82). Columns include both `mean_tanimoto` (raw) and `mean_tanimoto_adj` (subsampled).
- `per_protein_substrate_tanimoto.tsv` — one row per protein with ≥2 known substrates.
- `protein_cluster_substrate_boxplot.{pdf,png}` — box plot of raw pairwise Tanimoto for the 40 largest clusters.
- `protein_cluster_substrate_scatter.{pdf,png}` — substrate breadth vs chemical tightness.

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
`within_cluster_tanimoto_by_type_report.md` for the full rationale.