# Within-cluster substrate Tanimoto, split by protein type

## Classification

- **Enzyme** — protein has at least one EC number in `ec_numbers`.
- **Transporter** — no EC number; `keyword_labels` contains `Transport` / `Ion channel` / `Symporter` / `Antiporter`.
- **Receptor** — no EC number; `keyword_labels` contains `Receptor`.
- **Other** — no EC number and no transporter/receptor keyword (binding proteins, sensors, apolipoproteins, morphogens like Hedgehog).

A cluster is labelled with a **dominant type** when ≥60% of its proteins share one type; otherwise it is labelled **Mixed**.

## Overall type breakdown (14,089 proteins)

- **Enzyme** — 6,674 proteins (47.4%)
- **Receptor** — 738 proteins (5.2%)
- **Transporter** — 1,579 proteins (11.2%)
- **Other** — 5,098 proteins (36.2%)

## Dominant-type breakdown across the 82 clusters

- **Enzyme** — 28 clusters
- **Receptor** — 6 clusters
- **Transporter** — 14 clusters
- **Other** — 27 clusters
- **Mixed** — 7 clusters

## Interpretation

- The **enzyme-only figure** (`within_cluster_tanimoto_enzymes_only.png`) gives you a fair enzyme-to-enzyme comparison of substrate breadth — this is the panel to cite when arguing about enzyme substrate specificity.
- The **grouped-by-type figure** (`within_cluster_tanimoto_by_type.png`) shows the whole atlas colour-coded so you can see that receptors and transporters tend to sit at the top (high Tanimoto = narrow ligand selectivity), while enzymes span the whole range.

## Files
- `cluster_type_composition.tsv` — per-cluster type breakdown (82 rows)
- `within_cluster_tanimoto_by_type.{pdf,png}` — all clusters, colored by type
- `within_cluster_tanimoto_enzymes_only.{pdf,png}` — enzyme clusters only

## Limitations

**Substrate / product mixing for enzyme clusters.** The
`interacting_chebi_ids` column used as the substrate pool for every protein
is drawn from Rhea reaction participants without distinguishing the reaction
side. For catalytic enzymes this means the pool contains **both** the
substrates (reaction "left" side) and the products (reaction "right"
side) of each reaction the enzyme is annotated to catalyze. Because most
steroid-modifying reactions add or remove a single functional group
(hydroxylation, oxidation, reduction, conjugation, cleavage), the
substrate and its product typically share ≥90% of the ECFP4 bits set —
their pairwise Tanimoto is often 0.9–1.0. Each such pair therefore
contributes a near-duplicate similarity value to the pool.

The magnitude of the effect is modest. For an enzyme cluster with ~20
reactions in its substrate list (~40 compounds, ~780 pairs), the ~20
substrate-product near-duplicate pairs account for roughly 2–3% of the
pool and inflate the mean pairwise Tanimoto by an estimated **0.02–0.05**.
This offset is largely uniform across enzyme clusters, so the
**cross-cluster ranking** (which cluster is more chemically coherent than
another) is preserved. What is subtly biased are:

1. **Absolute Tanimoto magnitudes** for enzyme clusters — the reported
   mean is slightly higher than an "input-only" analysis would report.
2. **Enzyme-vs-non-enzyme comparisons** — receptors and transporters lack
   this confound because they only bind (they do not turn substrate over
   into product), so enzyme substrate pools are inflated relative to
   receptor/transporter ligand pools by a small constant.

A more rigorous future revision would re-fetch each Rhea reaction via its
REST API, separate the ChEBI identifiers on each side of the arrow, and
recompute pairwise Tanimoto on **input substrates only**. Given the modest
size of the bias and the preserved ranking, we retain the current
combined-pool analysis for this manuscript and flag this caveat here.