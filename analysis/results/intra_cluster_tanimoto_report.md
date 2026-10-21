# Intra-cluster Tanimoto — steroid molecule clusters

- Total molecules with valid SMILES: **676**
- Random-pair Tanimoto: mean **0.238**, median **0.205** (n = 19,956 sampled pairs)
- All intra-cluster pairs: mean **0.360**, median **0.333** (n = 45,650 pairs)

## Per-cluster

| Cluster | n | mean | median | std | min | max |
|---|---|---|---|---|---|---|
| C1 | 215 | 0.404 | 0.383 | 0.132 | 0.137 | 1.000 |
| C2 | 145 | 0.300 | 0.264 | 0.130 | 0.094 | 1.000 |
| C3 | 96 | 0.311 | 0.253 | 0.163 | 0.108 | 1.000 |
| C4 | 32 | 0.727 | 0.717 | 0.082 | 0.580 | 0.923 |
| C5 | 112 | 0.278 | 0.246 | 0.146 | 0.066 | 1.000 |
| C6 | 27 | 0.516 | 0.578 | 0.320 | 0.134 | 1.000 |
| C7 | 27 | 0.595 | 0.573 | 0.114 | 0.391 | 0.935 |
| C8 | 22 | 0.485 | 0.467 | 0.129 | 0.254 | 0.913 |

**Reading:** higher mean = the cluster is chemically tight (its members share substructure). Big gap between mean and min = the cluster contains some outliers RDKit sees as chemically distant despite sharing the UMAP neighborhood.

## Files
- `intra_cluster_tanimoto.tsv` — per-cluster summary
- `per_molecule_tanimoto.tsv` — per-molecule centrality inside its cluster
- `intra_cluster_tanimoto_distribution.{pdf,png}` — box plot
- `intra_vs_random_tanimoto.{pdf,png}` — KDE overlay vs null