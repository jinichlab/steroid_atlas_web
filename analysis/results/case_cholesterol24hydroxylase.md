# Case study — Cholesterol 24-hydroxylase (5 clusters)

## Motivation

Five separate UMAP clusters (C6, C18, C26, C48, C76 — total 1,353 proteins) all get labelled `Cholesterol 24-hydroxylase` by the dominant-name fingerprint. Why does k-means put them in five distinct clusters rather than one?

**Short answer:** ProtT5 embeddings encode the amino-acid sequence, so proteins that share a *recommended enzyme name* but belong to different sequence families / organisms / lengths land in different embedding-space clusters. The dominant-name label collapses that diversity into a single string, but the underlying biology is 5 distinct enzyme populations.

## Per-cluster composition

| Cluster | n proteins | median length (aa) | top kingdom | top organism |
|---|---|---|---|---|
| C6 | 713 | 501 (463–510) | Other/unknown (89%) | Oncorhynchus mykiss (Rainbow trout) (Salmo gairdneri) (33) |
| C18 | 182 | 500 (221–587) | Other/unknown (84%) | Macaca mulatta (Rhesus macaque) (4) |
| C26 | 200 | 394 (373–401) | Other/unknown (94%) | Thamnophis sirtalis (2) |
| C48 | 199 | 394 (394–531) | Other/unknown (81%) | Rattus norvegicus (Rat) (10) |
| C76 | 59 | 445 (237–542) | Other/unknown (59%) | Homo sapiens (Human) (8) |

## UMAP centroid distances (ProtT5 embedding proxy)

```
        C6    C18    C26    C48    C76
C6    0.00   4.31   5.91   7.68  18.94
C18   4.31   0.00   1.97   3.38  17.15
C26   5.91   1.97   0.00   2.13  18.05
C48   7.68   3.38   2.13   0.00  16.79
C76  18.94  17.15  18.05  16.79   0.00
```

Higher distance = more divergent embedding = more sequence-level difference. Clusters that share a name but sit far apart on the UMAP are the strongest evidence that the name is convergent, not reflecting a single sequence family.

## Substrate ChEBI overlap between clusters (Jaccard)

```
       C6   C18   C26   C48   C76
C6   1.00  0.43  0.77  0.33  0.30
C18  0.43  1.00  0.45  0.33  0.33
C26  0.77  0.45  1.00  0.30  0.31
C48  0.33  0.33  0.30  1.00  0.37
C76  0.30  0.33  0.31  0.37  1.00
```

## Mean pairwise substrate Tanimoto (within & cross)

```
        C6    C18    C26    C48    C76
C6   0.378  0.317  0.423  0.335  0.317
C18  0.317  0.272  0.334  0.279  0.271
C26  0.423  0.334  0.434  0.348  0.329
C48  0.335  0.279  0.348  0.301  0.288
C76  0.317  0.271  0.329  0.288  0.279
```

Diagonal = within-cluster mean pairwise Tanimoto. Off-diagonal = mean pairwise Tanimoto between substrates of cluster i and substrates of cluster j. If the diagonal >> off-diagonal for a pair, that pair of clusters really does act on chemically distinct substrate sets.

## Files
- `case_cholesterol24hydroxylase_composition.tsv` — the 5-row summary table with lengths, organisms, kingdoms, UMAP centroids.
- `case_cholesterol24hydroxylase.{pdf,png}` — 4-panel comparison figure (length distribution, kingdom composition, substrate Jaccard, substrate cross-Tanimoto).