# Cluster substrate-Tanimoto significance

## Method

For each cluster we take the observed mean pairwise Tanimoto over its unique substrate set, then draw **2,000 random substrate samples of the same size** from the full molecule pool (n = 591) and compute the same mean each time to build a null distribution. The empirical p-value is `(#null ≥ observed + 1) / (N_perm + 1)` (specialist direction) or `(#null ≤ observed + 1) / (N_perm + 1)` (promiscuous direction). We correct each direction separately by Benjamini-Hochberg FDR.

**Z-score sign convention:** positive = observed > null (specialist), negative = observed < null (promiscuous). Larger magnitudes mean the cluster is farther from what random sampling would produce.

## Summary

- Testable clusters (≥2 substrates): **77 of 82**
- Significantly SPECIALIST (q_specialist < 0.05, z > 0): **71**
- Significantly PROMISCUOUS (q_promiscuous < 0.05, z < 0): **0**

---

## Top 10 significantly SPECIALIST protein families
*(substrates are chemically MORE similar than random pools of the same size)*

| Protein family | # substrates | Observed | Null mean | Z | q |
|---|---|---|---|---|---|
| **Sterol O-acyltransferase 1** | 61 | 0.401 | 0.231 | +18.96 | 7.9e-04 |
| **Acyl-coenzyme A diphosphatase NUDT19** | 17 | 0.574 | 0.232 | +17.87 | 7.9e-04 |
| **Cytochrome P450 85A1** | 29 | 0.453 | 0.231 | +15.18 | 7.9e-04 |
| **Bile salt export pump** | 10 | 0.581 | 0.231 | +13.12 | 7.9e-04 |
| **Cholesterol 24-hydroxylase** | 20 | 0.434 | 0.231 | +11.78 | 7.9e-04 |
| **Sterol carrier protein 2** | 50 | 0.348 | 0.231 | +11.25 | 7.9e-04 |
| **Phosphatidylcholine-sterol acyltransferase** | 38 | 0.367 | 0.231 | +11.19 | 7.9e-04 |
| **Delta(24)-sterol reductase** | 58 | 0.327 | 0.231 | +10.39 | 7.9e-04 |
| **Cholesterol 24-hydroxylase** | 26 | 0.378 | 0.231 | +9.96 | 7.9e-04 |
| **G-protein coupled estrogen receptor 1** | 10 | 0.492 | 0.230 | +9.81 | 7.9e-04 |

## Top 10 significantly PROMISCUOUS protein families
*(substrates are chemically MORE diverse than random pools of the same size)*

| Protein family | # substrates | Observed | Null mean | Z | q |
|---|---|---|---|---|---|

---

## Files
- `cluster_tanimoto_significance.tsv` — full 82-row table with z-scores + q-values in both directions.
- `cluster_tanimoto_significance_volcano.{pdf,png}` — Z vs −log₁₀(q).