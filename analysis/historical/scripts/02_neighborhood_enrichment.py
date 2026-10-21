"""Neighborhood enrichment on the 2D Steroid Atlas UMAP.

For each protein p with a populated annotation label L_p on feature F, we:
  1. Find its k=20 nearest UMAP neighbors (Euclidean distance in 2D).
  2. Count how many share L_p    → observed_k
  3. Compute expected fraction   → atlas-wide baseline frequency of L_p
  4. Enrichment ratio            → (observed / k) / expected
  5. Binomial z-score            → statistical significance

This measures whether the exploratory experience a user has (clicking a point,
seeing its neighbors) surfaces biologically similar proteins more often than
by chance.

Inputs:
  ../analysis/annotations_matrix.tsv     (feature matrix from step 01)

Outputs:
  ../analysis/enrichment_per_protein.tsv         one row × feature × protein
  ../analysis/enrichment_summary_per_feature.tsv one row per feature
  ../analysis/enrichment_summary.txt             human-readable report
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = output("annotations_matrix.tsv", historical=True)
OUT_PER_PROTEIN = output("enrichment_per_protein.tsv", historical=True)
OUT_PER_FEATURE = output("enrichment_summary_per_feature.tsv", historical=True)
OUT_REPORT = output("enrichment_summary.txt", historical=True)

K_NEIGHBORS = 20  # excludes self

# Features to enrich against. For each, we require a populated non-blank label.
FEATURES = [
    ("ec_top_level", "EC top-level (1..7)"),
    ("ec_subclass", "EC subclass (e.g., 1.14)"),
    ("ec_sub_subclass", "EC sub-subclass (e.g., 1.14.14)"),
    ("protein_name_first_word", "Protein-name first word"),
    ("cyp_family", "CYP family (e.g., CYP17)"),
    ("cyp_subfamily", "CYP subfamily (e.g., CYP17A)"),
    ("genus", "Taxonomic genus"),
    ("species_binomial", "Species (binomial)"),
]

# Boolean flags — enriched separately (fraction of neighbors also flagged True)
BOOL_FEATURES = [
    ("is_p450", "Cytochrome P450"),
    ("is_bsh", "Bile salt hydrolase"),
    ("is_hsd", "Hydroxysteroid dehydrogenase"),
    ("is_reductase", "Reductase/dehydrogenase"),
    ("is_receptor", "Receptor"),
    ("is_transporter", "Transporter/carrier"),
]


def main() -> int:
    print(f"Loading {IN}...")
    df = pd.read_csv(IN, sep="\t", low_memory=False)
    n = len(df)
    print(f"  {n:,} proteins")

    # ─── Nearest-neighbor lookup on 2D UMAP ─────────────────────────────
    coords = df[["umap_1", "umap_2"]].to_numpy()
    print(f"\nFitting KDTree on 2D UMAP coordinates...")
    nn = NearestNeighbors(n_neighbors=K_NEIGHBORS + 1, algorithm="kd_tree").fit(coords)
    print(f"Querying {K_NEIGHBORS+1} nearest neighbors for each protein...")
    _, idx = nn.kneighbors(coords)  # idx[i, 0] == i (self)
    neighbors = idx[:, 1:]  # drop self column -> shape (n, K_NEIGHBORS)
    print(f"  neighbors matrix: {neighbors.shape}")

    per_protein_rows: list[dict] = []
    feature_summary_rows: list[dict] = []

    # ─── Categorical features ───────────────────────────────────────────
    for feat, label in FEATURES:
        col = df[feat].astype(str).replace({"nan": ""})
        populated_mask = col.astype(str).str.strip().ne("")
        n_pop = int(populated_mask.sum())
        if n_pop == 0:
            print(f"[SKIP] {feat}: no populated rows")
            continue
        # Atlas baseline frequency of each label (among populated rows)
        pop_labels = col[populated_mask]
        baseline_freq = pop_labels.value_counts(normalize=True).to_dict()
        n_labels = len(baseline_freq)

        # For each protein with a populated label, count how many of its k neighbors share it
        enrichments = []
        z_scores = []
        observed_fractions = []
        expected_fractions = []
        for i in np.where(populated_mask)[0]:
            my_label = col.iat[i]
            neigh_labels = col.iloc[neighbors[i]].values  # k labels (may include blanks)
            obs = int((neigh_labels == my_label).sum())
            expected_freq = baseline_freq.get(my_label, 0.0)
            observed_freq = obs / K_NEIGHBORS
            enrichment = (observed_freq / expected_freq) if expected_freq > 0 else float("nan")
            # binomial z-score
            var = expected_freq * (1 - expected_freq) / K_NEIGHBORS
            z = (observed_freq - expected_freq) / math.sqrt(var) if var > 0 else 0.0
            enrichments.append(enrichment)
            z_scores.append(z)
            observed_fractions.append(observed_freq)
            expected_fractions.append(expected_freq)
            per_protein_rows.append({
                "accession": df.iat[i, df.columns.get_loc("accession")],
                "feature": feat,
                "my_label": my_label,
                "observed_in_neighbors": obs,
                "expected_fraction": expected_freq,
                "observed_fraction": observed_freq,
                "enrichment_ratio": enrichment,
                "z_score": z,
            })

        enrichments = np.array(enrichments, dtype=float)
        z_scores = np.array(z_scores, dtype=float)
        obs_fr = np.array(observed_fractions)
        exp_fr = np.array(expected_fractions)

        # Filter valid enrichments (drop NaN / inf from labels with 0 baseline)
        valid = np.isfinite(enrichments) & (enrichments > 0)
        e_valid = enrichments[valid]
        z_valid = z_scores[valid]

        row = {
            "feature": feat,
            "label": label,
            "n_populated": n_pop,
            "n_unique_labels": n_labels,
            "median_enrichment": float(np.median(e_valid)) if len(e_valid) else float("nan"),
            "mean_enrichment": float(np.mean(e_valid)) if len(e_valid) else float("nan"),
            "p90_enrichment": float(np.percentile(e_valid, 90)) if len(e_valid) else float("nan"),
            "median_z": float(np.median(z_valid)) if len(z_valid) else float("nan"),
            "frac_enrichment_gt_2x": float((e_valid > 2).mean()) if len(e_valid) else float("nan"),
            "frac_enrichment_gt_5x": float((e_valid > 5).mean()) if len(e_valid) else float("nan"),
            "frac_z_gt_3": float((z_valid > 3).mean()) if len(z_valid) else float("nan"),
            "median_observed_fraction": float(np.median(obs_fr)) if len(obs_fr) else float("nan"),
            "median_expected_fraction": float(np.median(exp_fr)) if len(exp_fr) else float("nan"),
        }
        feature_summary_rows.append(row)
        print(f"[{feat:<24}] n={n_pop:>5}  labels={n_labels:>4}  "
              f"median enrichment={row['median_enrichment']:>7.2f}×  "
              f"median z={row['median_z']:>6.1f}  "
              f"frac(z>3)={row['frac_z_gt_3']*100:>4.1f}%")

    # ─── Boolean features ───────────────────────────────────────────────
    # For a boolean feature, "sharing" means both my flag and neighbor's flag are True.
    # We restrict to proteins where the flag is True.
    for feat, label in BOOL_FEATURES:
        col = df[feat].astype(bool).values
        pos_idx = np.where(col)[0]
        n_pop = int(len(pos_idx))
        if n_pop == 0:
            print(f"[SKIP] {feat}: no True rows")
            continue
        baseline_freq = col.mean()
        enrichments = []
        z_scores = []
        observed_fractions = []
        for i in pos_idx:
            neigh_flags = col[neighbors[i]]
            obs = int(neigh_flags.sum())
            observed_freq = obs / K_NEIGHBORS
            enrichment = (observed_freq / baseline_freq) if baseline_freq > 0 else float("nan")
            var = baseline_freq * (1 - baseline_freq) / K_NEIGHBORS
            z = (observed_freq - baseline_freq) / math.sqrt(var) if var > 0 else 0.0
            enrichments.append(enrichment)
            z_scores.append(z)
            observed_fractions.append(observed_freq)
            per_protein_rows.append({
                "accession": df.iat[i, df.columns.get_loc("accession")],
                "feature": feat,
                "my_label": "True",
                "observed_in_neighbors": obs,
                "expected_fraction": baseline_freq,
                "observed_fraction": observed_freq,
                "enrichment_ratio": enrichment,
                "z_score": z,
            })

        e_arr = np.array(enrichments, dtype=float)
        z_arr = np.array(z_scores, dtype=float)
        obs_fr = np.array(observed_fractions)
        valid = np.isfinite(e_arr) & (e_arr > 0)
        e_valid = e_arr[valid]
        z_valid = z_arr[valid]

        row = {
            "feature": feat,
            "label": label,
            "n_populated": n_pop,
            "n_unique_labels": 2,
            "median_enrichment": float(np.median(e_valid)) if len(e_valid) else float("nan"),
            "mean_enrichment": float(np.mean(e_valid)) if len(e_valid) else float("nan"),
            "p90_enrichment": float(np.percentile(e_valid, 90)) if len(e_valid) else float("nan"),
            "median_z": float(np.median(z_valid)) if len(z_valid) else float("nan"),
            "frac_enrichment_gt_2x": float((e_valid > 2).mean()) if len(e_valid) else float("nan"),
            "frac_enrichment_gt_5x": float((e_valid > 5).mean()) if len(e_valid) else float("nan"),
            "frac_z_gt_3": float((z_valid > 3).mean()) if len(z_valid) else float("nan"),
            "median_observed_fraction": float(np.median(obs_fr)) if len(obs_fr) else float("nan"),
            "median_expected_fraction": float(baseline_freq),
        }
        feature_summary_rows.append(row)
        print(f"[{feat:<24}] n={n_pop:>5}  base_freq={baseline_freq*100:>5.2f}%  "
              f"median enrichment={row['median_enrichment']:>7.2f}×  "
              f"median z={row['median_z']:>6.1f}  "
              f"frac(z>3)={row['frac_z_gt_3']*100:>4.1f}%")

    # ─── Write outputs ──────────────────────────────────────────────────
    per_protein_df = pd.DataFrame(per_protein_rows)
    per_protein_df.to_csv(OUT_PER_PROTEIN, sep="\t", index=False)
    print(f"\nWrote per-protein enrichment: {OUT_PER_PROTEIN.name} ({len(per_protein_df):,} rows)")

    feat_df = pd.DataFrame(feature_summary_rows)
    feat_df.to_csv(OUT_PER_FEATURE, sep="\t", index=False)
    print(f"Wrote per-feature summary:    {OUT_PER_FEATURE.name} ({len(feat_df)} rows)")

    # ─── Human-readable report ─────────────────────────────────────────
    lines = []
    lines.append("=== Neighborhood enrichment on 2D Steroid Atlas UMAP ===\n")
    lines.append(f"k = {K_NEIGHBORS} nearest neighbors (excluding self)")
    lines.append(f"n proteins = {n:,}")
    lines.append("\nInterpretation:")
    lines.append("  Enrichment ratio = (fraction of k neighbors sharing my label) / (atlas baseline frequency of my label)")
    lines.append("  z-score          = binomial standardized deviation of observed sharing from expected")
    lines.append("  Higher enrichment = neighbors are more biologically similar than random")
    lines.append("\n")
    lines.append(f"{'Feature':<24}{'n':>7}{'labels':>8}{'med enr':>10}{'p90 enr':>10}{'med z':>8}"
                 f"{'%>2×':>7}{'%>5×':>7}{'%z>3':>7}")
    lines.append("-" * 88)
    for r in feature_summary_rows:
        lines.append(f"{r['feature']:<24}{r['n_populated']:>7}{r['n_unique_labels']:>8}"
                     f"{r['median_enrichment']:>10.2f}{r['p90_enrichment']:>10.2f}"
                     f"{r['median_z']:>8.2f}"
                     f"{r['frac_enrichment_gt_2x']*100:>6.1f}%"
                     f"{r['frac_enrichment_gt_5x']*100:>6.1f}%"
                     f"{r['frac_z_gt_3']*100:>6.1f}%")

    OUT_REPORT.write_text("\n".join(lines) + "\n")
    print(f"Wrote report:                 {OUT_REPORT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
