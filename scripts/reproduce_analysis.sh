#!/bin/sh
# Re-run every downstream analysis from the published data/ CSVs, then rebuild
# the web app's data. No network access, no GPU. See analysis/README.md for
# what each step does and for the upstream steps this does NOT cover.
#
#   pip install -r requirements.txt
#   sh scripts/reproduce_analysis.sh            # everything (~10–20 min)
#   sh scripts/reproduce_analysis.sh --no-web   # analyses only
#
# Outputs are written in place (analysis/results/, analysis/figures/, web/public/atlas/),
# so `git diff` afterwards shows exactly what changed.
set -eu
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
export MPLBACKEND=Agg  # render figures without a display
export PYTHONHASHSEED=0  # fixed set/dict iteration order across runs
export SOURCE_DATE_EPOCH=0  # fixed PDF timestamps, so figures are byte-identical

run() {
  echo "── $1"
  start=$(date +%s)
  "$PY" "$@"
  echo "   done in $(( $(date +%s) - start ))s"
}

# 1. Data sync: dedupe molecules.csv (676) + rebuild the natural set and refit
#    the natural+synthetic UMAP. Idempotent on the published CSVs.
run analysis/scripts/18_sync_steroid_sets.py

# 2. Protein clustering: k-selection sweep + silhouette figure (k = 82 chosen
#    from this sweep; the k = 82 labels themselves are in data/proteins.csv).
run analysis/scripts/03_kmeans_kselection_extended.py
run analysis/scripts/08_silhouette_plot.py

# 3. Cluster identity: fingerprints (names used by the web app) + GO enrichment.
run analysis/scripts/07_cluster_fingerprints.py
# High-dimensional k-means (same k) for the web app's clustering toggle.
# Proteins need data/embeddings/prott5.h5; without it only steroids are done.
run analysis/scripts/19_highdim_kmeans.py
if [ -f analysis/results/protein_clusters_hd.tsv ]; then
  run analysis/scripts/07_cluster_fingerprints.py --labels cluster_hd
fi
run analysis/scripts/09_cluster_go_enrichment.py
run analysis/scripts/10_cluster_enrichment_figure.py
run analysis/scripts/12_organism_breakdown.py

# 4. Substrate chemistry (Tanimoto) per steroid cluster and per protein cluster.
run analysis/scripts/09_intra_cluster_tanimoto.py
run analysis/scripts/10_protein_cluster_substrate_tanimoto.py
run analysis/scripts/11_cluster_tanimoto_significance.py
run analysis/scripts/12_within_cluster_tanimoto_figure.py
run analysis/scripts/13_within_cluster_tanimoto_by_type.py

# 5. Case studies + paper figures.
run analysis/scripts/14_case_study_cholesterol24hydroxylase.py
run analysis/scripts/15_cyp_superfamily_landscape.py
run analysis/scripts/16_figureX_coherence_and_type.py
run analysis/scripts/16b_figureX_panels_separate.py
run analysis/scripts/17_type_stratification_by_ec_class.py

# 6. Web app data (web/public/atlas/*.json + structure PNGs).
if [ "${1:-}" != "--no-web" ]; then
  run scripts/build_atlas_data.py
  run scripts/build_structures.py
fi

echo "All steps finished."
