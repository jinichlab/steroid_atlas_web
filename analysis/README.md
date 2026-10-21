# analysis/

Every analysis behind the Steroid Atlas paper, plus its outputs.

```
analysis/
├── _paths.py        shared paths: scripts read data/, write results/ + figures/
├── scripts/         analyses that re-run from data/  (section 1)
├── results/         their tables (.tsv) and reports (.md / .txt)
├── figures/         their figures (.png + .pdf)
└── historical/      how data/proteins.csv was built — kept as a record,
    ├── scripts/     not re-runnable from this repo   (section 2)
    ├── results/
    └── figures/
```

**Script numbers are chronological (the order they were written), not a run
order** — several numbers are shared (two `09_`, two `12_`, …). The run order
is `scripts/reproduce_analysis.sh`.

## Reproduce

```bash
pip install -r requirements.txt
sh scripts/reproduce_analysis.sh            # all downstream analyses + web data
git diff --stat analysis/ web/public/atlas/     # what changed vs. the committed outputs
```

No network access or GPU needed. Randomized steps use fixed seeds
(`RANDOM_STATE = 42`, `RNG_SEED = 20260824`), so re-runs give the same numbers.

## 1. `scripts/` — reproducible from `data/` (run by `reproduce_analysis.sh`)

| Script | What it does | Main outputs |
|---|---|---|
| `18_sync_steroid_sets.py` | Dedupes `molecules.csv` (677 → 676), rebuilds the natural set of `natural_synthetic_steroids.csv` (2,976 rows), refits that UMAP (Morgan r=2/2048, Jaccard, n_neighbors=15, min_dist=0.1, seed 42). Idempotent. | `natsyn_umap_report.txt`, `natsyn_umap.joblib` (gitignored) |
| `03_kmeans_kselection_extended.py` | k-means sweep k = 10–200 on the 2D protein UMAP; silhouette + Davies-Bouldin. k = 82 was chosen from the [50, 95] window. | `kselection_extended.tsv`, `_report.txt` |
| `08_silhouette_plot.py` | Silhouette / D-B vs k figure. | `silhouette_vs_k.*` |
| `07_cluster_fingerprints.py` | Dominant identity per protein cluster (the cluster names in the web app). | `cluster_fingerprints.tsv`, `cluster_shortlist.txt` |
| `19_highdim_kmeans.py` | k-means in the full embedding with the same k (steroids: Morgan fingerprints; proteins: ProtT5 from `data/embeddings/prott5.h5`, if present), plus agreement (ARI) with the 2D clusters. Feeds the web app's clustering toggle; `07_cluster_fingerprints.py --labels cluster_hd` names the protein clusters. | `molecule_clusters_hd.tsv`, `protein_clusters_hd.tsv`, `highdim_kmeans_report.md` |
| `09_cluster_go_enrichment.py` | Per-cluster GO enrichment. | `cluster_go_enrichment.tsv`, `cluster_go_top_terms.tsv`, `cluster_go_shortlist.txt` |
| `10_cluster_enrichment_figure.py` | Atlas-wide GO enrichment figure. | `cluster_enrichment_figure.*` |
| `12_organism_breakdown.py` | Organism composition from the organism strings. | `organism_figure.*`, `organism_summary.txt` |
| `09_intra_cluster_tanimoto.py` | Intra- vs random-pair Tanimoto for the steroid (molecule) clusters. | `intra_cluster_tanimoto*.{tsv,md,png,pdf}`, `per_molecule_tanimoto.tsv` |
| `10_protein_cluster_substrate_tanimoto.py` | Substrate chemical diversity per protein cluster. | `protein_cluster_substrate_*` |
| `11_cluster_tanimoto_significance.py` | Permutation test (2,000 perms) of within-cluster substrate Tanimoto. | `cluster_tanimoto_significance.tsv`, `_report.md`, `_volcano.*` |
| `12_within_cluster_tanimoto_figure.py` | Within-cluster pairwise Tanimoto figure. | `within_cluster_tanimoto_{all,extremes,enzymes_only}.*` |
| `13_within_cluster_tanimoto_by_type.py` | Same, split by protein class (enzyme / transporter / …). | `cluster_type_composition.tsv`, `within_cluster_tanimoto_by_type*` |
| `14_case_study_cholesterol24hydroxylase.py` | Why 5 clusters are all named "cholesterol 24-hydroxylase". | `case_cholesterol24hydroxylase.*` |
| `15_cyp_superfamily_landscape.py` | CYP superfamily landscape across clusters. | `cyp_*` |
| `16_figureX_coherence_and_type.py`, `16b_figureX_panels_separate.py` | Figure X (cluster coherence + class stratification), combined and as separate panels. | `figure_X_*` |
| `17_type_stratification_by_ec_class.py` | Panel C refined with enzymes split by EC class. | `cluster_type_composition_ec.tsv`, `figure_X_C_type_stratification_ec.*` |

## 2. `historical/` — curation already applied to `data/proteins.csv`

These built the published `proteins.csv` from the earlier ~35k-entry corpus.
They are kept as the record of how the data was made, but **cannot be re-run
from this repo**: they need that pre-curation `proteins.csv` (not included),
call the live UniProt / PubMed APIs (results drift as UniProt is updated), and
several rewrite `data/proteins.csv` in place.

| Order | Script | Step |
|---|---|---|
| 1 | `04_validate_binder_annotations.py` → `05_apply_curation.py` | Binder audit: drop 231 vitamin-D-only / unannotated binder entries (`binder_audit.tsv`, `curation_diff_report.txt`). |
| 2 | `06_enrich_go_and_keywords.py` | Fetch GO + keyword annotations for every entry (network). |
| 3 | `13_rhea_reaudit.py` → `14_apply_rhea_curation.py` | Re-check Rhea links against current UniProt; drop stale ones (`rhea_audit.tsv`). |
| 4 | `15_annotation_score_audit.py` | Fetch annotation score + entry status (`score_audit.tsv`). |
| 5 | `16_fetch_pubmed_and_build_master.py` (driven by `run_final_master.sh`) | Keep score 5 + active + steroid-associated → 14,075; PubMed IDs; final master table (`master_table_report.txt`). The 14 literature-recruited entries were then restored → 14,089. |

Also historical (35k-era corpus, outputs kept for reference):
`01_feature_extraction.py` / `02_neighborhood_enrichment.py` (their
`annotations_matrix.tsv` input was pruned), `11_sequence_similarity.py`
(needs Biopython + the CD-HIT binary), `12_taxonomy_breakdown.py` (network).

## 3. Upstream — not in any repo yet

Corpus construction, embedding and clustering ran on the lab server
(`/home/adsiordia/steroid_core_classifier/`, see `docs/methods_corpus_construction.md`).
To make the atlas reproducible end-to-end these still need to be added:

- **Corpus build:** `20260511_comprehensive_steroid_uniprot.py` (sterane filter → Rhea → UniProt, incl. `modify_molecule`), `20260511_add_binders_and_old_only.py`, `20260511_make_fasta.py`, and the older `20251120_*` seed scripts.
- **Embedding:** `script_prot5.sh` + `prott5_embedder.py`, and the resulting ProtT5 matrix (`all_steroid_uniprot_comprehensive_v2.h5`, ~37k × 1,024 — too large for git; a Hugging Face dataset fits).
- **Protein UMAP + clusters:** `refit_full_umap.py` (UMAP parameters incl. metric and seed), `demo/_dedupe_sequences.py`, and the script that assigned the final k = 82 labels.
- **Steroid (molecule) UMAP + k = 8 clusters:** the fingerprint settings and UMAP parameters behind `molecules.csv` `umap_1/umap_2/cluster` are not recorded (k-means k = 8, `random_state=42`, `n_init=10` on the 2D UMAP reproduces the labels exactly).
- **Pre-curation `proteins.csv`** (~35k rows), to re-run section 2.
