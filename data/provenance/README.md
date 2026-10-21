# Literature audit pipeline

Everything under this directory is provenance for the 15 literature-recruited protein entries in `../proteins.csv` (rows where `is_literature_recruited=1`). Every script is runnable and every fetched artifact is preserved.

## Layout

```
data/provenance/
├── scripts/                              # runnable pipeline
├── supplementary/                        # SI files downloaded from papers
├── sequences/                            # fetched FASTAs organized per paper
├── embeddings/                           # ProtT5 embedding + projection reports
├── guzior2024_*.json / *.tsv             # per-paper audit artifacts
└── chebi_audit_report.tsv                # verified ChEBI IDs across the atlas
```

## The pipeline, in order

Each script is idempotent (safe to re-run) and writes backups before modifying any CSV.

```bash
# 1. Parse a paper's supplementary table (example: Guzior 2024)
python scripts/parse_guzior_table1.py

# 2. Fetch protein sequences from public databases
python scripts/fetch_guzior_sequences.py
python scripts/fetch_protein_sequences.py           # generic Tier-1..3 fetcher

# 3. Attach evidence tiers explaining WHY each entry qualifies
python scripts/annotate_guzior_evidence_tier.py

# 4. Append to the atlas
python scripts/append_cpbsht_only.py                # example: add just the
                                                    # biochemically-tested entry
python scripts/upgrade_rimal_bsh_entries.py         # example: fix annotations

# 5. Embed and project into UMAP (nearest-neighbor in ProtT5 space)
python scripts/embed_cpbsht_and_project.py

# 6. Curation policies
python scripts/strict_reaction_ecs_update.py        # only keep experimentally-
                                                    # confirmed ECs
python scripts/fix_chebi_ids.py                     # correct wrong ChEBI IDs
python scripts/fix_chebi_in_protein_csv.py          # ...in protein interacting list
python scripts/audit_all_chebi_ids.py               # verify all ChEBI IDs
                                                    # against ChEBI ontology
```

## Provenance JSON files

- `guzior2024_table1_normalized.json` — Guzior Supp Table 1 parsed to a flat list of 38 BSH/T records with type classification (RefSeq / GenBank / UniProt / JGI IMG numeric)
- `guzior2024_review_ready.json` — 19 fetched sequences with evidence-tier annotations, sequences, and proposed atlas rows (audit-ready)
- `guzior2024_dismissed_entries.json` — the 18 Guzior entries deliberately NOT added (no protein-level biochemical evidence), with reason documented
- `guzior2024_fetch_report.tsv` — per-sequence NCBI/UniProt fetch verification (genus match, length, keyword match)
- `chebi_audit_report.tsv` — every ChEBI ID in the atlas verified against the ChEBI ontology (594 MATCH, 1 stopword false-positive, 0 real mismatches)
- `embeddings/cpbsht_projection_report.json` — full nearest-neighbor projection provenance for the CpBSH/T UMAP position

## Embeddings

The h5 embedding files (`*.h5`) are **not tracked in git** (see `.gitignore`) — they're regenerable from the FASTAs via `scripts/embed_cpbsht_and_project.py` and the ProtT5 model from HuggingFace (`Rostlab/prot_t5_xl_half_uniref50-enc`, ~1.6 GB download).

If you want the pre-computed atlas embedding matrix (37,391 proteins × 1024 dims, ~150 MB), reach out — we can host it on Zenodo.
