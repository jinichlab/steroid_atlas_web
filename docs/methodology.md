# Methodology

## Corpus construction

1. **Seed set**: Query Rhea + UniProt for reactions involving steroid ChEBI IDs (the "Rhea steroid-reaction query"). Extract all UniProt accessions that catalyze at least one such reaction.
2. **Dedup**: Remove exact-sequence duplicates. Original set was 36,877 entries; dedup pass removed 1,528 duplicates → **35,349 unique protein sequences**.
3. **Small-molecule catalog**: For every ChEBI substrate/product from the Rhea query, keep the compound + SMILES + canonical name → 585 Rhea-derived molecules; adding 91 literature-recruited compounds gives **676 small molecules**.

## Embedding + UMAP

- **ProtT5 embeddings**: Compute per-residue embeddings via `Rostlab/prot_t5_xl_half_uniref50-enc`, mean-pool over the sequence → 1024-d per-protein vector.
- **UMAP**: Fit 2D UMAP (n_neighbors=15, min_dist=0.1) on the 36,877 × 1024 matrix. **The UMAP reducer was not pickled**, so new entries can't be `.transform()`ed; instead, they are placed via **cosine-nearest-neighbor** in ProtT5 space, inheriting the neighbor's UMAP coordinates plus a small Gaussian jitter (σ=0.15). This is a faithful substitute because UMAP preserves the local cosine k-NN structure it was built from.
- **HDBSCAN**: Cluster the 2D UMAP with default parameters → cluster IDs stored in the `cluster` column.

## Literature-recruited entries

The atlas grows by explicit paper audits. For each candidate paper:

1. **Ingest paper text** (PDF → RAG catalog for keyword-searchable chunks)
2. **Identify candidate proteins** (locus tags, UniProt accessions, GenBank IDs) mentioned in the paper's methods or supplementary tables
3. **Fetch sequences** from public databases (UniProt REST or NCBI E-utilities), verify by genus + length + keyword match against the paper's stated organism
4. **Attach evidence tier** based on what the paper *actually* demonstrates for each protein:
   - **Tier A**: Biochemistry + kinetics + mutagenesis on the purified enzyme (strongest)
   - **Tier B-strong**: Sole enzyme of that family in the species; species-culture activity unambiguously attributable
   - **Tier B-moderate/weak**: Multiple paralogs; canonical clade paralog inferred as likely carrier
   - **Tier C**: Group II/III paralog; activity not disambiguated
5. **Append to atlas** with full annotation (paper DOI, citation, evidence tier, sequence source)
6. **Embed + project** into UMAP via nearest-neighbor
7. **Update provenance table** `data/literature_recruited_proteins.csv`

## Curation policies

Enforced throughout:

- **No fabricated identifiers.** If a compound has no ChEBI entry, the `chebi_id` is blank. If a protein lacks a UniProt accession, its `accession` is the locus tag and `identifier_type=Locus_Tag`.
- **EC numbers reflect paper evidence.** For literature-recruited entries, `ec_numbers` lists only ECs directly demonstrated by the cited paper for that specific protein. Sequence-similarity-inferred ECs from UniProt are excluded.
- **Duplicate ChEBI IDs verified.** All 594 non-blank ChEBI IDs verified via the OLS/ChEBI ontology API before publication. See `data/provenance/chebi_audit_report.tsv`.
- **Dismissed entries documented.** Enzymes fetched but not added (e.g., 18 of 20 Guzior 2024 Supp Table 1 BSH/Ts without protein-level biochemical evidence) are preserved in `data/provenance/guzior2024_dismissed_entries.json` with the dismissal reason.
