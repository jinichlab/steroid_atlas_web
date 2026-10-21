# Data files

Three CSVs form the atlas. Every row is human-readable and traceable back to a paper or public database.

## `proteins.csv` — 35,349 rows × 21 columns

One row per protein sequence. Deduplicated (no identical sequences repeated).

| Column | Description |
|---|---|
| `accession` | Primary identifier (UniProt accession, RefSeq WP_ ID, or locus tag) |
| `entry_name` | UniProt-style entry name (e.g. `CBH_BIFL2`) |
| `protein_names` | Full protein name(s) with EC numbers as UniProt records them |
| `gene_names` | Gene name + locus tag(s) |
| `organism` | Species + strain designation |
| `length_aa` | Sequence length (amino acids) |
| `sequence` | Full protein sequence |
| `ec_numbers` | EC numbers **only** for reactions experimentally confirmed by the cited paper for this specific protein. Blank if no per-protein EC evidence. |
| `rhea_reactions` | Rhea reaction database IDs |
| `interacting_chebi_ids` | `;`-separated ChEBI IDs of steroid substrates/products |
| `interacting_compounds` | `;`-separated human-readable names matching the ChEBI IDs |
| `reaction_descriptions` | Free-text description of the reactions this protein catalyzes |
| `umap_1`, `umap_2` | 2D UMAP embedding of the ProtT5 sequence embedding |
| `cluster` | Cluster ID from HDBSCAN clustering of the UMAP |
| `prott5_cluster` | Original ProtT5-space cluster ID |
| `is_literature_recruited` | 1 if this row was added from a specific paper's audit; 0 if from the original Rhea/UniProt corpus |
| `paper_url` | DOI URL of the paper this entry was recruited from (only if `is_literature_recruited=1`) |
| `annotation` | Human-readable annotation citing the paper |
| `sequence_source` | Where the sequence was fetched from + evidence trail |
| `identifier_type` | UniProt / RefSeq_WP / Locus_Tag / NCBI_Genome |

## `molecules.csv` — 676 rows × 10 columns

One row per steroid/bile-acid small molecule.

| Column | Description |
|---|---|
| `compound_name` | Human-readable name (common or IUPAC) |
| `chebi_id` | ChEBI database ID (numeric, no `CHEBI:` prefix). Blank if no ChEBI entry known — never fabricated. |
| `smiles` | SMILES string (canonical or as-published) |
| `umap_1`, `umap_2` | 2D UMAP embedding of the RDKit Morgan fingerprint |
| `cluster` | Cluster ID from HDBSCAN clustering of the UMAP |
| `is_literature_recruited` | 1 if added from a specific paper; 0 if from original corpus |
| `paper_url` | Paper DOI (only if `is_literature_recruited=1`) |
| `interacting_protein_accessions` | List of protein accessions known to act on this compound |
| `chebi_url` | Direct link to the ChEBI page (blank if no ChEBI ID) |

## `literature_recruited_proteins.csv` — 15 rows × 12 columns

Provenance table: every row is a protein that was added to `proteins.csv` from a specific paper (i.e., `is_literature_recruited=1`).

| Column | Description |
|---|---|
| `Entry` | Same accession as in `proteins.csv` |
| `Paper_ShortRef` | Short citation (e.g. "Rimal 2024") |
| `Paper_Citation` | Full citation string |
| `Paper_DOI` | DOI |
| `Species` | Organism the enzyme is from |
| `Strain` | Strain designation |
| `Identifier_Type` | UniProt / RefSeq_WP / Locus_Tag / NCBI_Genome |
| `Sequence_Length` | aa |
| `Recruitment_Rationale` | Human-readable description of why this specific enzyme was added — what activity did the paper report |
| `Evidence_Level` | Biochemical / genetic / bioinformatic — nature of the evidence in the cited paper |
| `Sequence_Fetch_Source` | Which database the sequence came from + any special notes |
| `In_Atlas_UMAP` | "yes" if the entry has UMAP coordinates in `proteins.csv` |

A human-readable Markdown version is also provided at `literature_recruited_proteins.md`.

## Curation principles

1. **No fabricated identifiers.** If a compound has no known ChEBI ID, its `chebi_id` field is blank. Same for proteins that lack a UniProt accession — they use their locus tag instead of a made-up UniProt-like ID, and the `identifier_type` field tells you which.
2. **EC numbers reflect what the paper actually tested.** `ec_numbers` in `proteins.csv` lists only ECs directly demonstrated by the cited paper for that specific protein. Sequence-similarity inferences are excluded.
3. **Every literature-recruited row is auditable.** The `annotation` and `sequence_source` fields, plus the provenance table, form a complete chain from paper → gene → sequence → row.

## Optional: RAG AI Agent
Include your OpenAI API in the landing page of the atlas tool to use an AI agent that can answer questions about the data. The agent is trained on all available data on the small molecules visualized. It can provide answers for proteins based on general information available. 
