export interface Protein {
  accession: string;
  entry_name: string;
  protein_names: string;
  gene_names: string;
  organism: string;
  length_aa: number;
  ec_numbers: string;
  rhea_reactions: string;
  interacting_chebi_ids: string;
  interacting_compounds: string;
  go_ids: string;
  go_labels: string;
  keyword_ids: string;
  keyword_labels: string;
  pubmed_ids: string;
  pubmed_count: number;
  umap_1: number;
  umap_2: number;
  cluster: number;
  is_literature_recruited: number;
  paper_url: string;
  audit_decision: string;
}

export interface Molecule {
  compound_name: string;
  chebi_id: string;
  smiles: string;
  umap_1: number;
  umap_2: number;
  cluster: number;
  is_literature_recruited: number;
  paper_url: string;
  interacting_protein_accessions: string;
}

export interface NatsynEntry {
  compound_name: string;
  smiles: string;
  chebi_id: string;
  umap_1: number;
  umap_2: number;
  cluster: number | string;
  protein_entries: string;
}

export interface ClusterMeta {
  id: number;
  display_id: number;
  n: number;
  name: string;
  top_go: string;
}

export interface Summary {
  n_proteins: number;
  n_molecules: number;
  n_natsyn: number;
  n_clusters: number;
  n_new: number;
}

export type ViewKind = "protein" | "molecule" | "natsyn";
