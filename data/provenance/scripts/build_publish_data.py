"""Build the publish-ready CSVs in data/ from the audited source CSVs.

Reads:
  ../protein_sequence_embedding.DEDUP.csv   (audited protein atlas, 35,834 rows)
  ../small_molecule_centric.csv             (audited molecule atlas, 681 rows)
  ../literature_recruited_proteins.csv      (provenance table, 15 entries)

Writes into steroid-atlas/data/:
  proteins.csv                  — clean, publication-ready protein table
  molecules.csv                 — clean, publication-ready molecule table
  literature_recruited_proteins.csv  — provenance (unchanged)

Column renames follow the "simple names, snake_case" convention while
preserving all information a downstream user needs.
"""
from __future__ import annotations
from pathlib import Path

import pandas as pd

SRC = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer")
DST = Path("/home/adsiordia/marimo_visualizer/steroid-atlas/data")

# ─── PROTEINS ──────────────────────────────────────────────────────────────
print("Building proteins.csv ...")
prot = pd.read_csv(SRC / "protein_sequence_embedding.DEDUP.csv", low_memory=False,
                   dtype={"ChEBI ID": str, "Rhea ID": str, "SMILES": str})

prot_out = pd.DataFrame({
    "accession":                prot["Entry"],
    "entry_name":               prot.get("Entry Name", ""),
    "protein_names":            prot.get("Protein names", ""),
    "gene_names":               prot.get("Gene Names", ""),
    "organism":                 prot.get("Organism", ""),
    "length_aa":                prot.get("Length", ""),
    "sequence":                 prot.get("Sequence", ""),
    "ec_numbers":               prot.get("reaction_ecs", ""),
    "rhea_reactions":           prot.get("Rhea ID", ""),
    "interacting_chebi_ids":    prot.get("ChEBI ID", ""),
    "interacting_compounds":    prot.get("Compound Name", ""),
    "reaction_descriptions":    prot.get("reaction_descriptions", ""),
    "umap_1":                   prot.get("UMAP_1", ""),
    "umap_2":                   prot.get("UMAP_2", ""),
    "cluster":                  prot.get("clusters", ""),
    "prott5_cluster":           prot.get("prott5_cluster", ""),
    "is_literature_recruited":  prot.get("is_new", 0).fillna(0).astype(int),
    "paper_url":                prot.get("Paper", ""),
    "annotation":               prot.get("Annotation", ""),
    "sequence_source":          prot.get("Sequence_Source", ""),
    "identifier_type":          prot.get("Identifier_Type", ""),
    "uniprot_url":              prot["Entry"].apply(
        lambda e: f"https://www.uniprot.org/uniprotkb/{e}/entry" if isinstance(e, str) and e else ""),
    "alphafold_url":            prot["Entry"].apply(
        lambda e: f"https://alphafold.ebi.ac.uk/entry/{e}" if isinstance(e, str) and e else ""),
})
prot_out.to_csv(DST / "proteins.csv", index=False)
print(f"  {len(prot_out):,} rows · {len(prot_out.columns)} columns -> {DST/'proteins.csv'}")

# ─── MOLECULES ────────────────────────────────────────────────────────────
print("\nBuilding molecules.csv ...")
mol = pd.read_csv(SRC / "small_molecule_centric.csv", low_memory=False,
                  dtype={"ChEBI ID": str, "SMILES": str, "Compound Name": str})

mol_out = pd.DataFrame({
    "compound_name":            mol.get("Compound Name", ""),
    "chebi_id":                 mol.get("ChEBI ID", ""),
    "smiles":                   mol.get("SMILES", ""),
    "umap_1":                   mol.get("UMAP_1", ""),
    "umap_2":                   mol.get("UMAP_2", ""),
    "cluster":                  mol.get("clusters", ""),
    "is_literature_recruited":  mol.get("is_new", 0).fillna(0).astype(int),
    "paper_url":                mol.get("Paper", ""),
    "interacting_protein_accessions": mol.get("Entry", ""),
    "chebi_url":                mol.get("ChEBI ID", "").apply(
        lambda c: f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={c}" if str(c).strip() and str(c).lower() != 'nan' else ""),
})
mol_out.to_csv(DST / "molecules.csv", index=False)
print(f"  {len(mol_out):,} rows · {len(mol_out.columns)} columns -> {DST/'molecules.csv'}")

# ─── PROVENANCE (copy as-is) ───────────────────────────────────────────────
print("\nCopying literature_recruited_proteins.csv ...")
import shutil
shutil.copy2(SRC / "literature_recruited_proteins.csv", DST / "literature_recruited_proteins.csv")
print(f"  -> {DST/'literature_recruited_proteins.csv'}")

# Also copy the markdown version
shutil.copy2(SRC / "literature_recruited_proteins.md", DST / "literature_recruited_proteins.md")
print(f"  -> {DST/'literature_recruited_proteins.md'}")

print("\nDone.")
