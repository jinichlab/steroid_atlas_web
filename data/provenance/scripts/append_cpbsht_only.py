"""Append CpBSH/T (WP_243289361) to protein_sequence_embedding.DEDUP.csv.

This is the ONLY Guzior 2024 entry we add to the atlas — the single BSH/T with
direct biochemical evidence (Fig 1: kinetics + substrate scope + mutants).

The other 19 fetched Guzior entries are deliberately NOT added — no direct
protein-level activity data. A record of that decision is written to
literature/guzior2024_dismissed_entries.json for provenance.

UMAP coordinates left as NaN — user's memory: "Review literature before
embedding" — ProtT5 embedding + UMAP projection is a separate step, run only
after explicit review.

Run:
  LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \
  /home/adsiordia/miniconda3/bin/python \
  /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/literature/append_cpbsht_only.py
"""
from __future__ import annotations
import json
import shutil
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
CSV = HERE.parent / "protein_sequence_embedding.DEDUP.csv"
FASTA = HERE / "sequences" / "guzior" / "WP_243289361.fasta"
REVIEW_JSON = HERE / "guzior2024_review_ready.json"
DISMISSED = HERE / "guzior2024_dismissed_entries.json"

DOI = "https://doi.org/10.1038/s41586-024-07017-8"

# ─── Load fetched sequence ──────────────────────────────────────────────────
fasta_text = FASTA.read_text()
lines = fasta_text.strip().split("\n")
fetched_header = lines[0].lstrip(">").strip()
sequence = "".join(l.strip() for l in lines[1:]).replace("*", "")
print(f"Sequence:  {len(sequence)} aa")
print(f"Header:    {fetched_header[:100]}")

# ─── Rich annotation encoding Fig 1 biochemistry ────────────────────────────
protein_names = (
    "Bile salt hydrolase/transferase (BSH/T) "
    "(Bifunctional bile acid amine N-acyltransferase / choloylglycine hydrolase) "
    "(EC 2.3.1.-) (EC 3.5.1.24)"
)

annotation = (
    "Guzior et al., Nature 2024 (CpBSH/T from C. perfringens ATCC 13124 — "
    "the ONLY protein biochemically characterized in this paper: kinetics, "
    "substrate scope, N82Y and C2A mutagenesis, structural basis via PDB 2BJG) | "
    + DOI
)

reaction_descriptions = (
    "Bifunctional enzyme, single active site (Cys2 nucleophile). "
    "(1) Bile salt hydrolysis (EC 3.5.1.24): cleaves conjugated BAs (TCA/GCA -> CA + Tauro/Gly), "
    "active pH 3-7, 8 mM TCA saturates. "
    "(2) NOVEL Bile acid amine N-acyltransferase (EC 2.3.1.-): conjugates amino acids to CA "
    "forming microbially conjugated bile acids (MCBAs). "
    "pH optimum 5.3; at peak reaches 7% of hydrolysis rate (1 amino acid transferred per 15 TCA "
    "molecules hydrolyzed). Substrate scope: 16/20 amino acids from TCA + AA mix, 11/19 from GCA + "
    "AA mix, 12/20 from CA + AA mix. Never conjugated: Pro (secondary amine blocks nucleophilic "
    "attack), Asp. PheCA kinetics linear in [Phe] up to 5 mM. "
    "Unified mechanism: Cys2 nucleophilic attack -> enzyme-CA covalent intermediate, then water "
    "hydrolyzes (hydrolase pathway) OR amino acid attacks (acyltransferase pathway). "
    "Mutants: C2A abolishes BOTH activities (shared catalytic nucleophile); N82Y preserves "
    "activity but shifts amino-acid specificity (loses GluCA/LysCA/LeuCA, gains AlaCA - Asn82 "
    "shapes substrate preference). MCBA profile cluster 1, BSH/T phylogenetic group I."
)

sequence_source = (
    "RefSeq WP_243289361 [C. perfringens ATCC 13124]. Literature-recruited via Guzior et al., "
    "Nature 2024 Supp Table 1. This is the ONLY entry from Guzior 2024 added to the atlas "
    "because it is the only protein for which the paper provides direct biochemical evidence "
    "(purified enzyme + kinetics + N82Y/C2A mutants + PDB 2BJG structure). Species: "
    "Clostridium perfringens strain ATCC 13124 (Guzior 2024 strain source). MCBA profile "
    "cluster 1, BSH/T phylogenetic group I. Structure: PDB 2BJG (Kumar et al., refs 17,19)."
)

new_row = {
    "Entry": "WP_243289361",
    "Entry Name": "BSHT_CLOPE_ATCC13124",
    "Protein names": protein_names,
    "Gene Names": "bsh bsh/T",
    "Organism": "Clostridium perfringens (strain ATCC 13124)",
    "Length": len(sequence),
    "Sequence": sequence,
    "Annotation": annotation,
    "source": "Guzior 2024 Supp Table 1 (CpBSH/T)",
    "is_new": 1,
    "Paper": DOI,
    "reaction_ecs": "EC 3.5.1.24; EC 2.3.1.-",
    "reaction_descriptions": reaction_descriptions,
    "Sequence_Source": sequence_source,
    "Identifier_Type": "RefSeq_WP",
}

# ─── Load current CSV, add row, save ─────────────────────────────────────────
print(f"\nLoading {CSV.name}...")
df = pd.read_csv(CSV, low_memory=False)
print(f"  {len(df):,} rows currently")

# Duplicate check
existing = df["Entry"].astype(str).str.upper() == "WP_243289361"
if existing.any():
    print(f"ERROR: WP_243289361 already in CSV at row(s) {list(df.index[existing])}")
    return_code = 1
else:
    # Backup
    bak = CSV.with_suffix(CSV.suffix + ".bak2")
    shutil.copy2(CSV, bak)
    print(f"  backup -> {bak.name}")

    # Fill columns present in df but absent from new_row with NaN
    row_to_append = {c: new_row.get(c, pd.NA) for c in df.columns}
    df.loc[len(df)] = row_to_append

    df.to_csv(CSV, index=False)
    print(f"\n[OK] Appended CpBSH/T (WP_243289361)")
    print(f"  new row count: {len(df):,}")
    print(f"  is_new=1 rows: {int((df['is_new'] == 1).sum())} (was 15; +1 = 16 expected)")

    # Verify the row
    written = df[df["Entry"].astype(str) == "WP_243289361"].iloc[0]
    print(f"\nVerification (row read back):")
    for f in ("Entry", "Organism", "Length", "is_new", "Paper",
              "reaction_ecs", "Identifier_Type"):
        val = str(written.get(f, ""))[:100]
        print(f"  {f}: {val}")

    return_code = 0

# ─── Write dismissal record for the other 19 fetched entries ────────────────
review = json.loads(REVIEW_JSON.read_text())
dismissed = [r for r in review if r["accession_normalized"] != "WP_243289361"]
dismissal_record = {
    "decision": "Not added to atlas",
    "date_decided": "2026-07-23",
    "reason": (
        "No direct protein-level biochemical evidence in Guzior et al. 2024. "
        "Only CpBSH/T (WP_243289361) was purified, kinetically characterized, "
        "and mutagenized (N82Y, C2A). All other Table 1 entries carry only "
        "species-level MCBA-culture activity, which cannot be unambiguously "
        "attributed to a specific paralog protein. Kept in "
        "literature/guzior2024_review_ready.json for future reference if "
        "additional experimental evidence emerges."
    ),
    "n_dismissed": len(dismissed),
    "dismissed_accessions": [r["accession_normalized"] for r in dismissed],
    "paper_doi": DOI,
}
DISMISSED.write_text(json.dumps(dismissal_record, indent=2))
print(f"\nDismissal record -> {DISMISSED.name}  ({len(dismissed)} entries)")

raise SystemExit(return_code)
