"""Fix wrong ChEBI IDs in small_molecule_centric.csv.

The 7 CHEBI:-prefixed rows in the CSV were added later (all is_new=1) by a
separate pipeline that:
  - Formatted IDs with 'CHEBI:' prefix (unlike the other 589 numeric-only rows)
  - Assigned wrong ChEBI numbers for 6 of the 7 compounds

This script:
  1. Corrects the 6 wrong ChEBI IDs to the verified-correct ones (via OLS/ChEBI lookup)
  2. Strips 'CHEBI:' prefix to match the numeric-only format
     -> also fixes the 'CHEBI:CHEBI:...' doubling in the app display

Does NOT auto-merge with existing duplicate rows — those may be independent
literature entries with distinct SMILES/naming.

Backup: small_molecule_centric.csv.bak
"""
from __future__ import annotations
import shutil
from pathlib import Path

import pandas as pd

CSV = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/small_molecule_centric.csv")

FIXES = {
    # wrong CHEBI:-prefixed ID  ->  (correct numeric ID, compound display name for logging)
    "CHEBI:34964": ("11909", "isoallopregnanolone (3beta,5alpha-THP)"),
    "CHEBI:1156":  ("34461", "3alpha,5alpha-tetrahydrodeoxycorticosterone (THDOC)"),
    "CHEBI:27725": ("30154", "5beta-dihydroprogesterone"),
    "CHEBI:34958": ("16229", "epipregnanolone (3beta,5beta-THP)"),
    "CHEBI:16718": ("1712",  "pregnanolone (3alpha,5beta-THP)"),
    "CHEBI:34979": ("2150",  "5beta-dihydrotestosterone"),
    # This one was already the correct compound (brexanolone = allopregnanolone);
    # just strip the CHEBI: prefix for format consistency
    "CHEBI:50169": ("50169", "allopregnanolone (brexanolone) [correct compound, format-only fix]"),
}

def main() -> int:
    df = pd.read_csv(CSV, low_memory=False)
    print(f"Loaded {len(df)} molecule rows")

    bak = CSV.with_suffix(CSV.suffix + ".bak")
    shutil.copy2(CSV, bak)
    print(f"Backup -> {bak.name}")

    print("\n--- Applying fixes ---")
    for wrong_id, (correct_id, label) in FIXES.items():
        mask = df["ChEBI ID"].astype(str) == wrong_id
        n = mask.sum()
        if n == 0:
            print(f"  SKIP: {wrong_id!r} not found in CSV")
            continue
        df.loc[mask, "ChEBI ID"] = correct_id
        marker = "format-only" if wrong_id == "CHEBI:50169" else "WRONG->CORRECT"
        print(f"  [{marker}]  {wrong_id!r:<15} -> {correct_id!r:<8}  ({n} row) — {label}")

    df.to_csv(CSV, index=False)
    print(f"\n[OK] Wrote {CSV.name}")

    # Verify no CHEBI:-prefixed IDs remain
    remaining = df[df["ChEBI ID"].astype(str).str.startswith("CHEBI:")]
    print(f"\nRemaining CHEBI:-prefixed rows: {len(remaining)}  (should be 0)")
    if len(remaining):
        for _, r in remaining.iterrows():
            print(f"  {r['ChEBI ID']!r}  {r['Compound Name']!r}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
