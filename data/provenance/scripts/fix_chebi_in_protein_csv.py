"""Clean up the ChEBI ID column in protein_sequence_embedding.DEDUP.csv.

Two problems on this column:
  1. Some rows have `;`-separated lists of interacting metabolites where the
     entries use inconsistent formats: `17650.0` (float leak) alongside
     `CHEBI:50169;CHEBI:34964` (prefixed). The app does f'CHEBI:{id}' which
     produces `CHEBI:CHEBI:34964` for prefixed entries and `CHEBI:17650.0`
     (a broken URL) for float entries.
  2. Six of the CHEBI:-prefixed IDs point to the wrong compound (same audit
     as the molecule CSV: e.g. CHEBI:34964 -> guanyl-cysteine adduct instead
     of isoallopregnanolone).

This script:
  - Splits each ChEBI ID cell on ';'
  - For each token, strips 'CHEBI:' prefix
  - Rewrites the 6 wrong IDs to the correct ones
  - Strips trailing '.0' from float-leak tokens
  - Rejoins with ';'

The 7-entry ID mapping matches literature/fix_chebi_ids.py exactly.

Backup: protein_sequence_embedding.DEDUP.csv.bak6
"""
from __future__ import annotations
import re
import shutil
from pathlib import Path

import pandas as pd

CSV = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/protein_sequence_embedding.DEDUP.csv")

# Corrections (same 7 as literature/fix_chebi_ids.py)
WRONG_TO_CORRECT = {
    "34964": "11909",  # isoallopregnanolone
    "1156":  "34461",  # THDOC
    "27725": "30154",  # 5β-DHP
    "34958": "16229",  # epipregnanolone
    "16718": "1712",   # pregnanolone
    "34979": "2150",   # 5β-DHT
    "50169": "50169",  # brexanolone — right compound; format-only
}

def clean_token(tok: str) -> str:
    """Normalize one ChEBI-ID token: strip CHEBI: prefix, strip .0 float suffix,
    remap 6 wrong IDs to correct."""
    t = tok.strip()
    if not t or t.lower() == "nan":
        return ""
    # Strip CHEBI: prefix (case-insensitive)
    if t.lower().startswith("chebi:"):
        t = t.split(":", 1)[1].strip()
    # Strip trailing .0 float leak
    t = re.sub(r"\.0+$", "", t)
    # Apply corrections
    return WRONG_TO_CORRECT.get(t, t)

def clean_cell(cell: str) -> str:
    """Clean a full ';'-separated ChEBI ID cell."""
    if not isinstance(cell, str) or not cell.strip():
        return ""
    parts = [clean_token(t) for t in cell.split(";")]
    return ";".join(parts)

def main() -> int:
    df = pd.read_csv(CSV, low_memory=False, dtype={"ChEBI ID": str})
    print(f"Loaded {len(df):,} rows")

    bak = CSV.with_suffix(CSV.suffix + ".bak6")
    shutil.copy2(CSV, bak)
    print(f"Backup -> {bak.name}")

    # Track changes
    n_touched = 0
    n_corrected = 0
    changes: list[tuple[int, str, str]] = []
    for i, cell in df["ChEBI ID"].items():
        if not isinstance(cell, str) or not cell.strip():
            continue
        new = clean_cell(cell)
        if new != cell:
            n_touched += 1
            # Count corrections applied
            for w in WRONG_TO_CORRECT:
                if f"CHEBI:{w}" in cell and WRONG_TO_CORRECT[w] != w:
                    n_corrected += 1
            df.at[i, "ChEBI ID"] = new
            if len(changes) < 15:
                changes.append((i, cell[:200], new[:200]))

    print(f"\nRows touched:                  {n_touched}")
    print(f"Wrong-ID corrections applied:  {n_corrected}")
    print()
    print("--- Sample of changes ---")
    for i, old, new in changes:
        print(f"  row {i}:")
        print(f"    OLD: {old}")
        print(f"    NEW: {new}")
        print()

    df.to_csv(CSV, index=False)
    print(f"[OK] Wrote {CSV.name}")

    # Verify
    remaining_prefix = df[df["ChEBI ID"].astype(str).str.contains("CHEBI:", na=False)]
    print(f"\nRows still containing 'CHEBI:' after fix: {len(remaining_prefix)}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
