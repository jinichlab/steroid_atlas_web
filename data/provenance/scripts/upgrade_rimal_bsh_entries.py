"""Upgrade Rimal 2024 BSH annotations in protein_sequence_embedding.DEDUP.csv.

Changes:
  - P0DXD2 (BlBSH, B. longum NCTC 11818):
      is_new: 0 -> 1
      Paper: nan -> Rimal 2024 DOI
      Annotation: fills in Rimal citation with "primary purified enzyme" note
      Sequence_Source: notes literature recruitment via Rimal 2024

  - Q5LF84 (B. fragilis NCTC 9343 BSH):
      Clears the stale "UNVERIFIED / should likely be dropped" note in
      Sequence_Source now that the audit is complete and the paper
      attribution is confirmed correct (used as KO / complementation host
      in Rimal 2024).

Backup: writes <csv>.bak alongside before overwriting.

Run:
  LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \
  /home/adsiordia/miniconda3/bin/python \
  /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/literature/upgrade_rimal_bsh_entries.py
"""
from __future__ import annotations
import shutil
from pathlib import Path

import pandas as pd

CSV = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/protein_sequence_embedding.DEDUP.csv")

RIMAL_DOI = "https://doi.org/10.1038/s41586-023-06990-w"

P0DXD2_ANNOTATION = (
    "Rimal et al., Nature 2024 (BlBSH from B. longum NCTC 11818 — "
    "primary purified enzyme; novel amine N-acyltransferase activity "
    "forming BBAAs) | " + RIMAL_DOI
)
P0DXD2_SEQ_SRC = (
    "UniProt accession — B. longum NCTC 11818 BSH (P0DXD2, Swiss-Prot); "
    "literature-recruited via Rimal 2024 as the primary biochemically "
    "characterized BSH/T (amine N-acyltransferase, BBAA biosynthesis)."
)

Q5LF84_SEQ_SRC = (
    "UniProt accession — B. fragilis NCTC 9343 BSH; literature-recruited "
    "via Rimal 2024. Role in paper: KO / complementation host — deletion "
    "of this gene ablates BBAA production in vivo; complementation "
    "restores it. In vivo genetic evidence for the same novel amine "
    "N-acyltransferase activity characterized biochemically on BlBSH (P0DXD2)."
)

def main() -> int:
    print(f"Loading {CSV.name}...")
    df = pd.read_csv(CSV, low_memory=False)
    print(f"  {len(df):,} rows loaded")

    p0_mask = df["Entry"].astype(str).str.upper() == "P0DXD2"
    q5_mask = df["Entry"].astype(str).str.upper() == "Q5LF84"

    if p0_mask.sum() != 1:
        print(f"ERROR: expected 1 P0DXD2 row, got {p0_mask.sum()}")
        return 1
    if q5_mask.sum() != 1:
        print(f"ERROR: expected 1 Q5LF84 row, got {q5_mask.sum()}")
        return 1

    # Backup
    bak = CSV.with_suffix(CSV.suffix + ".bak")
    shutil.copy2(CSV, bak)
    print(f"  backup written -> {bak.name}")

    # --- P0DXD2 upgrade ---
    print("\nP0DXD2 (BlBSH) — before:")
    print(f"  is_new={df.loc[p0_mask, 'is_new'].iloc[0]}")
    print(f"  Paper={df.loc[p0_mask, 'Paper'].iloc[0]!r}")
    print(f"  Annotation={df.loc[p0_mask, 'Annotation'].iloc[0]!r}")
    print(f"  Sequence_Source={df.loc[p0_mask, 'Sequence_Source'].iloc[0]!r}")

    df.loc[p0_mask, "is_new"] = 1
    df.loc[p0_mask, "Paper"] = RIMAL_DOI
    df.loc[p0_mask, "Annotation"] = P0DXD2_ANNOTATION
    df.loc[p0_mask, "Sequence_Source"] = P0DXD2_SEQ_SRC

    print("\nP0DXD2 (BlBSH) — after:")
    print(f"  is_new={df.loc[p0_mask, 'is_new'].iloc[0]}")
    print(f"  Paper={df.loc[p0_mask, 'Paper'].iloc[0]!r}")
    print(f"  Annotation={df.loc[p0_mask, 'Annotation'].iloc[0]!r}")
    print(f"  Sequence_Source={df.loc[p0_mask, 'Sequence_Source'].iloc[0]!r}")

    # --- Q5LF84 clarification ---
    print("\nQ5LF84 (B. fragilis BSH) — before:")
    print(f"  Sequence_Source={df.loc[q5_mask, 'Sequence_Source'].iloc[0]!r}")

    df.loc[q5_mask, "Sequence_Source"] = Q5LF84_SEQ_SRC

    print("\nQ5LF84 (B. fragilis BSH) — after:")
    print(f"  Sequence_Source={df.loc[q5_mask, 'Sequence_Source'].iloc[0]!r}")

    # Save
    df.to_csv(CSV, index=False)
    print(f"\n✓ wrote {CSV.name}  ({len(df):,} rows)")

    n_new = int((df["is_new"] == 1).sum())
    print(f"\nTotal is_new=1 rows now: {n_new} (was 14; +1 = 15 expected)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
