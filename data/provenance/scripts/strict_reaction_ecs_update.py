"""Trim reaction_ecs to only ECs experimentally confirmed by the cited paper.

Policy: reaction_ecs should list only ECs that the paper cited in this row's
`Paper` field directly demonstrates for THIS specific protein. Sequence-
similarity inferences are excluded.

Changes:
  - P0DXD2:       3.5.1.24;3.5.1.74           -> 2.3.1.-;3.5.1.-
                  (Rimal 2024 experimentally confirmed via Swiss-Prot ECO:0000269;
                   the specific 3.5.1.24/74 were UniProt ECO:0000250 similarity)
  - Q5LF84:       nan                          -> 2.3.1.-
                  (Rimal 2024 Δbsh KO ablates BBAA production in vivo — genetic
                   evidence for the acyltransferase activity)
  - WP_243289361: "EC 3.5.1.24; EC 2.3.1.-"    -> 3.5.1.24;2.3.1.-
                  (Guzior 2024 direct kinetics for both; format normalized)

All other is_new=1 entries remain nan (no per-protein EC assignment in cited paper).

Writes backup to protein_sequence_embedding.DEDUP.csv.bak5
"""
from __future__ import annotations
import shutil
from pathlib import Path

import pandas as pd

CSV = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/protein_sequence_embedding.DEDUP.csv")

UPDATES = {
    "P0DXD2":       "2.3.1.-;3.5.1.-",
    "Q5LF84":       "2.3.1.-",
    "WP_243289361": "3.5.1.24;2.3.1.-",
}

def main() -> int:
    df = pd.read_csv(CSV, low_memory=False)
    print(f"Loaded {len(df):,} rows")

    bak = CSV.with_suffix(CSV.suffix + ".bak5")
    shutil.copy2(CSV, bak)
    print(f"Backup -> {bak.name}")

    print("\n--- reaction_ecs updates ---")
    for acc, new_val in UPDATES.items():
        mask = df["Entry"] == acc
        if mask.sum() != 1:
            print(f"  {acc}: SKIP — {mask.sum()} rows found")
            continue
        old = df.loc[mask, "reaction_ecs"].iloc[0]
        df.loc[mask, "reaction_ecs"] = new_val
        print(f"  {acc}:  {old!r:<40}  ->  {new_val!r}")

    df.to_csv(CSV, index=False)
    print(f"\n[OK] Wrote {CSV.name}")

    # Verify
    print("\n--- Post-update state (all is_new=1 rows) ---")
    for _, row in df[df["is_new"] == 1].iterrows():
        print(f"  {row['Entry']:<15}  reaction_ecs = {row['reaction_ecs']!r}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
