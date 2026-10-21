"""Apply the Rhea re-audit drops to data/proteins.csv.

Runs AFTER analysis/13_rhea_reaudit.py finishes. Reads the per-entry
audit table and removes every accession classified as DROP_STALE — the
individual entries whose current UniProt annotation no longer contains
any sterane-passing ChEBI, GO steroid vocabulary, or steroid keyword.

Backs up the pre-drop CSV so a rollback is trivial.

Reads:  ../data/proteins.csv, ../analysis/rhea_audit.tsv
Writes: ../data/proteins.csv                          (updated in place)
        ../data/proteins.csv.pre_rhea_audit_backup    (untouched copy)
        ../analysis/rhea_curation_diff_report.txt     before/after summary
"""
from __future__ import annotations

import shutil
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
PROTEINS = ROOT / "data" / "proteins.csv"
BACKUP = ROOT / "data" / "proteins.csv.pre_rhea_audit_backup"
AUDIT = output("rhea_audit.tsv", historical=True)
REPORT = output("rhea_curation_diff_report.txt", historical=True)


def main() -> int:
    print(f"Reading {PROTEINS.name}...")
    prot = pd.read_csv(PROTEINS, low_memory=False)
    n_before = len(prot)
    print(f"  {n_before:,} entries")

    print(f"Reading {AUDIT.name}...")
    audit = pd.read_csv(AUDIT, sep="\t", low_memory=False)
    print(f"  {len(audit):,} audit rows")
    print(f"  decisions:")
    for d, n in audit["decision"].value_counts().items():
        print(f"    {n:>6,}  {d}")

    if not BACKUP.exists():
        print(f"\nBacking up {PROTEINS.name} -> {BACKUP.name}")
        shutil.copy2(PROTEINS, BACKUP)
    else:
        print(f"\nBackup already exists at {BACKUP.name} (leaving in place)")

    drop_accs = set(audit.loc[audit["decision"] == "DROP_STALE", "accession"])
    n_dropped = int(prot["accession"].isin(drop_accs).sum())
    prot = prot[~prot["accession"].isin(drop_accs)].copy()
    print(f"\nDropped {n_dropped:,} DROP_STALE entries")
    print(f"  atlas now: {len(prot):,}")

    prot.to_csv(PROTEINS, index=False)
    print(f"\nWrote updated {PROTEINS.name}")

    # Diff report
    kept_source_counts = prot["sequence_source"].fillna("").astype(str).value_counts()
    lines = []
    lines.append("=" * 66)
    lines.append("Rhea re-audit curation diff")
    lines.append("=" * 66)
    lines.append(f"Before:   {n_before:,} entries")
    lines.append(f"Dropped:  {n_dropped:,} entries (DROP_STALE — no current steroid evidence)")
    lines.append(f"After:    {len(prot):,} entries")
    lines.append("")
    lines.append(f"sequence_source distribution (after drops):")
    for src, n in kept_source_counts.items():
        if src:
            lines.append(f"  {n:>6,}  {src}")
    lines.append("")
    lines.append(f"Backup:   {BACKUP.name}")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {REPORT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
