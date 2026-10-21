"""Apply the pragmatic curation to data/proteins.csv.

Actions:
1. Back up the original data/proteins.csv to data/proteins.csv.pre_audit_backup
2. Drop the 231 orphans flagged as DROP_VITD (208) or DROP_NONE (23) in binder_audit.tsv
3. For every retained entry, write accurate per-entry provenance:
   - `sequence_source`  fixed per entry (Rhea vs binder-search)
   - `binder_evidence`  new column — specific annotations that caught orphan entries
   - `audit_decision`   new column — KEEP / KEEP_BY_NAME / KEEP_WEAK / REVIEW / (empty for non-orphans)
   - `audit_reason`     new column — the human-readable evidence string
4. Write updated data/proteins.csv

Reads:  ../data/proteins.csv, ../analysis/binder_audit.tsv
Writes: ../data/proteins.csv                        (updated in place)
        ../data/proteins.csv.pre_audit_backup       (untouched copy)
        ../analysis/curation_diff_report.txt        (before/after summary)
"""
from __future__ import annotations

import shutil
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
PROTEINS = ROOT / "data" / "proteins.csv"
BACKUP = ROOT / "data" / "proteins.csv.pre_audit_backup"
AUDIT = output("binder_audit.tsv", historical=True)
REPORT = output("curation_diff_report.txt", historical=True)

DROP_DECISIONS = {"DROP_VITD", "DROP_NONE"}

# Human-readable descriptions for the specific GO terms we recognize
GO_LABELS = {
    "GO:0005496": "steroid binding (umbrella)",
    "GO:0005497": "androgen binding",
    "GO:0015485": "cholesterol binding",
    "GO:0032934": "sterol binding",
    "GO:0032810": "SREBP response element binding",
    "GO:0038181": "bile acid receptor activity",
    "GO:0038186": "bile acid binding",
    "GO:1902121": "lithocholate binding",
    "GO:0032052": "bile acid binding",
    "GO:0070330": "aromatase activity",
    "GO:0004769": "steroid delta-isomerase",
    "GO:1990239": "steroid hormone binding",
    "GO:0003707": "steroid hormone receptor activity",
    "GO:0008395": "steroid hydroxylase activity",
    "GO:0016125": "sterol metabolic process",
    "GO:0016126": "sterol biosynthetic process",
    "GO:0034185": "apolipoprotein binding",
    "GO:0005500": "sterol carrier activity",
}
KW_LABELS = {
    "KW-0754": "Steroid-binding",
    "KW-0675": "Steroid hormone receptor",
}
STEROID_GO_TERMS = set(GO_LABELS.keys())
STEROID_KEYWORDS = set(KW_LABELS.keys())


def summarize_binder_evidence(kws: str, gos: str, ligs: str) -> str:
    """Build a human-readable evidence string for a single orphan entry."""
    kw_hits = [k for k in (kws or "").split(";") if k.strip() in STEROID_KEYWORDS]
    go_hits = [g for g in (gos or "").split(";") if g.strip() in STEROID_GO_TERMS]
    lig_hits = [l for l in (ligs or "").split(";") if l.strip()]

    parts = []
    for k in kw_hits:
        parts.append(f"{k}={KW_LABELS.get(k, k)}")
    for g in go_hits:
        parts.append(f"{g}={GO_LABELS.get(g, g)}")
    if lig_hits:
        parts.append(f"ligand_chebis={','.join(lig_hits[:5])}"
                     + ("..." if len(lig_hits) > 5 else ""))
    return " | ".join(parts) if parts else ""


def main() -> int:
    print(f"Reading {PROTEINS.name}...")
    prot = pd.read_csv(PROTEINS, low_memory=False)
    n_before = len(prot)
    print(f"  {n_before:,} entries")

    print(f"Reading {AUDIT.name}...")
    audit = pd.read_csv(AUDIT, sep="\t", low_memory=False)
    audit["keywords"] = audit["keywords"].fillna("").astype(str)
    audit["go_ids"] = audit["go_ids"].fillna("").astype(str)
    audit["ligand_chebis"] = audit["ligand_chebis"].fillna("").astype(str)
    print(f"  {len(audit):,} audit rows")

    # Back up
    if not BACKUP.exists():
        print(f"Backing up {PROTEINS.name} → {BACKUP.name}")
        shutil.copy2(PROTEINS, BACKUP)
    else:
        print(f"Backup already exists at {BACKUP.name} (leaving in place)")

    # 1) Drop the flagged accessions
    drop_accs = set(audit.loc[audit["decision"].isin(DROP_DECISIONS), "accession"])
    n_dropped = int(prot["accession"].isin(drop_accs).sum())
    prot = prot[~prot["accession"].isin(drop_accs)].copy()
    print(f"\nDropped {n_dropped:,} entries ({', '.join(sorted(DROP_DECISIONS))})")
    print(f"  atlas now: {len(prot):,}")

    # 2) Build per-entry provenance
    # Map audit decisions and evidence
    audit_indexed = audit.set_index("accession")

    def evidence_for(acc):
        if acc not in audit_indexed.index:
            return "", "", "", ""   # non-orphan (Rhea)
        r = audit_indexed.loc[acc]
        decision = str(r["decision"])
        reason = str(r["reason"])
        binder = summarize_binder_evidence(r["keywords"], r["go_ids"], r["ligand_chebis"])
        # Correct sequence_source
        source = "UniProt binder search (KW-0754 / GO:0005496 hierarchy)"
        return source, binder, decision, reason

    # Vectorize by iterating the rows we need to update
    orphan_mask = prot["accession"].isin(audit_indexed.index)
    n_orphans_kept = int(orphan_mask.sum())
    n_rhea = len(prot) - n_orphans_kept
    print(f"  {n_rhea:,} Rhea-catalytic entries")
    print(f"  {n_orphans_kept:,} binder-search entries kept")

    # Prepare new columns
    src = prot["sequence_source"].astype(str).fillna("").tolist()
    binder_col = [""] * len(prot)
    decision_col = [""] * len(prot)
    reason_col = [""] * len(prot)

    for i, acc in enumerate(prot["accession"].tolist()):
        source, binder, decision, reason = evidence_for(acc)
        if source:
            src[i] = source
            binder_col[i] = binder
            decision_col[i] = decision
            reason_col[i] = reason
        else:
            # Non-orphan → correct label for Rhea catalytic path
            src[i] = "Rhea catalytic reaction (rhea2uniprot mapping)"

    prot["sequence_source"] = src
    prot["binder_evidence"] = binder_col
    prot["audit_decision"] = decision_col
    prot["audit_reason"] = reason_col

    # 3) Write updated proteins.csv
    prot.to_csv(PROTEINS, index=False)
    print(f"\nWrote updated {PROTEINS.name}  ({len(prot):,} rows, {len(prot.columns)} columns)")

    # 4) Diff report
    lines = []
    lines.append("=" * 60)
    lines.append("Curation diff — pragmatic audit applied")
    lines.append("=" * 60)
    lines.append(f"Before:   {n_before:,} entries")
    lines.append(f"Dropped:  {n_dropped:,} entries (DROP_VITD + DROP_NONE)")
    lines.append(f"After:    {len(prot):,} entries")
    lines.append("")
    lines.append("New/updated columns:")
    lines.append("  sequence_source   — corrected per entry (Rhea vs binder-search)")
    lines.append("  binder_evidence   — specific KW/GO/ligand annotations that caught each binder-search entry")
    lines.append("  audit_decision    — KEEP / KEEP_BY_NAME / KEEP_WEAK / REVIEW / (empty for Rhea entries)")
    lines.append("  audit_reason      — human-readable evidence string")
    lines.append("")
    lines.append("sequence_source distribution (after):")
    for src, n in prot["sequence_source"].value_counts().items():
        lines.append(f"  {n:>6,}  {src}")
    lines.append("")
    lines.append("audit_decision distribution (kept binder-search entries only):")
    for dec, n in prot["audit_decision"].value_counts().items():
        if not dec:
            continue
        lines.append(f"  {n:>6,}  {dec}")
    lines.append("")
    lines.append(f"Backup:   {BACKUP.name}")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"Wrote diff report → {REPORT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
