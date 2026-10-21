"""Fetch current UniProt annotation_score + entry status for every atlas entry.

The user's revised criterion: keep ONLY entries that (a) currently exist
in UniProtKB (not moved to UniParc / not deleted) AND (b) carry an
annotation_score of 5. This is much stricter than the original build
filter (score >= 4) and will drop most TrEMBL orthologs.

Reads:  ../data/proteins.csv
Writes: ../analysis/score_audit_cache.json    raw UniProt cache
        ../analysis/score_audit.tsv           per-entry decision table
        ../analysis/score_audit_summary.txt   human-readable summary
"""
from __future__ import annotations

import json
import time
from pathlib import Path
import sys

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
CACHE = output("score_audit_cache.json", historical=True)
OUT_TSV = output("score_audit.tsv", historical=True)
OUT_TXT = output("score_audit_summary.txt", historical=True)

BULK_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = "accession,annotation_score,reviewed"
BATCH_SIZE = 300


def fetch_bulk(accessions, retries=3):
    params = {"accessions": ",".join(accessions), "format": "json", "fields": FIELDS}
    for a in range(retries):
        try:
            r = requests.get(BULK_URL, params=params, timeout=60)
            if r.status_code == 200:
                return r.json().get("results", [])
            time.sleep(2 ** a)
        except requests.RequestException:
            time.sleep(2 ** a)
    return []


def fetch_single(acc, retries=2):
    """Direct entry endpoint: returns the entryType + annotationScore.

    For deleted entries, this returns entryType='Inactive' with an
    inactiveReason instead of the normal fields.
    """
    for a in range(retries):
        try:
            r = requests.get(f"{ENTRY_URL}/{acc}.json", timeout=15)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
            time.sleep(2 ** a)
        except requests.RequestException:
            time.sleep(2 ** a)
    return None


def extract(entry):
    acc = entry.get("primaryAccession", "")
    etype = entry.get("entryType", "")
    score = float(entry.get("annotationScore", 0)) if entry.get("annotationScore") else 0.0
    inactive = entry.get("inactiveReason", {}) or {}
    reason = inactive.get("deletedReason", "") or inactive.get("inactiveReasonType", "")
    return {
        "accession": acc,
        "entry_type": etype,
        "annotation_score": score,
        "inactive_reason": reason,
    }


def main() -> int:
    print(f"Reading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    accs = p["accession"].astype(str).tolist()
    n_total = len(accs)
    print(f"  {n_total:,} accessions")

    cache = {}
    if CACHE.exists():
        try:
            cache = json.loads(CACHE.read_text())
            print(f"  loaded {len(cache):,} cached")
        except Exception:
            pass

    remaining = [a for a in accs if a not in cache]
    print(f"\nBulk fetch {len(remaining):,} in batches of {BATCH_SIZE}...")
    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i:i + BATCH_SIZE]
        entries = fetch_bulk(batch)
        for e in entries:
            ann = extract(e)
            if ann["accession"]:
                cache[ann["accession"]] = ann
        if (i // BATCH_SIZE) % 10 == 9:
            CACHE.write_text(json.dumps(cache))
            print(f"  checkpoint: {len(cache):,} cached", flush=True)
        time.sleep(0.12)

    # Direct-entry fallback (TrEMBL entries the bulk endpoint skips)
    missing = [a for a in accs if a not in cache]
    if missing:
        print(f"\nDirect-endpoint retry for {len(missing):,}...")
        for j, acc in enumerate(missing):
            if j % 50 == 0:
                print(f"  {j}/{len(missing)}", flush=True)
            entry = fetch_single(acc)
            if entry is not None:
                ann = extract(entry)
                if ann["accession"]:
                    cache[ann["accession"]] = ann
                if ann["accession"] != acc:
                    cache[acc] = ann
            time.sleep(0.08)

    CACHE.write_text(json.dumps(cache))
    print(f"\nCache: {len(cache):,} / {n_total:,}")

    # Build per-entry decision
    print("\nClassifying...")
    rows = []
    for _, r in p.iterrows():
        acc = r["accession"]
        ann = cache.get(acc, {"entry_type": "MISSING", "annotation_score": 0.0,
                              "inactive_reason": "not fetched"})
        etype = ann.get("entry_type", "")
        score = float(ann.get("annotation_score", 0) or 0)
        reason = ann.get("inactive_reason", "") or ""

        if etype.startswith("Inactive") or "DELETED" in reason.upper() or "Not part" in reason:
            decision = "DROP_INACTIVE"
            note = f"UniProt inactive: {reason}"
        elif score < 5.0:
            decision = "DROP_LOWSCORE"
            note = f"annotation_score {score:g} < 5"
        elif score == 5.0:
            decision = "KEEP"
            note = "annotation_score 5 and active in UniProtKB"
        else:
            decision = "REVIEW"
            note = f"unexpected: entry_type={etype}, score={score}"

        rows.append({
            "accession": acc,
            "entry_name": r["entry_name"],
            "protein_names": str(r["protein_names"])[:100],
            "organism": r["organism"],
            "entry_type": etype,
            "annotation_score": score,
            "inactive_reason": reason,
            "decision": decision,
            "note": note,
        })

    audit = pd.DataFrame(rows)
    audit.to_csv(OUT_TSV, sep="\t", index=False)

    # Summary
    lines = []
    lines.append("=" * 76)
    lines.append(f"ANNOTATION-SCORE + STATUS AUDIT — all {n_total:,} atlas entries")
    lines.append("=" * 76)
    lines.append("")
    lines.append("Decision counts:")
    for d, n in audit["decision"].value_counts().items():
        lines.append(f"  {n:>7,}  {d}")
    lines.append("")
    lines.append("Annotation score distribution (all entries):")
    for s in sorted(audit["annotation_score"].unique()):
        n = int((audit["annotation_score"] == s).sum())
        lines.append(f"  score {s:>3.1f}:  {n:>7,}")
    lines.append("")
    lines.append("Entry type distribution:")
    for et, n in audit["entry_type"].value_counts().items():
        lines.append(f"  {n:>7,}  {et}")
    lines.append("")

    n_keep = int((audit["decision"] == "KEEP").sum())
    n_drop = int(audit["decision"].isin(["DROP_INACTIVE", "DROP_LOWSCORE"]).sum())
    lines.append(f"Recommended action (score-5-only + active): drop {n_drop:,} → atlas {n_total:,} → {n_keep:,}")
    OUT_TXT.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {OUT_TSV.name}, {OUT_TXT.name}")
    print("\n=== Decision summary ===")
    for d, n in audit["decision"].value_counts().items():
        print(f"  {n:>7,}  {d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
