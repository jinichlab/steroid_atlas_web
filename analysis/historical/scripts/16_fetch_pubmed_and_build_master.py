"""Build the final master proteins.csv.

Runs AFTER the score audit (15_annotation_score_audit.py) and the Rhea
re-audit (13_rhea_reaudit.py) have both finished. Computes the intersection
KEEP set — entries that pass EVERY criterion:
    1. annotation_score == 5
    2. entry is active in UniProtKB (not UniParc / not deleted)
    3. entry has some steroid signal per the Rhea audit
       (KEEP_RHEA, KEEP_LIGAND, or KEEP_ANNOTATION);
       or the entry is in the binder-search set with an audit_decision of
       KEEP, KEEP_BY_NAME, KEEP_WEAK, or REVIEW (from the earlier binder
       audit — those have already been filtered against secosteroids).

Then fetches PubMed IDs for the surviving entries and adds a pubmed_ids
column. Backs up the current proteins.csv and writes the final master.

Reads:
    ../data/proteins.csv
    ../analysis/score_audit.tsv
    ../analysis/rhea_audit.tsv
Writes:
    ../analysis/pubmed_cache.json                  raw pubmed fetch cache
    ../data/proteins.csv.pre_master_table_backup   backup of current CSV
    ../data/proteins.csv                           final master table
    ../analysis/master_table_report.txt            summary of the assembly
"""
from __future__ import annotations

import json
import shutil
import time
from pathlib import Path
import sys

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
PROTEINS = ROOT / "data" / "proteins.csv"
BACKUP = ROOT / "data" / "proteins.csv.pre_master_table_backup"
SCORE = output("score_audit.tsv", historical=True)
RHEA = output("rhea_audit.tsv", historical=True)
PUBMED_CACHE = output("pubmed_cache.json", historical=True)
REPORT = output("master_table_report.txt", historical=True)

BULK_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = "accession,lit_pubmed_id"
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
    for a in range(retries):
        try:
            r = requests.get(f"{ENTRY_URL}/{acc}.json",
                             params={"fields": FIELDS}, timeout=15)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
            time.sleep(2 ** a)
        except requests.RequestException:
            time.sleep(2 ** a)
    return None


def extract_pubmeds(entry) -> str:
    """Concatenate all PubMed IDs cited by this entry into a ';'-joined string."""
    ids = []
    for ref in entry.get("references", []):
        for xref in (ref.get("citation", {}) or {}).get("citationCrossReferences", []) or []:
            if xref.get("database") == "PubMed":
                pid = xref.get("id")
                if pid:
                    ids.append(str(pid))
    # dedup, keep order
    seen = set(); out = []
    for x in ids:
        if x not in seen:
            seen.add(x); out.append(x)
    return ";".join(out)


def do_pubmed_fetch(accessions):
    cache = {}
    if PUBMED_CACHE.exists():
        try:
            cache = json.loads(PUBMED_CACHE.read_text())
            print(f"  loaded {len(cache):,} pubmed-cached")
        except Exception:
            pass
    remaining = [a for a in accessions if a not in cache]
    print(f"  {len(remaining):,} to fetch (batch {BATCH_SIZE})")
    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i:i + BATCH_SIZE]
        entries = fetch_bulk(batch)
        for e in entries:
            acc = e.get("primaryAccession", "")
            if acc:
                cache[acc] = extract_pubmeds(e)
        if (i // BATCH_SIZE) % 10 == 9:
            PUBMED_CACHE.write_text(json.dumps(cache))
            print(f"    checkpoint: {len(cache):,} cached", flush=True)
        time.sleep(0.12)
    # Direct fallback for TrEMBL misses
    missing = [a for a in accessions if a not in cache]
    if missing:
        print(f"  direct-endpoint retry: {len(missing):,}")
        for j, acc in enumerate(missing):
            if j % 50 == 0:
                print(f"    {j}/{len(missing)}", flush=True)
            e = fetch_single(acc)
            if e is not None:
                cache[acc] = extract_pubmeds(e)
            time.sleep(0.08)
    PUBMED_CACHE.write_text(json.dumps(cache))
    return cache


def main() -> int:
    print(f"Reading {PROTEINS.name}, {SCORE.name}, {RHEA.name}...")
    prot = pd.read_csv(PROTEINS, low_memory=False)
    score = pd.read_csv(SCORE, sep="\t", low_memory=False)
    rhea = pd.read_csv(RHEA, sep="\t", low_memory=False)
    n_before = len(prot)
    print(f"  atlas: {n_before:,}   score audit: {len(score):,}   rhea audit: {len(rhea):,}")

    # ── Compute final KEEP set ─────────────────────────────────────────────
    score_pass = set(score.loc[score["decision"] == "KEEP", "accession"])
    rhea_keep_set = set(rhea.loc[
        rhea["decision"].isin(["KEEP_RHEA", "KEEP_LIGAND", "KEEP_ANNOTATION"]),
        "accession",
    ])
    # Binder-search entries were already curated in the earlier binder audit;
    # any survivor there counts as a steroid keep by prior evaluation.
    binder_kept = set(prot.loc[
        prot["audit_decision"].fillna("").isin(["KEEP", "KEEP_BY_NAME", "KEEP_WEAK", "REVIEW"]),
        "accession",
    ])
    steroid_kept = rhea_keep_set | binder_kept

    final_keep = score_pass & steroid_kept
    print(f"\nFinal KEEP set (intersection):")
    print(f"  score==5:                      {len(score_pass):,}")
    print(f"  steroid-associated:            {len(steroid_kept):,}")
    print(f"  INTERSECTION (final KEEP):     {len(final_keep):,}")

    to_drop = set(prot["accession"]) - final_keep
    print(f"  → drop {len(to_drop):,}, keep {len(final_keep):,}")

    # ── Fetch PubMed IDs for the survivors only ────────────────────────────
    kept_accs = sorted(final_keep)
    print(f"\nFetching PubMed IDs for the {len(kept_accs):,} survivors...")
    pubmed = do_pubmed_fetch(kept_accs)
    print(f"  pubmed cache: {len(pubmed):,} accessions ready")

    # ── Back up + apply drops + add columns ────────────────────────────────
    if not BACKUP.exists():
        print(f"\nBacking up {PROTEINS.name} → {BACKUP.name}")
        shutil.copy2(PROTEINS, BACKUP)

    master = prot[prot["accession"].isin(final_keep)].copy()
    # Attach annotation_score
    score_map = dict(zip(score["accession"], score["annotation_score"]))
    master["annotation_score"] = master["accession"].map(score_map)
    # Attach pubmed_ids
    master["pubmed_ids"] = master["accession"].map(lambda a: pubmed.get(a, ""))
    # Attach a pubmed_count for convenience
    master["pubmed_count"] = master["pubmed_ids"].apply(
        lambda s: 0 if not s else len(s.split(";"))
    )

    # Reorder columns so the master table reads left-to-right in a sensible way
    preferred_order = [
        "accession", "entry_name", "protein_names", "gene_names", "organism",
        "length_aa", "sequence",
        "annotation_score",
        "ec_numbers", "rhea_reactions", "interacting_chebi_ids",
        "interacting_compounds", "reaction_descriptions",
        "go_ids", "go_labels", "keyword_ids", "keyword_labels",
        "pubmed_ids", "pubmed_count",
        "umap_1", "umap_2", "cluster", "prott5_cluster",
        "is_literature_recruited", "paper_url",
        "sequence_source", "binder_evidence", "audit_decision", "audit_reason",
        "identifier_type", "annotation",
        "uniprot_url", "alphafold_url",
    ]
    ordered = [c for c in preferred_order if c in master.columns]
    other = [c for c in master.columns if c not in ordered]
    master = master[ordered + other]

    master.to_csv(PROTEINS, index=False)
    print(f"\nWrote master {PROTEINS.name}  ({len(master):,} rows × {len(master.columns)} cols)")

    # ── Report ────────────────────────────────────────────────────────────
    n_with_pm = int(master["pubmed_count"].fillna(0).gt(0).sum())
    median_pm = int(master["pubmed_count"].median())
    lines = []
    lines.append("=" * 74)
    lines.append("FINAL MASTER TABLE")
    lines.append("=" * 74)
    lines.append(f"Before:   {n_before:,} entries")
    lines.append(f"Filters applied (intersection):")
    lines.append(f"  1. annotation_score == 5             passed: {len(score_pass):,}")
    lines.append(f"  2. active in UniProtKB               (subset of #1)")
    lines.append(f"  3. steroid-associated per audits     passed: {len(steroid_kept):,}")
    lines.append(f"  → intersection kept:                 {len(final_keep):,}")
    lines.append(f"  → dropped:                           {n_before - len(final_keep):,}")
    lines.append("")
    lines.append(f"Master table:  {len(master):,} rows × {len(master.columns)} columns")
    lines.append(f"  columns:  {', '.join(master.columns)}")
    lines.append("")
    lines.append(f"PubMed coverage:")
    lines.append(f"  entries with ≥1 PubMed ID:  {n_with_pm:,} / {len(master):,}  "
                 f"({n_with_pm/len(master)*100:.1f}%)")
    lines.append(f"  median PubMed IDs per entry: {median_pm}")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {REPORT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
