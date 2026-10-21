"""Fetch GO IDs + labels + keyword IDs + keyword labels for every atlas entry
and merge them into data/proteins.csv as new columns.

Adds four columns:
  go_ids           semicolon-joined GO IDs (e.g. "GO:0005496;GO:0003707")
  go_labels        semicolon-joined GO labels (e.g. "steroid binding;steroid hormone receptor activity")
  keyword_ids      semicolon-joined KW-#### ids
  keyword_labels   semicolon-joined KW labels

The visualizer's protein table uses these columns to make GO/keyword text searchable
alongside protein name / gene / sequence.

Reads:  ../data/proteins.csv
Writes: ../data/proteins.csv                        (updated in place)
        ../analysis/all_annotations_cache.json      (raw UniProt cache — resumable)
        ../analysis/all_annotations.tsv             (accession → GO/keyword table for reference)
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
PROTEINS = ROOT / "data" / "proteins.csv"
CACHE = output("all_annotations_cache.json", historical=True)
OUT_TSV = output("all_annotations.tsv", historical=True)

BULK_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = "accession,keyword,go"          # minimal — GO xrefs carry labels, keywords carry labels
BATCH_SIZE = 300


def fetch_bulk(accessions, retries=3):
    params = {"accessions": ",".join(accessions), "format": "json", "fields": FIELDS}
    for attempt in range(retries):
        try:
            r = requests.get(BULK_URL, params=params, timeout=60)
            if r.status_code == 200:
                return r.json().get("results", [])
            print(f"  ! HTTP {r.status_code} on batch of {len(accessions)} (attempt {attempt+1})", flush=True)
            time.sleep(2 ** attempt)
        except requests.RequestException as e:
            print(f"  ! {e} (attempt {attempt+1})", flush=True)
            time.sleep(2 ** attempt)
    return []


def fetch_single(accession, retries=2):
    params = {"format": "json", "fields": FIELDS}
    for attempt in range(retries):
        try:
            r = requests.get(f"{ENTRY_URL}/{accession}.json", params=params, timeout=30)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
            time.sleep(2 ** attempt)
        except requests.RequestException:
            time.sleep(2 ** attempt)
    return None


def extract(entry):
    """Return {accession, go_ids, go_labels, keyword_ids, keyword_labels}."""
    acc = entry.get("primaryAccession", "")
    go_pairs = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        gid = xref.get("id", "")
        label = ""
        for p in xref.get("properties", []) or []:
            if p.get("key") == "GoTerm":
                # value looks like "F:steroid binding" — strip the F:/P:/C: prefix
                v = p.get("value", "")
                label = v.split(":", 1)[1] if ":" in v else v
                break
        if gid:
            go_pairs.append((gid, label))

    kw_pairs = []
    for k in entry.get("keywords", []):
        kid = k.get("id", "")
        klabel = k.get("name", "") or k.get("value", "")
        if kid:
            kw_pairs.append((kid, klabel))

    return {
        "accession": acc,
        "go_ids": ";".join(g for g, _ in go_pairs),
        "go_labels": ";".join(l for _, l in go_pairs),
        "keyword_ids": ";".join(k for k, _ in kw_pairs),
        "keyword_labels": ";".join(l for _, l in kw_pairs),
    }


def main():
    print(f"Reading {PROTEINS.name}...")
    prot = pd.read_csv(PROTEINS, low_memory=False)
    accs = prot["accession"].tolist()
    print(f"  {len(accs):,} accessions to enrich")

    # Load cache
    cache = {}
    if CACHE.exists():
        try:
            cache = json.loads(CACHE.read_text())
            print(f"  Loaded {len(cache):,} cached entries from {CACHE.name}")
        except Exception:
            pass

    remaining = [a for a in accs if a not in cache]
    print(f"  {len(remaining):,} entries to fetch")

    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i:i + BATCH_SIZE]
        bnum = i // BATCH_SIZE + 1
        btot = (len(remaining) + BATCH_SIZE - 1) // BATCH_SIZE
        print(f"  batch {bnum}/{btot}  ({len(batch)})", flush=True)
        entries = fetch_bulk(batch)
        for entry in entries:
            ann = extract(entry)
            if ann["accession"]:
                cache[ann["accession"]] = ann
        # checkpoint every 10 batches
        if bnum % 10 == 0:
            CACHE.write_text(json.dumps(cache))
            print(f"    checkpoint: {len(cache):,} cached", flush=True)
        time.sleep(0.2)

    # Retry any still-missing via direct endpoint
    missing = [a for a in accs if a not in cache]
    if missing:
        print(f"\nRetrying {len(missing)} missing via direct endpoint...")
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
            time.sleep(0.1)

    CACHE.write_text(json.dumps(cache))
    print(f"\nCache: {len(cache):,} entries")

    # Build the enrichment table
    tbl_rows = []
    for acc in accs:
        ann = cache.get(acc)
        if ann is None:
            tbl_rows.append({"accession": acc, "go_ids": "", "go_labels": "",
                             "keyword_ids": "", "keyword_labels": ""})
        else:
            tbl_rows.append(ann)
    tbl = pd.DataFrame(tbl_rows)
    tbl.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUT_TSV.name}")

    # Merge into proteins.csv (drop existing columns first to avoid duplicates)
    for c in ["go_ids", "go_labels", "keyword_ids", "keyword_labels"]:
        if c in prot.columns:
            prot = prot.drop(columns=[c])
    prot = prot.merge(tbl, on="accession", how="left")
    for c in ["go_ids", "go_labels", "keyword_ids", "keyword_labels"]:
        prot[c] = prot[c].fillna("").astype(str)
    prot.to_csv(PROTEINS, index=False)
    print(f"Updated {PROTEINS.name}: {len(prot):,} rows × {len(prot.columns)} columns")

    # Coverage report
    with_go = int((prot["go_ids"] != "").sum())
    with_kw = int((prot["keyword_ids"] != "").sum())
    print(f"  entries with any GO term:   {with_go:,}  ({with_go/len(prot)*100:.1f}%)")
    print(f"  entries with any keyword:   {with_kw:,}  ({with_kw/len(prot)*100:.1f}%)")


if __name__ == "__main__":
    main()
