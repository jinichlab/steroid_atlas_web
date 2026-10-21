"""Fetch UniProt lineage for every atlas entry and produce a taxonomy breakdown.

Every one of the 35,117 accessions is queried via the UniProt REST bulk
endpoint for its taxonomy lineage. Each protein is then categorized into
a broad lineage group (Mammals, Fish, Amphibians, Reptiles, Birds,
Invertebrates, Fungi, Plants, Bacteria, Archaea, Protists, Viruses, Other)
via keyword rules over the full lineage string. Three complementary
breakdowns are produced:

  A. Superkingdom  — Eukaryota / Bacteria / Archaea / Viruses
  B. Broad lineage — the categorized groups above
  C. Top-15 species — most represented individual organisms

Reads:  ../data/proteins.csv
Writes: ../analysis/taxonomy_cache.json          raw UniProt lineage cache
        ../analysis/taxonomy_per_entry.tsv       per-entry lineage table
        ../analysis/taxonomy_summary.txt         human-readable summary
        ../analysis/taxonomy_figure.png / .pdf   3-panel breakdown figure
"""
from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
CACHE = output("taxonomy_cache.json", historical=True)
OUT_TSV = output("taxonomy_per_entry.tsv", historical=True)
OUT_TXT = output("taxonomy_summary.txt", historical=True)
OUT_PNG = output("taxonomy_figure.png", historical=True)
OUT_PDF = output("taxonomy_figure.pdf", historical=True)

BULK_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = "accession,organism_name,organism_id,lineage,lineage_ids"
BATCH_SIZE = 300

# Ordered rules: the FIRST rule whose keyword appears in the lineage wins.
# Order matters — put more specific groups above broader ones so e.g. Mammals
# is checked before generic "Metazoa".
LINEAGE_RULES = [
    ("Mammals",       ["Mammalia"]),
    ("Fish",          ["Actinopterygii", "Actinopteri", "Sarcopterygii", "Chondrichthyes",
                       "Cyclostomata", "Agnatha", "Petromyzontidae", "Teleostei"]),
    ("Birds",         ["Aves"]),
    ("Reptiles",      ["Reptilia", "Squamata", "Testudines", "Crocodylia", "Lepidosauria"]),
    ("Amphibians",    ["Amphibia"]),
    ("Invertebrates", ["Arthropoda", "Nematoda", "Mollusca", "Annelida", "Echinodermata",
                       "Cnidaria", "Porifera", "Platyhelminthes", "Ecdysozoa", "Lophotrochozoa"]),
    ("Fungi",         ["Fungi"]),
    ("Plants",        ["Viridiplantae", "Streptophyta", "Chlorophyta"]),
    ("Protists",      ["Alveolata", "Stramenopiles", "Rhodophyta", "Euglenozoa", "Amoebozoa",
                       "Choanoflagellida", "Rhizaria", "Discoba"]),
    ("Bacteria",      ["Bacteria"]),
    ("Archaea",       ["Archaea"]),
    ("Viruses",       ["Viruses"]),
]


def fetch_bulk(accessions, retries=3):
    params = {"accessions": ",".join(accessions), "format": "json", "fields": FIELDS}
    for attempt in range(retries):
        try:
            r = requests.get(BULK_URL, params=params, timeout=60)
            if r.status_code == 200:
                return r.json().get("results", [])
            time.sleep(2 ** attempt)
        except requests.RequestException:
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
    org = entry.get("organism", {}) or {}
    lineage_names = org.get("lineage", []) or []
    return {
        "accession": entry.get("primaryAccession", ""),
        "organism_name": org.get("scientificName", ""),
        "organism_id": org.get("taxonId", ""),
        "lineage": "; ".join(lineage_names),
    }


def broad_group(lineage: str) -> str:
    lin = lineage or ""
    for group, keywords in LINEAGE_RULES:
        for kw in keywords:
            if kw in lin:
                return group
    return "Other / unclassified"


def superkingdom(lineage: str) -> str:
    lin = (lineage or "").lower()
    for k in ("bacteria", "archaea", "viruses", "eukaryota"):
        if k in lin:
            return k.capitalize()
    return "Unknown"


def do_fetch(accessions):
    cache = {}
    if CACHE.exists():
        try:
            cache = json.loads(CACHE.read_text())
            print(f"  loaded {len(cache):,} cached entries")
        except Exception:
            pass
    remaining = [a for a in accessions if a not in cache]
    print(f"  {len(remaining):,} to fetch (batch size {BATCH_SIZE})")
    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i:i + BATCH_SIZE]
        entries = fetch_bulk(batch)
        for e in entries:
            ann = extract(e)
            if ann["accession"]:
                cache[ann["accession"]] = ann
        if (i // BATCH_SIZE) % 10 == 9:
            CACHE.write_text(json.dumps(cache))
            print(f"    checkpoint: {len(cache):,} cached", flush=True)
        time.sleep(0.15)
    # Direct-entry fallback for what bulk skipped (usually TrEMBL)
    missing = [a for a in accessions if a not in cache]
    if missing:
        print(f"  Retrying {len(missing)} via direct endpoint...")
        for j, acc in enumerate(missing):
            if j % 50 == 0:
                print(f"    {j}/{len(missing)}", flush=True)
            entry = fetch_single(acc)
            if entry is not None:
                ann = extract(entry)
                if ann["accession"]:
                    cache[ann["accession"]] = ann
                if ann["accession"] != acc:
                    cache[acc] = ann
            time.sleep(0.08)
    CACHE.write_text(json.dumps(cache))
    print(f"  final cache: {len(cache):,} / {len(accessions):,}")
    return cache


def make_figure(df: pd.DataFrame, n_total: int):
    BAR = "#1F4B99"
    fig = plt.figure(figsize=(13.0, 4.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.2, 1.4])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])

    # ── Panel A — superkingdom breakdown (bar) ─────────────────────────
    sk = df["superkingdom"].value_counts()
    ypos = np.arange(len(sk))[::-1]
    axA.barh(ypos, sk.values, color=BAR, height=0.58,
             edgecolor="white", linewidth=0.4)
    for y, (name, count) in zip(ypos, sk.items()):
        pct = 100 * count / n_total
        axA.text(count * 1.03, y, f" {count:,} ({pct:.1f}%)",
                 va="center", ha="left", fontsize=9, color="#111827")
    axA.set_yticks(ypos)
    axA.set_yticklabels(sk.index, fontsize=9)
    axA.tick_params(axis="y", length=0, pad=4)
    axA.tick_params(axis="x", labelsize=8)
    axA.set_xlim(0, sk.max() * 1.35)
    axA.set_xlabel("Proteins", fontsize=9)
    axA.set_title("A. Superkingdom", fontsize=11, loc="left", pad=8)
    axA.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    for s in ("top", "right", "left"):
        axA.spines[s].set_visible(False)
    axA.grid(True, axis="x", alpha=0.25, linewidth=0.5)

    # ── Panel B — broad lineage groups ─────────────────────────────────
    order = [g for g, _ in LINEAGE_RULES] + ["Other / unclassified"]
    grp = df["broad_group"].value_counts().reindex(order).dropna()
    grp = grp[grp > 0].sort_values(ascending=True)
    ypos = np.arange(len(grp))
    axB.barh(ypos, grp.values, color=BAR, height=0.7,
             edgecolor="white", linewidth=0.4)
    for y, (name, count) in zip(ypos, grp.items()):
        pct = 100 * count / n_total
        axB.text(count * 1.03, y, f" {count:,} ({pct:.1f}%)",
                 va="center", ha="left", fontsize=8.5, color="#111827")
    axB.set_yticks(ypos)
    axB.set_yticklabels(grp.index, fontsize=9)
    axB.tick_params(axis="y", length=0, pad=4)
    axB.tick_params(axis="x", labelsize=8)
    axB.set_xlim(0, grp.max() * 1.4)
    axB.set_xlabel("Proteins", fontsize=9)
    axB.set_title("B. Broad lineage groups", fontsize=11, loc="left", pad=8)
    axB.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    for s in ("top", "right", "left"):
        axB.spines[s].set_visible(False)
    axB.grid(True, axis="x", alpha=0.25, linewidth=0.5)

    # ── Panel C — top-15 species ───────────────────────────────────────
    top_sp = df["organism_name"].value_counts().head(15)
    ypos = np.arange(len(top_sp))[::-1]
    axC.barh(ypos, top_sp.values, color=BAR, height=0.7,
             edgecolor="white", linewidth=0.4)
    for y, (name, count) in zip(ypos, top_sp.items()):
        pct = 100 * count / n_total
        axC.text(count * 1.02, y, f" {count:,}",
                 va="center", ha="left", fontsize=8.5, color="#111827")
    axC.set_yticks(ypos)
    axC.set_yticklabels(top_sp.index, fontsize=8.5)
    axC.tick_params(axis="y", length=0, pad=4)
    axC.tick_params(axis="x", labelsize=8)
    axC.set_xlim(0, top_sp.max() * 1.25)
    axC.set_xlabel("Proteins", fontsize=9)
    axC.set_title("C. Top 15 species", fontsize=11, loc="left", pad=8)
    axC.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    for s in ("top", "right", "left"):
        axC.spines[s].set_visible(False)
    axC.grid(True, axis="x", alpha=0.25, linewidth=0.5)

    fig.suptitle(f"Taxonomic composition of the atlas  (n = {n_total:,} proteins)",
                 fontsize=11.5, x=0.01, ha="left", y=1.04)

    for f in (OUT_PNG, OUT_PDF):
        fig.savefig(f, dpi=300, bbox_inches="tight")
        print(f"Wrote {f.name}")


def write_summary(df: pd.DataFrame, n_total: int):
    lines = []
    lines.append("=" * 72)
    lines.append(f"TAXONOMIC COMPOSITION — all {n_total:,} atlas proteins")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Coverage:")
    n_with_lineage = int((df["lineage"] != "").sum())
    lines.append(f"  proteins with UniProt lineage: {n_with_lineage:,} / {n_total:,} "
                 f"({n_with_lineage/n_total*100:.1f}%)")
    lines.append("")
    lines.append("[A] Superkingdom:")
    for k, v in df["superkingdom"].value_counts().items():
        lines.append(f"    {v:>7,}  ({v/n_total*100:5.1f}%)  {k}")
    lines.append("")
    lines.append("[B] Broad lineage group:")
    for k, v in df["broad_group"].value_counts().items():
        lines.append(f"    {v:>7,}  ({v/n_total*100:5.1f}%)  {k}")
    lines.append("")
    lines.append(f"[C] Distinct species: {df['organism_name'].nunique():,}")
    lines.append("    Top 20 species:")
    for k, v in df["organism_name"].value_counts().head(20).items():
        lines.append(f"      {v:>5}  {k}")
    OUT_TXT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_TXT.name}")


def main() -> int:
    print(f"Reading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    accessions = p["accession"].astype(str).tolist()
    n_total = len(accessions)
    print(f"  {n_total:,} accessions")

    print("\nFetching UniProt lineage...")
    cache = do_fetch(accessions)

    rows = []
    for acc in accessions:
        ann = cache.get(acc)
        if ann is None:
            rows.append({"accession": acc, "organism_name": "", "organism_id": "",
                         "lineage": "", "superkingdom": "Unknown",
                         "broad_group": "Other / unclassified"})
        else:
            rows.append({**ann,
                         "superkingdom": superkingdom(ann["lineage"]),
                         "broad_group": broad_group(ann["lineage"])})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUT_TSV.name}")

    write_summary(df, n_total)
    make_figure(df, n_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
