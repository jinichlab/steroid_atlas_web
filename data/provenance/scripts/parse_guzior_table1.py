"""Parse Guzior 2024 Supplementary Table 1 into a fetch-ready reference JSON.

Input:  literature/supplementary/Guzior2024_SupplTables.xlsx (sheet 'Table 1')
Output: literature/guzior2024_table1_normalized.json

Structure of Table 1:
  Each species-row has metadata (Species/Source/Strain/Family/Cluster/Genome)
  and a first BSH/T accession. Paralog rows follow with only an accession +
  BSH/T Group filled in.

We forward-fill species/strain/etc. across paralog rows so every accession has
its full context.

Output JSON is a list of records, each:
  {
    "species": ..., "source": ..., "strain": ...,
    "family": ..., "mcba_profile_cluster": ..., "genome_accession": ...,
    "protein_accession_raw": ...,
    "accession_type": "RefSeq_WP" | "GenBank" | "UniProt" | "IMG_numeric",
    "accession_normalized": ...,   # version-stripped, IMG# prefix removed
    "bsh_t_group": "I" | "II" | "III",
    "fetchable_public": True/False,
    "fetch_method": "ncbi_efetch_protein" | "uniprot_rest" | "jgi_img_deferred"
  }
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
XLSX = HERE / "supplementary" / "Guzior2024_SupplTables.xlsx"
OUT = HERE / "guzior2024_table1_normalized.json"

# UniProt accession pattern (Swiss-Prot + TrEMBL formats)
UNIPROT_RX = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
UNIPROT_ENTRY_NAME_RX = re.compile(r"^[A-Z0-9]+_[A-Z0-9]+$")
REFSEQ_WP_RX = re.compile(r"^WP_\d+(?:\.\d+)?$")
GENBANK_PROT_RX = re.compile(r"^[A-Z]{2,4}\d{5,}(?:\.\d+)?$")

def classify(acc: str) -> tuple[str, str, str]:
    """Returns (accession_type, accession_normalized, fetch_method)."""
    a = str(acc).strip()
    if a.startswith("IMG#"):
        return ("IMG_numeric", a.replace("IMG#", ""), "jgi_img_deferred")
    stripped = a.split(".")[0]
    if REFSEQ_WP_RX.match(a):
        return ("RefSeq_WP", stripped, "ncbi_efetch_protein")
    # UniProt entry name like D7VDZ4_LACPN
    if "_" in a and UNIPROT_ENTRY_NAME_RX.match(a):
        acc_only = a.split("_")[0]
        if UNIPROT_RX.match(acc_only):
            return ("UniProt", acc_only, "uniprot_rest")
    # Bare UniProt accession
    if UNIPROT_RX.match(a):
        return ("UniProt", a, "uniprot_rest")
    # GenBank-style protein ID (letters + digits + optional version)
    if GENBANK_PROT_RX.match(a):
        return ("GenBank", stripped, "ncbi_efetch_protein")
    return ("unknown", a, "unknown")

def main() -> int:
    df = pd.read_excel(XLSX, sheet_name="Table 1")

    # Forward-fill species-level metadata across paralog rows (NaN-species rows)
    ffill_cols = ["Species", "Source", "Strain", "Taxonomic Family",
                  "MCBA Profile Cluster", "Genome Accession#"]
    for c in ffill_cols:
        df[c] = df[c].ffill()

    records = []
    n_no_protein = 0
    for _, row in df.iterrows():
        prot = row.get("BSH/T Protein Accession#")
        if pd.isna(prot) or not str(prot).strip():
            n_no_protein += 1
            continue
        prot = str(prot).strip()
        acc_type, acc_norm, fetch_method = classify(prot)

        cluster = row["MCBA Profile Cluster"]
        # normalize cluster field
        if pd.isna(cluster):
            cluster_norm = "control"
        else:
            cluster_norm = str(cluster).strip()

        rec = {
            "species": str(row["Species"]).strip() if pd.notna(row["Species"]) else "",
            "source": str(row["Source"]).strip() if pd.notna(row["Source"]) else "",
            "strain": str(row["Strain"]).strip() if pd.notna(row["Strain"]) else "",
            "family": str(row["Taxonomic Family"]).strip() if pd.notna(row["Taxonomic Family"]) else "",
            "mcba_profile_cluster": cluster_norm,
            "genome_accession": str(row["Genome Accession#"]).strip() if pd.notna(row["Genome Accession#"]) else "",
            "protein_accession_raw": prot,
            "accession_type": acc_type,
            "accession_normalized": acc_norm,
            "bsh_t_group": str(row["BSH/T Group"]).strip() if pd.notna(row["BSH/T Group"]) else "",
            "fetchable_public": fetch_method != "jgi_img_deferred" and acc_type != "unknown",
            "fetch_method": fetch_method,
        }
        records.append(rec)

    OUT.write_text(json.dumps(records, indent=2))
    print(f"Wrote {OUT.name}  ({len(records)} enzyme rows)")
    print(f"  species-rows without any BSH/T accession: {n_no_protein}")

    # Summary
    print("\n--- by accession type ---")
    from collections import Counter
    types = Counter(r["accession_type"] for r in records)
    for t, n in types.most_common():
        print(f"  {t:15s} {n:3d}")

    print("\n--- by fetch method ---")
    methods = Counter(r["fetch_method"] for r in records)
    for m, n in methods.most_common():
        print(f"  {m:25s} {n:3d}")

    print("\n--- by MCBA profile cluster ---")
    clusters = Counter(r["mcba_profile_cluster"] for r in records)
    for c, n in sorted(clusters.items()):
        print(f"  cluster={c!r:6s}  {n:3d} enzymes")

    print("\n--- fetchable via public APIs (Tier 1-4) ---")
    fetchable = [r for r in records if r["fetchable_public"]]
    print(f"  count: {len(fetchable)}")
    for r in fetchable:
        print(f"    {r['accession_normalized']:15s} [{r['accession_type']:9s}] "
              f"cluster={r['mcba_profile_cluster']:>3s} grp={r['bsh_t_group']:>3s}  "
              f"{r['species']} {r['strain']}")

    print("\n--- deferred (JGI IMG numeric IDs) ---")
    deferred = [r for r in records if not r["fetchable_public"]]
    print(f"  count: {len(deferred)}")
    for r in deferred:
        print(f"    {r['accession_normalized']:12s}  cluster={r['mcba_profile_cluster']:>3s} grp={r['bsh_t_group']:>3s}  "
              f"{r['species']} {r['strain']}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
