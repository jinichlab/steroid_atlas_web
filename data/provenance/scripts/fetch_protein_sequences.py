"""Fetch FASTA for new_proteins_2024_2026.json with explicit verification.

Tiers (skips Tier 4 — those need a human decision before any fetch):

  Tier 1  UniProt accession set                -> rest.uniprot.org
  Tier 2  GenBank/RefSeq protein ID set        -> NCBI efetch (db=protein)
  Tier 3  Locus tag only                       -> NCBI esearch -> efetch
  Tier 4  Community-ortholog guess / no ID     -> SKIPPED, flagged

Writes:
  sequences/<safe_id>.fasta   one file per fetched protein
  sequences/all_new.fasta     combined
  sequence_fetch_report.tsv   per-row verification table

Run:
  LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \
  /home/adsiordia/miniconda3/bin/python \
  /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/literature/fetch_protein_sequences.py
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "new_proteins_2024_2026.json"
OUT_DIR = HERE / "sequences"
OUT_DIR.mkdir(exist_ok=True)
COMBINED = OUT_DIR / "all_new.fasta"
REPORT = HERE / "sequence_fetch_report.tsv"

UA = "literature-extraction (adsiordia@ucsd.edu)"
NCBI_DELAY = 0.4  # be polite to E-utilities

FUNCTION_KEYWORDS = {
    "bsh": ["hydrolase", "choloyl", "bile salt"],
    "bsh/T": ["hydrolase", "choloyl", "bile salt", "acyltransferase"],
    "OsrA": ["reductase", "oxidoreductase", "ene", "OYE"],
    "OsrB": ["reductase", "oxidoreductase", "ene", "OYE"],
    "OsrC": ["reductase", "oxidoreductase", "ene", "OYE", "dehydrogenase"],
    "ci2350": ["reductase", "oxidoreductase", "OYE", "ene"],
    "ci2349": ["reductase", "baiH", "oxidoreductase"],
    "Elen_2451": ["fdhD", "molybdo", "sulfurtransferase"],
    "Elen_2452": ["ferredoxin", "4Fe-4S"],
    "Elen_2453": ["molybdopterin", "oxidoreductase", "dehydrogenase"],
    "Elen_2454": ["SPFH", "band 7", "band-7", "stomatin", "prohibitin"],
}

def http_get(url: str, timeout: int = 45) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")

def fetch_uniprot(acc: str) -> str:
    return http_get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")

def fetch_ncbi_protein(acc: str) -> str:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={urllib.parse.quote(acc)}&rettype=fasta&retmode=text"
    )
    time.sleep(NCBI_DELAY)
    return http_get(url)

def ncbi_esearch_protein(query: str) -> list[str]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=protein&term={urllib.parse.quote(query)}&retmode=json"
    )
    time.sleep(NCBI_DELAY)
    payload = json.loads(http_get(url))
    return payload.get("esearchresult", {}).get("idlist", [])

def parse_fasta(text: str) -> tuple[str, str]:
    text = text.strip()
    if not text.startswith(">"):
        return "", ""
    lines = text.split("\n")
    header = lines[0].lstrip(">").strip()
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    seq = re.sub(r"[^A-Za-z*]", "", seq)
    return header, seq

def genus_match(our_org: str, header: str) -> bool:
    if not our_org or not header:
        return False
    genus = our_org.split()[0]
    return genus.lower() in header.lower()

def function_match(gene_name: str, header: str) -> bool:
    keys = FUNCTION_KEYWORDS.get(gene_name.split()[0], [])
    h = header.lower()
    return any(k.lower() in h for k in keys) if keys else False

def safe_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")

def classify(entry: dict) -> tuple[str, str | None]:
    """Returns (tier, identifier_to_use_or_None)."""
    uacc = (entry.get("uniprot_accession") or "").strip()
    ncbi = (entry.get("ncbi_or_refseq_locus") or "").strip()
    gene = (entry.get("gene_name") or "").strip()

    # Tier 4: explicit guesses / missing
    # Q5DWB7 was flagged as "community ortholog" in JSON.
    if "community ortholog" in uacc.lower():
        return "T4_guess", None
    if not uacc and not ncbi:
        return "T4_no_id", None
    # The 3beta-HSDH row has no specific locus
    if "PF01073" in entry.get("organism", "") or "(family-level" in ncbi:
        return "T4_family_only", None
    # Rimal BlBSH: UniProt empty, NCBI field is a genome ref + community ortholog
    if "community ortholog" in ncbi.lower():
        return "T4_guess", None

    # Tier 1: clean UniProt accession (covers both OPQ-prefix and other-prefix formats)
    UNIPROT_RX = re.compile(
        r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9])"      # e.g. Q5LF84, P12345
        r"|(?:[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"  # e.g. C8WL28, A0A024R161
    )
    if UNIPROT_RX.fullmatch(uacc):
        return "T1_uniprot", uacc

    # Tier 2: protein accession in NCBI field (RefSeq WP_/NP_, GenBank MFU/ACV/AAA-style with optional version)
    m = re.search(r"\b((?:WP_|NP_|YP_|XP_)?[A-Z]{2,}_?\d+(?:\.\d+)?)\b", ncbi)
    if m:
        return "T2_ncbi", m.group(1)

    # Tier 3: pure locus tag (e.g. ci2349, BF9343_1433)
    m = re.search(r"locus tag\s+([A-Za-z0-9_]+)", ncbi)
    if m:
        return "T3_locus", m.group(1)

    return "T4_unrecognized", None

def main() -> int:
    rows = json.loads(SRC.read_text())
    report_rows: list[dict] = []
    combined_fasta_parts: list[str] = []

    for entry in rows:
        gene = entry.get("gene_name", "")
        org = entry.get("organism", "")
        tier, ident = classify(entry)

        rec: dict = {
            "gene_name": gene,
            "organism": org,
            "tier": tier,
            "identifier_used": ident or "",
            "fetched_header": "",
            "length": "",
            "genus_match": "",
            "function_keyword_match": "",
            "status": "",
            "notes": "",
        }

        if tier.startswith("T4"):
            rec["status"] = "SKIPPED"
            rec["notes"] = {
                "T4_guess": "community-ortholog guess — needs user decision",
                "T4_no_id": "no identifier in source paper",
                "T4_family_only": "paper only specifies the protein family, not a locus",
                "T4_unrecognized": "identifier field present but pattern unrecognized",
            }[tier]
            report_rows.append(rec)
            continue

        try:
            if tier == "T1_uniprot":
                fasta = fetch_uniprot(ident)
            elif tier == "T2_ncbi":
                fasta = fetch_ncbi_protein(ident)
            elif tier == "T3_locus":
                ids = ncbi_esearch_protein(f"{ident}[Gene Name] AND \"{org.split('(')[0].strip()}\"[Organism]")
                if not ids:
                    ids = ncbi_esearch_protein(f"{ident} AND \"{org.split('(')[0].strip()}\"[Organism]")
                if not ids:
                    rec["status"] = "FAILED"
                    rec["notes"] = "esearch returned no hits"
                    report_rows.append(rec)
                    continue
                if len(ids) > 1:
                    rec["notes"] = f"esearch returned {len(ids)} hits, took first ({ids[0]})"
                rec["identifier_used"] = f"{ident} -> NCBI:{ids[0]}"
                fasta = fetch_ncbi_protein(ids[0])
            else:
                rec["status"] = "FAILED"
                rec["notes"] = f"unrecognized tier {tier}"
                report_rows.append(rec)
                continue

            header, seq = parse_fasta(fasta)
            if not seq:
                rec["status"] = "FAILED"
                rec["notes"] = "fetched body was not FASTA / empty"
                report_rows.append(rec)
                continue

            rec["fetched_header"] = header[:140]
            rec["length"] = str(len(seq))
            rec["genus_match"] = "yes" if genus_match(org, header) else "NO"
            rec["function_keyword_match"] = "yes" if function_match(gene, header) else "no"
            rec["status"] = "OK"

            out_name = safe_id(f"{gene}_{ident}") + ".fasta"
            (OUT_DIR / out_name).write_text(fasta if fasta.endswith("\n") else fasta + "\n")
            combined_fasta_parts.append(fasta.strip() + "\n")

        except Exception as e:
            rec["status"] = "ERROR"
            rec["notes"] = str(e)[:160]

        report_rows.append(rec)

    COMBINED.write_text("\n".join(combined_fasta_parts) + "\n")

    cols = ["gene_name", "organism", "tier", "identifier_used", "fetched_header",
            "length", "genus_match", "function_keyword_match", "status", "notes"]
    with REPORT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in report_rows:
            fh.write("\t".join(r.get(c, "") for c in cols) + "\n")

    print(f"\nReport: {REPORT}")
    print(f"FASTA dir: {OUT_DIR}")
    print(f"Combined: {COMBINED}")

    print("\n--- summary ---")
    ok = sum(1 for r in report_rows if r["status"] == "OK")
    skipped = sum(1 for r in report_rows if r["status"] == "SKIPPED")
    failed = sum(1 for r in report_rows if r["status"] in ("FAILED", "ERROR"))
    print(f"OK: {ok}   SKIPPED: {skipped}   FAILED/ERROR: {failed}")
    print("\n--- per-row ---")
    for r in report_rows:
        print(
            f"  [{r['status']:7s}] {r['gene_name']:40s} "
            f"len={r['length']:>5s} genus={r['genus_match']:>3s} func={r['function_keyword_match']:>3s}  "
            f"-> {r['identifier_used']}"
        )
        if r["notes"]:
            print(f"           note: {r['notes']}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
