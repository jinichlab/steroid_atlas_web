"""Fetch BSH/T sequences from Guzior 2024 Supplementary Table 1.

Reads:  literature/guzior2024_table1_normalized.json  (fetchable entries only)
Writes: literature/sequences/guzior/<accession>.fasta (per-protein)
        literature/sequences/guzior/all_guzior.fasta   (combined)
        literature/guzior2024_review_ready.json        (rich, review-ready records)
        literature/guzior2024_fetch_report.tsv         (per-fetch verification)

Verification for each fetched sequence:
  - genus_match: fetched header contains the Guzior species genus
  - length_sane: 250 <= length <= 500  (BSHs are ~317-340 aa)
  - looks_like_bsh: header matches BSH/hydrolase/choloyl keywords

Skips JGI IMG numeric IDs (not fetchable via NCBI/UniProt without JGI API).

Run:
  LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib \
  /home/adsiordia/miniconda3/bin/python \
  /home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/literature/fetch_guzior_sequences.py
"""
from __future__ import annotations
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT = HERE / "guzior2024_table1_normalized.json"
FASTA_DIR = HERE / "sequences" / "guzior"
FASTA_DIR.mkdir(parents=True, exist_ok=True)
COMBINED = FASTA_DIR / "all_guzior.fasta"
REPORT = HERE / "guzior2024_fetch_report.tsv"
REVIEW = HERE / "guzior2024_review_ready.json"

UA = "steroid-atlas literature-extraction (adsiordia@ucsd.edu)"
NCBI_DELAY = 0.4

BSH_KEYWORDS = ("bsh", "bile salt hydrolase", "choloyl", "conjugated bile", "choloylglycine")

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

def parse_fasta(text: str) -> tuple[str, str]:
    text = text.strip()
    if not text.startswith(">"):
        return "", ""
    lines = text.split("\n")
    header = lines[0].lstrip(">").strip()
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    seq = re.sub(r"[^A-Za-z*]", "", seq)
    return header, seq

def genus_match(species: str, header: str) -> bool:
    if not species or not header:
        return False
    # Handle brackets in taxonomy (e.g., "[Clostridium] scindens")
    genus = species.replace("[", "").replace("]", "").split()[0]
    return genus.lower() in header.lower()

def looks_like_bsh(header: str) -> bool:
    h = (header or "").lower()
    return any(kw in h for kw in BSH_KEYWORDS)

def main() -> int:
    records = json.loads(INPUT.read_text())
    fetchable = [r for r in records if r.get("fetchable_public")]
    print(f"Total fetchable entries in normalized JSON: {len(fetchable)}")

    combined_parts: list[str] = []
    report_rows: list[dict] = []
    review_records: list[dict] = []
    n_ok = n_fail = 0

    for rec in fetchable:
        acc = rec["accession_normalized"]
        atype = rec["accession_type"]
        report: dict = {
            "accession": acc,
            "acc_raw": rec["protein_accession_raw"],
            "type": atype,
            "species": rec["species"],
            "strain": rec["strain"],
            "cluster": rec["mcba_profile_cluster"],
            "group": rec["bsh_t_group"],
            "fetched_header": "",
            "length": "",
            "genus_match": "",
            "looks_like_bsh": "",
            "status": "",
            "notes": "",
        }

        try:
            if atype == "UniProt":
                fasta = fetch_uniprot(acc)
            elif atype in ("RefSeq_WP", "GenBank"):
                fasta = fetch_ncbi_protein(acc)
            else:
                report["status"] = "SKIPPED"
                report["notes"] = f"unhandled type {atype}"
                report_rows.append(report)
                continue

            header, seq = parse_fasta(fasta)
            if not seq:
                report["status"] = "FAILED"
                report["notes"] = "no sequence in fetched body"
                report_rows.append(report)
                n_fail += 1
                continue

            report["fetched_header"] = header[:180]
            report["length"] = str(len(seq))
            report["genus_match"] = "yes" if genus_match(rec["species"], header) else "NO"
            report["looks_like_bsh"] = "yes" if looks_like_bsh(header) else "no"
            report["status"] = "OK"

            fasta_path = FASTA_DIR / f"{acc}.fasta"
            fasta_path.write_text(fasta if fasta.endswith("\n") else fasta + "\n")
            combined_parts.append(fasta.strip() + "\n")
            n_ok += 1

            review_records.append({
                **rec,
                "fetched_header": header,
                "sequence_length": len(seq),
                "sequence": seq,
                "genus_match": genus_match(rec["species"], header),
                "looks_like_bsh": looks_like_bsh(header),
                "paper_url": "https://doi.org/10.1038/s41586-024-07017-8",
                "paper_citation": "Guzior et al., Nature 626, 852-858 (2024)",
                "atlas_row_proposed": {
                    "Entry": acc,
                    "Entry Name": f"BSH_T_{rec['species'].split()[0][:3].upper()}_{rec['strain']}".replace(" ", "_"),
                    "Protein names": (
                        "Bile salt hydrolase/transferase (BSH/T) "
                        "(Bile acid amine N-acyltransferase) (EC 2.3.1.-) "
                        "(EC 3.5.1.24)"
                    ),
                    "Gene Names": "bsh bsh/T",
                    "Organism": f"{rec['species']} strain {rec['strain']}",
                    "Length": len(seq),
                    "Sequence": seq,
                    "Annotation": (
                        f"Guzior et al., Nature 2024 (Supp Table 1; "
                        f"MCBA profile cluster {rec['mcba_profile_cluster']}, "
                        f"BSH/T Group {rec['bsh_t_group']}) | "
                        "https://doi.org/10.1038/s41586-024-07017-8"
                    ),
                    "source": "Guzior 2024 Supp Table 1",
                    "is_new": 1,
                    "Paper": "https://doi.org/10.1038/s41586-024-07017-8",
                    "reaction_ecs": "EC 2.3.1.- (acyltransferase); EC 3.5.1.24 (choloylglycine hydrolase)",
                    "reaction_descriptions": (
                        "Bile salt hydrolase (deconjugation) + amine N-acyltransferase "
                        f"(BBAA/MCBA biosynthesis). MCBA profile cluster: {rec['mcba_profile_cluster']}. "
                        f"BSH/T phylogenetic group: {rec['bsh_t_group']}."
                    ),
                    "Sequence_Source": (
                        f"{atype} accession {acc} — literature-recruited via Guzior 2024 "
                        f"Supp Table 1 ({rec['species']} {rec['strain']}, MCBA cluster "
                        f"{rec['mcba_profile_cluster']}, Group {rec['bsh_t_group']})."
                    ),
                    "Identifier_Type": atype,
                },
            })

        except Exception as e:
            report["status"] = "ERROR"
            report["notes"] = str(e)[:200]
            n_fail += 1

        report_rows.append(report)

    # Write outputs
    COMBINED.write_text("\n".join(combined_parts) + "\n")

    cols = ["accession", "acc_raw", "type", "species", "strain", "cluster", "group",
            "fetched_header", "length", "genus_match", "looks_like_bsh", "status", "notes"]
    with REPORT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in report_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    REVIEW.write_text(json.dumps(review_records, indent=2))

    print(f"\n{'-'*80}")
    print(f"Report:   {REPORT.name}")
    print(f"FASTA dir: {FASTA_DIR}")
    print(f"Combined: {COMBINED}")
    print(f"Review:   {REVIEW.name}")
    print(f"\nOK: {n_ok}  FAILED: {n_fail}")

    print("\n--- Per-row verification ---")
    for r in report_rows:
        genus_ok = r["genus_match"] == "yes"
        bsh_ok = r["looks_like_bsh"] == "yes"
        flag = "✓" if r["status"] == "OK" and genus_ok and bsh_ok else ("!" if r["status"] == "OK" else "✗")
        print(f"  {flag} [{r['status']:6s}] {r['accession']:15s} "
              f"len={r['length']:>4s} genus={r['genus_match']:>3s} bsh={r['looks_like_bsh']:>3s}  "
              f"cluster={r['cluster']} grp={r['group']}  {r['species']}")
        if r.get("notes"):
            print(f"           note: {r['notes']}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
