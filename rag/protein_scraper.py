#!/usr/bin/env python3
"""
protein_scraper.py

Scrapes UniProt "publications" data for a list of UniProt accessions and
builds a per-protein RAG training corpus (JSON + TXT) from the cited
literature (title, authors, journal, abstract, and full text when the
paper is open access via Europe PMC).

Usage:
    python protein_scraper.py \
        --csv data/proteins.csv \
        --out rag_data/rag_train/proteins_json \
        --email you@example.com

Notes:
    - --email is sent as a courtesy identifier in the User-Agent header to
      UniProt/Europe PMC (their guidelines ask for this; not authentication).
    - Re-running is safe: proteins with an existing {accession}.json in the
      output folder are skipped unless --overwrite is passed.
    - Publications come from rest.uniprot.org (the API behind the
      uniprot.org/uniprotkb/{acc}/publications page). Accessions that are
      genuinely unknown/invalid return HTTP 404/400 and are logged to
      {out}/_inactive.txt; an entry with simply zero publications is written
      out normally with an empty list.
"""

import argparse
import csv
import json
import logging
import re
import sys
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm

UNIPROT_ACCESSION_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)

PUBLICATIONS_URL = "https://rest.uniprot.org/uniprotkb/{accession}/publications"
EPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EPMC_FULLTEXT_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{source}/{ext_id}/fullTextXML"

_ROOT = Path(__file__).resolve().parents[1]  # repo root (rag/ → ..)
DEFAULT_CSV = _ROOT / "data" / "proteins.csv"
DEFAULT_OUT = _ROOT / "rag_data" / "rag_train" / "proteins_json"

CANDIDATE_ID_COLUMNS = [
    "uniprot_id", "uniprot", "uniprotkb", "uniprot_accession",
    "accession", "entry", "uniprotkb_id", "protein_id", "id",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("protein_scraper")


def make_session(email: str | None) -> requests.Session:
    session = requests.Session()
    ua = "steroid-atlas-protein-scraper/1.0"
    if email:
        ua += f" ({email})"
    session.headers.update({"Accept": "application/json", "User-Agent": ua})
    retry = Retry(
        total=5,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def detect_id_column(fieldnames: list[str]) -> str:
    lower_map = {f.lower(): f for f in fieldnames}
    for candidate in CANDIDATE_ID_COLUMNS:
        if candidate in lower_map:
            return lower_map[candidate]
    # fall back: first column whose values look like uniprot accessions
    return fieldnames[0]


def load_accessions(csv_path: Path, id_column: str | None) -> list[str]:
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise ValueError(f"No columns found in {csv_path}")
        col = id_column or detect_id_column(fieldnames)
        if col not in fieldnames:
            raise ValueError(
                f"Column '{col}' not found. Available columns: {fieldnames}"
            )
        log.info("Using CSV column '%s' for UniProt accessions", col)
        rows = list(reader)

    accessions = []
    seen = set()
    for row in rows:
        val = (row.get(col) or "").strip()
        if not val or val in seen:
            continue
        seen.add(val)
        accessions.append(val)
    log.info("Loaded %d unique accessions from %s", len(accessions), csv_path)
    return accessions


def fetch_publications(session: requests.Session, accession: str) -> list | None:
    """
    Fetch UniProt's publication list for an accession from rest.uniprot.org
    (the API behind the uniprot.org/uniprotkb/{acc}/publications page).

    Returns a list of result objects (following pagination), or None if the
    accession is not found / invalid (HTTP 404 or 400).
    """
    url = PUBLICATIONS_URL.format(accession=accession)
    params = {"size": 500}
    results = []
    while url:
        resp = session.get(url, params=params, timeout=30)
        if resp.status_code in (400, 404):
            return None
        resp.raise_for_status()
        params = None  # subsequent pages carry the cursor in the Link URL itself
        if resp.text.strip():
            data = resp.json()
            results.extend(data.get("results", []))
        # UniProt paginates via a Link header: <...cursor=...>; rel="next"
        url = resp.links.get("next", {}).get("url")
    return results


def is_inactive(payload) -> bool:
    """
    The publications sub-resource doesn't return the deleted-entry stub that
    the main entry endpoint does, so an empty result list here just means
    'zero publications', which is common for unreviewed entries and is NOT
    treated as inactive.
    """
    return False


def parse_publication_records(accession: str, payload) -> list[dict]:
    """
    Normalizes the rest.uniprot.org publications payload into a flat list of
    {title, authors, journal, year, pmid, doi, categories} records.

    Shape (per result object):
      {
        "citation": {
          "title": ..., "authors": [...], "journal": ...,
          "publicationDate": ...,
          "citationCrossReferences": [{"database": "PubMed", "id": ...},
                                      {"database": "DOI", "id": ...}]
        },
        "references": [{"sourceCategories": [...], "source": {...}}, ...],
        "statistics": {...}
      }
    """
    if not payload:
        return []

    # fetch_publications returns a list of result objects; be tolerant if a
    # raw {"results": [...]} dict is passed instead.
    if isinstance(payload, dict):
        results = payload.get("results", [])
    elif isinstance(payload, list):
        results = payload
    else:
        results = []

    out = []
    for result in results:
        if not isinstance(result, dict):
            continue
        citation = result.get("citation") or {}
        if not citation:
            continue

        title = citation.get("title") or ""
        authors = citation.get("authors") or []
        if isinstance(authors, str):
            authors = [authors]
        journal = citation.get("journal") or ""
        year = citation.get("publicationDate") or ""

        pmid, doi = None, None
        for xref in citation.get("citationCrossReferences") or []:
            db = (xref.get("database") or "").upper()
            if db == "PUBMED":
                pmid = xref.get("id")
            elif db == "DOI":
                doi = xref.get("id")

        # sourceCategories live inside each entry of the references array
        categories = []
        for ref in result.get("references") or []:
            for cat in ref.get("sourceCategories") or []:
                if cat not in categories:
                    categories.append(cat)

        out.append({
            "accession": accession,
            "title": title.strip() if isinstance(title, str) else "",
            "authors": authors,
            "journal": journal,
            "year": year,
            "pmid": pmid,
            "doi": doi,
            "categories": categories,
        })

    # de-dup by (pmid or title)
    seen = set()
    deduped = []
    for rec in out:
        key = rec["pmid"] or rec["title"]
        if key and key not in seen:
            seen.add(key)
            deduped.append(rec)
    return deduped


def fetch_europepmc(session: requests.Session, pmid: str) -> dict:
    """Fetch abstract (always) and full text XML (only if open access)."""
    result = {"abstract": "", "is_open_access": False, "fulltext": ""}
    try:
        resp = session.get(
            EPMC_SEARCH_URL,
            params={
                "query": f"EXT_ID:{pmid} AND SRC:MED",
                "format": "json",
                "resultType": "core",
            },
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        hits = data.get("resultList", {}).get("result", [])
        if not hits:
            return result
        hit = hits[0]
        result["abstract"] = hit.get("abstractText", "") or ""
        is_oa = (hit.get("isOpenAccess") or "N").upper() == "Y"
        result["is_open_access"] = is_oa
        if is_oa and hit.get("pmcid"):
            ft_resp = session.get(
                EPMC_FULLTEXT_URL.format(source="PMC", ext_id=hit["pmcid"]),
                timeout=30,
            )
            if ft_resp.status_code == 200 and ft_resp.text.strip():
                # Strip XML tags to plain text (lightweight, no extra deps)
                text = re.sub(r"<[^>]+>", " ", ft_resp.text)
                text = re.sub(r"\s+", " ", text).strip()
                result["fulltext"] = text
    except requests.RequestException as e:
        log.warning("Europe PMC fetch failed for PMID %s: %s", pmid, e)
    return result


def write_protein_outputs(out_dir: Path, accession: str, records: list[dict]) -> None:
    json_path = out_dir / f"{accession}.json"
    txt_path = out_dir / f"{accession}.txt"

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({"accession": accession, "publications": records}, f, indent=2, ensure_ascii=False)

    lines = [f"UniProt accession: {accession}", ""]
    for rec in records:
        lines.append(f"Title: {rec['title']}")
        if rec["authors"]:
            lines.append(f"Authors: {', '.join(rec['authors'])}")
        if rec["journal"] or rec["year"]:
            lines.append(f"Source: {rec['journal']} ({rec['year']})")
        if rec.get("categories"):
            lines.append(f"Cited for: {', '.join(rec['categories'])}")
        body = rec.get("fulltext") or rec.get("abstract") or ""
        if body:
            lines.append("")
            lines.append(body)
        lines.append("")
        lines.append("---")
        lines.append("")

    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description="Scrape UniProt publications into a RAG corpus.")
    parser.add_argument("--csv", default=DEFAULT_CSV, type=Path, help=f"Path to proteins.csv (default: {DEFAULT_CSV})")
    parser.add_argument("--out", default=DEFAULT_OUT, type=Path, help=f"Output directory for per-protein JSON/TXT (default: {DEFAULT_OUT})")
    parser.add_argument("--id-column", default=None, help="Column name holding UniProt accessions (auto-detected if omitted)")
    parser.add_argument("--email", default=None, help="Contact email sent in the User-Agent header (courtesy to EBI/Europe PMC)")
    parser.add_argument("--sleep", type=float, default=0.34, help="Seconds to sleep between accessions (default ~3 req/s)")
    parser.add_argument("--limit", type=int, default=None, help="Only process the first N accessions (for testing)")
    parser.add_argument("--overwrite", action="store_true", help="Re-scrape proteins that already have output files")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    session = make_session(args.email)

    accessions = load_accessions(args.csv, args.id_column)
    if args.limit:
        accessions = accessions[: args.limit]

    inactive_log = args.out / "_inactive.txt"
    error_log = args.out / "_errors.txt"

    n_done = n_skipped = n_inactive = n_errors = 0

    pbar = tqdm(accessions, desc="Scraping", unit="protein")
    for accession in pbar:
        pbar.set_postfix(done=n_done, skip=n_skipped, inactive=n_inactive, err=n_errors)
        json_path = args.out / f"{accession}.json"
        if json_path.exists() and not args.overwrite:
            n_skipped += 1
            continue

        try:
            payload = fetch_publications(session, accession)
        except requests.RequestException as e:
            tqdm.write(f"[ERROR] failed to fetch {accession}: {e}")
            with open(error_log, "a", encoding="utf-8") as f:
                f.write(f"{accession}\t{e}\n")
            n_errors += 1
            time.sleep(args.sleep)
            continue

        if payload is None:
            with open(inactive_log, "a", encoding="utf-8") as f:
                f.write(f"{accession}\tnot_found\n")
            n_inactive += 1
            time.sleep(args.sleep)
            continue

        if is_inactive(payload):
            with open(inactive_log, "a", encoding="utf-8") as f:
                f.write(f"{accession}\tinactive\n")
            n_inactive += 1
            time.sleep(args.sleep)
            continue

        records = parse_publication_records(accession, payload)

        payload_has_content = isinstance(payload, list) and len(payload) > 0
        if payload_has_content and not records:
            tqdm.write(
                f"[WARN] {accession}: payload had {len(payload)} result(s) but 0 "
                f"records parsed -- schema mismatch likely; raw JSON saved to _raw_unparsed/"
            )
            raw_dir = args.out / "_raw_unparsed"
            raw_dir.mkdir(exist_ok=True)
            with open(raw_dir / f"{accession}.json", "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, ensure_ascii=False)

        for rec in records:
            if rec["pmid"]:
                epmc = fetch_europepmc(session, rec["pmid"])
                rec["abstract"] = epmc["abstract"]
                rec["is_open_access"] = epmc["is_open_access"]
                rec["fulltext"] = epmc["fulltext"]
                time.sleep(0.1)

        write_protein_outputs(args.out, accession, records)
        n_done += 1
        time.sleep(args.sleep)

    log.info(
        "Done. scraped=%d skipped(existing)=%d inactive=%d errors=%d",
        n_done, n_skipped, n_inactive, n_errors,
    )


if __name__ == "__main__":
    main()