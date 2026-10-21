"""
small_molecule_scraper.py — end to end, from your CSV of ChEBI IDs to RAG documents.

Per compound it writes one .txt containing:

  1. Identity + chemical properties   (PubChem PUG-REST)
  2. ChEBI definition                 (ChEBI via OLS4, PubChem as fallback)
  3. Europe PMC abstracts             where a NAME VARIANT appears in
                                      TITLE or ABSTRACT
  4. Related-compound abstracts       only if step 3 is empty, clearly
                                      labelled so the LLM knows it's borrowed

Skips the junk that made the last corpus ~96% noise: no synonym dumps,
no co-occurrence table headers, no spectra, no cross-reference link lists.
Synonyms are fetched but used ONLY for name matching — never written out.

Why name variants matter: papers write "11-oxo cucurbitadienol", ChEBI calls
it "11-Oxocucurbitadienol". An exact-phrase search on one spelling finds
nothing. This builds spacing/hyphen variants plus short PubChem synonyms and
ORs them into a single query. Long IUPAC names are dropped — they never
appear in an abstract.

Set SINGLE to one ChEBI ID to preview before the full run.
Resumable — rerun after a stop and it picks up where it left off.

    python rag/small_molecule_scraper.py
"""

import csv
import json
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

# ============================ CONFIG ============================

ROOT = Path(__file__).resolve().parents[1]  # repo root (rag/ → ..)

CSV_PATH = ROOT / "data" / "molecules.csv"
CHEBI_COLUMN = "chebi_id"          # <-- your column name

SINGLE = None                      # e.g. "CHEBI:138973" to preview one
DEBUG = False                      # True prints variants, queries, parents

OUT_DIR = ROOT / "rag_data" / "rag_train" / "small_molecule"
RAW_DIR = OUT_DIR / "raw_json"
STATUS_CSV = OUT_DIR / "_status.csv"

MAX_OWN_PAPERS = 5
MAX_PARENT_PAPERS = 3
MAX_VARIANTS = 8                   # name spellings to OR together
MAX_NAME_CHARS = 60                # longer = IUPAC, never appears in abstracts
MIN_ABSTRACT_CHARS = 200           # skip stub abstracts
MAX_SYNONYMS = 25                  # scanned for matching, not written out
REQUEST_DELAY = 0.25               # PubChem allows max 5 req/sec
TIMEOUT = 30
RETRIES = 3

PROPERTIES = ("MolecularFormula,MolecularWeight,IUPACName,SMILES,"
              "InChI,InChIKey,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount")

PROP_LABELS = [
    ("IUPACName", "IUPAC Name"), ("MolecularFormula", "Molecular Formula"),
    ("MolecularWeight", "Molecular Weight"), ("SMILES", "SMILES"),
    ("InChI", "InChI"), ("InChIKey", "InChIKey"), ("XLogP", "XLogP"),
    ("TPSA", "TPSA"), ("HBondDonorCount", "H-Bond Donors"),
    ("HBondAcceptorCount", "H-Bond Acceptors"),
]

# ================================================================

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUGVIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"
OLS4 = "https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# ChEBI definitions state relationships in plain English — parse them.
DEF_PATTERNS = [
    (r"functionally related to (?:an?\s+)?([^.,;]+)", "is functionally related to"),
    (r"conjugate base of (?:an?\s+)?([^.,;]+)",       "is the conjugate base of"),
    (r"conjugate acid of (?:an?\s+)?([^.,;]+)",       "is the conjugate acid of"),
    (r"tautomer of (?:an?\s+)?([^.,;]+)",             "is a tautomer of"),
    (r"enantiomer of (?:an?\s+)?([^.,;]+)",           "is an enantiomer of"),
    (r"derives from (?:a hydride of )?(?:an?\s+)?([^.,;]+)", "derives from"),
]

# Vendor / registry codes that pollute PubChem synonym lists.
JUNK_SYNONYM = re.compile(
    r"(CHEBI|CHEMBL|SCHEMBL|DTXSID|DTXCID|AKOS|MFCD|CAS-|UNII|HMDB|LMBA|"
    r"NSC\d|BDBM|ZINC|EINECS|NCGC|SMR\d|MLS\d|Q\d{6,}|^\d[\d\-]+$)",
    re.IGNORECASE,
)

session = requests.Session()
session.headers.update({"User-Agent": "steroid-atlas-rag/1.0 (contact@example.com)"})


# ------------------------- HTTP helpers -------------------------

def get_json(url, params=None):
    """GET with backoff. Returns None on 404 or repeated failure."""
    for attempt in range(RETRIES):
        try:
            r = session.get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 404:
                return None
            if r.status_code in (429, 500, 502, 503):
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError):
            time.sleep(2 ** attempt)
    if DEBUG:
        print(f"    ! failed: {url}")
    return None


# ------------------------- ChEBI --------------------------------

def chebi_term(chebi_id):
    """Label + definition straight from ChEBI via OLS4."""
    num = chebi_id.replace("CHEBI:", "").replace("CHEBI_", "").strip()
    iri = f"http://purl.obolibrary.org/obo/CHEBI_{num}"
    data = get_json(OLS4, {"iri": iri})
    try:
        term = data["_embedded"]["terms"][0]
    except (TypeError, KeyError, IndexError):
        return None, None
    label = term.get("label")
    descs = term.get("description") or []
    return label, (descs[0].strip() if descs else None)


# ------------------------- PubChem ------------------------------

def chebi_to_cid(chebi_id, label):
    """ChEBI ID -> PubChem CID. Registry cross-ref first, name lookup second."""
    num = chebi_id.replace("CHEBI:", "").replace("CHEBI_", "").strip()
    data = get_json(f"{PUG}/compound/xref/RegistryID/CHEBI:{num}/cids/JSON")
    cids = (data or {}).get("IdentifierList", {}).get("CID", [])
    if cids:
        return cids[0]
    if label:
        time.sleep(REQUEST_DELAY)
        safe = requests.utils.quote(label, safe="")
        data = get_json(f"{PUG}/compound/name/{safe}/cids/JSON")
        cids = (data or {}).get("IdentifierList", {}).get("CID", [])
        if cids:
            return cids[0]
    return None


def pubchem_properties(cid):
    data = get_json(f"{PUG}/compound/cid/{cid}/property/{PROPERTIES}/JSON")
    try:
        return data["PropertyTable"]["Properties"][0]
    except (TypeError, KeyError, IndexError):
        return {}


def pubchem_synonyms(cid):
    data = get_json(f"{PUG}/compound/cid/{cid}/synonyms/JSON")
    try:
        syns = data["InformationList"]["Information"][0]["Synonym"]
    except (TypeError, KeyError, IndexError):
        return []
    return syns[:MAX_SYNONYMS]


def pubchem_description(cid):
    """Fallback definition if ChEBI/OLS4 came back empty."""
    data = get_json(f"{PUGVIEW}/{cid}/JSON", {"heading": "Record Description"})
    texts = []

    def walk(node):
        if isinstance(node, dict):
            for item in node.get("Information", []) or []:
                for s in (item.get("Value", {}).get("StringWithMarkup") or []):
                    if s.get("String"):
                        texts.append(s["String"].strip())
            for section in node.get("Section", []) or []:
                walk(section)
            if "Record" in node:
                walk(node["Record"])
        elif isinstance(node, list):
            for n in node:
                walk(n)

    walk(data or {})
    texts = [t for t in texts if len(t) > 80]
    return texts[0] if texts else None


# ------------------------- name variants ------------------------

def _spacing_variants(name):
    """11-Oxocucurbitadienol -> '11-oxo cucurbitadienol', '11 oxo ...', etc."""
    out = {name}
    out.add(name.replace("-", " "))
    out.add(name.replace("-", ""))
    # split a leading locant+substituent off a fused word: 11-Oxocucurbita... 
    m = re.match(r"^(\d+[a-zA-Z]?-)(oxo|hydroxy|keto|amino|methyl|deoxy)(\w{6,})$",
                 name, re.IGNORECASE)
    if m:
        out.add(f"{m.group(1)}{m.group(2)} {m.group(3)}")
        out.add(f"{m.group(1).rstrip('-')} {m.group(2)} {m.group(3)}")
    return out


def name_variants(primary, synonyms):
    """Short, plausible spellings — shortest first, capped at MAX_VARIANTS."""
    pool = set()
    for candidate in [primary, *synonyms]:
        if not candidate:
            continue
        candidate = candidate.strip()
        if len(candidate) > MAX_NAME_CHARS or len(candidate) < 4:
            continue
        if JUNK_SYNONYM.search(candidate) or '"' in candidate:
            continue
        pool |= _spacing_variants(candidate)

    seen, ordered = set(), []
    for v in sorted(pool, key=len):
        key = v.lower()
        if key in seen or len(v) < 4 or len(v) > MAX_NAME_CHARS:
            continue
        seen.add(key)
        ordered.append(v)
    return ordered[:MAX_VARIANTS]


# ------------------------- Europe PMC ---------------------------

def epmc_abstracts(variants, limit):
    """Abstracts where any variant is an exact phrase in TITLE or ABSTRACT."""
    if not variants:
        return [], ""
    clause = " OR ".join(f'(TITLE:"{v}" OR ABSTRACT:"{v}")' for v in variants)
    query = f"({clause}) AND (SRC:MED OR SRC:PMC) AND HAS_ABSTRACT:Y"
    if DEBUG:
        print(f"    query: {query[:220]}")

    data = get_json(EPMC, {
        "query": query, "format": "json",
        "resultType": "core", "pageSize": limit * 3,
    })
    results = (data or {}).get("resultList", {}).get("result", [])

    papers, seen = [], set()
    for r in results:
        abstract = (r.get("abstractText") or "").strip()
        abstract = re.sub(r"<[^>]+>", "", abstract)
        pmid = r.get("pmid") or r.get("id")
        if not abstract or len(abstract) < MIN_ABSTRACT_CHARS or pmid in seen:
            continue
        seen.add(pmid)
        papers.append({
            "pmid": pmid,
            "title": (r.get("title") or "").strip().rstrip("."),
            "journal": r.get("journalTitle") or "",
            "year": r.get("pubYear") or "",
            "abstract": abstract,
        })
        if len(papers) >= limit:
            break
    return papers, query


def parents_from_definition(definition):
    """['cucurbitadienol', ...] parsed out of the ChEBI definition text."""
    if not definition:
        return []
    found = []
    for pattern, relation in DEF_PATTERNS:
        for m in re.finditer(pattern, definition, re.IGNORECASE):
            term = m.group(1).strip().rstrip(".")
            if 3 < len(term) <= MAX_NAME_CHARS:
                found.append((term, relation))
    seen, out = set(), []
    for term, relation in found:
        if term.lower() not in seen:
            seen.add(term.lower())
            out.append((term, relation))
    return out


# ------------------------- document build -----------------------

def build_doc(chebi_id):
    label, definition = chebi_term(chebi_id)
    time.sleep(REQUEST_DELAY)

    cid = chebi_to_cid(chebi_id, label)
    props, synonyms = {}, []
    if cid:
        time.sleep(REQUEST_DELAY)
        props = pubchem_properties(cid)
        time.sleep(REQUEST_DELAY)
        synonyms = pubchem_synonyms(cid)
        if not definition:
            time.sleep(REQUEST_DELAY)
            definition = pubchem_description(cid)

    if not definition and not cid:
        raise ValueError("resolved to nothing in ChEBI or PubChem "
                         "(bad ID, or both services unreachable)")

    name = label or (synonyms[0] if synonyms else None) or chebi_id
    variants = name_variants(name, synonyms)
    if DEBUG:
        print(f"    variants: {variants}")

    papers, _ = epmc_abstracts(variants, MAX_OWN_PAPERS)
    status = "own_literature" if papers else None
    borrowed = []

    if not papers:
        for parent, relation in parents_from_definition(definition)[:2]:
            if DEBUG:
                print(f"    parent: {parent} ({relation})")
            time.sleep(REQUEST_DELAY)
            found, _ = epmc_abstracts(name_variants(parent, []), MAX_PARENT_PAPERS)
            if found:
                borrowed.append((parent, relation, found))
                break
        status = "parent_literature" if borrowed else "no_literature"

    # ---- render ----
    lines = [f"Compound: {name}", f"ChEBI ID: {chebi_id}"]
    if cid:
        lines.append(f"PubChem CID: {cid}")

    prop_lines = [f"{lab}: {props[key]}" for key, lab in PROP_LABELS
                  if props.get(key) not in (None, "")]
    if prop_lines:
        lines += ["", "== Chemical properties ==", *prop_lines]

    if definition:
        lines += ["", "== ChEBI definition ==", definition]

    if papers:
        lines += ["", f"== Literature on {name} =="]
        for p in papers:
            lines += ["", f"[PMID {p['pmid']}] {p['title']}",
                      f"{p['journal']} {p['year']}".strip(), p["abstract"]]

    for parent, relation, found in borrowed:
        lines += ["", f"--- RELATED COMPOUND: {parent} ---",
                  f"No literature found for {name} itself. {name} {relation} "
                  f"{parent}; the abstracts below are about {parent}."]
        for p in found:
            lines += ["", f"[PMID {p['pmid']}] {p['title']}",
                      f"{p['journal']} {p['year']}".strip(), p["abstract"]]

    doc = "\n".join(lines).strip() + "\n"
    meta = {
        "chebi_id": chebi_id, "name": name, "cid": cid, "status": status,
        "n_papers": len(papers) or sum(len(f) for _, _, f in borrowed),
        "chars": len(doc), "has_definition": bool(definition),
    }
    return doc, meta


# ------------------------- runner -------------------------------

def normalize_chebi(raw):
    """'1301.0', 'CHEBI:1301.0', 1301.0, 'chebi_1301' -> 'CHEBI:1301'.

    The int(float(...)) is load-bearing: pandas reads a numeric ChEBI column
    as float when it contains any NaN, and 'CHEBI:1301.0' matches nothing in
    ChEBI or PubChem — every downstream lookup then returns empty.
    """
    s = str(raw).strip().strip("[]'\" ")
    s = re.sub(r"(?i)^chebi[:_\s]*", "", s)
    if not s or s.lower() in ("nan", "none"):
        return None
    try:
        return f"CHEBI:{int(float(s))}"
    except ValueError:
        return None


def load_chebi_ids():
    if SINGLE:
        return [normalize_chebi(SINGLE)]
    df = pd.read_csv(os.path.expanduser(str(CSV_PATH)), low_memory=False,
                     dtype={CHEBI_COLUMN: str})     # read as text, never float
    if CHEBI_COLUMN not in df.columns:
        raise SystemExit(f"No '{CHEBI_COLUMN}' column. Found: {list(df.columns)[:20]}")
    ids = []
    for raw in df[CHEBI_COLUMN].dropna().astype(str):
        for part in re.split(r"[;,]", raw):
            cid = normalize_chebi(part)
            if cid:
                ids.append(cid)
    ids = sorted(set(ids), key=lambda x: int(x.split(":")[1]))
    print(f"sample: {ids[:5]}")                     # eyeball before a long run
    return ids


def main():
    (OUT_DIR / "data_txt").mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    chebi_ids = load_chebi_ids()
    print(f"{len(chebi_ids)} ChEBI IDs")

    rows, skipped, failures = [], 0, 0
    for chebi_id in tqdm(chebi_ids, desc="compounds"):
        slug = chebi_id.replace(":", "_")
        out_path = OUT_DIR / "data_txt" / f"{slug}.txt"
        if out_path.exists() and not SINGLE:
            skipped += 1
            continue
        if DEBUG or SINGLE:
            print(f"\n{chebi_id}")
        try:
            doc, meta = build_doc(chebi_id)
        except Exception as e:
            print(f"  ! {chebi_id}: {e}")
            failures += 1
            if failures >= 10 and not rows:
                raise SystemExit(
                    "First 10 compounds all failed — check the ChEBI IDs "
                    "printed above and your network before burning the rest."
                )
            continue
        out_path.write_text(doc, encoding="utf-8")
        (RAW_DIR / f"{slug}.json").write_text(json.dumps(meta, indent=2))
        rows.append(meta)
        if SINGLE:
            print("\n" + "=" * 70 + "\n" + doc)

    if rows:
        with open(STATUS_CSV, "a", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            if STATUS_CSV.stat().st_size == 0:
                w.writeheader()
            w.writerows(rows)

    total = len(rows) or 1
    print(f"\nwrote {len(rows)} · skipped {skipped} (already present) · "
          f"failed {failures}")
    for status in ("own_literature", "parent_literature", "no_literature"):
        n = sum(1 for r in rows if r["status"] == status)
        print(f"  {status:20s} {n:5d}  ({100 * n / total:.1f}%)")
    if rows:
        print(f"  median doc size: "
              f"{sorted(r['chars'] for r in rows)[len(rows) // 2]:,} chars")
    print(f"\noutput: {OUT_DIR}")


if __name__ == "__main__":
    main()