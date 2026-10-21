"""Full audit: for every ChEBI ID in small_molecule_centric.csv, look up the
canonical ChEBI label + synonyms via OLS and compare against the CSV's
Compound Name. Flag any mismatches.

Approach: cheap fuzzy match on lowercased alphanumerics. If ANY significant
token from the CSV compound name appears in the ChEBI label or a synonym
(or vice-versa), consider it a match. Otherwise flag it.

Output: literature/chebi_audit_report.tsv
         (columns: idx, csv_compound, chebi_id, chebi_label, top_synonym, verdict, note)

Then a summary of confirmed matches / mismatches / lookup errors.

Rate-limited: 0.15s between OLS calls. ~90 seconds for 589 unique IDs.
"""
from __future__ import annotations
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

CSV = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/small_molecule_centric.csv")
REPORT = Path("/home/adsiordia/marimo_visualizer/MarimoSteroidVisualizer/literature/chebi_audit_report.tsv")

UA = "adsiordia@ucsd.edu (steroid atlas audit)"
DELAY = 0.15

# tokens too generic to give a fuzzy-match signal
STOP_TOKENS = {
    "acid", "one", "ol", "diol", "triol", "steroid", "hydroxy", "keto",
    "alpha", "beta", "the", "and", "of", "in", "with", "an", "a"
}

def norm(s: str) -> str:
    """Lowercase, keep only alphanumerics."""
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())

def tokenize(s: str) -> list[str]:
    """Split into >2-char alphanumeric tokens (lowercase), filter stopwords."""
    parts = re.findall(r"[A-Za-z0-9]{3,}", (s or "").lower())
    return [p for p in parts if p not in STOP_TOKENS]

def ols_lookup(chebi_id: str) -> tuple[str, list[str]]:
    """Returns (label, synonyms) or ('', []) if not found."""
    iri = f"http://purl.obolibrary.org/obo/{chebi_id.replace(':', '_')}"
    url = f"https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms?iri={urllib.parse.quote(iri, safe=':/')}"
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=20) as r:
        data = json.load(r)
    terms = data.get("_embedded", {}).get("terms", [])
    if not terms:
        return "", []
    t = terms[0]
    return t.get("label", ""), list(t.get("synonyms", []) or [])

def looks_like_match(csv_name: str, chebi_label: str, synonyms: list[str]) -> tuple[bool, str]:
    """Fuzzy match. Returns (is_match, evidence)."""
    csv_tokens = set(tokenize(csv_name))
    chebi_bag = " ".join([chebi_label] + list(synonyms))
    chebi_tokens = set(tokenize(chebi_bag))
    if not csv_tokens or not chebi_tokens:
        return False, "empty tokens"
    overlap = csv_tokens & chebi_tokens
    csv_norm = norm(csv_name)
    label_norm = norm(chebi_label)
    syn_norms = [norm(s) for s in synonyms]
    # Substring match either direction: strong signal
    if len(csv_norm) > 5 and (csv_norm in label_norm or label_norm in csv_norm):
        return True, f"substring vs label ({chebi_label!r})"
    for sn in syn_norms:
        if len(csv_norm) > 5 and (csv_norm in sn or sn in csv_norm):
            return True, f"substring vs synonym"
    # Token overlap: need at least one meaningful token
    if overlap:
        return True, f"token overlap: {sorted(overlap)[:3]}"
    return False, f"no overlap (csv={sorted(csv_tokens)[:4]}, chebi={sorted(chebi_tokens)[:4]})"

def main() -> int:
    # Read ChEBI ID column as string so numeric IDs don't become "17263.0" floats
    df = pd.read_csv(CSV, low_memory=False, dtype={"ChEBI ID": str})
    # Only rows with a non-empty ChEBI ID and Compound Name
    df["_id"] = df["ChEBI ID"].astype(str).str.strip()
    # Strip stray trailing ".0" (in case anything else re-typed the column upstream)
    df["_id"] = df["_id"].str.replace(r"\.0$", "", regex=True)
    df["_name"] = df["Compound Name"].astype(str).str.strip()
    todo = df[(df["_id"] != "") & (df["_id"] != "nan") & (df["_name"] != "") & (df["_name"] != "nan")].copy()
    print(f"Rows to audit: {len(todo)}")

    unique_ids = sorted(set(todo["_id"].tolist()))
    print(f"Unique ChEBI IDs: {len(unique_ids)}")

    lookup_cache: dict[str, tuple[str, list[str]]] = {}
    errors: dict[str, str] = {}
    t0 = time.time()
    for i, cid in enumerate(unique_ids):
        chebi_id = f"CHEBI:{cid}" if not cid.startswith("CHEBI:") else cid
        try:
            lookup_cache[cid] = ols_lookup(chebi_id)
        except Exception as e:
            errors[cid] = str(e)[:120]
            lookup_cache[cid] = ("", [])
        if (i+1) % 50 == 0:
            print(f"  {i+1}/{len(unique_ids)}  ({(i+1)/(time.time()-t0):.1f} lookups/s)")
        time.sleep(DELAY)

    print(f"\nDone lookups in {time.time()-t0:.0f}s")
    print(f"  successful lookups: {len([v for v in lookup_cache.values() if v[0]])}")
    print(f"  lookup errors:      {len(errors)}")

    # Score each row
    rows_out = []
    for _, row in todo.iterrows():
        cid = row["_id"]
        name = row["_name"]
        label, synonyms = lookup_cache.get(cid, ("", []))
        if not label and cid in errors:
            verdict = "ERROR"
            evidence = errors[cid]
        elif not label:
            verdict = "NOT_FOUND"
            evidence = "ChEBI ID not in OLS"
        else:
            ok, why = looks_like_match(name, label, synonyms)
            verdict = "MATCH" if ok else "MISMATCH"
            evidence = why
        rows_out.append({
            "idx": row.name,
            "csv_compound": name,
            "chebi_id": cid,
            "chebi_label": label,
            "top_synonym": synonyms[0] if synonyms else "",
            "verdict": verdict,
            "note": evidence,
            "is_new": row.get("is_new", 0),
        })

    out = pd.DataFrame(rows_out)
    out.to_csv(REPORT, sep="\t", index=False)
    print(f"\nReport -> {REPORT.name}")

    counts = out["verdict"].value_counts()
    print("\nVerdict summary:")
    for v, n in counts.items():
        print(f"  {v:<10}  {n}")

    # Show all mismatches
    mism = out[out["verdict"] == "MISMATCH"]
    if len(mism):
        print(f"\n--- All {len(mism)} MISMATCHES ---")
        for _, r in mism.iterrows():
            marker = "[NEW]" if r["is_new"] == 1 else "[old]"
            print(f"  {marker} idx={r['idx']:<4}  ID={r['chebi_id']:<8}  csv={r['csv_compound']!r}")
            print(f"          chebi label={r['chebi_label']!r}")
            print(f"          note={r['note']}")
            print()

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
