"""Re-audit the 32,869 Rhea-recruited entries against current UniProt.

We previously trusted the rhea2uniprot mapping used at atlas construction
time — but UniProt curates annotations continuously, and TrEMBL entries in
particular can have their Rhea reactions refined or removed after build.
The A0A498N7D4 case exposed this: our atlas had it linked to RHEA:60128
(estradiol-glucuronide transport), but the current UniProt entry lists
only RHEA:85015 (2'->5' oligoadenylate transport, a non-steroid).

For every Rhea-source entry in the atlas, this script:
  1. Fetches the current UniProt annotation (CC CATALYTIC ACTIVITY + Rhea
     xrefs + FT BINDING features)
  2. Enumerates every ChEBI id that appears in the entry's current Rhea
     reactions or binding-site features
  3. Checks whether ANY of those ChEBI ids is in our 677-molecule
     sterane-passing set from molecules.csv
  4. Classifies:
       KEEP_RHEA         current Rhea reaction lists a sterane-passing ChEBI
       KEEP_LIGAND       no steroid Rhea, but a sterane-passing ligand ChEBI
                         appears in a binding-site feature
       DROP_STALE        neither — the recruitment reason has been removed
                         from UniProt since the atlas was built

Reads:  ../data/proteins.csv, ../data/molecules.csv
Writes: ../analysis/rhea_audit.tsv                  per-entry decision table
        ../analysis/rhea_audit_summary.txt          human-readable summary
        ../analysis/rhea_audit_cache.json           raw UniProt cache (resumable)
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
MOLECULES = ROOT / "data" / "molecules.csv"
CACHE = output("rhea_audit_cache.json", historical=True)
OUT_TSV = output("rhea_audit.tsv", historical=True)
OUT_TXT = output("rhea_audit_summary.txt", historical=True)

BULK_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = "accession,cc_catalytic_activity,rhea,ft_binding"
BATCH_SIZE = 300


def load_sterane_chebis() -> set[str]:
    """Same normalization as 04/05 — read molecules.csv and produce CHEBI:<n> strings.

    Also adds a small set of alternative ChEBI ids for compounds that molecules.csv
    covers under one id but UniProt now references under another (e.g. 5-alpha-DHT).
    These were found during the Rhea re-audit false-positive check.
    """
    mol = pd.read_csv(MOLECULES, low_memory=False)
    out = set()
    for x in mol["chebi_id"].dropna():
        try:
            out.add(f"CHEBI:{int(float(x))}")
        except (ValueError, TypeError):
            pass
    # Alternative ChEBI ids for compounds already in the sterane set under a
    # different id — added after Q9Y394 (DHRS7) false-positive was found.
    ALT_IDS = {
        "CHEBI:17336",   # 5-alpha-dihydrotestosterone (5-alpha-DHT); molecules.csv has 5-beta only (CHEBI:2150)
        "CHEBI:17898",   # 5-beta-dihydrotestosterone alternative id (same molecule as CHEBI:2150)
        "CHEBI:16330",   # 5-alpha-androstane-3,17-dione
        "CHEBI:36713",   # androgen (generic superclass)
    }
    out |= ALT_IDS
    return out


def norm_chebi(cid: str) -> str | None:
    s = str(cid).strip().upper()
    while s.startswith("CHEBI:"):
        s = s[len("CHEBI:"):]
    if s.isdigit():
        return f"CHEBI:{s}"
    return None


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


def extract(entry: dict) -> dict:
    """Return the sets of Rhea ids, reaction-ChEBIs, and ligand-ChEBIs currently listed."""
    acc = entry.get("primaryAccession", "")
    rhea_ids: set[str] = set()
    rxn_chebis: set[str] = set()
    ligand_chebis: set[str] = set()

    for c in entry.get("comments", []):
        if c.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        reaction = c.get("reaction", {}) or {}
        for xref in reaction.get("reactionCrossReferences", []) or []:
            db = (xref.get("database") or "").upper()
            xid = xref.get("id") or ""
            if db == "RHEA":
                # Rhea ids look like "RHEA:60128"
                if xid.upper().startswith("RHEA:"):
                    rhea_ids.add(xid)
                else:
                    rhea_ids.add(f"RHEA:{xid}")
            elif db == "CHEBI":
                normed = norm_chebi(xid)
                if normed:
                    rxn_chebis.add(normed)

    for feat in entry.get("features", []):
        if feat.get("type") != "Binding site":
            continue
        lig = feat.get("ligand", {}) or {}
        normed = norm_chebi(lig.get("id", ""))
        if normed:
            ligand_chebis.add(normed)

    return {
        "accession": acc,
        "rhea_ids": sorted(rhea_ids),
        "rxn_chebis": sorted(rxn_chebis),
        "ligand_chebis": sorted(ligand_chebis),
    }


STEROID_VOCAB = ("steroid", "sterol", "estrogen", "estradiol", "testosterone",
                 "androgen", "progesterone", "cortisol", "corticosteroid",
                 "aldosterone", "bile acid", "bile-acid", "cholesterol",
                 "pregnenolone", "ecdysone")
STEROID_KWS = ("KW-0754", "KW-0675")  # Steroid-binding, Steroid hormone receptor


def classify(ann: dict, sterane: set[str],
             go_labels: str = "", keyword_ids: str = "",
             protein_names: str = "") -> tuple[str, str]:
    """Per-entry classification with GO/keyword rescue.

    Priority:
      1. Current Rhea reaction lists a sterane-passing ChEBI          → KEEP_RHEA
      2. Current binding-site feature has sterane-passing ligand      → KEEP_LIGAND
      3. Atlas-stored GO label or keyword mentions steroid vocabulary → KEEP_ANNOTATION
         (catches entries whose Rhea was pruned but whose functional
         annotation still attests to steroid activity)
      4. Otherwise                                                    → DROP_STALE
    """
    rxn_hits = set(ann["rxn_chebis"]) & sterane
    lig_hits = set(ann["ligand_chebis"]) & sterane
    if rxn_hits:
        return ("KEEP_RHEA", f"current Rhea reaction has sterane-passing ChEBI: {sorted(rxn_hits)[:3]}")
    if lig_hits:
        return ("KEEP_LIGAND", f"binding-site feature has sterane-passing ligand: {sorted(lig_hits)[:3]}")

    # Rescue via atlas-stored GO / keyword / name annotations
    go_lo = (go_labels or "").lower()
    kw_str = keyword_ids or ""
    name_lo = (protein_names or "").lower()
    go_hits = [v for v in STEROID_VOCAB if v in go_lo]
    kw_hits = [k for k in STEROID_KWS if k in kw_str]
    name_hits = [v for v in STEROID_VOCAB if v in name_lo]
    if go_hits or kw_hits or name_hits:
        reasons = []
        if go_hits:   reasons.append(f"GO mentions {go_hits[:2]}")
        if kw_hits:   reasons.append(f"keyword {kw_hits}")
        if name_hits: reasons.append(f"name mentions {name_hits[:2]}")
        return ("KEEP_ANNOTATION", "; ".join(reasons))

    if not ann["rxn_chebis"] and not ann["ligand_chebis"]:
        return ("DROP_STALE", "no current Rhea reactions or ligand annotations at all")
    return ("DROP_STALE",
            f"current Rhea reactions/ligands do not contain any sterane-passing ChEBI "
            f"(has {len(ann['rxn_chebis'])} rxn ChEBIs, {len(ann['ligand_chebis'])} ligand ChEBIs, none pass)")


def main() -> int:
    print(f"Loading {PROTEINS.name}...")
    p = pd.read_csv(PROTEINS, low_memory=False)
    p["sequence_source"] = p["sequence_source"].fillna("").astype(str)
    rhea_entries = p[p["sequence_source"].str.startswith("Rhea catalytic")].copy()
    print(f"  {len(rhea_entries):,} Rhea-source entries to re-audit")
    accs = rhea_entries["accession"].astype(str).tolist()

    print(f"Loading sterane-passing ChEBI set from {MOLECULES.name}...")
    sterane = load_sterane_chebis()
    print(f"  {len(sterane):,} sterane-passing ChEBI ids")

    cache = {}
    if CACHE.exists():
        try:
            cache = json.loads(CACHE.read_text())
            print(f"  loaded {len(cache):,} cached entries")
        except Exception:
            pass

    remaining = [a for a in accs if a not in cache]
    print(f"\nFetching {len(remaining):,} entries in batches of {BATCH_SIZE}...")
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

    # Direct-entry fallback for TrEMBL entries the bulk endpoint skipped
    missing = [a for a in accs if a not in cache]
    if missing:
        print(f"\nRetrying {len(missing)} via direct endpoint...")
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
    print(f"\nCache size: {len(cache):,} / {len(accs):,}")

    # Classify — with GO/keyword rescue using atlas-stored annotations
    print("\nClassifying (with GO/keyword rescue)...")
    for c in ["go_labels", "keyword_ids", "protein_names"]:
        if c in rhea_entries.columns:
            rhea_entries[c] = rhea_entries[c].fillna("").astype(str)
        else:
            rhea_entries[c] = ""

    rows = []
    missing_ct = 0
    for _, r in rhea_entries.iterrows():
        acc = r["accession"]
        ann = cache.get(acc)
        if ann is None:
            missing_ct += 1
            rows.append({
                "accession": acc,
                "entry_name": r["entry_name"],
                "protein_names": str(r["protein_names"])[:120],
                "organism": r["organism"],
                "decision": "MISSING",
                "reason": "UniProt fetch failed",
                "current_rhea_ids": "",
                "current_rxn_chebis": "",
                "current_ligand_chebis": "",
            })
            continue
        decision, reason = classify(
            ann, sterane,
            go_labels=r["go_labels"],
            keyword_ids=r["keyword_ids"],
            protein_names=r["protein_names"],
        )
        rows.append({
            "accession": acc,
            "entry_name": r["entry_name"],
            "protein_names": str(r["protein_names"])[:120],
            "organism": r["organism"],
            "decision": decision,
            "reason": reason,
            "current_rhea_ids": ";".join(ann["rhea_ids"]),
            "current_rxn_chebis": ";".join(ann["rxn_chebis"]),
            "current_ligand_chebis": ";".join(ann["ligand_chebis"]),
        })

    audit = pd.DataFrame(rows)
    audit.to_csv(OUT_TSV, sep="\t", index=False)

    # Summary
    counts = audit["decision"].value_counts()
    lines = []
    lines.append("=" * 76)
    lines.append(f"RHEA RE-AUDIT — {len(audit):,} Rhea-source entries")
    lines.append("=" * 76)
    lines.append("")
    lines.append("Decision counts:")
    for d, n in counts.items():
        lines.append(f"  {n:>6,}  {d}")
    lines.append("")
    n_drop = int((audit["decision"] == "DROP_STALE").sum())
    n_keep = int(audit["decision"].isin(["KEEP_RHEA", "KEEP_LIGAND"]).sum())
    lines.append(f"Recommended action:  drop {n_drop:,} stale entries → "
                 f"Rhea corpus shrinks from {len(audit):,} → {len(audit)-n_drop:,}")
    lines.append("")
    for d in ("KEEP_RHEA", "KEEP_LIGAND", "DROP_STALE", "MISSING"):
        sub = audit[audit["decision"] == d]
        if not len(sub):
            continue
        lines.append(f"--- Top protein-name substrings in {d} ({len(sub):,}) ---")
        first = sub["protein_names"].astype(str).str.split("(").str[0].str.strip().str[:70]
        for pn, n in first.value_counts().head(12).items():
            lines.append(f"    {n:>5}  {pn}")
        lines.append("")

    OUT_TXT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_TSV.name}, {OUT_TXT.name}")

    print("\n=== Decision summary ===")
    for d, n in counts.items():
        print(f"  {n:>6,}  {d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
