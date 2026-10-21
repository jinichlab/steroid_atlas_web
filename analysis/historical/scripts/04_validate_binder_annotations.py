"""Audit the 2,480 orphan-looking entries in the atlas.

An "orphan" is an atlas entry with no Rhea reactions, no interacting ChEBI ids,
no interacting compounds, and no EC numbers. These slipped in via the binder
search (KW-0754 or GO:0005496 including its GO-hierarchy children such as
GO:0005499 vitamin-D binding).

Our operational steroid definition requires a sterane 4-ring core. Vitamin D
(secosteroid, ring B open) fails that test, so entries whose only steroid-family
association is vitamin-D binding should be dropped from the atlas.

This script fetches the actual UniProt annotations for each orphan (GO terms,
keywords, ligand ChEBI ids) and classifies each entry:

  KEEP        - has a direct sterane-passing ChEBI ligand (in molecules.csv)
                OR has KW-0754 / GO:0005496 with at least one non-vitD steroid GO child
  DROP_VITD   - only steroid-family term is GO:0005499 (vitamin-D binding)
  DROP_NONE   - no steroid keyword or GO term present at all (likely stale annotation)
  REVIEW      - has generic GO:0005496 but no specific ligand named

Reads:   ../data/proteins.csv, ../data/molecules.csv
Writes:  ../analysis/binder_audit.tsv           per-entry decision table
         ../analysis/binder_audit_summary.txt   human-readable summary
         ../analysis/binder_audit_uniprot.json  raw UniProt response cache
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
OUT_TSV = output("binder_audit.tsv", historical=True)
OUT_REPORT = output("binder_audit_summary.txt", historical=True)
OUT_CACHE = output("binder_audit_uniprot.json", historical=True)

UNIPROT_ACCESSIONS_URL = "https://rest.uniprot.org/uniprotkb/accessions"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
BATCH_SIZE = 300
FIELDS = ",".join([
    "accession",
    "protein_name",
    "keyword",
    "go_id",
    "cc_function",
    "cc_catalytic_activity",
    "ft_binding",
    "cc_activity_regulation",
    "cc_pathway",
])

# GO term references
GO_STEROID_BINDING     = "GO:0005496"   # umbrella "steroid binding"
GO_VITAMIN_D_BINDING   = "GO:0005499"   # vitamin D binding (secosteroid — DROP)
# Non-vitD children of GO:0005496 that are unambiguously sterane
GO_STERANE_CHILDREN = {
    "GO:0005497",   # androgen binding
    "GO:0005500",   # sterol carrier activity → cholesterol/oxysterols
    "GO:0005539",   # glycosaminoglycan binding (not relevant, exclude)
    "GO:0005543",   # phospholipid binding (not relevant)
    "GO:0008144",   # drug binding (too generic)
    "GO:0015485",   # cholesterol binding
    "GO:0032934",   # sterol binding
    "GO:0032810",   # sterol response element binding
    "GO:0034185",   # apolipoprotein binding
    "GO:0034186",   # apolipoprotein A-I binding
    "GO:0038181",   # bile acid receptor activity
    "GO:0038186",   # bile acid binding
    "GO:0043178",   # alcohol binding (too generic)
    "GO:0050544",   # arachidonate ω-hydroxylase (not steroid)
    "GO:0070330",   # aromatase activity
    "GO:0032052",   # bile acid binding
    "GO:0004769",   # steroid delta-isomerase
    "GO:0033218",   # amide binding (not steroid)
    "GO:1990239",   # steroid hormone binding
    "GO:0003707",   # steroid hormone receptor activity
    "GO:0005502",   # 11-cis retinal binding (NOT steroid — exclude)
    "GO:0005501",   # retinoid binding (NOT steroid — exclude)
    "GO:0016918",   # retinal binding (NOT steroid — exclude)
    "GO:0043404",   # corticotropin-releasing hormone (not steroid)
    "GO:1902121",   # lithocholate binding
    "GO:0043325",   # phosphatidylinositol-3,4-bisphosphate binding
    "GO:0008395",   # steroid hydroxylase activity
    "GO:0016125",   # sterol metabolic process
    "GO:0016126",   # sterol biosynthetic process
}
# Whitelist of sterane-associated GO terms (curated subset of the above)
GO_STERANE_WHITELIST = {
    "GO:0005497",   # androgen binding
    "GO:0015485",   # cholesterol binding
    "GO:0032934",   # sterol binding
    "GO:0032810",   # SREBP response element binding
    "GO:0038181",   # bile acid receptor
    "GO:0038186",   # bile acid binding
    "GO:1902121",   # lithocholate binding
    "GO:0032052",   # bile acid binding
    "GO:0070330",   # aromatase (metabolizes androgens/estrogens)
    "GO:0004769",   # steroid delta-isomerase
    "GO:1990239",   # steroid hormone binding
    "GO:0003707",   # steroid hormone receptor
    "GO:0008395",   # steroid hydroxylase
    "GO:0016125",   # sterol metabolic process
    "GO:0016126",   # sterol biosynthetic process
    "GO:0034185",   # apolipoprotein binding
    "GO:0005500",   # sterol carrier activity
}

KW_STEROID_BINDING     = "KW-0754"      # UniProt keyword: steroid-binding
KW_STEROID_HORMONE_RECEPTOR = "KW-0675" # steroid hormone receptor


def load_orphans() -> pd.DataFrame:
    df = pd.read_csv(PROTEINS, low_memory=False)
    orphan = df[
        (df['rhea_reactions'].astype(str).isin(['nan', ''])) &
        (df['interacting_chebi_ids'].astype(str).isin(['nan', ''])) &
        (df['interacting_compounds'].astype(str).isin(['nan', ''])) &
        (df['ec_numbers'].astype(str).isin(['nan', '']))
    ][['accession', 'entry_name', 'protein_names', 'gene_names', 'organism',
       'sequence_source']].copy()
    return orphan.reset_index(drop=True)


def load_sterane_chebi_ids() -> set[str]:
    """molecules.csv stores chebi_id as float64 (e.g. 17263.0). Normalize to 'CHEBI:17263'."""
    mol = pd.read_csv(MOLECULES, low_memory=False)
    ids = set()
    for x in mol['chebi_id'].dropna():
        try:
            n = int(float(x))
            ids.add(f'CHEBI:{n}')
        except (ValueError, TypeError):
            s = str(x).strip().upper()
            if s.startswith('CHEBI:'):
                ids.add(s)
    return ids


# Known vitamin-D / secosteroid ChEBI ids (ring B open — fail sterane test).
# If an entry's ONLY named ligand comes from this set, it should be dropped.
KNOWN_VITD_CHEBIS = {
    'CHEBI:27300',   # vitamin D
    'CHEBI:28940',   # cholecalciferol (vitamin D3)
    'CHEBI:28935',   # ergocalciferol (vitamin D2)
    'CHEBI:17933',   # 25-hydroxycholecalciferol / calcidiol
    'CHEBI:17823',   # 1α,25-dihydroxycholecalciferol / calcitriol
    'CHEBI:73558',   # 1α,25-dihydroxyvitamin D3 (alt id)
    'CHEBI:41927',   # 1α-hydroxycholecalciferol / alfacalcidol
    'CHEBI:34026',   # 24,25-dihydroxycholecalciferol
    'CHEBI:71305',   # 25-hydroxyvitamin D2
    'CHEBI:29131',   # 1α,25-dihydroxyergocalciferol
    'CHEBI:73559',   # calcipotriol
    'CHEBI:63099',   # paricalcitol
    'CHEBI:34918',   # secalciferol
    'CHEBI:73561',   # doxercalciferol
    'CHEBI:33487',   # 1α-hydroxyvitamin D3 alt
    'CHEBI:78272',   # 3-epi-1α,25-dihydroxyvitamin D3
}


def fetch_batch(accessions: list[str], retries: int = 3) -> list[dict]:
    params = {"accessions": ",".join(accessions), "format": "json", "fields": FIELDS}
    for attempt in range(retries):
        try:
            r = requests.get(UNIPROT_ACCESSIONS_URL, params=params, timeout=60)
            if r.status_code == 200:
                return r.json().get("results", [])
            print(f"  ! HTTP {r.status_code} on batch of {len(accessions)} (attempt {attempt+1})")
            time.sleep(2 ** attempt)
        except requests.RequestException as e:
            print(f"  ! {e} (attempt {attempt+1})")
            time.sleep(2 ** attempt)
    return []


def fetch_single(accession: str, retries: int = 2) -> dict | None:
    """Direct entry endpoint — works for TrEMBL entries the bulk endpoint skips."""
    params = {"format": "json", "fields": FIELDS}
    url = f"{UNIPROT_ENTRY_URL}/{accession}.json"
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
            time.sleep(2 ** attempt)
        except requests.RequestException:
            time.sleep(2 ** attempt)
    return None


def extract_annotations(entry: dict) -> dict:
    """Extract the fields we need from one UniProt JSON entry."""
    acc = entry.get("primaryAccession", "")
    keywords = [k.get("id", "") for k in entry.get("keywords", []) if k.get("id")]
    go_ids = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GO":
            go_ids.append(xref.get("id", ""))

    # Ligand ChEBI ids come from FT BINDING features and CC CATALYTIC ACTIVITY xrefs.
    def _norm_chebi(cid: str) -> str | None:
        """Normalize any ChEBI id form to 'CHEBI:<n>' (fixes double-prefix like 'CHEBI:CHEBI:17823')."""
        s = str(cid).strip().upper()
        while s.startswith("CHEBI:"):
            s = s[len("CHEBI:"):]
        if s.isdigit():
            return f"CHEBI:{s}"
        return None

    ligand_chebis = []
    for feat in entry.get("features", []):
        if feat.get("type") != "Binding site":
            continue
        lig = feat.get("ligand", {})
        norm = _norm_chebi(lig.get("id", ""))
        if norm:
            ligand_chebis.append(norm)

    ft_text_bits = []
    for c in entry.get("comments", []):
        if c.get("commentType") in {"FUNCTION", "CATALYTIC ACTIVITY", "PATHWAY"}:
            for t in c.get("texts", []):
                ft_text_bits.append(t.get("value", "") or "")
            reaction = c.get("reaction", {})
            for xref in reaction.get("reactionCrossReferences", []) or []:
                if xref.get("database", "").upper() == "CHEBI":
                    norm = _norm_chebi(xref.get("id", ""))
                    if norm:
                        ligand_chebis.append(norm)

    return {
        "accession": acc,
        "keywords": keywords,
        "go_ids": go_ids,
        "ligand_chebis": sorted(set(ligand_chebis)),
        "free_text": " || ".join(ft_text_bits)[:2000],
    }


def classify(ann: dict, sterane_chebis: set[str], protein_name: str = "") -> tuple[str, str]:
    """Return (decision, reason).

    Ligand ChEBI checks take precedence over keyword/GO checks. A protein whose
    only named ligand is a secosteroid (vitamin D family) is dropped even if it
    carries KW-0675 or KW-0754, because those keywords are family-level curation
    and don't imply sterane-passing ligand.
    """
    kws = set(ann["keywords"])
    gos = set(ann["go_ids"])
    ligs = set(ann["ligand_chebis"])
    text = (ann["free_text"] or "").lower()
    pname_lower = (protein_name or "").lower()

    sterane_ligs = ligs & sterane_chebis
    vitd_ligs = ligs & KNOWN_VITD_CHEBIS

    # HARD KEEP: any named ligand is a sterane-passing steroid
    if sterane_ligs:
        return ("KEEP", f"binds sterane-passing ChEBI: {sorted(sterane_ligs)[:3]}")

    # HARD DROP: has vitamin-D-family ligand(s) and NO sterane-passing ligand
    if vitd_ligs and not sterane_ligs:
        return ("DROP_VITD", f"only vitamin-D-family ligands: {sorted(vitd_ligs)[:3]}")

    # Name-based hard drop for well-known secosteroid-only proteins.
    # The parenthesized VDR form ("Vitamin D (1,25-...) receptor") won't match the plain
    # substrings, so also detect names that start with "vitamin d" and contain "receptor".
    vdr_name_patterns = [
        "vitamin d3 receptor", "vitamin d receptor", "1,25-dihydroxyvitamin d3 receptor",
    ]
    is_vdr_variant = (
        any(p in pname_lower for p in vdr_name_patterns)
        or (pname_lower.startswith("vitamin d") and "receptor" in pname_lower)
    )
    if is_vdr_variant and not sterane_ligs:
        return ("DROP_VITD", "VDR — vitamin D nuclear receptor (secosteroid-only)")
    if "vitamin d 25-hydroxylase" in pname_lower and not sterane_ligs:
        return ("DROP_VITD", "vitamin D 25-hydroxylase acts on secosteroid substrate")
    if "vitamin d-binding protein" in pname_lower and not sterane_ligs:
        return ("DROP_VITD", "vitamin-D-binding protein (secosteroid carrier)")
    if "calbindin" in pname_lower and not sterane_ligs:
        return ("DROP_VITD", "calbindin — vitamin-D-induced Ca-binding protein, not a steroid binder")

    has_kw_steroid = (KW_STEROID_BINDING in kws) or (KW_STEROID_HORMONE_RECEPTOR in kws)
    has_go_steroid_umbrella = GO_STEROID_BINDING in gos
    has_go_vitd = GO_VITAMIN_D_BINDING in gos
    has_go_sterane_child = bool(gos & GO_STERANE_WHITELIST)

    # DROP: only steroid annotation is vitamin D
    if has_go_vitd and not (has_go_sterane_child or has_kw_steroid or has_go_steroid_umbrella):
        return ("DROP_VITD", "only steroid-family GO term is GO:0005499 (vitamin D binding)")
    if has_go_vitd and has_go_steroid_umbrella and not has_go_sterane_child and not has_kw_steroid:
        return ("DROP_VITD", "steroid GO umbrella + vitamin D child, no sterane-specific term")

    # KEEP: has explicit sterane-family GO child
    if has_go_sterane_child:
        matches = sorted(gos & GO_STERANE_WHITELIST)
        return ("KEEP", f"GO sterane-specific: {matches[:3]}")

    # KEEP: has UniProt steroid-hormone-receptor keyword (KW-0675) — specific
    if KW_STEROID_HORMONE_RECEPTOR in kws:
        return ("KEEP", "KW-0675 steroid hormone receptor")

    # KEEP but weak: has KW-0754 (steroid-binding keyword) — curated
    if KW_STEROID_BINDING in kws:
        return ("KEEP_WEAK", "KW-0754 steroid-binding keyword (curator-assigned)")

    # KEEP but weak: umbrella GO only, no vitD, no sterane child
    if has_go_steroid_umbrella and not has_go_vitd:
        return ("KEEP_WEAK", "GO:0005496 umbrella only (no vitD, no specific child)")

    # Nothing steroid at all — free-text fallback
    if not (has_kw_steroid or has_go_steroid_umbrella or has_go_vitd):
        vocab = ["steroid", "sterol", "estrogen", "estradiol", "testosterone",
                 "androgen", "progesterone", "cortisol", "corticoster", "aldosterone",
                 "bile acid", "bile-acid", "cholesterol", "pregnenolone", "ecdysone",
                 "cardiac glycoside", "ouabain", "bufalin", "digoxin"]
        if "vitamin d" in text and not any(v in text for v in vocab if v != "sterol"):
            return ("DROP_VITD", "free-text mentions vitamin D but no other steroid vocabulary")
        if any(v in text for v in vocab):
            return ("REVIEW", "free-text mentions steroid vocabulary but no formal annotation")

        # Rescue: TrEMBL orthologs with zero annotations but a clearly steroid-family name.
        # Only trust names that unambiguously identify sterane-binding proteins.
        name_rescues = [
            ("estrogen receptor",          "sterane — estradiol binds sterane"),
            ("androgen receptor",          "sterane — testosterone/DHT bind sterane"),
            ("glucocorticoid receptor",    "sterane — cortisol binds sterane"),
            ("progesterone receptor",      "sterane — progesterone binds sterane"),
            ("mineralocorticoid receptor", "sterane — aldosterone binds sterane"),
            ("g-protein coupled estrogen receptor", "sterane — GPER1 binds estradiol"),
            ("g protein-coupled estrogen receptor", "sterane — GPER1 binds estradiol"),
            ("bile acid receptor",         "sterane — bile acids have sterane core"),
            ("farnesoid",                  "sterane — FXR binds bile acids"),
            ("smoothened",                 "sterane — cholesterol-modulated"),
            ("patched",                    "sterane — cholesterol-modulated"),
            ("caveolin",                   "sterane — caveolin binds cholesterol"),
            ("apolipoprotein a-i", "sterane — APOA1 shuttles cholesterol"),
            ("apolipoprotein a1",  "sterane — APOA1 shuttles cholesterol"),
            ("ecdysone receptor",          "sterane — ecdysone has sterane core"),
            ("rar related orphan receptor a", "sterane — RORα binds cholesterol sulfate"),
            ("ror alpha",                  "sterane — RORα binds cholesterol sulfate"),
            ("insulin-induced gene",       "sterane — INSIG binds cholesterol"),
            ("oxysterol-binding",          "sterane — OSBP binds oxysterols"),
            ("srebp cleavage-activating",  "sterane — SCAP is a cholesterol sensor"),
            ("sterol regulatory element-binding protein cleavage-activating",
                                            "sterane — SCAP is a cholesterol sensor"),
            ("hedgehog",                   "sterane — cholesterol-modified ligand"),
            ("sodium/potassium-transporting atpase",
                                            "sterane — Na/K-ATPase α binds cardiac glycosides"),
            ("nfe2l1",                     "sterane — ER cholesterol/lipid sensor"),
            ("endoplasmic reticulum membrane sensor nfe2l1",
                                            "sterane — NFE2L1 cholesterol sensor"),
        ]
        for pat, why in name_rescues:
            if pat in pname_lower:
                return ("KEEP_BY_NAME", why)

        return ("DROP_NONE", "no steroid keyword, GO term, or ligand annotation")

    return ("REVIEW", "residual case — inspect manually")


def main() -> int:
    print(f"Loading orphans from {PROTEINS.name}...")
    orphans = load_orphans()
    print(f"  {len(orphans):,} orphan entries")

    print(f"Loading sterane-passing ChEBI set from {MOLECULES.name}...")
    sterane_chebis = load_sterane_chebi_ids()
    print(f"  {len(sterane_chebis):,} sterane-passing ChEBI ids")

    # Try to load cache
    all_anns: dict[str, dict] = {}
    if OUT_CACHE.exists():
        try:
            all_anns = json.loads(OUT_CACHE.read_text())
            print(f"  Loaded {len(all_anns):,} cached UniProt annotations from {OUT_CACHE.name}")
        except Exception:
            pass

    accs = orphans["accession"].tolist()
    remaining = [a for a in accs if a not in all_anns]
    print(f"\nFetching UniProt annotations for {len(remaining):,} orphans in batches of {BATCH_SIZE}...")

    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i:i + BATCH_SIZE]
        print(f"  batch {i//BATCH_SIZE + 1}/{(len(remaining)+BATCH_SIZE-1)//BATCH_SIZE}  ({len(batch)} accs)", flush=True)
        entries = fetch_batch(batch)
        for entry in entries:
            ann = extract_annotations(entry)
            if ann["accession"]:
                all_anns[ann["accession"]] = ann
        if (i // BATCH_SIZE) % 5 == 4:
            OUT_CACHE.write_text(json.dumps(all_anns))
        time.sleep(0.3)

    # Retry missing ones via direct entry endpoint (bulk /accessions skips TrEMBL)
    missing = [a for a in accs if a not in all_anns]
    if missing:
        print(f"\nRetrying {len(missing)} missing accessions via direct entry endpoint...")
        for j, acc in enumerate(missing):
            if j % 20 == 0:
                print(f"  {j}/{len(missing)}", flush=True)
            entry = fetch_single(acc)
            if entry is not None:
                ann = extract_annotations(entry)
                if ann["accession"]:
                    all_anns[ann["accession"]] = ann
                if ann["accession"] != acc:
                    all_anns[acc] = ann
            time.sleep(0.15)

    OUT_CACHE.write_text(json.dumps(all_anns))
    print(f"\n  Cached {len(all_anns):,} entries → {OUT_CACHE.name}")

    # Now classify
    print("\nClassifying...")
    rows = []
    missing = 0
    for _, r in orphans.iterrows():
        acc = r["accession"]
        ann = all_anns.get(acc)
        if ann is None:
            missing += 1
            rows.append({
                **r.to_dict(),
                "keywords": "", "go_ids": "", "ligand_chebis": "",
                "decision": "MISSING", "reason": "UniProt fetch failed",
            })
            continue
        decision, reason = classify(ann, sterane_chebis, r["protein_names"])
        rows.append({
            **r.to_dict(),
            "keywords": ";".join(ann["keywords"]),
            "go_ids": ";".join(ann["go_ids"]),
            "ligand_chebis": ";".join(ann["ligand_chebis"]),
            "decision": decision,
            "reason": reason,
        })
    if missing:
        print(f"  ! {missing} entries could not be fetched — marked MISSING")

    audit = pd.DataFrame(rows)
    audit.to_csv(OUT_TSV, sep="\t", index=False)

    # Report
    counts = audit["decision"].value_counts()
    lines = []
    lines.append(f"=== Orphan-entry audit ({len(audit):,} entries) ===\n")
    lines.append("An 'orphan' has no Rhea, no ChEBI, no compounds, and no EC in the atlas.")
    lines.append("Classification based on UniProt keywords, GO terms, and ligand ChEBI ids.\n")
    lines.append("Decision counts:")
    for d, n in counts.items():
        lines.append(f"  {n:>5}  {d}")
    lines.append("")

    # Top protein names per decision
    for d in ["KEEP", "KEEP_WEAK", "REVIEW", "DROP_VITD", "DROP_NONE", "MISSING"]:
        sub = audit[audit["decision"] == d]
        if not len(sub):
            continue
        lines.append(f"--- Top protein-name substrings in {d} ({len(sub)} entries) ---")
        first = sub["protein_names"].astype(str).str.split("(").str[0].str.strip().str[:70]
        for pn, n in first.value_counts().head(12).items():
            lines.append(f"  {n:>4}  {pn}")
        lines.append("")

    OUT_REPORT.write_text("\n".join(lines) + "\n")

    print("\n=== Decision summary ===")
    for d, n in counts.items():
        print(f"  {n:>5}  {d}")
    print(f"\nWrote {OUT_TSV.name}, {OUT_REPORT.name}, {OUT_CACHE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
