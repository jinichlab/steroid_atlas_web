"""Convert Steroid Atlas CSVs into optimized JSON for the Next.js app.

Reads:
    ../steroid_atlas/data/proteins.csv
    ../steroid_atlas/data/molecules.csv
    ../steroid_atlas/data/natural_synthetic_steroids.csv
    ../steroid_atlas/analysis/cluster_fingerprints.tsv
    ../steroid_atlas/analysis/cluster_go_top_terms.tsv (optional)

Writes:
    public/atlas/proteins.json          slim per-protein rows, sorted by cluster
    public/atlas/molecules.json         per-molecule rows with SMILES
    public/atlas/natsyn.json            natural + synthetic entries
    public/atlas/clusters.json          per-cluster metadata (name, GO, size)
    public/atlas/summary.json           top-level stats for the landing page

The JSON is chunked and slim (numeric fields as short floats, no huge free-text
where possible) so the initial page load stays snappy. Structure images are
pre-rendered separately if RDKit is installed (via build_structure_images.py).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT.parent / "steroid_atlas"
OUT = ROOT / "public" / "atlas"
OUT.mkdir(parents=True, exist_ok=True)

# ── PROTEINS ─────────────────────────────────────────────────────────────
print("Reading proteins.csv...")
prot = pd.read_csv(SRC / "data" / "proteins.csv", low_memory=False)
print(f"  {len(prot):,} proteins")

_PROT_COLS = [
    "accession", "entry_name", "protein_names", "gene_names", "organism",
    "length_aa", "ec_numbers", "rhea_reactions",
    "interacting_chebi_ids", "interacting_compounds",
    "go_ids", "go_labels", "keyword_ids", "keyword_labels",
    "pubmed_ids", "pubmed_count",
    "umap_1", "umap_2", "cluster",
    "is_literature_recruited", "paper_url",
    "audit_decision",
]
_kept = [c for c in _PROT_COLS if c in prot.columns]
prot = prot[_kept].fillna("")
# Round UMAP to 3 decimals — nothing else needs more than that on screen.
for _c in ("umap_1", "umap_2"):
    if _c in prot.columns:
        prot[_c] = prot[_c].astype(float).round(3)
prot["cluster"] = prot["cluster"].astype(int)
prot["is_literature_recruited"] = prot["is_literature_recruited"].astype(int)
prot["length_aa"] = prot["length_aa"].astype(int)

(OUT / "proteins.json").write_text(
    prot.to_json(orient="records", double_precision=3), encoding="utf-8"
)
print(f"  wrote proteins.json ({(OUT / 'proteins.json').stat().st_size / 1e6:.1f} MB)")

# ── MOLECULES ─────────────────────────────────────────────────────────────
print("Reading molecules.csv...")
mol = pd.read_csv(SRC / "data" / "molecules.csv", low_memory=False)
print(f"  {len(mol):,} molecules")

_MOL_COLS = [
    "compound_name", "chebi_id", "smiles",
    "umap_1", "umap_2", "cluster",
    "is_literature_recruited", "paper_url",
    "interacting_protein_accessions",
]
mol = mol[[c for c in _MOL_COLS if c in mol.columns]].fillna("")
for _c in ("umap_1", "umap_2"):
    if _c in mol.columns:
        mol[_c] = mol[_c].astype(float).round(3)
if "cluster" in mol.columns:
    mol["cluster"] = mol["cluster"].astype(int)
if "is_literature_recruited" in mol.columns:
    mol["is_literature_recruited"] = mol["is_literature_recruited"].astype(int)

(OUT / "molecules.json").write_text(
    mol.to_json(orient="records", double_precision=3), encoding="utf-8"
)
print(f"  wrote molecules.json ({(OUT / 'molecules.json').stat().st_size / 1e6:.1f} MB)")

# ── NATURAL + SYNTHETIC ───────────────────────────────────────────────────
print("Reading natural_synthetic_steroids.csv...")
natsyn = pd.read_csv(SRC / "data" / "natural_synthetic_steroids.csv", low_memory=False)
print(f"  {len(natsyn):,} entries")
if "UMAP_1" in natsyn.columns:
    natsyn = natsyn.rename(columns={
        "UMAP_1": "umap_1", "UMAP_2": "umap_2",
        "SMILES": "smiles", "Compound Name": "compound_name",
        "ChEBI ID": "chebi_id", "Entry": "protein_entries",
        "clusters": "cluster",
    })
_NS_COLS = ["compound_name", "smiles", "chebi_id", "umap_1", "umap_2",
            "cluster", "protein_entries"]
natsyn = natsyn[[c for c in _NS_COLS if c in natsyn.columns]].fillna("")
for _c in ("umap_1", "umap_2"):
    if _c in natsyn.columns:
        natsyn[_c] = natsyn[_c].astype(float).round(3)

(OUT / "natsyn.json").write_text(
    natsyn.to_json(orient="records", double_precision=3), encoding="utf-8"
)
print(f"  wrote natsyn.json ({(OUT / 'natsyn.json').stat().st_size / 1e6:.1f} MB)")

# ── CLUSTERS: PROTEIN ─────────────────────────────────────────────────────
print("Reading cluster_fingerprints.tsv (protein clusters)...")
fp = pd.read_csv(SRC / "analysis" / "cluster_fingerprints.tsv",
                 sep="\t", low_memory=False)
clusters = fp[["cluster", "size", "dominant_stem", "dominant_go"]].copy()
clusters.columns = ["id", "n", "name", "top_go"]
clusters["id"] = clusters["id"].astype(int)
clusters["n"] = clusters["n"].astype(int)
clusters["display_id"] = clusters["id"] + 1  # 1-based for display
clusters = clusters.sort_values("id")

(OUT / "protein_clusters.json").write_text(
    clusters.to_json(orient="records"), encoding="utf-8"
)
print(f"  wrote protein_clusters.json ({(OUT / 'protein_clusters.json').stat().st_size / 1e3:.1f} KB)")

# ── CLUSTERS: MOLECULE (steroid centric) ──────────────────────────────────
# Molecules have their own clusters (many fewer than proteins). Label each
# with a family stem inferred from the top compound names using steroid
# vocabulary — same idea as the protein stem, but keyword-based.
print("Building molecule cluster metadata (steroid-family classifier)...")

_FAMILY_PATTERNS = [
    ("Bile acids",       r"(chol|chenodeoxy|deoxychol|lithochol|ursodeoxy|hyodeoxy|muricho|hyocho)"),
    ("Estrogens",        r"(estradiol|estrone|estriol|estetrol|\bestr)"),
    ("Androgens",        r"(testoster|androst|dihydrotestos|dehydroepi|androsterone)"),
    ("Progestins",       r"(progest|\bpregn|pregnan)"),
    ("Corticosteroids",  r"(cort|dexameth|prednis|aldos|betameth)"),
    ("Sterols",          r"(sterol|cholesterol|ergost|sitoster|lanoster|desmoster|campester|stigmaster)"),
    ("Ecdysteroids",     r"(ecdys|ponaster)"),
    ("Brassinosteroids", r"(brassin|teasteron|castasteron|typhaster|katasteron|cathasteron)"),
    ("Vitamin D",        r"(cholecalcif|calcitriol|calciol|ergocalcif|vitamin d)"),
]
import re as _re

def _family_of(name: str) -> str:
    _low = str(name).lower()
    for _lbl, _pat in _FAMILY_PATTERNS:
        if _re.search(_pat, _low):
            return _lbl
    return ""

mol_clusters_out = []
if "cluster" in mol.columns:
    for _cid, _sub in mol.groupby("cluster"):
        _cid_int = int(_cid)
        _n = len(_sub)
        _fams = {}
        for _nm in _sub["compound_name"].astype(str):
            _f = _family_of(_nm)
            if _f:
                _fams[_f] = _fams.get(_f, 0) + 1
        if _fams:
            _top_fam = sorted(_fams, key=lambda k: -_fams[k])[0]
            _cov = round(_fams[_top_fam] / _n * 100)
            _name = f"{_top_fam} · {_cov}%"
        else:
            _name = "Mixed steroids"
        mol_clusters_out.append({
            "id": _cid_int,
            "display_id": _cid_int + 1,
            "n": _n,
            "name": _name,
            "top_go": "",
        })

(OUT / "molecule_clusters.json").write_text(
    json.dumps(sorted(mol_clusters_out, key=lambda r: r["id"])),
    encoding="utf-8",
)
print(f"  wrote molecule_clusters.json "
      f"({len(mol_clusters_out)} clusters, "
      f"{(OUT / 'molecule_clusters.json').stat().st_size / 1e3:.1f} KB)")

# ── CLUSTERS: NATSYN (2 groups only) ──────────────────────────────────────
natsyn_clusters = []
if "cluster" in natsyn.columns:
    for _i, (_cid, _sub) in enumerate(natsyn.groupby("cluster")):
        natsyn_clusters.append({
            "id": _i,
            "display_id": _i + 1,
            "n": int(len(_sub)),
            "name": str(_cid),
            "top_go": "",
        })
(OUT / "natsyn_clusters.json").write_text(
    json.dumps(natsyn_clusters), encoding="utf-8"
)
print(f"  wrote natsyn_clusters.json ({len(natsyn_clusters)} clusters)")

# Keep the old clusters.json as an alias for the protein clusters
# (existing pages still import it).
(OUT / "clusters.json").write_text(
    clusters.to_json(orient="records"), encoding="utf-8"
)

# ── SUMMARY ───────────────────────────────────────────────────────────────
summary = {
    "n_proteins": int(len(prot)),
    "n_molecules": int(len(mol)),
    "n_natsyn": int(len(natsyn)),
    "n_clusters": int(len(clusters)),
    "n_new": int((prot["is_literature_recruited"] == 1).sum()),
}
(OUT / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
print(f"summary: {summary}")

print("\nDone — outputs in", OUT)
