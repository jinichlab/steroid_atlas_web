"""Pre-render 2D structure PNGs from SMILES with RDKit.

Reads:   ../steroid_atlas/data/molecules.csv (677 rows w/ SMILES + ChEBI)
         ../steroid_atlas/data/natural_synthetic_steroids.csv
Writes:  public/atlas/structures/<safe_id>.png

`<safe_id>` = urllib-safe slug of the compound name (lower-cased). Also
writes an alias index `public/atlas/structures/index.json` that maps
both compound-name (lowercased) and `chebi:<id>` → PNG filename, so the
web app can look up a structure by either key exactly like the marimo
`structure_cache`.

Requires: rdkit (via miniconda), pandas. Run under:
    LD_LIBRARY_PATH=~/miniconda3/lib ~/miniconda3/bin/python3 scripts/build_structures.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError as e:
    raise SystemExit("rdkit not installed. Run under the conda python.") from e

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT.parent / "steroid_atlas"
OUT = ROOT / "public" / "atlas" / "structures"
OUT.mkdir(parents=True, exist_ok=True)

_SLUG_RE = re.compile(r"[^a-z0-9]+")


def slugify(s: str) -> str:
    s = _SLUG_RE.sub("-", str(s).lower()).strip("-")
    return s[:80] or "unnamed"


def render(smiles: str, out_path: Path, size: int = 260) -> bool:
    if not smiles:
        return False
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    img = Draw.MolToImage(m, size=(size, size))
    img.save(out_path, "PNG", optimize=True)
    return True


def main() -> int:
    index: dict[str, str] = {}
    n_ok = 0
    n_skip = 0

    for csv, cols in (
        (SRC / "data" / "molecules.csv",
         {"name": "compound_name", "chebi": "chebi_id", "smiles": "smiles"}),
        (SRC / "data" / "natural_synthetic_steroids.csv",
         {"name": "Compound Name", "chebi": "ChEBI ID", "smiles": "SMILES"}),
    ):
        if not csv.exists():
            continue
        print(f"Reading {csv.name}...")
        df = pd.read_csv(csv, low_memory=False)
        for _, row in df.iterrows():
            name = str(row.get(cols["name"], "") or "").strip()
            chebi = str(row.get(cols["chebi"], "") or "").strip()
            smiles = str(row.get(cols["smiles"], "") or "").strip()
            if not name and not chebi:
                n_skip += 1; continue
            slug = slugify(name or f"chebi{chebi}")
            fname = f"{slug}.png"
            out_path = OUT / fname
            if not out_path.exists():
                if not render(smiles, out_path):
                    n_skip += 1; continue
            n_ok += 1
            if name:
                index[name.lower()] = fname
            if chebi:
                index[f"chebi:{chebi}"] = fname
        print(f"  running total: {n_ok:,} rendered, {n_skip:,} skipped")

    idx_path = OUT.parent / "structures_index.json"
    idx_path.write_text(json.dumps(index), encoding="utf-8")
    print(f"\nDone. Rendered {n_ok:,} structures, {n_skip:,} skipped.")
    print(f"  images:  {OUT}  ({sum(1 for _ in OUT.iterdir())} files)")
    print(f"  index :  {idx_path}  ({len(index):,} keys)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
