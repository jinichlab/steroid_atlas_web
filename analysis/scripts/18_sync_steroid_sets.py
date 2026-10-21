"""Sync the natural steroid set with molecules.csv and refit the
natural + synthetic UMAP.

Two inconsistencies fixed here:

  1. molecules.csv listed 5beta-dihydrotestosterone (ChEBI 2150) twice — once
     from Rhea (with its interacting proteins) and once as a literature
     recruitment (no proteins). The pair is merged into the Rhea row, which
     inherits the literature flag + paper_url.  677 → 676 rows.

  2. The "natural" rows of natural_synthetic_steroids.csv were built from an
     older molecule list: they lacked the 91 distinct literature-recruited
     compounds and still held the 4 benzo[c]phenanthridine alkaloids that
     v0.2.1 dropped from molecules.csv.  They are rebuilt from the 676
     molecules; the 2,300 synthetic rows are kept unchanged.  2,889 → 2,976.

Because the original natural + synthetic UMAP reducer was never saved, the
91 new compounds can't be projected into it — the whole map is refit:

  - RDKit Morgan fingerprints, radius 2, 2048 bits (ECFP4)
  - umap.UMAP(n_neighbors=15, min_dist=0.1, metric="jaccard",
              random_state=42), standard initialization
    (1 − Jaccard on bit vectors = Tanimoto similarity)

The fitted reducer is saved so future compounds can be added with
.transform() instead of another refit.

Inputs are read from the *.bak originals (created on first run), so the
script is idempotent: re-running it reproduces the same outputs.

Reads:  ../data/molecules.csv, ../data/natural_synthetic_steroids.csv
        (their *.bak copies once those exist)
Writes: ../data/molecules.csv                     676 rows
        ../data/natural_synthetic_steroids.csv    2,976 rows, new UMAP coords
        ../analysis/natsyn_umap.joblib            fitted UMAP reducer
        ../analysis/natsyn_umap_report.txt        human-readable summary
"""
from __future__ import annotations

import shutil
import sys
from importlib.metadata import version
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

# umap's __init__ tries to import TensorFlow for ParametricUMAP (unused here);
# on some macOS installs that import deadlocks. Block it so umap falls back.
sys.modules["tensorflow"] = None
import umap  # noqa: E402

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
MOLECULES = ROOT / "data" / "molecules.csv"
NATSYN = ROOT / "data" / "natural_synthetic_steroids.csv"
OUT_MODEL = output("natsyn_umap.joblib")
OUT_REPORT = output("natsyn_umap_report.txt")

ALKALOID_CHEBIS = {"17183", "17209", "17031", "17789"}  # dropped in v0.2.1
FP_RADIUS, FP_BITS = 2, 2048
UMAP_PARAMS = dict(n_neighbors=15, min_dist=0.1, metric="jaccard", random_state=42)


def read_original(path: Path) -> pd.DataFrame:
    """Back up `path` once, then always read the pristine backup.

    dtype=str + keep_default_na=False round-trips both CSVs byte-for-byte, so
    rows we don't touch are written back unchanged.
    """
    bak = path.with_suffix(path.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(path, bak)
    return pd.read_csv(bak, dtype=str, keep_default_na=False)


def merge_duplicate_molecules(mol: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Collapse rows sharing (name, ChEBI) into the row that has proteins."""
    key = mol["compound_name"].str.strip().str.lower() + "|" + mol["chebi_id"]
    notes, drop = [], []
    for _, grp in mol.groupby(key, sort=False):
        if len(grp) < 2:
            continue
        # Keeper: the Rhea row (non-empty protein list); ties → first row.
        has_prot = grp["interacting_protein_accessions"].str.strip() != ""
        keep = grp[has_prot].index[0] if has_prot.any() else grp.index[0]
        if (grp["is_literature_recruited"] == "1").any():
            mol.at[keep, "is_literature_recruited"] = "1"
            urls = [u for u in grp["paper_url"] if u]
            mol.at[keep, "paper_url"] = " | ".join(dict.fromkeys(urls))
        drop += [i for i in grp.index if i != keep]
        notes.append(
            f"{grp['compound_name'].iloc[0].strip()} (ChEBI {grp['chebi_id'].iloc[0]}): "
            f"merged {len(grp)} rows into row {keep}"
        )
    return mol.drop(index=drop).reset_index(drop=True), notes


def natural_rows_from(mol: pd.DataFrame) -> pd.DataFrame:
    """natural_synthetic_steroids.csv rows (natural group) for every molecule."""
    return pd.DataFrame({
        "SMILES": mol["smiles"],
        "Compound Name": mol["compound_name"],
        "UMAP_1": "",
        "UMAP_2": "",
        "clusters": "natural",
        "ChEBI ID": [f"['{c}']" if c else "[]" for c in mol["chebi_id"]],
        "Entry": [e if e.strip() else "[]" for e in mol["interacting_protein_accessions"]],
    })


def canonical(smi: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


def main() -> int:
    lines: list[str] = []
    log = lambda s="": (print(s), lines.append(s))  # noqa: E731

    # ── 1. molecules.csv ──────────────────────────────────────────────────
    mol = read_original(MOLECULES)
    n_mol_before = len(mol)
    mol, merge_notes = merge_duplicate_molecules(mol)
    assert len(mol) == 676, f"expected 676 molecules, got {len(mol)}"
    key = mol["compound_name"].str.strip().str.lower() + "|" + mol["chebi_id"]
    assert not key.duplicated().any(), "duplicate (name, ChEBI) rows remain"

    # ── 2. natural + synthetic ───────────────────────────────────────────
    ns = read_original(NATSYN)
    old_nat = ns[ns["clusters"] == "natural"]
    syn = ns[ns["clusters"] != "natural"].reset_index(drop=True)
    nat = natural_rows_from(mol)
    ns = pd.concat([syn, nat], ignore_index=True)  # synthetic first, as before

    old_chebis = {c.strip("[]'\" ") for c in old_nat["ChEBI ID"]}
    new_chebis = {c.strip("[]'\" ") for c in nat["ChEBI ID"]}
    dropped_alkaloids = sorted(ALKALOID_CHEBIS & old_chebis - new_chebis)
    assert not (ALKALOID_CHEBIS & new_chebis), "alkaloid ChEBIs still present"
    assert (ns["clusters"] == "natural").sum() == 676
    assert len(syn) == 2300

    # Report-only overlap check: same structure in both groups.
    nat_can = {canonical(s) for s in nat["SMILES"]}
    overlap = syn[[canonical(s) in nat_can for s in syn["SMILES"]]]

    # ── 3. refit UMAP ────────────────────────────────────────────────────
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    fps = np.array(
        [gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in ns["SMILES"]],
        dtype=bool,
    )
    reducer = umap.UMAP(**UMAP_PARAMS)
    coords = reducer.fit_transform(fps)
    assert np.isfinite(coords).all(), "non-finite UMAP coordinates"
    ns["UMAP_1"] = [repr(float(x)) for x in coords[:, 0]]
    ns["UMAP_2"] = [repr(float(y)) for y in coords[:, 1]]

    # ── 4. write ─────────────────────────────────────────────────────────
    mol.to_csv(MOLECULES, index=False)
    ns.to_csv(NATSYN, index=False)
    joblib.dump(reducer, OUT_MODEL)

    log("=" * 74)
    log("STEROID SET SYNC + NATURAL/SYNTHETIC UMAP REFIT")
    log("=" * 74)
    log(f"molecules.csv:   {n_mol_before} → {len(mol)} rows")
    for n in merge_notes:
        log(f"  {n}")
    log()
    log(f"natural_synthetic_steroids.csv:   {len(old_nat) + len(syn):,} → {len(ns):,} rows")
    log(f"  natural:    {len(old_nat)} → {len(nat)}  (rebuilt from molecules.csv)")
    log(f"  synthetic:  {len(syn):,} (unchanged)")
    log(f"  alkaloid ChEBIs removed from natural: {', '.join(dropped_alkaloids) or 'none'}")
    log(f"  natural compounds added: {len(new_chebis - old_chebis)} distinct ChEBI ids "
        f"(+ compounds without a ChEBI id)")
    log(f"  synthetic rows whose canonical SMILES equals a natural one: {len(overlap)}")
    for _, r in overlap.head(20).iterrows():
        log(f"    {r['Compound Name'].strip()}")
    log()
    log("UMAP refit (all rows, standard initialization):")
    log(f"  fingerprints:  Morgan radius {FP_RADIUS}, {FP_BITS} bits")
    log(f"  params:        {UMAP_PARAMS}")
    log(f"  model:         {OUT_MODEL.name}")
    log("  versions:      " + ", ".join(
        f"{p} {version(p)}" for p in ("umap-learn", "rdkit", "numpy", "scikit-learn")
    ))
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
