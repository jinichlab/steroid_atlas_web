"""High-dimensional k-means — the "full embedding" clustering behind the web
app's cluster toggle.

The atlas's published clusters are k-means on the 2D UMAP coordinates
(k = 82 proteins, k = 8 steroids). This script clusters the same entries in
the space the UMAP was computed FROM, with the same k, so the two can be
compared directly:

  - steroids: RDKit Morgan fingerprints (radius 2, 2048 bits) of
    data/molecules.csv SMILES — same generator as 18_sync_steroid_sets.py.
  - proteins: ProtT5 per-protein vectors (1,024-d) from
    data/embeddings/prott5.h5, L2-normalised so Euclidean k-means ranks
    neighbours like cosine similarity. Skipped (with a message) when the file
    is absent; proteins without a vector get cluster_hd = -1.

KMeans(k, random_state=42, n_init=10) everywhere. Agreement with the 2D
clusters is reported as the adjusted Rand index (ARI).

Reads:  ../../data/molecules.csv, ../../data/proteins.csv,
        ../../data/embeddings/prott5.h5 (optional)
Writes: ../results/molecule_clusters_hd.tsv   row-aligned with molecules.csv → cluster_hd
        ../results/protein_clusters_hd.tsv    accession → cluster_hd  (if .h5 present)
        ../results/highdim_kmeans_report.md
"""
from __future__ import annotations

import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402

RDLogger.DisableLog("rdApp.*")

MOLECULES = ROOT / "data" / "molecules.csv"
PROTEINS = ROOT / "data" / "proteins.csv"
EMBEDDINGS = ROOT / "data" / "embeddings" / "prott5.h5"
K_MOLECULES, K_PROTEINS = 8, 82
FP_RADIUS, FP_BITS = 2, 2048
RANDOM_STATE, N_INIT = 42, 10


def kmeans(X: np.ndarray, k: int) -> np.ndarray:
    return KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=N_INIT).fit_predict(X)


def crosstab_lines(umap_labels: pd.Series, hd_labels: pd.Series, top: int = 10) -> list[str]:
    """For the largest high-dim clusters: which 2D clusters they draw from."""
    lines = ["| High-dim cluster | Size | Main 2D cluster (share) | 2D clusters spanned |",
             "|---:|---:|---|---:|"]
    ct = pd.crosstab(hd_labels, umap_labels)
    for hd in ct.sum(axis=1).sort_values(ascending=False).index[:top]:
        row = ct.loc[hd]
        main = int(row.idxmax())
        lines.append(f"| {int(hd) + 1} | {int(row.sum())} | {main + 1} "
                     f"({row.max() / row.sum():.0%}) | {int((row > 0).sum())} |")
    return lines


def load_embeddings(accessions: list[str]) -> tuple[np.ndarray, list[str]]:
    """Vectors for the accessions present in the .h5. Keys are FASTA headers:
    "ACCESSION ENTRY_NAME | Organism" in the lab's file (first whitespace
    token = accession, as build_protein_centric_v2.py reads it), or the UniProt
    style "sp|P12345|NAME ..." (middle field)."""
    import h5py

    vecs: dict[str, np.ndarray] = {}
    def take(name: str, obj) -> None:
        # Walk every dataset: a "/" inside a header (e.g. an organism strain)
        # makes HDF5 nest it in groups, so top-level keys alone can miss some.
        if isinstance(obj, h5py.Dataset):
            first = name.split()[0]
            acc = first.split("|")[1] if first.count("|") >= 2 else first
            vecs[acc] = np.asarray(obj, dtype=np.float32).reshape(-1)

    with h5py.File(EMBEDDINGS, "r") as f:
        f.visititems(take)
    found = [a for a in accessions if a in vecs]
    return np.stack([vecs[a] for a in found]), found


def main() -> int:
    lines = ["# High-dimensional k-means vs. 2D UMAP k-means", "",
             "Same k as the published clusters; only the space differs. "
             "ARI = adjusted Rand index (1 = identical partitions, 0 = chance).", ""]

    # ── steroids: Morgan fingerprints ────────────────────────────────────
    mol = pd.read_csv(MOLECULES, low_memory=False)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    X = np.array([gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in mol["smiles"]],
                 dtype=np.float32)
    mol["cluster_hd"] = kmeans(X, K_MOLECULES)
    # Same row order as molecules.csv (one compound has no name, so names
    # can't key it); smiles is kept so consumers can check the alignment.
    mol[["compound_name", "smiles", "cluster_hd"]].to_csv(
        output("molecule_clusters_hd.tsv"), sep="\t", index=False)
    ari = adjusted_rand_score(mol["cluster"], mol["cluster_hd"])
    print(f"steroids: {len(mol)} molecules, k={K_MOLECULES}, ARI vs 2D = {ari:.3f}")
    lines += [f"## Steroids ({len(mol)} molecules, k = {K_MOLECULES})", "",
              f"- Space: Morgan fingerprints, radius {FP_RADIUS}, {FP_BITS} bits",
              f"- ARI vs 2D UMAP clusters: **{ari:.3f}**", "",
              *crosstab_lines(mol["cluster"], mol["cluster_hd"]), ""]

    # ── proteins: ProtT5 embeddings ──────────────────────────────────────
    prot = pd.read_csv(PROTEINS, usecols=["accession", "cluster"], low_memory=False)
    if not EMBEDDINGS.exists():
        msg = (f"proteins skipped: {EMBEDDINGS.relative_to(ROOT)} not found "
               "(ProtT5 vectors from the lab server — see analysis/README.md)")
        print(msg)
        lines += [f"## Proteins (k = {K_PROTEINS})", "", f"_{msg}_", ""]
    else:
        V, found = load_embeddings(prot["accession"].tolist())
        V /= np.linalg.norm(V, axis=1, keepdims=True)  # cosine geometry
        hd = pd.Series(kmeans(V, K_PROTEINS), index=found)
        prot["cluster_hd"] = prot["accession"].map(hd).fillna(-1).astype(int)
        prot[["accession", "cluster_hd"]].to_csv(
            output("protein_clusters_hd.tsv"), sep="\t", index=False)
        have = prot[prot["cluster_hd"] >= 0]
        missing = prot.loc[prot["cluster_hd"] < 0, "accession"].tolist()
        ari = adjusted_rand_score(have["cluster"], have["cluster_hd"])
        print(f"proteins: {len(have)}/{len(prot)} embedded, k={K_PROTEINS}, ARI vs 2D = {ari:.3f}")
        lines += [f"## Proteins ({len(have):,} of {len(prot):,} embedded, k = {K_PROTEINS})", "",
                  f"- Space: ProtT5 (`{EMBEDDINGS.name}`, {V.shape[1]}-d, L2-normalised)",
                  f"- ARI vs 2D UMAP clusters: **{ari:.3f}**",
                  f"- No vector (cluster_hd = -1): {len(missing)}"
                  + (f" — {', '.join(missing[:20])}" if missing else ""), "",
                  *crosstab_lines(have["cluster"], have["cluster_hd"]), ""]

    lines += ["## Settings", "",
              f"- `KMeans(random_state={RANDOM_STATE}, n_init={N_INIT})`",
              "- " + ", ".join(f"{p} {version(p)}" for p in ("scikit-learn", "rdkit", "numpy"))]
    output("highdim_kmeans_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote highdim_kmeans_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
