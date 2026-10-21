"""Compute ProtT5 embedding for CpBSH/T (WP_243289361), find its cosine nearest
neighbor in the existing 37,391-protein embedding matrix, and inherit the
neighbor's UMAP_1/UMAP_2 coordinates (with small Gaussian jitter).

Same approach as literature/project_new_into_umap.py — the original UMAP reducer
was never pickled, so nearest-neighbor projection is the faithful substitute
(UMAP preserves local cosine structure, so a top-cosine neighbor lies in the
same local manifold patch).

Reads:
  literature/sequences/guzior/WP_243289361.fasta   input sequence
  /home/adsiordia/steroid_core_classifier/embeddings/all_steroid_uniprot_comprehensive_v2.h5
  protein_sequence_embedding.DEDUP.csv

Writes:
  literature/embeddings/cpbsht.h5                  new ProtT5 embedding
  protein_sequence_embedding.DEDUP.csv             UMAP_1/UMAP_2/clusters/prott5_cluster updated for WP_243289361
  literature/embeddings/cpbsht_projection_report.json  full provenance
"""
from __future__ import annotations
import json
import re
import shutil
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from transformers import T5EncoderModel, T5Tokenizer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FASTA = HERE / "sequences" / "guzior" / "WP_243289361.fasta"
DEDUP = ROOT / "protein_sequence_embedding.DEDUP.csv"
H5_OLD = Path("/home/adsiordia/steroid_core_classifier/embeddings/all_steroid_uniprot_comprehensive_v2.h5")
EMB_DIR = HERE / "embeddings"
EMB_DIR.mkdir(exist_ok=True)
H5_NEW = EMB_DIR / "cpbsht.h5"
REPORT = EMB_DIR / "cpbsht_projection_report.json"

TARGET_ACC = "WP_243289361"

# ─── 1. Read sequence from FASTA ───────────────────────────────────────────
lines = FASTA.read_text().strip().split("\n")
header = lines[0].lstrip(">").strip()
sequence = "".join(l.strip() for l in lines[1:]).replace("*", "").upper()
# Replace non-standard AAs the same way ProtT5 tokenizer expects
sequence_spaced = " ".join(list(re.sub(r"[UZOB]", "X", sequence)))
print(f"Sequence:  {len(sequence)} aa")
print(f"Header:    {header[:80]}")

# ─── 2. Load ProtT5 model ──────────────────────────────────────────────────
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"\nDevice: {device}")
print(f"Loading Rostlab/prot_t5_xl_half_uniref50-enc ...")
t0 = time.time()
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
if device.type == "cpu":
    model = model.to(torch.float32)
model = model.to(device).eval()
tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False)
print(f"  loaded in {time.time()-t0:.1f}s")

# ─── 3. Compute embedding ─────────────────────────────────────────────────
print(f"\nEmbedding CpBSH/T ({len(sequence)} aa) ...")
t0 = time.time()
with torch.no_grad():
    token_encoding = tokenizer(
        [sequence_spaced], add_special_tokens=True, padding=True, return_tensors="pt"
    ).to(device)
    embedding = model(
        input_ids=token_encoding["input_ids"],
        attention_mask=token_encoding["attention_mask"],
    )
    # embedding.last_hidden_state: (1, seq_len_with_special, 1024)
    # Trim last token (EOS) and mean-pool
    per_res = embedding.last_hidden_state[0, :len(sequence)].cpu().numpy().astype(np.float32)
    per_protein = per_res.mean(axis=0)  # (1024,)
print(f"  done in {time.time()-t0:.1f}s")
print(f"  per-protein vector shape: {per_protein.shape}")
print(f"  vector norm: {np.linalg.norm(per_protein):.3f}")

# ─── 4. Save new h5 ────────────────────────────────────────────────────────
with h5py.File(H5_NEW, "w") as f:
    f.create_dataset(TARGET_ACC, data=per_protein)
print(f"\nWrote {H5_NEW}")

# ─── 5. Load old h5 for nearest-neighbor lookup ────────────────────────────
print(f"\nLoading old embedding matrix from {H5_OLD.name} ...")
t0 = time.time()
with h5py.File(H5_OLD, "r") as f:
    old_keys = list(f.keys())
    old_mat = np.stack([f[k][:] for k in old_keys]).astype(np.float32)
print(f"  shape: {old_mat.shape}   loaded in {time.time()-t0:.1f}s")
old_acc = [k.split()[0] for k in old_keys]

# ─── 6. Cosine nearest neighbor ───────────────────────────────────────────
old_norm = old_mat / (np.linalg.norm(old_mat, axis=1, keepdims=True) + 1e-9)
new_norm = per_protein / (np.linalg.norm(per_protein) + 1e-9)
sims = old_norm @ new_norm
top_k_idx = np.argsort(-sims)[:5]

print("\nTop 5 nearest neighbors (cosine):")
for rank, idx in enumerate(top_k_idx):
    print(f"  {rank+1}. {old_acc[idx]:<15} sim={sims[idx]:.4f}")

nn_idx = int(top_k_idx[0])
nn_acc = old_acc[nn_idx]
nn_sim = float(sims[nn_idx])
print(f"\nSelected NN: {nn_acc} (cosine sim = {nn_sim:.4f})")

# ─── 7. Look up NN's UMAP coords in DEDUP CSV ─────────────────────────────
df = pd.read_csv(DEDUP, low_memory=False)
nn_rows = df[df["Entry"] == nn_acc]
if len(nn_rows) == 0:
    # Fall back: nn was removed by dedup. Try next-best neighbors.
    print(f"WARN: NN {nn_acc} not in DEDUP CSV (probably dropped in dedup). Trying next best...")
    for rank in range(1, 5):
        candidate = old_acc[int(top_k_idx[rank])]
        if len(df[df["Entry"] == candidate]) > 0:
            nn_idx = int(top_k_idx[rank])
            nn_acc = candidate
            nn_sim = float(sims[nn_idx])
            nn_rows = df[df["Entry"] == nn_acc]
            print(f"  fell back to {nn_acc} (sim={nn_sim:.4f})")
            break
    if len(nn_rows) == 0:
        raise SystemExit("ERROR: no top-5 neighbor found in DEDUP CSV")

nn_row = nn_rows.iloc[0]
print(f"\nNN row in DEDUP CSV:")
print(f"  Entry: {nn_row['Entry']}")
print(f"  Organism: {str(nn_row['Organism'])[:80]}")
print(f"  Protein: {str(nn_row['Protein names'])[:80]}")
print(f"  UMAP_1: {nn_row['UMAP_1']}   UMAP_2: {nn_row['UMAP_2']}")
print(f"  cluster: {nn_row['clusters']}")

# ─── 8. Update target row with jittered coords ─────────────────────────────
rng = np.random.default_rng(seed=42)
jitter_x, jitter_y = rng.normal(0, 0.15, size=2)
target_umap_1 = float(nn_row["UMAP_1"]) + float(jitter_x)
target_umap_2 = float(nn_row["UMAP_2"]) + float(jitter_y)
target_cluster = int(nn_row["clusters"]) if pd.notna(nn_row["clusters"]) else -1
target_prott5_cluster = nn_row.get("prott5_cluster", None)

target_mask = df["Entry"] == TARGET_ACC
if target_mask.sum() != 1:
    raise SystemExit(f"ERROR: expected 1 target row, got {target_mask.sum()}")

# Backup
bak = DEDUP.with_suffix(DEDUP.suffix + ".bak4")
shutil.copy2(DEDUP, bak)
print(f"\nBackup -> {bak.name}")

df.loc[target_mask, "UMAP_1"] = target_umap_1
df.loc[target_mask, "UMAP_2"] = target_umap_2
df.loc[target_mask, "clusters"] = target_cluster
if target_prott5_cluster is not None and pd.notna(target_prott5_cluster):
    df.loc[target_mask, "prott5_cluster"] = target_prott5_cluster

df.to_csv(DEDUP, index=False)

# Verify
verify = df[target_mask].iloc[0]
print(f"\n[OK] Updated row for {TARGET_ACC}:")
print(f"  UMAP_1: {verify['UMAP_1']}   UMAP_2: {verify['UMAP_2']}")
print(f"  clusters: {verify['clusters']}   prott5_cluster: {verify.get('prott5_cluster', 'n/a')}")

# ─── 9. Provenance report ─────────────────────────────────────────────────
report = {
    "target_accession": TARGET_ACC,
    "target_sequence_length": len(sequence),
    "target_species": "Clostridium perfringens ATCC 13124",
    "paper": "Guzior et al., Nature 2024",
    "paper_doi": "10.1038/s41586-024-07017-8",
    "prott5_model": "Rostlab/prot_t5_xl_half_uniref50-enc",
    "prott5_embedding_dim": int(per_protein.shape[0]),
    "prott5_embedding_l2_norm": float(np.linalg.norm(per_protein)),
    "projection_method": "cosine nearest neighbor in ProtT5 space + Gaussian jitter (sigma=0.15)",
    "projection_rationale": (
        "Original UMAP reducer was not pickled; nearest-neighbor projection is "
        "the faithful substitute since UMAP preserves local cosine structure "
        "(k=15 neighbor graph). Same method used to add the previous 17 entries."
    ),
    "nearest_neighbor": {
        "accession": nn_acc,
        "cosine_similarity": nn_sim,
        "organism": str(nn_row["Organism"]),
        "protein_names": str(nn_row["Protein names"]),
        "umap_1": float(nn_row["UMAP_1"]),
        "umap_2": float(nn_row["UMAP_2"]),
        "cluster": int(nn_row["clusters"]) if pd.notna(nn_row["clusters"]) else None,
    },
    "top5_neighbors": [
        {"accession": old_acc[int(i)], "cosine_similarity": float(sims[int(i)])}
        for i in top_k_idx
    ],
    "assigned_coordinates": {
        "umap_1": target_umap_1,
        "umap_2": target_umap_2,
        "jitter_sigma": 0.15,
        "jitter_x_applied": float(jitter_x),
        "jitter_y_applied": float(jitter_y),
        "cluster": target_cluster,
    },
    "backup_file": str(bak.relative_to(ROOT)),
}
REPORT.write_text(json.dumps(report, indent=2))
print(f"\nProvenance report -> {REPORT.relative_to(ROOT)}")
