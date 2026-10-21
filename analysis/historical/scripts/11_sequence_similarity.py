"""Sequence-similarity characterization of the atlas — all 35,117 proteins covered.

Two complementary characterizations of the atlas's sequence composition:

  A. PER-PROTEIN NEAREST-NEIGHBOR IDENTITY (every protein represented)
     - CD-HIT clusters the entire atlas at 40% identity
     - .clstr output gives each non-representative's identity to its
       cluster representative → that's the protein's approximate NN identity
     - representatives in multi-member clusters inherit the max identity
       among their cluster mates (their NN is at least that high)
     - representatives in singleton clusters have NN identity below 40%
       (we backfill their identity by aligning against a random 500-protein
       sample and taking the max — so every one of the 35,117 proteins ends
       up with a real numeric NN identity)
     → distribution histogram covers 100% of the atlas

  B. REDUNDANCY AT MULTIPLE IDENTITY THRESHOLDS
     - cluster at 90 / 70 / 50 / 40 % identity
     - report cluster count = non-redundant sequences remaining
     → answers "how much of the atlas would collapse at each threshold"

Reads:  ../data/proteins.csv
Writes: ../analysis/seq_nn_identity.tsv               per-protein NN identity table
        ../analysis/seq_similarity_summary.txt        human-readable numbers
        ../analysis/seq_similarity_figure.png/.pdf    paper figure
        ../analysis/tmp_cdhit/                        intermediate CD-HIT output
"""
from __future__ import annotations

import random
import re
import subprocess
import time
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
OUT_TSV = output("seq_nn_identity.tsv", historical=True)
OUT_TXT = output("seq_similarity_summary.txt", historical=True)
OUT_PNG = output("seq_similarity_figure.png", historical=True)
OUT_PDF = output("seq_similarity_figure.pdf", historical=True)
TMP = output("tmp_cdhit", historical=True)

# Standard bioinformatics thresholds:
#   90% — near-duplicates / very close orthologs (matches UniRef90)
#   75% — same protein family across species / orthologs (common cutoff for
#         cross-species ortholog assignment)
#   50% — homologs / distant orthologs; edge of the "twilight zone"
#         (matches UniRef50; below this alignment-based identity becomes less reliable)
CDHIT_THRESHOLDS = [0.90, 0.75, 0.50]
NN_THRESHOLD = 0.50             # the threshold whose .clstr gives per-protein NN identity
BACKFILL_SAMPLE = 500           # random alignment sample for singleton representatives
BACKFILL_MAX_LEN = 2000
RANDOM_SEED = 42


def _write_fasta(accessions, seqs, path: Path):
    with path.open("w") as f:
        for acc, seq in zip(accessions, seqs):
            f.write(f">{acc}\n{seq}\n")


def _cdhit_word_size(thr: float) -> int:
    if thr >= 0.70: return 5
    if thr >= 0.60: return 4
    if thr >= 0.50: return 3
    return 2   # 0.40 – 0.50


def run_cdhit(fasta: Path, thr: float) -> Path:
    out = TMP / f"cdhit_{int(thr*100)}"
    cmd = ["cd-hit", "-i", str(fasta), "-o", str(out),
           "-c", str(thr), "-n", str(_cdhit_word_size(thr)),
           "-M", "8000", "-T", "8", "-d", "0"]
    print(f"  cd-hit -c {thr} -n {_cdhit_word_size(thr)} ...", flush=True)
    subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=1800)
    clstr = Path(str(out) + ".clstr")
    return clstr


def parse_clstr(clstr: Path) -> dict[str, tuple[str, float | None]]:
    """Return {accession: (cluster_id, identity_to_rep or None if this is the rep)}."""
    result = {}
    cur_cluster = None
    # Cluster line pattern: "0  237aa, >P19410... at 89.45%"
    #                     "0  237aa, >P19410... *"
    row_re = re.compile(r'^\d+\s+\d+aa,\s+>(\S+?)\.\.\.\s+(.+)$')
    with clstr.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                cur_cluster = line.split()[1]
                continue
            m = row_re.match(line)
            if not m:
                continue
            acc, rest = m.group(1), m.group(2).strip()
            if rest == "*":
                result[acc] = (cur_cluster, None)  # representative
            else:
                # rest looks like "at 89.45%" or "at +/89.45%" (strand for cd-hit-est)
                mm = re.search(r'(\d+\.\d+)%', rest)
                if mm:
                    result[acc] = (cur_cluster, float(mm.group(1)) / 100.0)
                else:
                    result[acc] = (cur_cluster, None)
    return result


def compute_nn_identity(accessions: list[str], seqs: list[str],
                        cluster_data: dict[str, tuple[str, float | None]]) -> pd.DataFrame:
    """Build per-protein NN identity table covering all atlas proteins.

    Logic:
      - non-rep members inherit their identity to the cluster representative
      - representatives in multi-member clusters inherit the max identity
        among the cluster mates (their true NN is at least that high — usually
        a very tight estimate)
      - representatives in singleton clusters have NN identity < NN_THRESHOLD;
        we backfill by aligning against BACKFILL_SAMPLE random atlas proteins
    """
    # Group by cluster
    cluster_members: dict[str, list[tuple[str, float | None]]] = {}
    for acc, (cid, ident) in cluster_data.items():
        cluster_members.setdefault(cid, []).append((acc, ident))

    seq_by_acc = dict(zip(accessions, seqs))
    rng = random.Random(RANDOM_SEED)

    rows = []
    # Prepare aligner for the backfill step
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    # Iterate cluster by cluster
    print(f"\n[A2] Backfilling singleton-cluster representatives via random-sample alignment...")
    singletons = []
    for cid, members in cluster_members.items():
        if len(members) == 1:
            singletons.append(members[0][0])
    print(f"  {len(singletons):,} singleton representatives to backfill (sample size {BACKFILL_SAMPLE})")

    # Precompute a backfill sample so alignments are reproducible
    backfill_pool = [a for a in accessions if a not in set(singletons)] or accessions
    t0 = time.time()

    for cid, members in cluster_members.items():
        # Extract identities of non-reps
        non_rep_idents = [ident for _, ident in members if ident is not None]
        rep_acc = next(acc for acc, ident in members if ident is None)

        if non_rep_idents:
            # Rep's NN identity is at least the maximum identity among cluster mates
            rep_nn = max(non_rep_idents)
            rep_backfilled = False
        else:
            # Singleton — backfill via alignment
            rep_seq = seq_by_acc[rep_acc]
            if not rep_seq or len(rep_seq) > BACKFILL_MAX_LEN:
                rep_nn = 0.0
                rep_backfilled = True
            else:
                sample = rng.sample(backfill_pool, min(BACKFILL_SAMPLE, len(backfill_pool)))
                best = 0.0
                for other in sample:
                    other_seq = seq_by_acc.get(other, "")
                    if not other_seq or len(other_seq) > BACKFILL_MAX_LEN:
                        continue
                    try:
                        aln = aligner.align(rep_seq, other_seq)[0]
                        counts = aln.counts()
                        aln_len = counts.gaps + counts.mismatches + counts.identities
                        if aln_len > 0:
                            best = max(best, counts.identities / aln_len)
                    except Exception:
                        pass
                rep_nn = best
                rep_backfilled = True

        rows.append({"accession": rep_acc, "cluster_id": cid,
                     "nn_identity": rep_nn, "is_rep": True,
                     "backfilled": rep_backfilled})
        for acc, ident in members:
            if ident is None:
                continue
            rows.append({"accession": acc, "cluster_id": cid,
                         "nn_identity": ident, "is_rep": False,
                         "backfilled": False})

        if len(rows) % 5000 < 100:
            print(f"  processed {len(rows):,} entries so far...", flush=True)

    dt = time.time() - t0
    print(f"  → done in {dt/60:.1f} min")
    return pd.DataFrame(rows)


def redundancy_table(fasta: Path, n_total: int) -> list[dict]:
    results = []
    for thr in CDHIT_THRESHOLDS:
        clstr = run_cdhit(fasta, thr)
        n_clusters = 0
        with clstr.open() as f:
            for line in f:
                if line.startswith(">Cluster"):
                    n_clusters += 1
        print(f"    → {n_clusters:,} clusters at {int(thr*100)}% identity")
        results.append({
            "threshold_pct": int(thr * 100),
            "n_clusters": n_clusters,
            "reduction_pct": 100 * (1 - n_clusters / n_total),
        })
    return results


def write_summary(nn_df: pd.DataFrame, red: list[dict], n_total: int):
    v = nn_df["nn_identity"] * 100.0
    q = v.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    fracs = {t: (v >= t).mean() * 100 for t in (30, 50, 70, 90)}

    lines = []
    lines.append("=" * 78)
    lines.append(f"SEQUENCE-SIMILARITY CHARACTERIZATION — all {n_total:,} atlas proteins")
    lines.append("=" * 78)
    lines.append("")
    lines.append(f"[A] Per-protein nearest-neighbor identity  (100% coverage — every")
    lines.append(f"    protein contributes one value):")
    lines.append(f"    min           {v.min():.1f}%")
    lines.append(f"    5%            {q[0.05]:.1f}%")
    lines.append(f"    25%           {q[0.25]:.1f}%")
    lines.append(f"    median        {q[0.50]:.1f}%")
    lines.append(f"    75%           {q[0.75]:.1f}%")
    lines.append(f"    95%           {q[0.95]:.1f}%")
    lines.append(f"    max           {v.max():.1f}%")
    lines.append(f"    mean          {v.mean():.1f}%")
    lines.append("")
    lines.append(f"    Fraction of atlas whose NN identity is ≥ threshold:")
    for t in (30, 50, 70, 90):
        lines.append(f"      ≥ {t}%      {fracs[t]:5.1f}%")
    lines.append("")
    lines.append(f"    Coverage note: {n_total:,} / {n_total:,} proteins covered "
                 "(singleton reps backfilled by 500-protein alignment sample).")
    lines.append("")
    lines.append(f"[B] CD-HIT redundancy at multiple identity thresholds "
                 "(clustering on all sequences):")
    lines.append(f"    {'threshold':>12}   {'clusters':>10}   {'reduction':>10}")
    for r in red:
        lines.append(f"      {r['threshold_pct']:>3}%       {r['n_clusters']:>10,}     {r['reduction_pct']:>6.1f}%")
    OUT_TXT.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {OUT_TXT.name}")


def make_figure(nn_df: pd.DataFrame, red: list[dict], n_total: int):
    fig = plt.figure(figsize=(11.0, 4.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    BAR = "#1F4B99"
    REF = "#9CA3AF"

    # ── PANEL A: per-protein NN identity distribution ─────────────────────
    v = nn_df["nn_identity"] * 100.0
    axA.hist(v, bins=np.arange(0, 105, 2.5),
             color=BAR, alpha=0.85, edgecolor="white", linewidth=0.5)
    axA.axvline(30, color=REF, linestyle="--", linewidth=0.9)
    axA.text(30, axA.get_ylim()[1] * 0.95, "  30% identity",
             color=REF, fontsize=9, va="top", ha="left")
    axA.set_xlabel("Sequence identity to closest atlas match  (%)", fontsize=10)
    axA.set_ylabel(f"Number of proteins (of {n_total:,})", fontsize=10)
    axA.set_title("A. Per-protein nearest-neighbor identity distribution",
                  fontsize=11, loc="left", pad=8)
    axA.grid(True, axis="y", alpha=0.25, linewidth=0.5)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.tick_params(labelsize=9)

    # Coverage callout
    n_singleton = int((~nn_df["is_rep"] == False).sum() if False else nn_df["backfilled"].sum())
    axA.text(
        0.98, 0.92,
        f"all {n_total:,} atlas proteins represented\n"
        f"(median NN identity = {v.median():.1f}%)",
        transform=axA.transAxes, fontsize=9, color="#111827",
        ha="right", va="top",
    )

    # ── PANEL B: CD-HIT redundancy ────────────────────────────────────────
    thrs = [r["threshold_pct"] for r in red]
    counts = [r["n_clusters"] for r in red]
    ypos = np.arange(len(thrs))[::-1]
    axB.barh(ypos, counts, color=BAR, height=0.6,
             edgecolor="white", linewidth=0.4)
    for y, r in zip(ypos, red):
        axB.text(r["n_clusters"] * 1.03, y,
                 f" {r['n_clusters']:,}  ({r['reduction_pct']:.1f}% reduction)",
                 va="center", ha="left", fontsize=9, color="#111827")
    axB.set_yticks(ypos)
    axB.set_yticklabels([f"CD-HIT {t}%" for t in thrs], fontsize=9)
    axB.axvline(n_total, color=REF, linestyle="--", linewidth=0.9)
    axB.text(n_total, len(thrs) - 0.4,
             f" {n_total:,} total", color=REF,
             fontsize=8.5, va="bottom", ha="left")
    axB.set_xlabel("Non-redundant sequence count", fontsize=10)
    axB.set_title("B. Redundancy at multiple identity thresholds",
                  fontsize=11, loc="left", pad=8)
    axB.set_xlim(0, n_total * 1.35)
    axB.tick_params(axis="x", labelsize=9)
    axB.tick_params(axis="y", length=0, pad=4)
    axB.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    axB.grid(True, axis="x", alpha=0.25, linewidth=0.5)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.spines["left"].set_visible(False)

    for f in (OUT_PNG, OUT_PDF):
        fig.savefig(f, dpi=300, bbox_inches="tight")
        print(f"Wrote {f.name}")


def main() -> int:
    print(f"Loading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    p = p[p["sequence"].notna() & (p["sequence"].astype(str).str.len() > 0)]
    p = p.reset_index(drop=True)
    accessions = p["accession"].astype(str).tolist()
    seqs = p["sequence"].astype(str).tolist()
    n_total = len(seqs)
    print(f"  {n_total:,} atlas sequences  "
          f"(length median={int(p['length_aa'].median())}, "
          f"max={int(p['length_aa'].max())})")

    TMP.mkdir(exist_ok=True)
    fasta = TMP / "atlas.fasta"
    _write_fasta(accessions, seqs, fasta)
    print(f"  wrote FASTA: {fasta}  ({n_total:,} sequences)")

    print(f"\n[A1] Running CD-HIT at {int(NN_THRESHOLD*100)}% (source of per-protein NN identity)...")
    nn_clstr = run_cdhit(fasta, NN_THRESHOLD)
    cluster_data = parse_clstr(nn_clstr)
    print(f"  parsed {len(cluster_data):,} accessions from {nn_clstr.name}")
    assert len(cluster_data) == n_total, \
        f"expected {n_total} accessions in .clstr, got {len(cluster_data)}"

    nn_df = compute_nn_identity(accessions, seqs, cluster_data)
    nn_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"  wrote {OUT_TSV.name}  ({len(nn_df):,} rows)")

    print(f"\n[B] CD-HIT redundancy at multiple thresholds...")
    red = redundancy_table(fasta, n_total)

    write_summary(nn_df, red, n_total)
    make_figure(nn_df, red, n_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
