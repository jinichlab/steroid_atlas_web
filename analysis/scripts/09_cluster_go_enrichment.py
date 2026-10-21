"""Per-cluster GO enrichment — find clusters with clear vs. unclear functional identity.

For every (cluster, GO term) pair we compute a hypergeometric enrichment
statistic:

    background:      how often the GO term appears across the whole atlas
    cluster:         how often it appears in this cluster
    enrichment_ratio (cluster freq / background freq)
    log10_p_value    -log10 of the hypergeometric survival function

Then we rank clusters by the strength of their top enriched GO term. Clusters
with a strong dominant term (enrichment > 5×, p < 1e-20) have a "clear"
functional identity and are excellent case-study candidates. Clusters whose
top term is only mildly enriched (enrichment < 2×) are the "mixed/unclear"
group — worth discussing separately in the paper as either sub-family
composites or embedding-space artifacts.

Housekeeping GO terms (cytoplasm, membrane, protein binding, etc.) are
filtered so we don't just rediscover that "most steroid enzymes are in the
ER" for every ER-resident cluster.

Reads:  ../data/proteins.csv
Writes: ../analysis/cluster_go_enrichment.tsv   long form — one row per (cluster, GO term)
        ../analysis/cluster_go_top_terms.tsv    wide — top 5 enriched terms per cluster
        ../analysis/cluster_go_shortlist.txt    human-readable ranked shortlist
"""
from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
OUT_LONG = output("cluster_go_enrichment.tsv")
OUT_TOP = output("cluster_go_top_terms.tsv")
OUT_SHORT = output("cluster_go_shortlist.txt")

# Only consider GO terms that appear at least this often in the atlas.
# Very rare terms produce noisy enrichment ratios (1 hit in a 300-member
# cluster is 100+× if the atlas only has 3 hits total).
MIN_BG_COUNT = 15

# Report top-N enriched GO terms per cluster
TOP_N_PER_CLUSTER = 5

# For the shortlist file: show top-20 clusters ranked by enrichment strength
TOP_CLUSTERS_TO_SHOW = 20

GO_HOUSEKEEPING = {
    "cytoplasm", "cytosol", "nucleus", "nucleoplasm",
    "membrane", "plasma membrane", "integral component of membrane",
    "endoplasmic reticulum", "endoplasmic reticulum membrane",
    "mitochondrion", "mitochondrial matrix", "mitochondrial inner membrane",
    "mitochondrial outer membrane", "extracellular region", "extracellular space",
    "extracellular exosome", "extracellular vesicle", "vesicle",
    "protein binding", "metal ion binding", "atp binding", "zinc ion binding",
    "identical protein binding", "dna binding", "rna binding",
    "hydrolase activity", "transferase activity", "oxidoreductase activity",
    "catalytic activity", "iron ion binding",
    "cell surface", "perinuclear region of cytoplasm",
    "golgi apparatus", "golgi membrane", "lysosome",
    "reference proteome",
    "cellular response to xenobiotic stimulus",
    "atp hydrolysis activity",
    "nucleotide binding",
}


def _canon_cluster(x):
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return ""
    try:
        return str(int(float(s)))
    except (ValueError, TypeError):
        return s


def main() -> int:
    print(f"Reading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    for c in ["cluster", "go_labels", "protein_names"]:
        if c not in p.columns:
            p[c] = ""
        p[c] = p[c].fillna("").astype(str)
    p["cluster"] = p["cluster"].apply(_canon_cluster)
    p = p[p["cluster"] != ""].reset_index(drop=True)
    total_n = len(p)
    print(f"  {total_n:,} entries · {p['cluster'].nunique()} clusters")

    # Build a wide indicator: for each protein, the set of GO labels it carries
    print("Parsing GO labels per protein...")
    p["_go_set"] = p["go_labels"].apply(
        lambda s: {g.strip() for g in s.split(";") if g.strip()
                   and g.strip().lower() not in GO_HOUSEKEEPING}
    )

    # Background frequency of each GO term across the atlas
    bg_ctr = Counter()
    for s in p["_go_set"]:
        bg_ctr.update(s)
    # Filter to terms that meet the minimum background count
    considered_terms = {t for t, c in bg_ctr.items() if c >= MIN_BG_COUNT}
    print(f"  {len(bg_ctr):,} unique GO terms · {len(considered_terms):,} with count ≥ {MIN_BG_COUNT}")

    # Precompute cluster membership vectors
    cluster_ids = sorted(p["cluster"].unique(), key=lambda x: int(x))
    print(f"\nComputing hypergeometric enrichment for {len(cluster_ids)} clusters × "
          f"{len(considered_terms)} GO terms ...")

    records = []
    N = total_n
    for cid in cluster_ids:
        cluster_mask = (p["cluster"] == cid).values
        n_cluster = int(cluster_mask.sum())
        # Count GO term occurrences within this cluster
        cluster_ctr = Counter()
        for s in p.loc[cluster_mask, "_go_set"]:
            cluster_ctr.update(s)

        # sorted: set order varies between runs (string-hash randomization)
        for term in sorted(considered_terms):
            k = cluster_ctr.get(term, 0)
            if k == 0:
                continue
            K = bg_ctr[term]         # background hits
            # Hypergeometric survival function P(X >= k)
            pval = float(hypergeom.sf(k - 1, N, K, n_cluster))
            if pval <= 0:
                log10_p = 300.0
            else:
                log10_p = min(300.0, -math.log10(pval))
            cluster_freq = k / n_cluster
            bg_freq = K / N
            enrichment = cluster_freq / bg_freq if bg_freq > 0 else np.nan
            records.append({
                "cluster": cid,
                "cluster_size": n_cluster,
                "go_term": term,
                "k_in_cluster": k,
                "K_in_atlas": K,
                "cluster_freq": round(cluster_freq, 4),
                "bg_freq": round(bg_freq, 4),
                "enrichment_ratio": round(enrichment, 2),
                "log10_p": round(log10_p, 2),
            })

    long_df = pd.DataFrame(records)
    long_df.to_csv(OUT_LONG, sep="\t", index=False)
    print(f"  wrote {OUT_LONG.name}  ({len(long_df):,} rows)")

    # Top-N enriched terms per cluster (ranked by log10_p, then enrichment)
    top_rows = []
    for cid, sub in long_df.groupby("cluster"):
        top = (sub[sub["k_in_cluster"] >= 3]
               .sort_values(["log10_p", "enrichment_ratio"], ascending=[False, False])
               .head(TOP_N_PER_CLUSTER))
        for rank, (_, r) in enumerate(top.iterrows(), 1):
            top_rows.append({**r.to_dict(), "rank": rank})
    top_df = pd.DataFrame(top_rows)
    top_df.to_csv(OUT_TOP, sep="\t", index=False)
    print(f"  wrote {OUT_TOP.name}  (top {TOP_N_PER_CLUSTER} enriched terms per cluster)")

    # ─── Shortlist: rank clusters by their TOP enrichment strength ───────────
    # Case-study "clear-story" clusters have a strong top term (high log10_p + high ratio).
    # "Unclear" clusters have a weak top term (small log10_p or low ratio).
    per_cluster_summary = []
    for cid in cluster_ids:
        rows = top_df[top_df["cluster"] == cid]
        if len(rows) == 0:
            per_cluster_summary.append({
                "cluster": cid,
                "cluster_size": int((p["cluster"] == cid).sum()),
                "top_go": "",
                "top_enrichment": 0.0,
                "top_log10_p": 0.0,
                "n_terms_10x": 0,
                "clarity": "no_go",
            })
            continue
        first = rows.iloc[0]
        # Number of terms enriched at least 10x
        n_10x = int((rows["enrichment_ratio"] >= 10.0).sum())
        # Clarity label
        if first["enrichment_ratio"] >= 10 and first["log10_p"] >= 20:
            clarity = "clear"
        elif first["enrichment_ratio"] >= 3 and first["log10_p"] >= 10:
            clarity = "moderate"
        else:
            clarity = "unclear"
        per_cluster_summary.append({
            "cluster": cid,
            "cluster_size": int(first["cluster_size"]),
            "top_go": first["go_term"],
            "top_enrichment": float(first["enrichment_ratio"]),
            "top_log10_p": float(first["log10_p"]),
            "n_terms_10x": n_10x,
            "clarity": clarity,
        })
    summary = pd.DataFrame(per_cluster_summary)
    # Ranking: clear > moderate > unclear; within clarity, higher log10_p first
    _clarity_order = {"clear": 0, "moderate": 1, "unclear": 2, "no_go": 3}
    summary["_ord"] = summary["clarity"].map(_clarity_order)
    summary = summary.sort_values(["_ord", "top_log10_p"], ascending=[True, False])
    summary = summary.drop(columns=["_ord"])

    # Human-readable shortlist
    lines = []
    lines.append("=" * 82)
    lines.append(f"PER-CLUSTER GO ENRICHMENT — {total_n:,} proteins · {len(cluster_ids)} clusters")
    lines.append("=" * 82)
    lines.append("")
    lines.append("Clarity buckets:")
    lines.append(f"  clear    : top GO enrichment ≥ 10×  AND  -log10(p) ≥ 20")
    lines.append(f"  moderate : top GO enrichment ≥ 3×   AND  -log10(p) ≥ 10")
    lines.append(f"  unclear  : below both thresholds")
    lines.append(f"  no_go    : no annotated members with GO terms")
    lines.append("")
    for label in ("clear", "moderate", "unclear", "no_go"):
        n = int((summary["clarity"] == label).sum())
        lines.append(f"  {n:>3}  {label}")
    lines.append("")
    lines.append(f"--- TOP {TOP_CLUSTERS_TO_SHOW} CLUSTERS (strongest to weakest enrichment story) ---")
    lines.append("")
    for _, r in summary.head(TOP_CLUSTERS_TO_SHOW).iterrows():
        lines.append(f"cluster {r['cluster']}  size={r['cluster_size']}  [{r['clarity']}]")
        lines.append(f"    top: {r['top_go']}")
        lines.append(f"         enrichment = {r['top_enrichment']:.1f}×    "
                     f"-log10(p) = {r['top_log10_p']:.1f}    "
                     f"terms ≥10×: {r['n_terms_10x']}")
        # Also list next 2-3 enriched terms if present
        others = top_df[top_df["cluster"] == r["cluster"]].iloc[1:4]
        for _, o in others.iterrows():
            lines.append(f"    also: {o['go_term']}  "
                         f"({o['enrichment_ratio']:.1f}×, -log10p={o['log10_p']:.1f})")
        lines.append("")
    lines.append("")
    lines.append(f"--- BOTTOM 10 CLUSTERS (weakest — potential 'mixed/unclear' case study) ---")
    lines.append("")
    for _, r in summary.tail(10).iterrows():
        lines.append(f"cluster {r['cluster']}  size={r['cluster_size']}  [{r['clarity']}]"
                     f"  top: {r['top_go'] or '(no annotated GO)'}  "
                     f"({r['top_enrichment']:.1f}×, -log10p={r['top_log10_p']:.1f})")
    lines.append("")
    OUT_SHORT.write_text("\n".join(lines) + "\n")
    print(f"  wrote {OUT_SHORT.name}")
    print()
    print("Clarity buckets:")
    print(summary["clarity"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
