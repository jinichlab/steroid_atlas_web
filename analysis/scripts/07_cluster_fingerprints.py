"""Build a per-cluster fingerprint of the atlas to identify case-study candidates.

For each cluster, compute:
  - size
  - dominant_stem      most common protein-name first token
  - stem_coherence     % of cluster sharing that stem
  - dominant_go        most-informative GO label (housekeeping terms filtered)
  - go_coherence       % of cluster carrying that GO label
  - dominant_kw        most-informative UniProt keyword (housekeeping terms filtered)
  - kw_coherence       % of cluster carrying that keyword
  - dominant_ec        most common EC top-level (or none)
  - ec_coherence       % of cluster with that EC top-level (of those with any EC)
  - outliers           accession list of members whose stem doesn't match the dominant
  - case_study_score   heuristic: prefer high coherence, medium size, small but nonzero outliers

Reads:  ../data/proteins.csv
Writes: ../analysis/cluster_fingerprints.tsv     one row per cluster, ranked
        ../analysis/cluster_shortlist.txt        top-N clusters with a scannable summary
"""
from __future__ import annotations

import argparse

import re
from collections import Counter
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"

TOP_N_SHORTLIST = 15

# GO terms that appear in most eukaryotic proteins and add no discriminative power
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
}

# UniProt keywords that are curation infrastructure, not functional discriminators
KW_HOUSEKEEPING = {
    "Reference proteome", "3D-structure", "Direct protein sequencing",
    "Signal", "Transmembrane", "Transmembrane helix", "Membrane",
    "Cell membrane", "Cytoplasm", "Nucleus", "Mitochondrion",
    "Endoplasmic reticulum", "Alternative splicing", "Phosphoprotein",
    "Glycoprotein", "Disulfide bond", "Metal-binding", "Zinc",
    "Isopeptide bond", "Ubl conjugation", "Acetylation",
    "ATP-binding", "Nucleotide-binding",
}

STEM_STOPWORDS = {"protein", "putative", "probable", "uncharacterized",
                  "hypothetical", "sub", "unnamed", "fragment"}


def _clean_first_name(name: str) -> str:
    """Strip parenthesized aliases and drop trailing EC / accession noise so we
    keep just the primary protein name (still in original casing).

    Aliases-in-parens are recognized only when there's whitespace before the
    open-paren — that keeps identifiers like ``Delta(24)-sterol reductase``
    intact while still stripping ``(EC 1.2.3)`` and other trailing aliases.
    """
    if not name:
        return ""
    s = str(name).strip()
    s = re.sub(r"\s+\(.*", "", s)          # drop " (alias) (…)" only when the paren stands apart
    s = re.split(r"\s*[;]\s*", s, 1)[0]    # drop everything past a ";"
    s = re.sub(r"\s+EC\s*[\d.]+\s*$", "", s, flags=re.IGNORECASE)  # trailing " EC 1.2.3"
    # Drop trailing "/17,20-lyase" style fragments that indicate joined names
    s = re.sub(r"/\d.*$", "", s).strip()
    return s


def name_key(name: str) -> str:
    """Coarse grouping key: first two meaningful tokens, lowercased.
    Used to compute stem coherence — many families share word 1 but differ
    on word 2 (e.g. "Cytochrome P450" vs "Cytochrome b5")."""
    s = _clean_first_name(name)
    tokens = [t for t in re.split(r"[\s,;]+", s) if t and t.lower() not in STEM_STOPWORDS]
    if not tokens:
        return ""
    return " ".join(tokens[:2]).lower()


def display_name(name: str) -> str:
    """Full, casing-preserved display name — used as the cluster label
    everywhere the user sees it (picker, legend, on-plot centroid text).
    Keeps the primary protein name; drops parenthesized aliases and EC noise."""
    return _clean_first_name(name)


# Kept for downstream imports that referenced the old name.
name_stem = name_key


def top_informative(items, housekeeping) -> tuple[str, int]:
    """Return (label, count) of the most frequent item not in the housekeeping set."""
    if not items:
        return "", 0
    ctr = Counter(x for x in items if x)
    for label, count in ctr.most_common():
        if label.lower() not in {h.lower() for h in housekeeping}:
            return label, count
    return "", 0


def semicolon_flat(series) -> list[list[str]]:
    """Turn a Series of 'a;b;c' strings into a list-of-lists."""
    return [
        [x.strip() for x in str(s).split(";") if str(x).strip()]
        for s in series.fillna("")
    ]


def main() -> int:
    # --labels cluster_hd names the high-dimensional k-means clusters
    # (19_highdim_kmeans.py) instead of the published 2D ones.
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", choices=["cluster", "cluster_hd"], default="cluster")
    labels = ap.parse_args().labels
    suffix = "_hd" if labels == "cluster_hd" else ""
    out_tsv = output(f"cluster_fingerprints{suffix}.tsv")
    out_txt = output(f"cluster_shortlist{suffix}.txt")

    print(f"Reading {IN.name}...")
    p = pd.read_csv(IN, low_memory=False)
    if labels == "cluster_hd":
        hd = pd.read_csv(output("protein_clusters_hd.tsv"), sep="\t")
        p = p.drop(columns="cluster").merge(
            hd.rename(columns={"cluster_hd": "cluster"}), on="accession")
        p = p[p["cluster"] >= 0]  # proteins without an embedding
    for c in ["cluster", "protein_names", "gene_names", "organism",
              "ec_numbers", "go_labels", "keyword_labels"]:
        if c not in p.columns:
            p[c] = ""
        p[c] = p[c].fillna("").astype(str)

    # Two columns: `_stem_key` for coherence grouping (lowercased first two
    # meaningful tokens), and `_display` for the human-readable label
    # (preserved casing, full primary name minus parenthesized aliases).
    p["_stem_key"] = p["protein_names"].apply(name_key)
    p["_display"] = p["protein_names"].apply(display_name)

    print(f"  {len(p):,} entries · {p['cluster'].nunique():,} clusters")

    rows = []
    for cid, sub in p.groupby("cluster"):
        n = len(sub)
        stem_ctr = Counter(s for s in sub["_stem_key"].tolist() if s)
        dom_stem_key, dom_stem_n = ("", 0)
        for s, k in stem_ctr.most_common():
            if s not in STEM_STOPWORDS:
                dom_stem_key, dom_stem_n = s, k
                break
        stem_coh = dom_stem_n / n if n else 0

        # Pick the most-common FULL display name from members whose stem-key
        # matches the dominant. This yields readable, correctly-cased labels
        # like "Estrogen receptor" or "Cholesterol 24-hydroxylase" instead of
        # the terse "estrogen receptor" / "cholesterol 24-hydroxylase" pair.
        if dom_stem_key:
            _dom_group = sub[sub["_stem_key"] == dom_stem_key]
            _display_ctr = Counter(x for x in _dom_group["_display"].tolist() if x)
            dom_stem = _display_ctr.most_common(1)[0][0] if _display_ctr else dom_stem_key
        else:
            dom_stem = ""

        # GO labels — flatten across members
        go_flat = [g for row in semicolon_flat(sub["go_labels"]) for g in row]
        dom_go, dom_go_n = top_informative(go_flat, GO_HOUSEKEEPING)
        # go count is per-mention across members; convert to per-member coverage
        go_holders = sum(1 for row in semicolon_flat(sub["go_labels"])
                         if any(g == dom_go for g in row))
        go_coh = go_holders / n if n else 0

        kw_flat = [k for row in semicolon_flat(sub["keyword_labels"]) for k in row]
        dom_kw, _ = top_informative(kw_flat, KW_HOUSEKEEPING)
        kw_holders = sum(1 for row in semicolon_flat(sub["keyword_labels"])
                         if any(k == dom_kw for k in row))
        kw_coh = kw_holders / n if n else 0

        # EC top-level
        ec_top = []
        for ec_str in sub["ec_numbers"]:
            for e in str(ec_str).split(";"):
                e = e.strip()
                if e and e[0].isdigit():
                    ec_top.append(e.split(".")[0])
        ec_ctr = Counter(ec_top)
        dom_ec = ec_ctr.most_common(1)[0][0] if ec_ctr else ""
        ec_coh = ec_ctr.get(dom_ec, 0) / len(ec_top) if ec_top else 0

        outliers = sub[sub["_stem_key"] != dom_stem_key] if dom_stem_key else sub.iloc[0:0]

        # Case-study score: reward coherent clusters with SOME outliers of medium size
        n_outliers = len(outliers)
        # size sweet spot 15–200
        size_penalty = 0.0 if 15 <= n <= 200 else min(1.0, abs(n - 100) / 500.0)
        # want at least 3, no more than 30 outliers
        outlier_reward = min(n_outliers, 20) / 20 if 3 <= n_outliers <= 40 else 0.0
        score = stem_coh * 0.55 + go_coh * 0.15 + outlier_reward * 0.25 - size_penalty * 0.2

        rows.append({
            "cluster": cid,
            "size": n,
            "dominant_stem": dom_stem,
            "stem_coherence": round(stem_coh, 3),
            "dominant_go": dom_go,
            "go_coherence": round(go_coh, 3),
            "dominant_kw": dom_kw,
            "kw_coherence": round(kw_coh, 3),
            "dominant_ec_top": dom_ec,
            "ec_coherence": round(ec_coh, 3),
            "n_outliers": n_outliers,
            "case_study_score": round(score, 3),
            "outlier_accessions": ";".join(outliers["accession"].tolist()[:30]),
            "outlier_protein_names": ";".join(
                outliers["protein_names"].str[:60].tolist()[:30]
            ),
        })

    out = pd.DataFrame(rows).sort_values("case_study_score", ascending=False)
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv.name}  ({len(out)} clusters ranked)")

    # Human-readable shortlist
    lines = []
    lines.append("=" * 78)
    lines.append(f"CLUSTER CASE-STUDY SHORTLIST — top {TOP_N_SHORTLIST} of {len(out)} clusters")
    lines.append("=" * 78)
    lines.append("Ranking heuristic:")
    lines.append("  score = 0.55·stem_coherence + 0.15·go_coherence + 0.25·outlier_reward - 0.2·size_penalty")
    lines.append("  outlier_reward peaks when a cluster has 3–40 members whose protein-name")
    lines.append("  first token differs from the dominant stem (interesting to write about).")
    lines.append("")
    for _, r in out.head(TOP_N_SHORTLIST).iterrows():
        lines.append(f"--- CLUSTER {r['cluster']}  (size {r['size']}, score {r['case_study_score']}) ---")
        lines.append(f"  Dominant identity : {r['dominant_stem']}  "
                     f"({int(r['stem_coherence']*r['size'])}/{r['size']} = {r['stem_coherence']*100:.0f}%)")
        if r["dominant_go"]:
            lines.append(f"  Dominant GO term  : {r['dominant_go']}  ({r['go_coherence']*100:.0f}%)")
        if r["dominant_kw"]:
            lines.append(f"  Dominant keyword  : {r['dominant_kw']}  ({r['kw_coherence']*100:.0f}%)")
        if r["dominant_ec_top"]:
            lines.append(f"  Dominant EC       : EC {r['dominant_ec_top']}.-.-.-  ({r['ec_coherence']*100:.0f}% of EC-carriers)")
        lines.append(f"  Outliers ({r['n_outliers']}):")
        if r["outlier_protein_names"]:
            for pn in r["outlier_protein_names"].split(";")[:8]:
                if pn.strip():
                    lines.append(f"      · {pn}")
            if r["n_outliers"] > 8:
                lines.append(f"      · ... (+{r['n_outliers']-8} more)")
        lines.append("")

    out_txt.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_txt.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
