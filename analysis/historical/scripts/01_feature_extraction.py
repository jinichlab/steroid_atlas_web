"""Feature extraction for the Steroid Atlas neighborhood-enrichment analysis.

Reads:  ../data/proteins.csv (35,349 rows × 23 cols)
Writes: ../analysis/annotations_matrix.tsv (one row per protein × extracted features)
        ../analysis/sequence_diversity_summary.txt (aggregate diversity stats)
        ../analysis/feature_coverage.tsv (per-feature % populated)

Features extracted per protein:
  accession
  cluster              (existing k-means k=95)
  umap_1, umap_2       (existing coordinates)

  # Function annotations
  ec_full              (raw EC field, first entry only if multiple)
  ec_top_level         (first digit of EC — 1..7)
  ec_subclass          (first two levels, e.g., "1.14")
  ec_sub_subclass      (first three levels, e.g., "1.14.14")

  # Name-derived function bucket (fallback for unpopulated EC)
  protein_name_first_word  (first meaningful noun-phrase word)
  is_p450              (regex CYP\d, P450, cytochrome)
  cyp_family           (CYP + family digit, e.g., "CYP17")
  cyp_subfamily        (e.g., "CYP17A")
  is_bsh               (bile salt hydrolase / choloylglycine hydrolase)
  is_hsd               (hydroxysteroid dehydrogenase)
  is_reductase         (reductase / dehydrogenase)
  is_receptor          (nuclear/membrane receptor)
  is_transporter       (transporter / carrier / ABC / SLC)

  # Taxonomy
  genus                (first word of organism)
  species_binomial     (first two words of organism)

  # Sequence-level properties
  length_aa            (existing)
  shannon_entropy      (per-sequence amino-acid Shannon entropy — proxy for compositional diversity)

  # Provenance
  is_literature_recruited
  has_ec               (boolean)
  has_rhea             (boolean)
"""
from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # analysis/
from _paths import ROOT, output  # noqa: E402
IN = ROOT / "data" / "proteins.csv"
OUT_MATRIX = output("annotations_matrix.tsv", historical=True)
OUT_DIVERSITY = output("sequence_diversity_summary.txt", historical=True)
OUT_COVERAGE = output("feature_coverage.tsv", historical=True)

# ─── Regex patterns ──────────────────────────────────────────────────────────
CYP_GENE_RE = re.compile(r"\b(?:cyp|CYP)(\d+)([A-Z]?)(\d*)\b", re.IGNORECASE)
CYP_NAME_RE = re.compile(r"\bP\s*[-]?\s*450\b|cytochrome\s*P.?450|Cytochrome\s+P450", re.IGNORECASE)
BSH_RE = re.compile(r"bile\s+salt\s+hydrolase|choloylglycine\s+hydrolase|BSH/T?", re.IGNORECASE)
HSD_RE = re.compile(r"hydroxysteroid\s+dehydrogenase|HSD\d|17-?beta-?HSD|3-?(alpha|beta)-?HSD", re.IGNORECASE)
REDUCTASE_RE = re.compile(r"\breductase\b|\bdehydrogenase\b|\boxidoreductase\b", re.IGNORECASE)
RECEPTOR_RE = re.compile(r"receptor|nuclear\s+receptor|steroid[- ]binding", re.IGNORECASE)
TRANSPORTER_RE = re.compile(r"transporter|carrier|ABC|SLC\d|permease|efflux", re.IGNORECASE)
EC_RE = re.compile(r"(\d+)\.(\d+)\.([\d\-]+)\.(\d+|-)")

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Words to skip when finding the first "meaningful" word in a protein name
NAME_STOPWORDS = {
    "the", "a", "an", "of", "for", "from", "protein", "putative",
    "probable", "uncharacterized", "hypothetical", "predicted", "family",
    "domain", "subunit", "type", "class", "isoform", "member",
}

def parse_ec(ec_field: str) -> tuple[str, str, str, str]:
    """Return (full, top, sub, subsub) — empty strings if no EC found."""
    if not isinstance(ec_field, str) or not ec_field.strip():
        return "", "", "", ""
    m = EC_RE.search(ec_field)
    if not m:
        return "", "", "", ""
    a, b, c, d = m.group(1), m.group(2), m.group(3), m.group(4)
    return f"{a}.{b}.{c}.{d}", a, f"{a}.{b}", f"{a}.{b}.{c}"

def cyp_designation(gene_field: str, name_field: str) -> tuple[str, str, bool]:
    """Extract (cyp_family, cyp_subfamily, is_p450)."""
    is_p450 = False
    family, subfamily = "", ""
    for f in (gene_field, name_field):
        if not isinstance(f, str):
            continue
        m = CYP_GENE_RE.search(f)
        if m:
            family = f"CYP{m.group(1)}"
            subfamily = f"CYP{m.group(1)}{m.group(2)}" if m.group(2) else family
            is_p450 = True
            break
    if not is_p450 and isinstance(name_field, str) and CYP_NAME_RE.search(name_field):
        is_p450 = True
    return family, subfamily, is_p450

def protein_name_first_word(name_field: str) -> str:
    if not isinstance(name_field, str) or not name_field.strip():
        return ""
    # Take part before the first opening parenthesis or bracket
    core = re.split(r"[\(\[]", name_field, maxsplit=1)[0].strip()
    tokens = re.findall(r"[A-Za-z][A-Za-z0-9\-]+", core)
    for tok in tokens:
        if tok.lower() not in NAME_STOPWORDS and len(tok) >= 3:
            return tok
    return tokens[0] if tokens else ""

def genus_of(org_field: str) -> tuple[str, str]:
    if not isinstance(org_field, str) or not org_field.strip():
        return "", ""
    # Strip brackets from taxonomy like [Clostridium]
    org = org_field.replace("[", "").replace("]", "").strip()
    parts = org.split()
    if not parts:
        return "", ""
    genus = parts[0]
    species = " ".join(parts[:2]) if len(parts) >= 2 else genus
    return genus, species

def shannon_entropy_aa(seq: str) -> float:
    if not isinstance(seq, str) or not seq:
        return 0.0
    seq_up = seq.upper()
    counts = Counter(c for c in seq_up if c in AA_ALPHABET)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    H = 0.0
    for n in counts.values():
        p = n / total
        H -= p * math.log2(p)
    return H

def bool_match(pattern: re.Pattern, s: str) -> bool:
    if not isinstance(s, str):
        return False
    return bool(pattern.search(s))

# ─── Main ───────────────────────────────────────────────────────────────────
def main() -> int:
    print(f"Loading {IN}...")
    df = pd.read_csv(IN, low_memory=False)
    n = len(df)
    print(f"  {n:,} proteins")

    rows: list[dict] = []
    for _, r in df.iterrows():
        acc = r.get("accession", "")
        ec_full, ec_top, ec_sub, ec_subsub = parse_ec(str(r.get("ec_numbers", "")))
        cyp_fam, cyp_sub, is_p450 = cyp_designation(
            str(r.get("gene_names", "")), str(r.get("protein_names", ""))
        )
        genus, sp = genus_of(str(r.get("organism", "")))
        name_first = protein_name_first_word(str(r.get("protein_names", "")))
        seq = str(r.get("sequence", ""))
        rows.append({
            "accession": acc,
            "cluster": r.get("cluster"),
            "umap_1": r.get("umap_1"),
            "umap_2": r.get("umap_2"),

            "ec_full": ec_full,
            "ec_top_level": ec_top,
            "ec_subclass": ec_sub,
            "ec_sub_subclass": ec_subsub,

            "protein_name_first_word": name_first,
            "is_p450": is_p450,
            "cyp_family": cyp_fam,
            "cyp_subfamily": cyp_sub,
            "is_bsh": bool_match(BSH_RE, str(r.get("protein_names", ""))),
            "is_hsd": bool_match(HSD_RE, str(r.get("protein_names", ""))),
            "is_reductase": bool_match(REDUCTASE_RE, str(r.get("protein_names", ""))),
            "is_receptor": bool_match(RECEPTOR_RE, str(r.get("protein_names", ""))),
            "is_transporter": bool_match(TRANSPORTER_RE, str(r.get("protein_names", ""))),

            "genus": genus,
            "species_binomial": sp,

            "length_aa": r.get("length_aa"),
            "shannon_entropy": round(shannon_entropy_aa(seq), 4) if seq else 0.0,

            "is_literature_recruited": r.get("is_literature_recruited"),
            "has_ec": bool(ec_full),
            "has_rhea": bool(
                isinstance(r.get("rhea_reactions"), str) and r.get("rhea_reactions").strip()
            ),
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_MATRIX, sep="\t", index=False)
    print(f"\nWrote annotations matrix: {OUT_MATRIX.name}  ({len(out):,} rows × {len(out.columns)} cols)")

    # ─── Feature coverage summary ───────────────────────────────────────────
    cov_records = []
    for c in out.columns:
        if c in ("accession", "cluster", "umap_1", "umap_2"):
            continue
        v = out[c]
        if v.dtype == bool:
            n_true = int(v.sum())
            n_pop = len(v)
            pct = n_true / n_pop * 100
            cov_records.append({
                "feature": c,
                "n_populated": n_pop,
                "n_true_or_nonblank": n_true,
                "pct_true_or_nonblank": round(pct, 2),
                "unique_values": 2,
            })
        else:
            v_str = v.astype(str).str.strip()
            n_pop = int((v_str != "").sum() - (v_str == "nan").sum())
            n_unique = v_str[v_str != ""].nunique()
            cov_records.append({
                "feature": c,
                "n_populated": n_pop,
                "n_true_or_nonblank": n_pop,
                "pct_true_or_nonblank": round(n_pop / len(v) * 100, 2),
                "unique_values": n_unique,
            })
    cov_df = pd.DataFrame(cov_records)
    cov_df.to_csv(OUT_COVERAGE, sep="\t", index=False)
    print(f"Wrote feature-coverage report: {OUT_COVERAGE.name}")
    print()
    print("=== Feature coverage ===")
    print(cov_df.to_string(index=False))

    # ─── Sequence diversity ─────────────────────────────────────────────────
    lengths = out["length_aa"].dropna()
    entropies = out["shannon_entropy"].dropna()
    entropies = entropies[entropies > 0]

    diversity_txt = []
    diversity_txt.append("=== Sequence diversity summary ===\n")
    diversity_txt.append(f"Total proteins: {n:,}\n")

    diversity_txt.append("\n--- Length (aa) ---")
    diversity_txt.append(f"  min:    {int(lengths.min())}")
    diversity_txt.append(f"  25%:    {int(lengths.quantile(0.25))}")
    diversity_txt.append(f"  median: {int(lengths.median())}")
    diversity_txt.append(f"  mean:   {lengths.mean():.1f}")
    diversity_txt.append(f"  75%:    {int(lengths.quantile(0.75))}")
    diversity_txt.append(f"  max:    {int(lengths.max())}")

    diversity_txt.append("\n--- Amino-acid Shannon entropy per sequence ---")
    diversity_txt.append("  (theoretical max for 20 AAs = 4.32 bits)")
    diversity_txt.append(f"  min:    {entropies.min():.3f}")
    diversity_txt.append(f"  25%:    {entropies.quantile(0.25):.3f}")
    diversity_txt.append(f"  median: {entropies.median():.3f}")
    diversity_txt.append(f"  mean:   {entropies.mean():.3f}")
    diversity_txt.append(f"  75%:    {entropies.quantile(0.75):.3f}")
    diversity_txt.append(f"  max:    {entropies.max():.3f}")

    n_genera = out["genus"].astype(str).replace("", pd.NA).dropna().nunique()
    n_species = out["species_binomial"].astype(str).replace("", pd.NA).dropna().nunique()
    diversity_txt.append("\n--- Taxonomic diversity ---")
    diversity_txt.append(f"  Distinct genera:   {n_genera:,}")
    diversity_txt.append(f"  Distinct species:  {n_species:,}")

    diversity_txt.append("\n  Top 15 genera by count:")
    top_g = out["genus"].astype(str).replace("", pd.NA).dropna().value_counts().head(15)
    for g, c in top_g.items():
        diversity_txt.append(f"    {c:>5,}  {g}")

    diversity_txt.append("\n--- Functional annotation coverage ---")
    for flag in ("is_p450", "is_bsh", "is_hsd", "is_reductase", "is_receptor", "is_transporter"):
        n_true = int(out[flag].sum())
        diversity_txt.append(f"  {flag:<18}  {n_true:>5,}  ({n_true/n*100:.1f}%)")

    diversity_txt.append("\n--- EC number availability ---")
    diversity_txt.append(f"  With any EC:         {int(out['has_ec'].sum()):>5,}  ({out['has_ec'].mean()*100:.1f}%)")
    diversity_txt.append(f"  With Rhea reactions: {int(out['has_rhea'].sum()):>5,}  ({out['has_rhea'].mean()*100:.1f}%)")

    diversity_txt.append("\n--- EC top-level distribution (for proteins with EC) ---")
    ec_top = out["ec_top_level"].astype(str).replace("", pd.NA).dropna().value_counts().sort_index()
    ec_labels = {
        "1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases",
        "4": "Lyases", "5": "Isomerases", "6": "Ligases", "7": "Translocases",
    }
    for k, v in ec_top.items():
        diversity_txt.append(f"    EC {k}.-.-.- {ec_labels.get(k, '?'):<20}  {v:>5,}")

    diversity_txt.append("\n--- Cluster size distribution ---")
    csz = out["cluster"].value_counts()
    diversity_txt.append(f"  n_clusters: {len(csz)}")
    diversity_txt.append(f"  min size:   {int(csz.min())}")
    diversity_txt.append(f"  median size:{int(csz.median())}")
    diversity_txt.append(f"  mean size:  {csz.mean():.1f}")
    diversity_txt.append(f"  max size:   {int(csz.max())}")

    OUT_DIVERSITY.write_text("\n".join(diversity_txt) + "\n")
    print(f"\nWrote sequence-diversity summary: {OUT_DIVERSITY.name}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
