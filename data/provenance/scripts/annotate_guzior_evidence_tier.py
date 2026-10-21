"""Annotate Guzior fetched records with evidence tier + full traceability strings.

Reads:  literature/guzior2024_review_ready.json
        literature/guzior2024_table1_normalized.json (for full-species paralog count)
Writes: back to literature/guzior2024_review_ready.json (in place)

Tiers:
  A         — biochemically characterized (purified enzyme + mutagenesis)
  B_strong  — sole BSH/T in species → species-culture activity unambiguously assigned
  B_moderate— sole Group I paralog when species has multiple paralogs → likely canonical
  B_weak    — one of ≥2 Group I paralogs in species
  C         — Group II or III paralog; activity not disambiguated from species-level fingerprint

Also adds:
  evidence_tier_short           — e.g. "A", "B_strong"
  evidence_tier_explanation     — the human-readable rationale (per row)
  traceability_chain            — plain-english provenance chain
"""
from __future__ import annotations
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NORM = HERE / "guzior2024_table1_normalized.json"
REVIEW = HERE / "guzior2024_review_ready.json"

DOI = "https://doi.org/10.1038/s41586-024-07017-8"

def classify_tier(acc: str, group: str, species: str, all_recs_by_species: dict) -> tuple[str, str]:
    if acc == "WP_243289361":
        return "A", (
            "Purified CpBSH/T from C. perfringens was assayed biochemically "
            "(kinetics: Km, Vmax; pH optimum 5.3; substrate scope) and heterologously "
            "expressed in E. coli DH5α with N82Y and C2A mutants (Fig 1, Fig 2, ED Fig 3). "
            "This is the ONLY protein in the paper with direct biochemical activity data."
        )
    species_recs = all_recs_by_species.get(species, [])
    n_all = len(species_recs)
    group_I_recs = [r for r in species_recs if r.get("bsh_t_group") == "I"]
    n_group_I = len(group_I_recs)

    if n_all == 1:
        return "B_strong", (
            f"This is the ONLY BSH/T detected in {species}'s genome. The species-culture "
            "MCBA profile is unambiguously attributable to this single enzyme."
        )
    if group == "I" and n_group_I == 1:
        return "B_moderate", (
            f"{species} has {n_all} BSH/T paralog(s), of which this is the sole Group I "
            "(canonical BSH/T clade). Group I sequences are the most conserved and best-"
            "characterized BSH/T family; likely the primary carrier of activity, but "
            f"contribution from the {n_all - 1} non-Group-I paralog(s) is not ruled out."
        )
    if group == "I":
        return "B_weak", (
            f"{species} has {n_group_I} Group I paralogs. The species's MCBA profile "
            "cannot be attributed to this specific paralog without biochemical isolation."
        )
    return "C", (
        f"This is a Group {group} paralog in {species}. Guzior 2024 defines Group I as "
        f"the canonical/well-characterized BSH/T clade; Group {group} sequences are more "
        "divergent and may or may not contribute to the observed species-culture activity. "
        "In vivo contribution is not disambiguated from the species-level fingerprint."
    )

def main() -> int:
    all_recs = json.loads(NORM.read_text())
    review = json.loads(REVIEW.read_text())

    # Group all records by species (using species field from all_recs)
    by_species = defaultdict(list)
    for r in all_recs:
        if r.get("species"):
            by_species[r["species"]].append(r)

    tier_counts = defaultdict(int)
    for rec in review:
        tier, explanation = classify_tier(
            rec["accession_normalized"],
            rec.get("bsh_t_group") or "",
            rec.get("species") or "",
            by_species,
        )
        rec["evidence_tier_short"] = tier
        rec["evidence_tier_explanation"] = explanation
        rec["traceability_chain"] = (
            f"Bacterial strain: {rec['species']} strain {rec['strain']} "
            f"(source: {rec['source']}). "
            f"Species cultured with bile acids + amino acids; LC-MS gave species-level "
            f"MCBA profile cluster '{rec['mcba_profile_cluster']}'. "
            f"Genome mining of assembly {rec['genome_accession']} identified BSH/T "
            f"paralog {rec['protein_accession_raw']} (this row), assigned to "
            f"phylogenetic Group {rec['bsh_t_group']}. "
            f"Sequence source database: {rec['accession_type']}. "
            f"Evidence tier for protein-level activity attribution: {tier}."
        )
        # Enrich the proposed atlas row too
        rec["atlas_row_proposed"]["Evidence_Tier"] = tier
        rec["atlas_row_proposed"]["Sequence_Source"] = (
            f"{rec['accession_type']} accession {rec['accession_normalized']} — "
            f"literature-recruited via Guzior 2024 Supp Table 1. "
            f"Species-culture MCBA cluster: {rec['mcba_profile_cluster']}. "
            f"BSH/T phylogenetic group: {rec['bsh_t_group']}. "
            f"Evidence tier: {tier} — {explanation}"
        )
        tier_counts[tier] += 1

    REVIEW.write_text(json.dumps(review, indent=2))
    print(f"Wrote {REVIEW.name}  ({len(review)} enriched entries)")
    print("\nEvidence tier distribution:")
    for t, n in sorted(tier_counts.items()):
        print(f"  {t:12s}  {n:2d}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
