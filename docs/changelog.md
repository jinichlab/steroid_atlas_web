# Changelog

## v0.2.2 — 2026-10-12 — Steroid set sync + natural/synthetic UMAP refit

Applied by `analysis/scripts/18_sync_steroid_sets.py` (summary in `analysis/results/natsyn_umap_report.txt`).

- **Duplicate removed from `molecules.csv`:** 5β-dihydrotestosterone (`CHEBI:2150`) was listed twice, once from Rhea and once as a literature recruitment. Merged into the Rhea row (which carries the interacting proteins), now flagged `is_literature_recruited=1` with the PubChem source link. Small-molecule atlas: 677 → **676**.
- **Natural set rebuilt from `molecules.csv`:** the `natural` rows of `natural_synthetic_steroids.csv` (589) came from an older molecule list. They lacked the 91 distinct literature-recruited compounds and still held the 4 alkaloids dropped in v0.2.1. Natural is now exactly the 676 molecules; the 2,300 synthetic rows are unchanged. Total: 2,889 → **2,976**.
- **Natural + synthetic UMAP refit:** the original reducer was never saved, so the new compounds couldn't be projected into it. All 2,976 compounds were re-embedded from RDKit Morgan fingerprints (radius 2, 2048 bits) with `umap.UMAP(n_neighbors=15, min_dist=0.1, metric="jaccard", random_state=42)`. The layout differs from earlier figures. The fitted reducer is saved to `analysis/results/natsyn_umap.joblib` (gitignored; the script regenerates it deterministically) for future `.transform()` projections.

## v0.2.1 — 2026-07-24 — Alkaloid cleanup

Dropped 4 benzo[c]phenanthridine alkaloids (sanguinarine, dihydrosanguinarine, chelirubine, dihydrochelirubine — `CHEBI:17183`, `CHEBI:17209`, `CHEBI:17031`, `CHEBI:17789`) from `molecules.csv`. These plant secondary metabolites were captured by the permissive sterane substructure filter (their pentacyclic ring system matches sterane after N→C substitution) but are not steroids by any standard biochemical definition. None had associated proteins.

Small-molecule atlas: 681 → **677 unique molecules**.

## v0.2.0 — 2026-07-24 — Scope tightening

- **Dropped 485 legacy bulk-curated entries** that were off-scope for the atlas's steroid-metabolism focus. These came from three earlier bulk-curation sources: 374 ChEBI cross-references (mostly animal transporters and vertebrate P450s), 65 PDB structural cross-references (mostly Human signaling proteins with a bound sterol in their crystal), 40 manual curation entries (mostly Human nuclear receptors), and 6 composite-source entries. None had specific paper-level provenance and all diluted the metabolic/enzymatic focus of the atlas.
- Atlas now: **35,349 unique protein sequences** (was 35,834), still with **15 literature-recruited STAR entries** fully preserved.
- Methods narrative simplified to two automated evidence sources (Rhea-catalytic + UniProt steroid-binding) plus the literature audit.
- Data dictionary in `data/README.md` reflects the current 21-column schema (dropped `uniprot_url` / `alphafold_url`, which are trivially derivable from the accession).

## v0.1.0 — 2026-07-23 — Initial public release

### Literature audit pass
- Verified paper attribution for every `is_literature_recruited=1` entry against source-paper text via RAG search
- Removed `A0A8F5DVT9` (wrong Guzior 2024 attribution for *C. minuta* — paper doesn't cover this species)
- Upgraded `P0DXD2` (Rimal 2024 BlBSH) annotation from garbage `"5"` to full citation + role description
- Cleared stale "should be dropped" note from `Q5LF84` (Rimal 2024 B. fragilis BSH) — verified as legitimate KO/complementation host
- Added `WP_243289361` (Guzior 2024 CpBSH/T) — the only Guzior 2024 protein with direct biochemical evidence; ProtT5-embedded and UMAP-projected via cosine-nearest-neighbor (nearest = P54965, cosine sim = 0.9994)

### Strict EC-number curation policy applied
- `P0DXD2`: `3.5.1.24;3.5.1.74` → `2.3.1.-;3.5.1.-` (only Rimal-experimentally-confirmed ECs)
- `Q5LF84`: blank → `2.3.1.-` (KO-evidence-supported)
- `WP_243289361`: `EC 3.5.1.24; EC 2.3.1.-` → `3.5.1.24;2.3.1.-` (Guzior direct kinetics; format normalized)

### ChEBI ID cleanup
- Fixed 6 wrong ChEBI IDs in `molecules.csv`:
  - `CHEBI:34964` → `11909` (isoallopregnanolone; was S-guanyl-cysteine adduct)
  - `CHEBI:1156` → `34461` (THDOC; was 2-hydroxyestrone)
  - `CHEBI:27725` → `30154` (5β-DHP; was butin)
  - `CHEBI:34958` → `16229` (epipregnanolone; was glutathione derivative)
  - `CHEBI:16718` → `1712` (pregnanolone; was reticuline alkaloid)
  - `CHEBI:34979` → `2150` (5β-DHT; was stanolone benzoate)
- Normalized format across atlas: no `CHEBI:` prefix, no `.0` float leak — pure numeric IDs
- Full audit: 594 ChEBI IDs verified against the ChEBI ontology, 0 real mismatches remaining

### Repo structure for public release
- Clean data/ folder with snake_case column names + data dictionary
- app/ with launcher script + README
- literature/ with full audit pipeline preserved
- docs/ with methodology + changelog
- MIT LICENSE, CITATION.cff, requirements.txt
