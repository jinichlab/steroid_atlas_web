# Literature-recruited proteins in Nature's Steroid Atlas

This table lists every protein added to the atlas beyond the initial 
Rhea-based UniProt corpus. Each row documents the paper it came from, 
the exact evidence for including it, and where the sequence was fetched from.

Total: **15 entries** across 6 papers.

## Rimal 2024

**Rimal, B. et al. Bile salt hydrolase catalyses formation of amine-conjugated bile acids. Nature 626, 859–863 (2024)**  
DOI: [10.1038/s41586-023-06990-w](https://doi.org/10.1038/s41586-023-06990-w)

### `P0DXD2` — Bifidobacterium longum subsp. longum (NCTC 11818 (= ATCC 15707 / DSM 20219 / JCM 1217 / E194b))

- **Identifier type:** UniProt  
- **Length:** 317 aa  
- **Evidence:** Biochemical (purified enzyme + kinetics + mutagenesis)  
- **Sequence source:** UniProt REST (Swiss-Prot curated)  
- **In atlas UMAP:** yes  
- **Why recruited:** Primary purified enzyme (BlBSH). Biochemistry demonstrated the previously-unknown amine N-acyltransferase activity of BSH — conjugates amino acids to bile acids to form BBAAs. Structure-guided mutagenesis (C1A, R17A, N172A, R225A, T2A, W21F, T171A, K264A, D267A) mapped catalytic residues.

### `Q5LF84` — Bacteroides fragilis (NCTC 9343 (= ATCC 25285 / DSM 2151 / CCUG 4856 / JCM 11019 / LMG 10263 / VPI 2553))

- **Identifier type:** UniProt  
- **Length:** 331 aa  
- **Evidence:** Genetic (KO / complementation in vivo)  
- **Sequence source:** UniProt REST  
- **In atlas UMAP:** yes  
- **Why recruited:** In vivo genetic evidence for BBAA biosynthesis. B. fragilis Δbsh strain loses BBAA production; complementation with wild-type bsh restores it. Shows the amine N-acyltransferase activity operates in a live commensal gut bacterium.

## Guzior 2024

**Guzior, D.V. et al. Bile salt hydrolase acyltransferase activity expands bile acid diversity. Nature 626, 852–858 (2024)**  
DOI: [10.1038/s41586-024-07017-8](https://doi.org/10.1038/s41586-024-07017-8)

### `WP_243289361` — Clostridium perfringens (ATCC 13124)

- **Identifier type:** RefSeq_WP  
- **Length:** 329 aa  
- **Evidence:** Biochemical (kinetics + mutagenesis + structure)  
- **Sequence source:** NCBI efetch (RefSeq)  
- **In atlas UMAP:** PENDING (embed + project)  
- **Why recruited:** CpBSH/T — the only Guzior 2024 protein with direct biochemical evidence. Bifunctional: bile salt hydrolase (EC 3.5.1.24; pH 3-7) AND novel amine N-acyltransferase (EC 2.3.1.-; pH optimum 5.3). Peak acyl transfer = 7% of hydrolysis rate. Substrate scope: 16/20 amino acids from TCA, 12/20 from CA; Pro & Asp never conjugated. Mutants: C2A abolishes both activities (shared Cys2 nucleophile); N82Y preserves activity but shifts specificity (loses GluCA/LysCA/LeuCA, gains AlaCA). Structure: PDB 2BJG. Other Guzior 2024 Table 1 entries deliberately NOT added — no protein-level biochemical evidence in the paper.

## McCurry 2024

**McCurry, M.D. et al. Gut bacteria convert glucocorticoids into progestins in the presence of hydrogen gas. Cell 187, 2949-2963.e19 (2024)**  
DOI: [10.1016/j.cell.2024.05.005](https://doi.org/10.1016/j.cell.2024.05.005)

### `C8WL28` — Eggerthella lenta (DSM 2243 (= ATCC 25559 / CCUG 17323 / JCM 9979 / KCTC 3265 / NCTC 11813 / VPI 0255))

- **Identifier type:** UniProt  
- **Length:** 281 aa  
- **Evidence:** Genetic (cluster identification + heterologous activity)  
- **Sequence source:** UniProt REST (from 2009 E. lenta genome deposit, DOI 10.4056/sigs.33592)  
- **In atlas UMAP:** yes  
- **Why recruited:** Elen_2451 — one of four proteins in the elen_2451-2454 cluster responsible for 21-dehydroxylation of glucocorticoids. Cluster identified via ATc-inducible expression + LC-MS activity assay in E. lenta.

### `C8WL29` — Eggerthella lenta (DSM 2243)

- **Identifier type:** UniProt  
- **Length:** 203 aa  
- **Evidence:** Genetic (cluster identification)  
- **Sequence source:** UniProt REST (2009 E. lenta genome deposit)  
- **In atlas UMAP:** yes  
- **Why recruited:** Elen_2452 — ferredoxin (4Fe-4S) component of the elen_2451-2454 21-dehydroxylation cluster.

### `C8WL30` — Eggerthella lenta (DSM 2243)

- **Identifier type:** UniProt  
- **Length:** 909 aa  
- **Evidence:** Genetic (cluster identification)  
- **Sequence source:** UniProt REST (2009 E. lenta genome deposit)  
- **In atlas UMAP:** yes  
- **Why recruited:** Elen_2453 — molybdopterin oxidoreductase/dehydrogenase in the elen_2451-2454 21-dehydroxylation cluster.

### `C8WL31` — Eggerthella lenta (DSM 2243)

- **Identifier type:** UniProt  
- **Length:** 307 aa  
- **Evidence:** Genetic (cluster identification)  
- **Sequence source:** UniProt REST (2009 E. lenta genome deposit)  
- **In atlas UMAP:** yes  
- **Why recruited:** Elen_2454 — SPFH / band-7 / stomatin-prohibitin family membrane protein in the elen_2451-2454 21-dehydroxylation cluster.

## Jacoby 2025

**Jacoby, C. et al. (Cell Host & Microbe, 2025)**  
DOI: [10.1016/j.chom.2025.09.014](https://doi.org/10.1016/j.chom.2025.09.014)

### `MFU7515415` — Clostridium steroidoreducens (HCS.1)

- **Identifier type:** NCBI_Genome  
- **Length:** 254 aa  
- **Evidence:** Genetic (pathway identification)  
- **Sequence source:** NCBI genome (from paper's supplementary genome)  
- **In atlas UMAP:** yes  
- **Why recruited:** OsrA — one of three genes (OsrABC) in the oxidative steroid reduction pathway identified in C. steroidoreducens HCS.1.

### `MFU7516964` — Clostridium steroidoreducens (HCS.1)

- **Identifier type:** NCBI_Genome  
- **Length:** 620 aa  
- **Evidence:** Genetic (pathway identification)  
- **Sequence source:** NCBI genome  
- **In atlas UMAP:** yes  
- **Why recruited:** OsrB — component of the OsrABC pathway from C. steroidoreducens HCS.1.

### `MFU7517346` — Clostridium steroidoreducens (HCS.1)

- **Identifier type:** NCBI_Genome  
- **Length:** 645 aa  
- **Evidence:** Genetic (pathway identification)  
- **Sequence source:** NCBI genome  
- **In atlas UMAP:** yes  
- **Why recruited:** OsrC — component of the OsrABC pathway from C. steroidoreducens HCS.1.

## Arp 2025

**Arp, G. et al. (Nature Communications, 2025)**  
DOI: [10.1038/s41467-025-61425-6](https://doi.org/10.1038/s41467-025-61425-6)

### `dw0526` — Dysosmobacter welbionis (J115)

- **Identifier type:** Locus_Tag  
- **Length:** 652 aa  
- **Evidence:** Biochemical (activity in Arp 2025 Supp Data 1)  
- **Sequence source:** Arp 2025 Supp Data 1 (GCA_005121165.3 genome)  
- **In atlas UMAP:** yes  
- **Why recruited:** Δ4-3-ketosteroid 5β-reductase (652 aa). Bile acid ring reduction activity characterized in Arp 2025 Supplementary Data 1.

### `cp1309` — Clostridium paraputrificum (NCTC 11833)

- **Identifier type:** Locus_Tag  
- **Length:** 645 aa  
- **Evidence:** Biochemical (Arp 2025)  
- **Sequence source:** Arp 2025 Supp Data 1 (GCA_900447045.1 genome)  
- **In atlas UMAP:** yes  
- **Why recruited:** Δ4-3-ketosteroid 5β-reductase (645 aa) from Arp 2025 Supp Data 1.

### `mf2052` — Mediterraneibacter faecis (JMC 15917)

- **Identifier type:** Locus_Tag  
- **Length:** 1155 aa  
- **Evidence:** Biochemical (Arp 2025)  
- **Sequence source:** Arp 2025 Supp Data 2 (GCA_001312505.1 genome)  
- **In atlas UMAP:** yes  
- **Why recruited:** Fused 3β-HSDH / Δ5-4 isomerase (1155 aa fused enzyme) from Arp 2025 Supp Data 2.

### `mf0519` — Mediterraneibacter faecis (ATCC BAA-2716)

- **Identifier type:** Locus_Tag  
- **Length:** 641 aa  
- **Evidence:** Biochemical (Arp 2025)  
- **Sequence source:** Arp 2025 Supp Data 2 (GCA_000153925.1 genome)  
- **In atlas UMAP:** yes  
- **Why recruited:** Δ6-3-ketosteroid reductase (641 aa) from Arp 2025 Supp Data 2.

## Yao 2018 / see note

**Referenced in downstream MCBA literature (paper attribution to be re-verified)**  
DOI: [10.1038/s41467-026-68556-4](https://doi.org/10.1038/s41467-026-68556-4)

### `Q8A6H3` — Bacteroides thetaiotaomicron (VPI-5482)

- **Identifier type:** UniProt  
- **Length:** 259 aa  
- **Evidence:** Bioinformatic (needs re-audit)  
- **Sequence source:** UniProt REST  
- **In atlas UMAP:** yes  
- **Why recruited:** B. theta BSH — canonical gut Bacteroides bile salt hydrolase. Paper attribution flagged for re-verification in a future audit.
