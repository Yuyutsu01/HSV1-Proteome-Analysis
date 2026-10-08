# Project Status: HSV_Computational_analysis

## Current Status
- **CURRENT PHASE**: **Phase 2 COMPLETE (Biological Preprocessing & Dataset Freeze)**
- **PHASE 2.5**: **Sequence Authenticity, Provenance, and Final Dataset Integrity Validation: COMPLETE**
- **PHASE 3**: **NOT STARTED (Locked pending Phase 2 sign-off)**

---

## Dataset Invariants & Accounting Summary

- **TOTAL CANONICAL SEQUENCES**: `22,689` (100% hash-deduplicated from 22,710 raw records)
- **SUPERVISED TEMPORAL DATASET**: `16,657` (100% evidence-grounded Immediate-Early, Early, or Late)
  - **IMMEDIATE_EARLY**: `1,552`
  - **EARLY**: `4,140`
  - **LATE**: `10,965`
- **UNKNOWN_ELIGIBLE (Curated Corpus)**: `5,356` (Biologically retained uncharacterized viral proteins)
- **CONFLICTING (Quarantined Tier C)**: `5` (Cross-file taxonomic conflicts)
- **EXCLUDED (Tier D)**: `671` (Non-HSV hosts [411], technical PDB [147], recombinant [50], synthetic [49], unassigned [14])
- **UNACCOUNTED SEQUENCES**: `0`

$$\text{Total Canonical } (22,689) = \text{Supervised } (16,657) + \text{Unknown Eligible } (5,356) + \text{Conflicting } (5) + \text{Excluded } (671)$$

---

## Core Scientific Constraints
1. **Raw FASTA Immutability**: All raw FASTA files in `./raw/` are strictly read-only and preserved with exact SHA-256 checksums.
2. **Retired 74-Protein Dataset**: The historical 74-protein HSV-1 Strain 17 dataset is retired and isolated in `./archive/legacy_74_protein_project/`.
3. **No Inference in Preprocessing**: Zero temporal labels were computationally inferred from sequence similarity, embeddings, or ML models.
4. **Scope Boundary**: No feature engineering, ProtBERT/ESM embeddings, or ML classifier training has occurred in Phases 1-2.5.
