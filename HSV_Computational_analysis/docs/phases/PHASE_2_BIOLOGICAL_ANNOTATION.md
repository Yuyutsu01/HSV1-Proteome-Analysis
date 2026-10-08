# Phase 2: Biological Annotation & Temporal Expression Curation

## 1. Objective
Establish an evidence-grounded biological annotation layer for all 22,689 unique HSV sequences without relying on sequence similarity, machine learning, ProtBERT embeddings, or downstream computational inferences.

## 2. Evidence Hierarchy
- **Level 1 (Direct Experimental Evidence)**: Kinetic transcriptional inhibition assays (cycloheximide / phosphonoacetic acid) and replication kinetics.
- **Level 2 (Authoritative Curated Database)**: Curated UniProtKB/Swiss-Prot and RefSeq annotations.
- **Level 3 (Primary Literature)**: Peer-reviewed primary virology studies with explicit kinetic characterization.
- **Level 4 (Secondary Literature)**: Review articles and structural summaries.
- **Level 5 (Insufficient Evidence)**: Uncharacterized/hypothetical viral proteins where temporal kinetics are not established (`UNKNOWN`).

## 3. Dataset Tiers
- **Tier A (High-Confidence Labeled)**: `10,482` sequences. Primary eligible full-length HSV proteins with Level 1 experimental temporal evidence.
- **Tier B (Verified Unlabeled / Sequence-Only)**: `11,531` sequences. Valid HSV biological sequences (fragments or uncharacterized proteins).
- **Tier C (Conflicting / Uncertain)**: `5` sequences. Quarantined cross-file taxonomic conflicts.
- **Tier D (Excluded)**: `671` sequences. Excluded host, PDB, synthetic, and recombinant entities.

## 4. Ground-Truth Traceability
Every labeled sequence references its exact literature evidence record in `results/tables/phase2_annotation_evidence.csv`.
