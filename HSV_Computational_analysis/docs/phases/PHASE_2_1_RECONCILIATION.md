# Phase 2.1: Biological Annotation & Ground-Truth Reconciliation

## 1. Executive Summary & Objective
Phase 2.1 formally decouples two critical attributes that were previously conflated:
1. **Temporal Annotation Status**: Whether a viral sequence has an evidence-backed biological temporal classification (`LABELED` vs `UNKNOWN` vs `CONFLICTING`).
2. **Ground-Truth Modeling Eligibility**: Whether that sequence satisfies all criteria required for the primary supervised ground-truth cohort (`ELIGIBLE_GROUND_TRUTH` vs `NOT_ELIGIBLE_GROUND_TRUTH` vs `CONFLICTED`).

> **Core Principle**: Temporal annotation and supervised ground-truth eligibility are treated as separate attributes.

## 2. Revised Dataset Tiers
- **Tier A (`TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH`)**: `10,482` sequences. Primary supervised dataset consisting of full-length, biologically verified HSV proteins with high-confidence Level 1 experimental evidence.
- **Tier B (`TIER_B_VERIFIED_EXTENDED_BIOLOGICAL_CORPUS`)**: `11,531` sequences. Verified HSV biological corpus containing both labeled sequences (`6,175` sequences) and uncharacterized/unknown sequences (`5,356` sequences).
- **Tier C (`TIER_C_CONFLICTING_QUARANTINED`)**: `5` sequences. Quarantined cross-file taxonomic conflicts.
- **Tier D (`TIER_D_EXCLUDED`)**: `671` sequences. Excluded non-HSV hosts, PDB chains, recombinant vectors, and synthetic constructs.

## 3. Evidence Provenance & Propagation
- **Direct Sequence Evidence**: `12` sequences.
- **Gene/Protein-Level Evidence**: `16,645` sequences.
- **Insufficient Evidence (`UNKNOWN`)**: `6,027` sequences.
- **Conflicting Evidence (`CONFLICTING`)**: `5` sequences.
