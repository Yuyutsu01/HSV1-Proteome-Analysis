# Phase 2 Final Dataset Freeze: Temporal Class Extraction, Cleaning, and Dataset Freeze

## 1. Overview & Objective

This document formalizes the **Final Dataset Freeze** of the biological preprocessing, eligibility curation, and temporal-expression annotation pipeline for the **canonical HSV protein sequence corpus** ($N=22,689$).

The objective of this phase is to produce the final clean, frozen, and reproducible temporal-class supervised dataset and the curated extended biological corpus from the researcher-supplied FASTA files.

> [!IMPORTANT]
> The supervised dataset contains only FASTA-derived HSV protein sequences with biologically supported **Immediate-Early**, **Early**, or **Late** temporal annotations.
> 
> The 22,689-sequence corpus is referred to as the **curated HSV protein sequence dataset** or **canonical HSV protein sequence corpus**, reflecting its diverse collection of laboratory reference strains, clinical isolates, and sequence fragments.

---

## 2. End-to-End Pipeline Summary

### 2.1. Raw Ingestion & Deduplication (Phase 1A)
- **Raw FASTA Files**: `raw/HSV2_non_redundant (2).fasta` ($10,540$ records) and `raw/non_redundant.fasta` ($12,170$ records).
- **Total Raw Records**: `22,710`.
- **Exact Duplicate Collapsing**: $21$ cross-file redundant sequence instances collapsed deterministically into canonical SHA-256 identifiers (`CANONICAL_<hash>`).
- **Canonical Sequence Population**: `22,689` unique sequences.

### 2.2. Biological Eligibility & Conflict Curation (Phases 1B & 1B.1)
- **Biologically Eligible Sequences**: `22,013` sequences ($13,321$ full-length primary + $8,692$ sequence-analysis-only fragments).
- **Excluded Non-HSV / Technical Sequences**: `671` sequences ($411$ non-herpesvirus hosts, $147$ PDB technical crystallization chains, $50$ recombinant/engineered mutants, $49$ synthetic patent constructs, $14$ unassigned).
- **Quarantined Cross-File Conflicts**: `5` sequences preserved in manual review quarantine (Tier C).

### 2.3. Temporal Expression & Provenance Auditing (Phases 2, 2.1, 2.1.1, 2.1.2A)
- Temporal kinetic expression classes were mapped from primary experimental herpesvirus literature (e.g., *Honess & Roizman 1974*, *Conley et al. 1981*, *McGeoch et al. 1988*, *Dolan et al. 1998*).
- Distinct evidentiary standards were enforced:
  1. Biological identity sources (NCBI GenBank, RefSeq, UniProtKB) verify sequence identity.
  2. Primary experimental time-course literature verifies kinetic expression timing.
  3. Target-species validation was verified across 28 core viral genes with independent HSV-2 experimental time-course publications.

---

## 3. Final Dataset Structure

### A. Final Supervised Temporal Dataset ($N=16,657$)
- **File Path**: [final_temporal_supervised_dataset.csv](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV_Computational_analysis/data/processed/final_temporal_supervised_dataset.csv)
- **Scope**: Contains only biologically verified sequences with unambiguous temporal classes (`IMMEDIATE_EARLY`, `EARLY`, `LATE`).
- **Columns**: `canonical_id`, `sequence`, `length`, `species`, `strain`, `gene`, `protein`, `temporal_class`, `temporal_evidence_source`, `temporal_evidence_type`, `temporal_confidence`.
- **Zero Inclusions**: 0 UNKNOWN, 0 CONFLICTING, 0 Excluded.

### B. Final Curated HSV Protein Corpus ($N=22,013$)
- **File Path**: [final_curated_hsv_protein_corpus.csv](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV_Computational_analysis/data/processed/final_curated_hsv_protein_corpus.csv)
- **Scope**: Contains all biologically retained sequences surviving Phase 1 curation ($16,657$ labeled + $5,356$ uncharacterized sequences).
- **Columns**: `canonical_id`, `sequence`, `length`, `species`, `strain`, `gene`, `protein`, `identity_status`, `completeness`, `temporal_class`, `temporal_evidence_source`, `temporal_evidence_type`, `temporal_confidence`, `temporal_annotation_status`, `biological_tier`, `biological_eligibility`.

---

## 4. Exact Final Counts & Cross-Tabulations

### 4.1. Temporal Class Distribution ($N=22,689$)

| Temporal Class | Total Sequences | HSV-1 Count | HSV-2 Count | Unique Genes | Unique Proteins |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **IMMEDIATE_EARLY** | `1,552` | `986` | `566` | `10` | `43` |
| **EARLY** | `4,140` | `2,344` | `1,796` | `16` | `223` |
| **LATE** | `10,965` | `5,670` | `5,295` | `66` | `1,045` |
| **UNKNOWN** | `6,027` | `2,802` | `2,606` | `1` | `47` |
| **CONFLICTING** | `5` | `0` | `0` | `5` | `5` |
| **Total** | **`22,689`** | **`11,802`** | **`10,263`** | **`92`** | **`1,311`** |

### 4.2. Species $\times$ Temporal Class Breakdown

| Species | Immediate-Early | Early | Late | Unknown | Conflicting | Total |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **HSV-1** | `986` | `2,344` | `5,670` | `2,802` | `0` | **`11,802`** |
| **HSV-2** | `566` | `1,796` | `5,295` | `2,606` | `0` | **`10,263`** |
| **Other / Non-HSV** | `0` | `0` | `0` | `619` | `5` | **`624`** |
| **Total** | **`1,552`** | **`4,140`** | **`10,965`** | **`6,027`** | **`5`** | **`22,689`** |

### 4.3. Exhaustive Status Reconciliation

| Final Status | Sequence Count | Description |
| :--- | :--- | :--- |
| **`SUPERVISED_TEMPORAL`** | `16,657` | Biologically retained with supported IE, Early, or Late labels |
| **`UNKNOWN_TEMPORAL`** | `5,356` | Biologically retained uncharacterized viral proteins (Tier B) |
| **`CONFLICTING_TEMPORAL`** | `5` | Quarantined cross-file taxonomic conflicts (Tier C) |
| **`EXCLUDED`** | `671` | Non-HSV host / technical / synthetic / engineered (Tier D) |
| **Total Partition** | **`22,689`** | **Zero orphan sequences; 100% data closure** |

---

## 5. Verification & Test Suite Summary

- **Total Test Cases**: **89 passed** across all 9 test modules:
  - `test_final_temporal_dataset.py` (14/14 passed)
  - `test_phase1.py` (8/8 passed)
  - `test_phase1b_biological_curation.py` (11/11 passed)
  - `test_phase1b1_conflict_resolution.py` (3/3 passed)
  - `test_phase2_annotation.py` (8/8 passed)
  - `test_phase2_1_reconciliation.py` (12/12 passed)
  - `test_phase2_1_1_provenance.py` (12/12 passed)
  - `test_phase2_1_2_cross_species.py` (12/12 passed)
  - `test_phase2_1_2A_target_species.py` (9/9 passed)

---

## 6. Stop Condition

> [!IMPORTANT]
> **Biological Preprocessing is Complete and Frozen.**
> 
> Phase 3 (Sequence & Representation Engineering: physicochemical profiling, ProtBERT / ESM embeddings, homology partitioning, and supervised ML model training) will commence in the next standalone execution upon user instruction.
