# Computational Representation and Classification Analysis of the HSV-1 Proteome

## Central Documentation Archive (Phases 1–10)

This directory provides the authoritative, permanent, human-readable record of all completed datasets, methodologies, validation audits, and exploratory analyses for the HSV-1 proteome study.

---

### Key Project Specifications

| Attribute | Value |
| :--- | :--- |
| **Project Title** | Computational Representation and Classification Analysis of the HSV-1 Proteome |
| **Target Organism** | Human alphaherpesvirus 1 (HSV-1) |
| **Reference Genome** | HSV-1 strain 17 (`NC_001806.2`) |
| **Dataset Size** | **74 unique proteins** (from 77 raw translated CDS records; 3 exact diploid duplicates removed) |
| **Temporal Classes** | **Immediate-Early (IE)**: 5 \| **Early (E)**: 15 \| **Late (L)**: 54 |
| **Current Completed Phase** | **Phase 10: Exploratory Representation Analysis** |
| **Next Planned Phase** | **Phase 11: Unsupervised Clustering & Validation** |

---

### Documentation Map

#### Core Documents
1. [Project Overview](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/PROJECT_OVERVIEW.md) — Scientific rationale, comparative framework, and explicit research boundaries.
2. [Dataset & Annotation](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/DATASET_AND_ANNOTATION.md) — NCBI ingestion, exact deduplication mappings, length profiles, and ground-truth temporal classes.
3. [Methodology](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/METHODOLOGY.md) — Complete end-to-end computational biology workflow and zero-leakage protocols.
4. [Reproducibility](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/REPRODUCIBILITY.md) — Deterministic seeds, runtime environments, and hash verification logs.

#### Phase-by-Phase Reports (`phases/`)
- [Phase 01: NCBI Data Acquisition](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_01_NCBI_DATA_ACQUISITION.md)
- [Phase 02: Protein Extraction](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_02_PROTEIN_EXTRACTION.md)
- [Phase 03: Quality Control](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_03_QUALITY_CONTROL.md)
- [Phase 04: Exact Deduplication](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_04_EXACT_DEDUPLICATION.md)
- [Phase 05: Annotation Template](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_05_ANNOTATION_TEMPLATE.md)
- [Phase 06: Temporal Annotation Finalization](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_06_TEMPORAL_ANNOTATION.md)
- [Phase 07: Physicochemical Feature Extraction](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_07_PHYSICOCHEMICAL_FEATURES.md)
- [Phase 08: ProtBERT Embedding Generation](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_08_PROTBERT_EMBEDDINGS.md)
- [Phase 09: Representation Matrix Construction](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_09_REPRESENTATION_MATRICES.md)
- [Phase 10: Exploratory Dimensionality Analysis](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/phases/PHASE_10_EXPLORATORY_ANALYSIS.md)

#### Decision Logs (`decisions/`)
- [Scientific Decisions](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/decisions/SCIENTIFIC_DECISIONS.md) — Foundational assumptions, strain choice, and comparative experimental principles.
- [Data Processing Decisions](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/decisions/DATA_PROCESSING_DECISIONS.md) — Handling repeats, sliding-window chunking, and pooling logic.
- [Modeling Decisions](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/decisions/MODELING_DECISIONS.md) — Unsupervised separation, future cross-validation safeguards, and class imbalance awareness.

#### Verification & Quality Control (`validation/`)
- [Validation Summary](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/validation/VALIDATION_SUMMARY.md) — Comprehensive master audit table across all phases.
- [Data Integrity](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/validation/DATA_INTEGRITY.md) — Bijective row alignments, NaN/Inf absence, and boundary checks.
- [Reproducibility Checks](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/validation/REPRODUCIBILITY_CHECKS.md) — Deterministic re-run verification logs.

#### Planning & Open Issues (`future/`)
- [Phase 11 Plan (Unsupervised Clustering)](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/future/PHASE_11_PLAN.md) — Planned K-Means ($K=2..10$), Hierarchical Clustering, stability metrics, and external validation protocol.
- [Open Methodological Issues](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/future/OPEN_ISSUES.md) — Unresolved decisions before Phase 11 execution (block weighting, PCA-before-clustering, stability protocols).

#### Audit Status
- [Documentation Audit](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/docs/DOCUMENTATION_AUDIT.md) — Verification of documentation completeness and integrity.

---

### Phase Status Overview

| Phase | Phase Title | Status | Primary Output Files |
| :--- | :--- | :---: | :--- |
| **Phase 1** | NCBI Data Acquisition | **COMPLETE** | `data/raw/ncbi/NC_001806.2.gbk`, `download_manifest.json` |
| **Phase 2** | Protein Extraction | **COMPLETE** | `data/intermediate/raw_proteins.csv`, `data/raw/fasta/raw_proteins.fasta` |
| **Phase 3** | Quality Control | **COMPLETE** | `data/intermediate/qc_passed_proteins.csv`, `qc_report.csv`, `qc_summary.json` |
| **Phase 4** | Exact Deduplication | **COMPLETE** | `data/processed/unique_proteins.csv`, `unique_proteins.fasta`, `duplicates.csv` |
| **Phase 5** | Annotation Template | **COMPLETE** | `data/annotations/temporal_annotations.csv` |
| **Phase 6** | Temporal Annotation | **COMPLETE** | `data/annotations/temporal_annotations_final.csv`, `annotation_audit.csv` |
| **Phase 7** | Physicochemical Features | **COMPLETE** | `data/processed/physicochemical_features.csv`, `feature_metadata.csv` |
| **Phase 8** | ProtBERT Embeddings | **COMPLETE** | `data/processed/protbert_embeddings.npy`, `protbert_metadata.csv` |
| **Phase 9** | Representation Matrices | **COMPLETE** | `X_physicochemical.npy`, `X_protbert.npy`, `X_combined_raw.npy`, `master_index.csv` |
| **Phase 10** | Exploratory Analysis | **COMPLETE** | PCA, UMAP, t-SNE tables & 24 publication figures in `results/figures/phase10/` |
| **Phase 11** | Unsupervised Clustering | **PLANNED** | *Not yet performed* |

---

### Current Status

- **Phase 1–10:** COMPLETED & VALIDATED
- **Phase 11:** PLANNED — NOT EXECUTED
- **Raw data:** FROZEN
- **Temporal annotations:** FROZEN
- **Phase 7 physicochemical features:** FROZEN
- **Phase 8 ProtBERT embeddings:** FROZEN
- **Phase 9 representation matrices:** FROZEN
- **Phase 10 exploratory outputs:** COMPLETED
- **Next scientific task:** Finalize Phase 11 clustering methodology before execution.
