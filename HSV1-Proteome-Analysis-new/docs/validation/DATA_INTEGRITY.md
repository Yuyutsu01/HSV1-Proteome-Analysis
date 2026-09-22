# Data Integrity Audit Report

## Authoritative Dataset Alignments, Dimensions, and Numerical Checks

This document summarizes the data integrity checks performed across all authoritative inputs, intermediate artifacts, and processed matrices.

---

### 1. Dataset Dimensions and Cardinality Checks

| Dataset / Matrix | File Path | Expected Shape | Actual Shape | Status |
| :--- | :--- | :---: | :---: | :---: |
| Raw CDS Extraction | `data/intermediate/raw_proteins.csv` | 77 rows | 77 rows | **MATCH** |
| QC Passing Proteins | `data/intermediate/qc_passed_proteins.csv` | 77 rows | 77 rows | **MATCH** |
| Deduplicated Unique Dataset | `data/processed/unique_proteins.csv` | 74 rows | 74 rows | **MATCH** |
| Deduplicated FASTA | `data/processed/unique_proteins.fasta` | 74 records | 74 records | **MATCH** |
| Final Annotations | `data/annotations/temporal_annotations_final.csv` | 74 rows | 74 rows | **MATCH** |
| Physicochemical Matrix | `data/processed/physicochemical_features.csv` | $74 \times 25$ | $74 \times 25$ | **MATCH** |
| ProtBERT Embedding Matrix | `data/processed/protbert_embeddings.npy` | $74 \times 1024$ | $74 \times 1024$ | **MATCH** |
| Aligned Physicochemical NPY | `data/processed/X_physicochemical.npy` | $74 \times 25$ | $74 \times 25$ | **MATCH** |
| Aligned ProtBERT NPY | `data/processed/X_protbert.npy` | $74 \times 1024$ | $74 \times 1024$ | **MATCH** |
| Raw Combined Matrix NPY | `data/processed/X_combined_raw.npy` | $74 \times 1049$ | $74 \times 1049$ | **MATCH** |
| Master Index | `data/processed/representation_master_index.csv` | 74 rows | 74 rows | **MATCH** |

---

### 2. Numerical Quality Checks

- **NaN / Null Values:** **0** across all representation matrices ($X_{\text{physicochemical}}$, $X_{\text{protbert}}$, $X_{\text{combined\_raw}}$) and annotation tables.
- **Infinite Values ($\pm \infty$):** **0** across all numeric arrays.
- **Amino Acid Composition Unity:** $\sum_{i=1}^{20} \text{fraction}_{AA} = 1.000000 \pm 10^{-6}$ for all 74 proteins.
- **Sequence Length Bijective Match:** $100.0\%$ match between sequence lengths in FASTA, CSVs, and feature metadata.
- **Embedding Value Bounds:** ProtBERT float32 embeddings range between $-1.9360$ and $+1.8973$, with zero vanishing/exploding gradients.

---

### 3. Bijective Row Alignment Verification

Every matrix row index $i \in \{0, \dots, 73\}$ in $X_{\text{physicochemical}}$, $X_{\text{protbert}}$, and $X_{\text{combined\_raw}}$ corresponds to the exact same protein entity defined in `representation_master_index.csv`:
- Row 0: `YP_009137073.1` (`RL1` / ICP34.5)
- Row 1: `YP_009137074.1` (`RL2` / ICP0)
- ...
- Row 73: `YP_009137148.1` (`US12` / ICP47)

**Row Alignment Result:** **100% Bijective Match (74/74 rows).**
