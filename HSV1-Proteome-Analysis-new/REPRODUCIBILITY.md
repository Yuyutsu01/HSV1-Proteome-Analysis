# Computational Reproducibility Guide

## 1. Computational Environment
- **Operating System:** Windows 11 / Linux (Cross-platform compatible)
- **Python Version:** Python 3.11.8
- **Core Dependencies:**
  - `numpy >= 1.26.0`
  - `pandas >= 2.1.0`
  - `scikit-learn >= 1.4.0`
  - `xgboost >= 2.0.0`
  - `biopython >= 1.83`
  - `torch >= 2.2.0`
  - `transformers >= 4.38.0`
  - `matplotlib >= 3.8.0`
  - `seaborn >= 0.13.0`
  - `pytest >= 8.0.0`

## 2. Reference Genome & Dataset Accession
- **Reference Genome:** Human herpesvirus 1 strain 17
- **NCBI RefSeq Accession:** `NC_001806.2`
- **Authoritative Protein Count:** N = 74 unique canonical proteins (77 raw records minus 3 exact duplicates)
- **Deduplicated Identical Terminal Repeats:** `RL2` (YP_009137133.1), `RL1` (YP_009137134.1), `RS1` (YP_009137149.1)

## 3. Predefined Pseudo-Random Seeds
- **Primary Cross-Validation Seeds (Phase 12):** `[42, 123, 456, 789, 2026]`
- **Independent Sensitivity Seeds (Phase 14):** `[7, 17, 37, 73, 97]`
- **Permutation Null Iterations (Phase 11):** N = 1,000 label permutations with Davison-Hinkley `(r+1)/(N+1)` finite-sample correction.

## 4. Pipeline Execution Order
1. `scripts/01_ncbi_download.py` — Download RefSeq NC_001806.2
2. `scripts/02_extract_proteins.py` — Extract protein FASTA and metadata
3. `scripts/03_quality_control.py` — Sequence integrity and amino acid alphabet QC
4. `scripts/04_deduplication.py` — Exact sequence deduplication (77 -> 74 proteins)
5. `scripts/05_temporal_annotations.py` — Literature annotation curation (IE=5, Early=15, Late=54)
6. `scripts/06_annotation_audit.py` — Independent audit of temporal literature evidence
7. `scripts/07_feature_extraction.py` — Compute exact 25 physicochemical features
8. `scripts/08_protbert_embeddings.py` — Generate 1024-dim ProtBERT embeddings (Rostlab/prot_bert)
9. `scripts/09_representation_matrices.py` — Build X_physicochemical, X_protbert, X_combined_raw
10. `scripts/10_exploratory_analysis.py` — Exploratory PCA, UMAP, t-SNE dimensionality analysis
11. `scripts/11_unsupervised_clustering.py` — Unsupervised clustering evaluation across K=2..10
12. `scripts/12_supervised_classification.py` — Leakage-controlled repeated Stratified CV classification
13. `scripts/13_error_analysis.py` — Protein-level error, margin, and disagreement interpretation
14. `scripts/14_robustness_analysis.py` — Systematic sensitivity suite (seeds, weighting, PCA, LOO, outliers)
15. `scripts/15_biological_context.py` — Literature-anchored biological contextualization
16. `scripts/final_synthesis.py` — Master synthesis, manuscript, tables, and claim audit

## 5. Automated Verification
Run the full unit test suite:
```bash
pytest tests/ -v
```
