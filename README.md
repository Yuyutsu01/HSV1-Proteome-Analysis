# Computational Analysis of the Herpes Simplex Virus Type 1 (HSV-1) Proteome
## Physicochemical Profiling, Dimensionality Reduction, Unsupervised Clustering, and Deep Embeddings

[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![NCBI RefSeq](https://img.shields.io/badge/NCBI_RefSeq-NC__001806.2-orange.svg)](https://www.ncbi.nlm.nih.gov/nuccore/NC_001806.2)
[![ProtBERT](https://img.shields.io/badge/Model-Rostlab%2Fprot__bert-purple.svg)](https://huggingface.co/Rostlab/prot_bert)

This repository contains the complete, reproducible computational biology codebase for the characterization and comparative representation benchmarking of the Herpes Simplex Virus Type 1 (HSV-1, strain 17) proteome.

---

## 🔬 Research Overview & Central Experimental Question

> **Central Question:** *"Do protein-language-model representations such as ProtBERT provide biologically meaningful information beyond conventional physicochemical descriptors for HSV-1 proteome characterization?"*

The pipeline systematically benchmarks:
1. **Classical Physicochemical Profiling (25-dimensional)**: Sequence length, molecular weight, theoretical isoelectric point ($pI$), Guruprasad instability index ($II$), aromaticity ($Ar$), and 20 canonical amino acid molar compositions ($f_A \dots f_Y$).
2. **Deep Transformer Embeddings (1024-dimensional)**: Residue-level representations pooled from `Rostlab/prot_bert` trained on 216M sequences (UniRef100).
3. **Combined Multimodal Fusion (1049-dimensional)**: Concatenation of standardized physicochemical and transformer embeddings.

Evaluation is conducted across **unsupervised clustering** ($K=2 \dots 10$), **manifold projections** (PCA, t-SNE, UMAP), **repeated stratified cross-validation classification** (Logistic Regression, SVM, Random Forest, XGBoost), **systematic feature ablation**, and **candidate outlier / disagreement profiling**.

---

## 📂 Project Architecture

```
HSV1-Proteome-Analysis/
├── data/
│   ├── NC_001806.2.gbk              # Raw NCBI RefSeq GenBank record
│   ├── raw_proteins.fasta           # All 77 raw extracted CDS proteins
│   └── curated_annotations.csv      # Curated temporal and functional categories
├── processed_data/
│   ├── non_redundant.fasta          # 74 deduplicated unique protein sequences
│   ├── dataset_statistics.csv       # Sequence lengths, duplicate removal, counts
│   └── physicochemical_features.csv # 25 physicochemical feature matrix
├── embeddings/
│   ├── protbert_embeddings.npy      # 74 x 1024 ProtBERT embedding matrix
│   ├── protbert_metadata.csv        # Protein IDs, lengths, and embedding metadata
│   └── combined_features.npy        # 74 x 1049 standardized multimodal representation
├── results/
│   ├── exploratory_statistics.csv   # Mean, median, std, IQR, skewness per feature
│   ├── feature_correlation_matrix.csv # Full Pearson correlation matrix
│   ├── pca_variance_explained.csv   # Explained & cumulative variance for all PCs
│   ├── pca_loadings.csv             # Feature loadings across principal components
│   ├── pca_results.csv              # Transformed PC coordinates per protein
│   ├── cluster_evaluation.csv       # K-Means metrics (WCSS, Silhouette, CH, DB)
│   ├── clustering_results.csv       # K-Means and Hierarchical cluster assignments
│   ├── biological_validation_metrics.csv # Temporal/Functional ARI, NMI, Purity
│   ├── protbert_clustering_evaluation.csv # ProtBERT clustering evaluation across K
│   ├── protbert_clustering_results.csv # ProtBERT cluster assignments
│   ├── representation_comparison.csv # Comprehensive head-to-head comparison table
│   ├── classification_results.csv   # Repeated 5-Fold CV metrics (Accuracy, F1)
│   ├── ablation_results.csv         # Feature ablation performance & deltas
│   └── candidate_outliers.csv       # Disagreement & low-confidence candidate proteins
├── figures/                         # 15 publication-grade figures (300 DPI)
├── scripts/                         # Modular Python execution scripts (01 - 12)
├── notebooks/
│   └── HSV1_Proteome_Complete_Pipeline.ipynb # Interactive notebook walkthrough
├── RESULTS_REPORT.md                # Comprehensive empirical results report
└── README.md                        # Documentation and reproduction guide
```

---

## ⚙️ Environment & Dependencies

- **Operating System**: Windows / Linux / macOS
- **Python Version**: Python 3.10+ (tested on Python 3.11)
- **Core Dependencies**:
  - `biopython >= 1.84`
  - `torch >= 2.0.0`
  - `transformers >= 4.35.0`
  - `scikit-learn >= 1.3.0`
  - `umap-learn >= 0.5.5`
  - `xgboost >= 2.0.0`
  - `matplotlib >= 3.8.0`
  - `seaborn >= 0.13.0`
  - `pandas >= 2.0.0`
  - `numpy >= 1.24.0`

### Installation
```bash
pip install biopython torch transformers scikit-learn umap-learn xgboost matplotlib seaborn pandas numpy
```

---

## 🚀 Reproduction Instructions

To execute the entire computational pipeline from scratch and reproduce all empirical tables, embeddings, and figures:

```bash
# 1. Fetch reference GenBank genome & deduplicate sequences
python scripts/01_fetch_and_deduplicate.py

# 2. Extract 25-dimensional physicochemical descriptors
python scripts/02_physicochemical_profiling.py

# 3. Exploratory statistics & correlation analysis
python scripts/03_exploratory_and_correlation.py

# 4. Dimensionality reduction via PCA
python scripts/04_dimensionality_reduction_pca.py

# 5. t-SNE and UMAP manifold visualization
python scripts/05_tsne_umap_visualization.py

# 6. Unsupervised clustering evaluation across K=2..10
python scripts/06_unsupervised_clustering.py

# 7. Post-hoc biological validation against curated ground truth
python scripts/07_biological_validation.py

# 8. Generate 1024-dimensional ProtBERT transformer embeddings
python scripts/08_protbert_embedding_generation.py

# 9. Evaluate ProtBERT manifold projections and clustering
python scripts/09_protbert_evaluation.py

# 10. Head-to-head representation comparison and multimodal fusion
python scripts/10_comparative_and_combined_analysis.py

# 11. Supervised machine learning classification benchmarks (Repeated 5-Fold CV)
python scripts/11_supervised_classification.py

# 12. Feature ablation study & candidate outlier profiling
python scripts/12_ablation_and_outliers.py

# 13. Render all 15 publication figures
python scripts/generate_all_figures.py
```

All random seeds are fixed to `42` (`random_state=42`) across all dimensionality reduction, clustering, and classification algorithms.

---

## 📄 License
This project is open source under the MIT License.
