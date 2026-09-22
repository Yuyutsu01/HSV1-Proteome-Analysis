# HSV-1 Proteome Analysis (Revised & Enhanced Pipeline)

This directory contains the modernized, modular, and fully reproducible computational biology pipeline for the characterization and comparative representation benchmarking of the Herpes Simplex Virus Type 1 (HSV-1) proteome.

---

## 📁 Repository Directory Structure

```
HSV1-Proteome-Analysis-new/
│
├── data/
│   ├── raw/
│   │   ├── genomes/        # Complete FASTA genome sequences (e.g., NC_001806.2)
│   │   ├── genbank/        # Full GenBank flat files (.gbk) with CDS annotations
│   │   └── proteins/       # Raw translated CDS FASTA files
│   │
│   ├── curated/
│   │   ├── proteins/       # Deduplicated non-redundant protein sequences
│   │   ├── annotations/    # Curated functional, structural, and gene annotations
│   │   └── labels/         # Ground-truth temporal (α/β/γ) and functional target labels
│   │
│   ├── homology/           # CD-HIT / BLAST sequence clustering and redundancy filtering
│   └── splits/             # Train/test and cross-validation split indices
│
├── features/
│   ├── physicochemical/    # 25-dimensional classical physicochemical descriptors
│   └── protbert/           # 1024-dimensional ProtBERT transformer embeddings (.npy)
│
├── results/                # Quantitative evaluation metrics, CSV tables, PCA loadings
├── figures/                # 300 DPI publication-quality visualizations
├── notebooks/              # Interactive Jupyter walkthroughs
├── scripts/                # Modular, reproducible Python execution scripts
└── README.md               # Project overview and reproduction guide
```

---

## 🔬 Pipeline Workflow Overview

1. **Data Acquisition & Curation** (`data/`):
   - Fetch reference `NC_001806.2` GenBank records.
   - Programmatic CDS extraction, exact deduplication, and homology-based filtering.
2. **Feature Representation Engineering** (`features/`):
   - 25-dimensional physicochemical feature extraction (Length, MW, pI, Instability, Aromaticity, 20 AA fractions).
   - ProtBERT (`Rostlab/prot_bert`) deep transformer embeddings with sliding-window pooling for large proteins.
3. **Exploratory & Dimensionality Analysis** (`results/`, `figures/`):
   - Linear variance decomposition via PCA.
   - Non-linear neighborhood projections via t-SNE and UMAP.
4. **Clustering & Biological Validation** (`results/`):
   - Unsupervised clustering ($K=2 \dots 10$) using K-Means and Hierarchical Agglomerative clustering.
   - External validation against curated temporal classes ($\alpha, \beta, \gamma$) via ARI, NMI, and cluster purity.
5. **Supervised Machine Learning & Ablation** (`results/`):
   - Repeated Stratified 5-Fold Cross-Validation with strict in-fold pipeline scaling.
   - 5-subset feature ablation to quantify the marginal utility of deep language models.
