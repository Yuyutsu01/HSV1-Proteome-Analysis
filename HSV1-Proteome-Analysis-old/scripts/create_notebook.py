"""
Generate comprehensive, self-contained Jupyter Notebook for the entire HSV-1 Proteome Analysis pipeline.

Author: Computational Biology Pipeline
"""

import json
import os

def create_pipeline_notebook(output_path="notebooks/HSV1_Proteome_Complete_Pipeline.ipynb"):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# Computational Analysis of the Herpes Simplex Virus Type 1 (HSV-1) Proteome\n",
                "## Physicochemical Profiling, Dimensionality Reduction, Unsupervised Clustering, and Deep Embeddings\n",
                "\n",
                "**Authors**: Bioinformatics Research Team  \n",
                "**Reference Accession**: NCBI RefSeq `NC_001806.2` (HSV-1 Strain 17)  \n",
                "\n",
                "### Central Research Question:\n",
                "> *\"Do protein-language-model representations such as ProtBERT provide biologically meaningful information beyond conventional physicochemical descriptors for HSV-1 proteome characterization?\"*\n",
                "\n",
                "---"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 1. Biological Background and Dataset Reconstruction\n",
                "HSV-1 has a 152 kbp double-stranded DNA genome encoding ~74 unique non-redundant proteins. Gene expression proceeds through an ordered temporal cascade:\n",
                "- **Immediate-Early ($\\alpha$)**: Regulatory transactivators and immune evasion factors (e.g. ICP0, ICP4, ICP22, ICP27, ICP47).\n",
                "- **Early ($\\beta$)**: DNA replication replisome machinery and nucleotide metabolism (e.g. UL30 Pol, UL42 processivity, UL5/8/52 helicase-primase, UL23 TK).\n",
                "- **Late ($\\gamma$)**: Structural capsid, tegument, and envelope glycoproteins ($gB, gC, gD, gH, gL$, etc.).\n",
                "\n",
                "Duplicate proteins in the terminal/internal repeats ($TR_L/IR_L, TR_S/IR_S$) are deduplicated to prevent sampling bias."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "import os, sys\n",
                "import pandas as pd\n",
                "import numpy as np\n",
                "import matplotlib.pyplot as plt\n",
                "import seaborn as sns\n",
                "\n",
                "# View Dataset Statistics\n",
                "stats_df = pd.read_csv('../processed_data/dataset_statistics.csv')\n",
                "print('HSV-1 Proteome Dataset Statistics:')\n",
                "display(stats_df)"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 2. Physicochemical Feature Descriptors (25-Dimensional Representation)\n",
                "We extract 25 numerical descriptors using BioPython ProtParam:\n",
                "- **Global Descriptors (5)**: Sequence Length ($L$), Molecular Weight ($MW$), Isoelectric Point ($pI$), Guruprasad Instability Index ($II$), Aromaticity ($Ar$).\n",
                "- **Amino Acid Composition (20)**: Molar fraction of all 20 canonical amino acids ($f_A \\dots f_Y$).\n",
                "\n",
                "*Dimensionality Note*: While earlier manuscript drafts referred to 26 features, mathematical enumeration of 5 scalar descriptors + 20 amino acid frequencies confirms a 25-dimensional feature space."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "features_df = pd.read_csv('../processed_data/physicochemical_features.csv')\n",
                "print(f'Physicochemical Feature Matrix: {features_df.shape[0]} proteins x {features_df.shape[1]} columns')\n",
                "display(features_df.head(10))"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 3. Exploratory Data Analysis & Feature Correlation\n",
                "Distributional parameters (mean, median, standard deviation, skewness) and Pearson correlations across physicochemical features."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "exploratory_stats = pd.read_csv('../results/exploratory_statistics.csv')\n",
                "display(exploratory_stats)\n",
                "\n",
                "# Display correlation heatmap\n",
                "plt.figure(figsize=(10, 8))\n",
                "corr_img = plt.imread('../figures/fig05_correlation_heatmap.png')\n",
                "plt.imshow(corr_img)\n",
                "plt.axis('off')\n",
                "plt.show()"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 4. Principal Component Analysis (PCA)\n",
                "Variance decomposition and component loadings on standardized physicochemical features."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "pca_var = pd.read_csv('../results/pca_variance_explained.csv')\n",
                "print('Top 10 Principal Components Variance:')\n",
                "display(pca_var.head(10))\n",
                "\n",
                "plt.figure(figsize=(14, 5.5))\n",
                "pca_img = plt.imread('../figures/fig06_pca_scree_and_biplot.png')\n",
                "plt.imshow(pca_img)\n",
                "plt.axis('off')\n",
                "plt.show()"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 5. Unsupervised Clustering & Quantitative K Evaluation\n",
                "Testing K-Means across $K=2 \\dots 10$ using WCSS (Inertia), Silhouette Score, Calinski-Harabasz, and Davies-Bouldin indices."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "cluster_eval = pd.read_csv('../results/cluster_evaluation.csv')\n",
                "display(cluster_eval)\n",
                "\n",
                "plt.figure(figsize=(12, 9))\n",
                "clust_img = plt.imread('../figures/fig08_kmeans_evaluation_curves.png')\n",
                "plt.imshow(clust_img)\n",
                "plt.axis('off')\n",
                "plt.show()"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 6. Biological Validation & Ground Truth Alignment\n",
                "Post-hoc evaluation comparing unsupervised clusters against independently curated viral annotations using Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), and Purity."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "val_df = pd.read_csv('../results/biological_validation_metrics.csv')\n",
                "display(val_df)"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 7. ProtBERT Deep Embeddings & Multimodal Comparison\n",
                "Direct benchmarking of 25-dim Physicochemical Descriptors vs 1024-dim ProtBERT Deep Language Model Embeddings vs 1049-dim Combined Representation."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "comp_df = pd.read_csv('../results/representation_comparison.csv')\n",
                "display(comp_df)\n",
                "\n",
                "plt.figure(figsize=(14, 5.5))\n",
                "comp_img = plt.imread('../figures/fig11_protbert_vs_physicochem_clustering.png')\n",
                "plt.imshow(comp_img)\n",
                "plt.axis('off')\n",
                "plt.show()"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 8. Supervised Classification & Ablation Studies\n",
                "Repeated Stratified 5-Fold Cross Validation comparing Logistic Regression, SVM, Random Forest, and XGBoost across representations."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "clf_results = pd.read_csv('../results/classification_results.csv')\n",
                "display(clf_results.sort_values(by='Macro_F1_Mean', ascending=False))\n",
                "\n",
                "abl_results = pd.read_csv('../results/ablation_results.csv')\n",
                "display(abl_results)\n",
                "\n",
                "plt.figure(figsize=(14, 5.5))\n",
                "abl_img = plt.imread('../figures/fig13_ablation_delta_analysis.png')\n",
                "plt.imshow(abl_img)\n",
                "plt.axis('off')\n",
                "plt.show()"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "### 9. Candidate Outlier & Disagreement Analysis\n",
                "Candidate proteins exhibiting boundary characteristics, classifier low confidence, or multimodal representation conflict."
            ]
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "outliers_df = pd.read_csv('../results/candidate_outliers.csv')\n",
                "print(f'Identified {len(outliers_df)} candidate outlier/disagreement proteins:')\n",
                "display(outliers_df[['gene', 'true_temporal_class', 'phys_pred_class', 'bert_pred_class', 'disagreement_flag', 'instability_index']].head(15))"
            ]
        }
    ]
    
    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.11"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 5
    }
    
    with open(output_path, 'w') as f:
        json.dump(nb, f, indent=2)
    print(f"Created notebook at {output_path}")

if __name__ == "__main__":
    create_pipeline_notebook()
