# Phase 10: Exploratory Representation Analysis

## 1. Objective
Conduct an unsupervised exploratory dimensionality analysis across the three representation spaces ($X_{\text{physicochemical}}$, $X_{\text{protbert}}$, and block-standardized $X_{\text{combined}}$) using Principal Component Analysis (PCA), Uniform Manifold Approximation and Projection (UMAP), and t-Distributed Stochastic Neighbor Embedding (t-SNE).

## 2. Inputs
- `data/processed/X_physicochemical.npy` ($74 \times 25$)
- `data/processed/X_protbert.npy` ($74 \times 1024$)
- `data/processed/X_combined_raw.npy` ($74 \times 1049$)
- `data/processed/representation_master_index.csv` (74 proteins with temporal labels)

## 3. Processing
- Executed `scripts/10_exploratory_representation_analysis.py`.
- **Dimensionality Reduction Methods:**
  - **PCA:** Full SVD deterministic decomposition. Standardized features for physicochemical and combined blocks.
  - **UMAP:** Primary $n\_neighbors = 10$, min_dist = 0.1, metric = 'euclidean' / 'cosine', `random_state=42`. Sensitivity sweeps at $n\_neighbors \in \{5, 20\}$.
  - **t-SNE:** Primary $\text{perplexity} = 15$, metric = 'euclidean' / 'cosine', `random_state=42`. Sensitivity sweeps at $\text{perplexity} \in \{5, 25\}$.
- Identified geometric outliers in projection spaces and recorded coordinate profiles.
- Generated 24 publication figures (scree, cumulative variance, 2D/3D PCA biplots, UMAP/t-SNE projections, and sensitivity grids) in `results/figures/phase10/`.
- Executed `tests/test_phase10_exploration.py` (6 automated test cases).

## 4. Outputs
- `results/tables/pca_explained_variance.csv` (Variance ratios per PC)
- `results/tables/pca_variance_thresholds.csv` (PCs required for 50%, 75%, 90%, 95% variance)
- `results/tables/pca_scores.csv` (Principal component scores per protein)
- `results/tables/umap_coordinates.csv` (UMAP coordinates across parameter sets)
- `results/tables/umap_parameters.csv`
- `results/tables/tsne_coordinates.csv` (t-SNE coordinates across parameter sets)
- `results/tables/tsne_parameters.csv`
- `results/tables/phase10_extreme_observations.csv` (Outlier observations)
- `results/tables/phase10_representation_comparison.csv`
- 24 high-resolution figures in `results/figures/phase10/`
- `results/logs/phase10_validation_report.txt`

## 5. Validation and Quantitative Findings

### Intrinsic Dimensionality (PCA Variance Thresholds)

| Representation Space | Input Dimensions | PC1 Variance | PC1+PC2 Cumulative | 50% Variance | 75% Variance | 90% Variance | 95% Variance |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Physicochemical** | 25 | **19.33%** | **32.56%** | 4 PCs | 9 PCs | 14 PCs | 16 PCs |
| **ProtBERT** | 1024 | **63.49%** | **70.79%** | 1 PC | 4 PCs | 11 PCs | 19 PCs |
| **Combined (Block-Std)** | 1049 | **20.34%** | **30.77%** | 5 PCs | 13 PCs | 26 PCs | 37 PCs |

### Key Exploratory Observations
- **Spectral Concentration:** ProtBERT exhibits strong spectral concentration along its leading principal axis (PC1 captures 63.49% variance), reflecting a dominant sequence length and compositional alignment axis across protein language models.
- **Multivariate Spread:** Physicochemical descriptors exhibit a more dispersed variance spectrum, requiring 16 components to capture 95% total variance.
- **Outlier Observations Identified:**
  - `RL1` (ICP34.5): Marked compositional skew (high Ala/Pro/Arg fractions).
  - `US11`: Marked basicity (Arg-rich RNA-binding domain, pI = 11.80).
  - `UL36`: Extreme sequence length (3,139 aa) and molecular weight (333.4 kDa).
  - `US5` (gJ): Small hydrophobic envelope glycoprotein.
  - `UL23` (Thymidine kinase): Distinct positioning in non-linear manifold embeddings.
- **Validation Suite:** `tests/test_phase10_exploration.py` (**6/6 tests passed**).

## 6. Scientific Decisions & Interpretation Rules
- **Strict Zero-Leakage Protocol:** Temporal labels were completely excluded during PCA, UMAP, and t-SNE fitting. Class labels were used strictly post-hoc for neutral visualization.
- **Cautious Interpretation Rule:** Visual clustering in 2D UMAP/t-SNE projections does not constitute statistical evidence of biological separability. Hypotheses generated in this phase must be formally tested via clustering metrics and cross-validated supervised classifiers.
- **Comparative Neutrality:** No representation is characterized as "superior" or "better" based on exploratory visuals alone.

## 7. Important Limitations
- Non-linear projections (UMAP and t-SNE) are sensitive to hyperparameters (perplexity, n_neighbors, distance metric) and random initialization; distance between disconnected clusters does not represent true metric space distance.

## 8. Files Generated
- [`results/tables/pca_explained_variance.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/tables/pca_explained_variance.csv)
- [`results/tables/pca_variance_thresholds.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/tables/pca_variance_thresholds.csv)
- [`results/tables/phase10_representation_comparison.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/tables/phase10_representation_comparison.csv)
- [`tests/test_phase10_exploration.py`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/tests/test_phase10_exploration.py)

## 9. Status
**COMPLETE & VALIDATED**
