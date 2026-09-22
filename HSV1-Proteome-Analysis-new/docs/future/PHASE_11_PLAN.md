# Phase 11: Unsupervised Clustering and Biological Concordance Plan

> [!IMPORTANT]
> **PLANNED PHASE — NOT YET EXECUTED**
> This document details the experimental design, algorithms, metrics, and validation protocols for Phase 11. No clustering or validation runs have been executed.

---

### 1. Scientific Objective
Systematically evaluate the intrinsic clusterability of the HSV-1 proteome across the three representation spaces ($X_{\text{physicochemical}}$, $X_{\text{protbert}}$, and $X_{\text{combined}}$) using multiple clustering algorithms, and assess whether unsupervised structural clusters correlate with known biological temporal expression classes.

---

### 2. Experimental Design and Scope

#### Evaluated Representation Spaces
1. **Physicochemical Matrix:** $X_{\text{physicochemical}} \in \mathbb{R}^{74 \times 25}$ (standardized features).
2. **ProtBERT Embedding Matrix:** $X_{\text{protbert}} \in \mathbb{R}^{74 \times 1024}$ (standardized / cosine normalized).
3. **Combined Multi-Modal Matrix:** $X_{\text{combined}} \in \mathbb{R}^{74 \times 1049}$ (block-standardized / modality-weighted).

#### Clustering Algorithms
- **K-Means Clustering:**
  - Evaluated across $K \in \{2, 3, 4, 5, 6, 7, 8, 9, 10\}$.
  - Multiple initializations (`n_init=100`, `k-means++`, predefined seeds).
- **Agglomerative Hierarchical Clustering:**
  - Evaluated across predefined linkage/metric combinations:
    1. Ward's linkage with Euclidean distance.
    2. Complete linkage with Cosine distance / Euclidean distance.
    3. Average linkage with Cosine distance / Euclidean distance.

---

### 3. Evaluation Metrics and Validation Protocol

```
Clustering Evaluation Architecture
 ├── Internal Cluster Quality (Unsupervised)
 │    ├── Silhouette Coefficient (SC)
 │    ├── Calinski-Harabasz Index (CH)
 │    └── Davies-Bouldin Index (DB)
 │
 ├── Cluster Stability & Robustness
 │    ├── Multi-Seed Jaccard / Adjusted Rand Stability
 │    └── Subsampling / Bootstrap Consensus Stability
 │
 └── External Biological Concordance (Post-Hoc Ground Truth Comparison)
      ├── Adjusted Rand Index (ARI)
      ├── Normalized Mutual Information (NMI)
      ├── Cluster Purity
      └── Full Contingency Confusion Matrices (IE / Early / Late)
```

---

### 4. Critical Methodological Safeguards

1. **Clustering Input Space:**
   - Primary clustering must be executed directly in the **full-dimensional representation spaces** ($25$-D, $1024$-D, $1049$-D), **never on 2D UMAP or t-SNE projection coordinates**.
2. **Internal K-Selection Rule:**
   - The optimal number of clusters $K^*$ must be selected based **exclusively on internal metrics (Silhouette, CH, DB) and stability analysis**.
   - **Never select $K$ using external ground-truth metrics (ARI / NMI)** to avoid circular optimization.
3. **Class Labels for External Evaluation Only:**
   - Temporal class labels are withheld completely during clustering. They are applied solely for post-hoc biological concordance reporting.
4. **Class Imbalance Awareness:**
   - Because Late proteins constitute 72.97% ($54/74$) of the proteome, unadjusted accuracy is misleading. Evaluation must emphasize Adjusted Rand Index (ARI), Normalized Mutual Information (NMI), and explicit contingency matrices.
