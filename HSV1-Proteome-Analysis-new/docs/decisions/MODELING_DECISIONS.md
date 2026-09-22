# Modeling Decision Record

## Machine Learning, Clustering, and Validation Protocol

This document outlines the modeling protocols, architectural constraints, and the strict distinction between completed exploratory phases and planned analytical phases.

---

### 1. Status of Computational Tasks

| Computational Task | Status in Current Pipeline | Notes / Implementation Strategy |
| :--- | :---: | :--- |
| **Exploratory Dimensionality Reduction (PCA, UMAP, t-SNE)** | **COMPLETED** | Fully unsupervised; class labels used strictly for post-hoc visualization. |
| **Unsupervised Clustering (K-Means, Hierarchical)** | **PLANNED** | Phase 11; planned over $K \in [2, 10]$ across representation spaces. |
| **Cluster Stability & Internal Metrics** | **PLANNED** | Phase 11; Silhouette, Calinski-Harabasz, Davies-Bouldin, bootstrap stability. |
| **External Biological Validation** | **PLANNED** | Phase 11; Adjusted Rand Index (ARI), Normalized Mutual Info (NMI), Purity. |
| **Supervised Classification Benchmarking** | **PLANNED** | Later phase; repeated stratified nested CV, class imbalance handling. |
| **Feature Ablation Studies** | **PLANNED** | Later phase; systematically evaluate modality contributions. |

---

### 2. Core Methodological Constraints for Future Modeling

#### Rule 1: No Label Leakage into Unsupervised Clustering
- Cluster algorithms ($K$-Means, Agglomerative Hierarchical) must be fitted strictly on feature matrices without access to temporal class labels.
- Optimal $K$ selection must rely exclusively on internal geometric metrics (Silhouette, Calinski-Harabasz, Davies-Bouldin, consensus stability), **never** on ARI or NMI with ground truth labels.

#### Rule 2: In-Fold Scaling for Supervised Learning
- In all planned supervised classification benchmarks, feature scaling (e.g., `StandardScaler`, `MinMaxScaler`) and dimensionality reduction (e.g., PCA) must be fitted strictly on the training partition of each cross-validation fold.
- Pre-scaling the complete dataset before cross-validation is strictly forbidden.

#### Rule 3: Handling Class Imbalance
- The HSV-1 temporal distribution is heavily imbalanced:
  - Immediate-Early: $N=5$ (6.76%)
  - Early: $N=15$ (20.27%)
  - Late: $N=54$ (72.97%)
- Modeling protocols must prioritize balanced evaluation metrics:
  - Macro-averaged F1 score, Balanced Accuracy, Cohen's Kappa, and contingency confusion matrices.
  - Class-weighted loss functions or balanced sampling inside training splits.

#### Rule 4: Regularization and Small-Sample Safeguards ($N=74$)
- Classifiers must incorporate strict regularization (e.g., $L_2$ ridge penalties, low tree depths, linear kernels) to prevent overfitting in high-dimensional embedding spaces ($D=1024$ and $D=1049$ with $N=74$).
