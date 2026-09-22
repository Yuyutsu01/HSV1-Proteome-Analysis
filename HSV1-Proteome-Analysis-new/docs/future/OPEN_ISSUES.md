# Open Methodological Issues and Pre-Phase 11 Decision Register

This document records the specific open methodological questions, priority rankings, and action items that must be explicitly aligned and finalized before commencing Phase 11 execution.

---

### 1. Pre-Phase 11 Open Issues Register

| Issue | Priority | Required Action / Decision |
| :--- | :---: | :--- |
| **Combined block weighting** | 🔴 **Critical** | Define **equal-feature** vs **equal-block** weighting strategy. Feature standardization assigns equal variance to each feature ($25$ vs $1024$ features $\rightarrow 97.6\%$ ProtBERT variance dominance). Equal-block scaling scales each matrix block by its Frobenius norm or $\sqrt{D}$ to grant $50\%/50\%$ total variance contribution. |
| **Clustering input space** | 🔴 **Critical** | Execute primary clustering directly in defined **high-dimensional representation spaces** ($25$-D, $1024$-D, $1049$-D). Low-dimensional 2-D UMAP/t-SNE coordinates are strictly prohibited as primary clustering inputs. |
| **K-selection protocol** | 🔴 **Critical** | Predefine $K \in [2, 10]$ search range. Formally define mathematical selection logic based exclusively on internal metrics (max Silhouette, max Calinski-Harabasz, min Davies-Bouldin) without reference to biological ground truth. |
| **Class labels** | 🔴 **Critical** | Strict isolation: temporal class labels must be utilized **exclusively for post-hoc external validation** (ARI, NMI, Purity, contingency matrices). |
| **Class imbalance** | 🔴 **Critical** | Given heavy class imbalance (IE: 5, Early: 15, Late: 54), emphasize chance-corrected metrics (**Adjusted Rand Index, Normalized Mutual Information**) and complete contingency confusion matrices over unadjusted purity. |
| **K-Means seeds** | 🟠 **High** | Define reproducibility and stability protocol across multiple predefined random seeds (`n_init=100`, fixed seed list `[42, 123, 456, 789, 999]`) and consensus tracking. |
| **Hierarchical linkage** | 🟠 **High** | Predefine exact linkage/metric combinations: Ward (Euclidean), Complete (Cosine, Euclidean), Average (Cosine, Euclidean). |
| **PCA-before-clustering** | 🟠 **High** | Formally decide whether dimensionality-reduced representations (e.g., retaining $95\%$ PCA variance) serve as the primary input or as a sensitivity analysis against full-dimensional spaces. |
| **Cluster stability** | 🟠 **High** | Add repeated-seed stability and bootstrap/subsampling perturbation analysis to test cluster robustness. |
| **Phase 10 outlier wording** | 🟡 **Moderate** | Ensure all documentation removes unsupported biological claims regarding exploratory outliers (`RL1`, `US11`, `UL36`, `US5`, `UL23`), restricting descriptions to empirical sequence/coordinate observations. |
| **UMAP/t-SNE wording** | 🟡 **Moderate** | Use cautious, scientifically rigorous robustness language acknowledging hyperparameter sensitivity and non-metric distance distortions. |

---

### 2. Methodological Decision Details

#### Combined Block Weighting Rationale
When concatenating physicochemical features ($D_1 = 25$) and ProtBERT embeddings ($D_2 = 1024$), standard individual z-score normalization yields:
$$\text{Total Variance}_{\text{phys}} = 25, \quad \text{Total Variance}_{\text{protbert}} = 1024$$
ProtBERT features contribute $\frac{1024}{1049} \approx 97.6\%$ of the Euclidean distance calculation.
An alternative block-weighting strategy normalizes each block by $\frac{1}{\sqrt{D_k}}$:
$$\tilde{X}_{\text{combined}} = \left[ \frac{1}{\sqrt{25}} Z_{\text{phys}} \;\Big\|\; \frac{1}{\sqrt{1024}} Z_{\text{protbert}} \right]$$
granting equal $50\%/50\%$ aggregate variance contribution to each biological modality. This must be predefined before Phase 11 execution.

#### PCA-Before-Clustering Rationale
High-dimensional ProtBERT space ($1024$-D) on small sample size ($N=74$) exhibits distance concentration. Evaluating clustering on top principal components capturing $95\%$ variance ($19$ PCs for ProtBERT, $16$ PCs for Physicochemical, $37$ PCs for Combined) provides a well-conditioned sensitivity benchmark alongside raw high-dimensional clustering.
