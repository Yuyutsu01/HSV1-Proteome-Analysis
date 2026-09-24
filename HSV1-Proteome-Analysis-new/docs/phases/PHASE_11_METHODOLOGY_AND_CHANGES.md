# Phase 11 Methodology and Changes

## 1. Phase 11 Objective
The primary scientific objective of Phase 11 is:
> *"Determine whether physicochemical, ProtBERT, and combined protein representations contain reproducible unsupervised structure that corresponds to independently assigned HSV-1 temporal classes (Immediate-Early=5, Early=15, Late=54)."*

This phase evaluates unsupervised clustering in strict isolation from temporal class labels to understand the geometric organization of the HSV-1 proteome across representation spaces.

---

## 2. Baseline Before Phase 11
Prior to Phase 11, the pipeline established:
- **Phase 9 Multi-Modal Matrices:**
  - $X_{\text{physicochemical}} \in \mathbb{R}^{74 \times 25}$
  - $X_{\text{protbert}} \in \mathbb{R}^{74 \times 1024}$
  - $X_{\text{combined\_raw}} \in \mathbb{R}^{74 \times 1049}$ (Raw unscaled concatenation, permanently frozen).
- **Phase 10 Exploratory Manifolds:**
  - Unsupervised Principal Component Analysis (PCA), Uniform Manifold Approximation and Projection (UMAP), and t-SNE projections with post-hoc coloring.
  - Identification of extreme observations (`RL1`, `US11`, `UL36`, `US5`, `UL23`).

---

## 3. Change 1 — Combined Representation Weighting

### The Dimensionality Disparity Problem
In previous exploratory concatenation, individual feature standardization was applied:
$$X_{\text{combined\_feat\_std}} = [Z_{\text{physicochemical}} \,\|\, Z_{\text{ProtBERT}}]$$
While each single feature has unit variance ($\sigma^2 = 1.0$), the aggregate variance contributed by the modalities is highly asymmetric:
$$\text{Var}_{\text{total}}(\text{Physicochemical}) = 25, \quad \text{Var}_{\text{total}}(\text{ProtBERT}) = 1024$$
Consequently, ProtBERT features contribute $\frac{1024}{1049} \approx 97.6\%$ to pairwise Euclidean distance calculations, dominating clustering geometry.

### Equal-Block Weighting Formulation
To evaluate multi-modal integration where both feature modalities exert balanced geometric influence, an equal-block-weighted representation was formulated:
1. Standardize each block independently:
   $$Z_{\text{phys}} = \text{StandardScaler}(X_{\text{physicochemical}})$$
   $$Z_{\text{pb}} = \text{StandardScaler}(X_{\text{protbert}})$$
2. Scale each block by the inverse square root of its dimensionality:
   $$Z_{\text{phys\_balanced}} = \frac{Z_{\text{phys}}}{\sqrt{25}}$$
   $$Z_{\text{pb\_balanced}} = \frac{Z_{\text{pb}}}{\sqrt{1024}}$$
3. Concatenate the balanced blocks:
   $$X_{\text{combined\_equal\_block}} = [Z_{\text{phys\_balanced}} \,\|\, Z_{\text{pb\_balanced}}]$$

This enforces that:
$$\sum_{j=1}^{25} \text{Var}(Z_{\text{phys\_balanced}, j}) = 1.0, \quad \sum_{k=1}^{1024} \text{Var}(Z_{\text{pb\_balanced}, k}) = 1.0$$
Both combined representations are evaluated in parallel as predefined sensitivity modalities.

---

## 4. Change 2 — Clustering Space Definition
- **Previous Potential Pitfall:** Clustering directly on 2-D UMAP or t-SNE projections introduces severe stochastic distortion, as non-linear projections warp inter-cluster distances and crowd data points.
- **Revised Method:** All primary clustering is executed directly in the **full high-dimensional representation spaces**:
  1. Physicochemical ($74 \times 25$ standardized)
  2. ProtBERT ($74 \times 1024$ native)
  3. Combined Feature-Standardized ($74 \times 1049$)
  4. Combined Equal-Block-Weighted ($74 \times 1049$)
- Low-dimensional PCA projections are reserved exclusively for post-hoc visualization.

---

## 5. Change 3 — Unbiased K Selection ($K \in [2, 10]$)
- **Previous Potential Pitfall:** Assuming $K=4$ or selecting $K$ based on maximum biological class overlap creates circular logic and inflates biological concordance claims.
- **Revised Method:** All integer values $K \in \{2, 3, 4, 5, 6, 7, 8, 9, 10\}$ are evaluated systematically.
- **Independence Rule:** Internal cluster validation metrics (Silhouette, Calinski-Harabasz, Davies-Bouldin) and stability metrics guide structural evaluation without reference to temporal annotations.

---

## 6. Change 4 — Multiple Clustering Algorithm Families
- **K-Means Clustering:** Centroid-based partitioning evaluated across 5 predefined random seeds with `n_init=100` restarts to prevent local minima entrapment.
- **Agglomerative Hierarchical Clustering:** Deterministic hierarchical tree construction evaluated across three mathematically rigorous linkage schemes:
  - Ward's minimum variance linkage (Euclidean metric)
  - Average linkage (Euclidean metric)
  - Complete linkage (Euclidean metric)

---

## 7. Change 5 — Multi-Seed and Subsampling Stability Analysis
Small sample size ($N=74$) makes unsupervised clustering sensitive to sample perturbation and seed initialization:
1. **Multi-Seed Stability:** Pairwise Adjusted Rand Index (ARI) computed across runs with seeds `[42, 123, 456, 789, 2026]`.
2. **Subsampling Co-Clustering:** 100 iterations of 80% subsampling (59/74 proteins sampled without replacement) tracking pairwise co-clustering frequencies to detect stable core clusters.

---

## 8. Change 6 — Post-Hoc External Biological Validation
- Temporal annotations (Immediate-Early: 5, Early: 15, Late: 54) are applied **strictly after cluster assignments are fixed**.
- Metrics computed:
  - **Adjusted Rand Index (ARI):** Measures chance-corrected cluster-class agreement.
  - **Normalized Mutual Information (NMI):** Quantifies shared information between cluster labels and temporal classes.
  - **Cluster Purity:** Proportion of dominant class members per cluster.
  - **Contingency Matrices:** Full $3 \times K$ cross-tabulations saved in `results/tables/contingency/`.

---

## 9. Change 7 — Empirical Permutation Testing
To establish whether observed ARI and NMI exceed what could arise by chance given the marginal class distributions:
- Cluster assignments are held fixed while temporal labels are randomly shuffled 1,000 times (`random_state=42`).
- Empirical null means, standard deviations, 95th percentiles, and empirical $p$-values are computed:
  $$p_{\text{empirical}} = \frac{1}{1000} \sum_{p=1}^{1000} \mathbb{I}(\text{Metric}_{\text{null}, p} \ge \text{Metric}_{\text{observed}})$$

---

## 10. Change 8 — Secondary PCA-Reduced Sensitivity Analysis
To test whether high-dimensional distance concentration affects clustering in ProtBERT ($D=1024$) and Combined spaces ($D=1049$):
- PCA is fitted independently on each representation.
- Principal components retaining $\ge 95\%$ cumulative variance are extracted:
  - **Physicochemical:** 16 PCs ($95.01\%$ cumulative variance)
  - **ProtBERT (Standardized):** 36 PCs ($95.11\%$ cumulative variance)
  - **Combined Equal-Block:** 33 PCs ($95.07\%$ cumulative variance)
- K-Means and Hierarchical clustering are re-evaluated and reported separately as a secondary sensitivity analysis.

---

## 11. Class Imbalance Considerations
The HSV-1 proteome is heavily skewed toward Late structural proteins:
- Immediate-Early: $5$ ($6.76\%$)
- Early: $15$ ($20.27\%$)
- Late: $54$ ($72.97\%$)

Because trivial majority-class assignment yields high purity ($\approx 73\%$), purity alone is never interpreted as evidence of biological partitioning. Primary emphasis is placed on chance-corrected metrics (**ARI, NMI**) and full contingency tables.

---

## 12. Strict Data Leakage Prevention
- Zero labels passed into `StandardScaler`, `PCA`, `KMeans`, `AgglomerativeClustering`, or subsampling algorithms.
- Complete algorithmic separation confirmed by automated unit tests in `tests/test_phase11_clustering.py`.

---

## 13. Files Modified
- *None.* No existing data files, upstream scripts, or authoritative Phase 1–10 outputs were modified.

---

## 14. Files Created
1. `scripts/11_unsupervised_clustering.py`
2. `tests/test_phase11_clustering.py`
3. `results/tables/phase11_internal_metrics.csv`
4. `results/tables/phase11_external_validation.csv`
5. `results/tables/phase11_kmeans_stability.csv`
6. `results/tables/phase11_stability_summary.csv`
7. `results/tables/phase11_coclustering_matrix_*.csv` (4 representation matrices)
8. `results/tables/phase11_permutation_tests.csv`
9. `results/tables/phase11_pca_clustering_metrics.csv`
10. `results/tables/phase11_representation_comparison.csv`
11. `results/tables/phase11_summary_for_manuscript.csv`
12. `results/tables/contingency/` (144 individual contingency CSVs)
13. `results/figures/phase11/` (PCA cluster biplots, dendrograms, co-clustering heatmaps, evaluation curves)
14. `results/logs/phase11_validation_report.txt`
15. `results/logs/phase11_interpretation_notes.txt`
16. `results/logs/phase11_environment.txt`
17. `docs/phases/PHASE_11_METHODOLOGY_AND_CHANGES.md`

---

## 15. Reproducibility Specifications
- **Python Version:** 3.11+
- **Primary Random Seed:** `42`
- **Multi-Seed Stability List:** `[42, 123, 456, 789, 2026]`
- **K-Means Initializations:** `n_init = 100`, `init = 'k-means++'`
- **Permutation Iterations:** `1,000`
- **Subsampling Iterations:** `100` (fraction: `0.80`)

---

## 16. Scientific Consequences
- Eliminates circularity and false claims of biological separation.
- Provides objective mathematical benchmarks comparing classical physicochemical properties against contextual transformer embeddings.
- Rigorously distinguishes between intrinsic geometric clusterability and biological class correspondence.

---

## 17. Limitations
- **Sample Size ($N=74$):** Small sample constraints increase sensitivity to boundary observations.
- **Single Reference Strain:** Limited to HSV-1 strain 17; strain-specific polymorphisms not evaluated.
- **Discrete Classes vs Kinetic Continuum:** Viral transcription is a continuous temporal cascade that is discretized into three nominal classes.

---

## 18. What Was NOT Changed
- Raw protein dataset `unique_proteins.csv` and `unique_proteins.fasta` remained untouched.
- Temporal annotations `temporal_annotations_final.csv` remained frozen.
- Phase 7 physicochemical features `physicochemical_features.csv` remained frozen.
- Phase 8 ProtBERT embeddings `protbert_embeddings.npy` remained frozen.
- Phase 9 raw multi-modal matrix `X_combined_raw.npy` remained frozen.

---

## 19. Final Phase 11 Audit Corrections

### A. PCA Dimension Consistency
The exact retained PCA dimensions capturing $\ge 95\%$ cumulative variance were derived from `results/tables/phase11_pca_clustering_metrics.csv`:
- **Physicochemical:** 16 PCs ($95.0112\%$ variance)
- **ProtBERT (Standardized):** 36 PCs ($95.1140\%$ variance)
- **Combined Equal-Block:** 33 PCs ($95.0713\%$ variance)

All documentation, validation logs, and summary notes reflect these exact numerical values.

### B. Singleton-Dominated Hierarchical Solutions
Audit of cluster size vectors revealed that Hierarchical Average Linkage consistently produces pathological singleton-peeling artifacts:
- Physicochemical $K=2$: `[73, 1]` (Silhouette = $0.347$)
- Physicochemical $K=3$: `[71, 1, 2]` (Silhouette = $0.323$)
- Physicochemical $K=4$: `[70, 1, 2, 1]` (Silhouette = $0.301$)
- Combined Equal-Block $K=2$: `[73, 1]` (Silhouette = $0.289$)

These solutions achieve deceptively high silhouette scores simply because isolated outliers (`UL36`, `RL1`, `US11`) are separated into singleton clusters while the remaining 70+ proteins form a single undifferentiated mass. In response, a formal `cluster_balance_warning` field (`SINGLETON_DOMINATED`, `SMALL_CLUSTER`, `NONE`) was added to all metric tables. Internal metrics must always be evaluated in tandem with cluster-size distributions.

### C. Permutation Test Formulation and Interpretation
1. **Finite-Sample Formulation:**
   In accordance with standard Monte Carlo permutation test theory (Davison & Hinkley 1997; North et al. 2002), empirical p-values are computed using the unbiased finite-sample estimator:
   $$p_{\text{empirical}} = \frac{\sum_{b=1}^B \mathbb{I}(T_b \ge T_{\text{obs}}) + 1}{B + 1}$$
   where $B = 1,000$. This ensures $p$-values are bounded in $(0, 1]$ (minimum achievable $p$-value $= \frac{1}{1001} \approx 0.000999001$). Zero $p$-values resulting from $0$ extreme permutations in the original $\frac{\text{count}}{B}$ estimator were updated to $\frac{1}{1001}$.
2. **Multiple Testing & Exploratory Framing:**
   Across the 144 evaluated clustering configurations, permutation $p$-values provide exploratory evidence of association relative to the label-permutation null; they are **not** interpreted as independent confirmatory hypothesis tests.
3. **Scientific Language Rules:**
   Claims such as *"confirmed biological validity"*, *"proved biological organization"*, or *"confirmed temporal regulation"* are strictly excluded. The permutation test evaluates whether observed agreement exceeds a randomized label null distribution; it does not demonstrate biological mechanisms.

### D. Limitations of Purity and Class Imbalance
With Late proteins representing $72.97\%$ ($54/74$) of the proteome:
- Assigning all samples to a single cluster trivially yields purity $= 0.7297$.
- Fragmenting data into $K=10$ clusters inflates purity ($>0.80$) without improving the genuine recovery of the three biological classes.
- Chance-corrected metrics (**ARI, NMI**) and full contingency matrices are required for meaningful biological assessment.

### E. Internal Geometry vs. Biological Concordance
No clustering configuration cleanly partitions the proteome into the three biological temporal classes. Intrinsic geometric clusterability is distinct from biological temporal regulation.

