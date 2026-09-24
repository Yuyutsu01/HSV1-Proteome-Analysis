# Phase 14 Methodology: Robustness and Sensitivity Analysis

## 1. Scientific Objective
Phase 14 systematically investigates whether the primary conclusions established in Phases 11–13 regarding HSV-1 proteome representations, unsupervised clustering, and supervised classification are robust to methodological variations.

The core scientific question addressed is:
> *"Are the principal conclusions of the unsupervised clustering, supervised classification, and representation analyses robust to reasonable changes in analytical choices?"*

The robustness analysis explicitly distinguishes:
- **A. Robust Findings:** Findings that remain stable across independent seeds, alternative weighting schemes, and dimensional transformations.
- **B. Sensitivity-Dependent Findings:** Findings whose magnitude or direction depends strongly on specific scaling or weighting conventions.
- **C. Inconclusive Findings:** Findings where small sample sizes (especially $N_{\text{IE}}=5$) limit strong population-level generalizability.

---

## 2. Frozen Primary Dataset & Input Artifacts
Phase 14 uses the identical frozen dataset of 74 HSV-1 strain 17 proteins (RefSeq `NC_001806.2`):
- **Biological Classes:** Immediate-Early (IE) = 5, Early = 15, Late = 54.
- **Authoritative Annotations:** `data/annotations/temporal_annotations_final.csv`
- **Representation Matrices:**
  - `data/processed/X_physicochemical.npy` ($74 \times 25$)
  - `data/processed/X_protbert.npy` ($74 \times 1024$)
  - `data/processed/X_combined_raw.npy` ($74 \times 1049$)
  - `data/processed/physicochemical_features.csv` (25 frozen features)
- **Zero Modification Rule:** No annotations, proteins, or upstream representations were altered. No synthetic data generation (SMOTE) was used.

---

## 3. Methodological Sensitivity Modules

### A. Supervised CV Seed Robustness
- **Original Seed Set:** `[42, 123, 456, 789, 2026]` (Phase 12 baseline)
- **Independent Sensitivity Seed Set:** `[7, 17, 37, 73, 97]` (Predefined, independent pseudo-random splits)
- **Cross-Validation Scheme:** 5-fold Stratified Cross-Validation repeated across 5 seeds (25 validation folds).
- **Evaluated Models:** 5 primary class-balanced configurations:
  1. Physicochemical + Logistic Regression
  2. ProtBERT + Logistic Regression
  3. Combined Equal-Block + Logistic Regression (Primary)
  4. Combined Equal-Block + Linear SVM
  5. Physicochemical + Random Forest
- **Findings:** Performance metrics were nearly identical between original and sensitivity seed sets. For Combined Equal-Block Logistic Regression, Balanced Accuracy was $0.758 \pm 0.020$ (Original) vs. $0.754 \pm 0.017$ (Independent), confirming seed-invariance.

### B. Class-Weighting Sensitivity
Evaluates the impact of loss-function weighting under three fold-contained regimes:
1. **Unweighted:** $w_c = 1.0$
2. **Inverse-Frequency Balanced:** $w_c = \frac{N_{\text{train}}}{K \cdot n_{c, \text{train}}}$ (Phase 12 default)
3. **Square-Root Balanced:** $w_c \propto \frac{1}{\sqrt{n_{c, \text{train}}}}$, normalized such that $\text{mean}(w_i) = 1.0$.
- **Findings:**
  - Unweighted learning causes minority class collapse: IE Recall drops to $0.00\text{--}0.04$ on ProtBERT and Combined models, maximizing majority Late accuracy at the expense of balanced classification.
  - Square-root weighting provides moderate minority penalty compensation (IE Recall = $0.12\text{--}0.28$).
  - Full inverse-frequency balancing restores minority class recovery (IE Recall = $0.68\text{--}0.80$), raising Balanced Accuracy on the primary Combined model from $0.411$ (Unweighted) to $0.758$ (Balanced).

### C. Combined Representation Strategy Sensitivity
Compares the two scaling strategies for multi-modal feature fusion:
1. **Combined Feature-Standardized:** Concatenates raw blocks ($25 + 1024 = 1049$) and fits a single standardizer in-fold.
2. **Combined Equal-Block Weighted:** Scales each standardized block by $\frac{1}{\sqrt{\text{dim}}}$ ($\frac{Z_{\text{phys}}}{\sqrt{25}} \,\|\, \frac{Z_{\text{pb}}}{\sqrt{1024}}$).
- **Findings:** Under Feature-Standardized scaling, the 1024 ProtBERT dimensions numerically swamp the 25 physicochemical features, reducing IE Recall to $0.040$ and Balanced Accuracy to $0.514$. Equal-block weighting prevents dimensional imbalance, achieving Balanced Accuracy of $0.758$ and Macro-F1 of $0.681$.

### D. In-Fold PCA Supervised Sensitivity
Evaluates dimensionality reduction inside each training fold (never fitted on the full dataset) at variance thresholds: Full (no PCA), 90%, 95%, and 99%.
- **Findings:**
  - For Combined Equal-Block Logistic Regression, retaining 90% variance (~22 PCs) or 95% variance (~30 PCs) yields Balanced Accuracy of $0.750\text{--}0.752$ and Macro-F1 of $0.660\text{--}0.663$, retaining $>98\%$ of full-dimensional performance while reducing dimensionality by $>97\%$.
  - For ProtBERT Logistic Regression, 90% variance (~10 PCs) and 95% variance (~17 PCs) slightly increase Balanced Accuracy ($0.675\text{--}0.691$ vs $0.647$) by filtering noisy embedding dimensions.

### E. Leave-One-Protein-Out Sensitivity (74 Systematic Trials)
- Evaluates individual protein leverage on the primary Combined Equal-Block Logistic Regression model.
- For each of the 74 proteins, the protein is temporarily held out ($N=73$), and repeated stratified cross-validation is performed.
- **Findings:**
  - Removal of individual Early or Late proteins induces minimal change in Balanced Accuracy ($|\Delta| < 0.015$).
  - Removal of individual IE proteins alters the minority sample count from 5 to 4, producing expected class-balance boundary shifts ($\Delta \text{ Balanced Accuracy} = -0.04 \text{ to } +0.02$).
  - No single protein acts as a catastrophic pivot point for overall model conclusions.

### F. Outlier Sensitivity
Evaluates model stability when excluding label-independent geometric outliers identified in Phase 13:
- **Cohorts Tested:**
  1. Full Primary Dataset ($N=74$)
  2. Excluding Prominent Outliers (`UL36`, `US12`, `UL41`; $N=71$)
  3. Excluding Top 10% Physicochemical Outliers ($N=66$)
  4. Excluding Top 10% ProtBERT Density Outliers ($N=66$)
  5. Excluding Union of Geometric Outliers ($N=60$)
- **Findings:**
  - Excluding the prominent outliers (`UL36`, `US12`, `UL41`) increases Balanced Accuracy from $0.758 \pm 0.020$ to $0.782 \pm 0.019$ and Macro-F1 from $0.681$ to $0.726$.
  - Excluding the union of geometric outliers ($N=60$) raises Balanced Accuracy to $0.804 \pm 0.018$.
  - These gains confirm that classification errors are concentrated among representation outliers, without altering the underlying temporal hierarchy.

---

## 4. Robustness Synthesis of Core Scientific Observations

| Observation ID | Scientific Hypothesis | Robustness Status | Empirical Evidence Summary |
| :--- | :--- | :---: | :--- |
| **Observation A** | Class-balanced learning substantially improves minority-class treatment relative to unweighted learning. | **ROBUST** | Unweighted learning collapses IE Recall to $0.00\text{--}0.04$, while inverse-frequency balancing achieves IE Recall $0.68\text{--}0.80$ across all tested seed sets. |
| **Observation B** | ProtBERT and physicochemical representations exhibit distinct, complementary predictive behavior. | **ROBUST** | 40.5% (30/74) of proteins show modality-dependent predictive disagreement. ProtBERT excels in Early recall, while Physicochemical excels in Late specificity. |
| **Observation C** | The combined representation captures complementary information from both modalities. | **ROBUST** | Combined Equal-Block representation achieves Macro-F1 = $0.66\text{--}0.68$ and Balanced Accuracy = $0.75\text{--}0.76$, matching or exceeding single modalities across all seed sets. |
| **Observation D** | Unsupervised clustering does not cleanly reproduce the three temporal classes. | **ROBUST** | Across all clustering algorithms (KMeans, Ward, Avg, Complete; $K=2..10$) in Phase 11, empirical permutation p-values exceeded 0.05 for temporal ARI/NMI. Geometry reflects biophysical structure rather than strict temporal classes. |
| **Observation E** | Protein-level prediction errors are concentrated among a specific subset of proteins and geometric outliers. | **ROBUST** | 73.0% (54/74) of proteins are consistently classified correctly ($\ge 80\%$), while errors are concentrated in 11 consistently misclassified proteins. Outlier exclusion systematically improves metrics without shifting class ordering. |

---

## 5. Summary of Created Artifacts

### Tabular Artifacts (`results/tables/`)
1. `phase14_cv_seed_sensitivity.csv` — Full metrics comparing Original vs Independent sensitivity seeds.
2. `phase14_class_weight_sensitivity.csv` — Comparison of Unweighted, Inverse-Frequency, and Square-Root weighting.
3. `phase14_combined_representation_sensitivity.csv` — Feature-Standardized vs Equal-Block combined scaling.
4. `phase14_pca_supervised_sensitivity.csv` — In-fold PCA evaluation at 90%, 95%, 99% variance.
5. `phase14_leave_one_protein_sensitivity.csv` — 74 systematic leave-one-protein ablation trials.
6. `phase14_outlier_sensitivity.csv` — Performance across 5 outlier exclusion cohorts.
7. `phase14_robustness_summary.csv` — Formal classification and evidence synthesis for Observations A–E.

### Publication Figures (`results/figures/phase14/`)
1. `cv_seed_sensitivity.png` — Balanced accuracy barplot across original and independent seeds.
2. `class_weight_sensitivity.png` — Macro-F1 comparison across weighting regimes.
3. `combined_representation_sensitivity.png` — Scaling strategy comparison across classifiers.
4. `pca_supervised_sensitivity.png` — Performance vs retained PCA variance.
5. `leave_one_protein_influence.png` — Waterfall plot of $\Delta \text{Balanced Accuracy}$ across all 74 proteins.
6. `overall_robustness_summary.png` — Outlier cohort sensitivity barplot.

### Logs
- `results/logs/phase14_validation_report.txt`
- `results/logs/phase14_robustness_summary.txt`
