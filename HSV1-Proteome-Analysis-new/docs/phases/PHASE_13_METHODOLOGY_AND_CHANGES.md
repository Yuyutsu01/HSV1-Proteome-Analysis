# Phase 13 Methodology: Supervised Model Error Analysis and Representation Interpretation

## 1. Scientific Objective
Phase 13 performs a systematic, individual-protein error and disagreement analysis on the supervised classification models evaluated in Phase 12.

The core scientific questions addressed are:
1. **Consistency:** Which HSV-1 proteins are consistently classified correctly versus consistently misclassified across repeated cross-validation folds?
2. **Transitions:** Which temporal-class transitions dominate classification confusion (e.g., IE $\rightarrow$ Early vs Early $\rightarrow$ Late)?
3. **Representation Differences:** Do physicochemical descriptors and ProtBERT embeddings produce distinct, representation-dependent predictive behaviors?
4. **Predictive Confidence:** Are model errors characterized by high-confidence mistakes or low-margin ambiguous predictions?
5. **Biochemical Profiles:** Are classification errors associated with unusual physicochemical properties (e.g., extreme length, charge, or hydrophobicity)?
6. **Cross-Phase Concordance:** How do supervised classification outcomes correspond to unsupervised clustering structures identified in Phase 11?

---

## 2. Frozen Upstream Inputs
Phase 13 strictly analyzes the existing outputs generated during Phase 12 and frozen representations:
- **Out-of-Fold Predictions:** `results/tables/phase12_out_of_fold_predictions.csv` ($14,800$ validation prediction records across 5 CV repeats).
- **Authoritative Annotations:** `data/annotations/temporal_annotations_final.csv` ($N=74$, IE=5, Early=15, Late=54).
- **Physicochemical Matrix:** `data/processed/physicochemical_features.csv` ($74 \times 25$).
- **Representation Matrices:** `X_physicochemical.npy`, `X_protbert.npy`, `X_combined_raw.npy`.
- **Phase 11 Clustering Outputs:** `results/tables/phase11_coclustering_matrix_combined_equal_block.csv`.

> [!IMPORTANT]
> **Zero Label Modification Rule:**
> Under no circumstances are temporal class annotations altered or reclassified based on model predictions. A computational misclassification is treated strictly as an informative observation about representation space geometry, not as biological evidence of misannotation.

---

## 3. Primary Models Analyzed

Five predefined class-balanced baseline configurations were selected for in-depth protein-level error profiling:
1. **Physicochemical + Logistic Regression** (Class-Balanced)
2. **ProtBERT + Logistic Regression** (Class-Balanced)
3. **Combined Equal-Block + Logistic Regression** (Class-Balanced Primary)
4. **Combined Equal-Block + SVM Linear** (Class-Balanced)
5. **Physicochemical + Random Forest** (Class-Balanced)

---

## 4. Analytical Modules and Methodological Definitions

### A. Protein Prediction Consistency
For each protein, out-of-fold validation predictions across all 5 cross-validation repeats ($5$ predictions per protein) were classified into four objective consistency tiers:
- **Consistently Correct:** Out-of-fold accuracy $\ge 80.0\%$ ($4/5$ or $5/5$ repeats correct).
- **Frequently Correct:** Out-of-fold accuracy $= 60.0\%$ ($3/5$ repeats correct).
- **Frequently Misclassified:** Out-of-fold accuracy $20.0\%\text{--}40.0\%$ ($1/5$ or $2/5$ repeats correct).
- **Consistently Misclassified:** Out-of-fold accuracy $= 0.0\%$ ($0/5$ repeats correct).

### B. Class-Confusion Transitions
Tracks the directional flow of predictions between biological classes ($3 \times 3 = 9$ transitions):
$$\text{Transition Fraction} = \frac{\text{Count}(\text{True Class } i \rightarrow \text{Predicted Class } j)}{\text{Total Count}(\text{True Class } i)}$$

### C. Representation Disagreement Matrix
Compares majority predictions across the three core modalities (Physicochemical, ProtBERT, Combined Equal-Block under Logistic Regression):
- `ALL_AGREE`: Unanimous prediction across all three modalities.
- `PHYS_PB_DISAGREE_COMBINED_MATCHES_PB`: Physicochemical differs from ProtBERT; Combined model aligns with ProtBERT.
- `PHYS_PB_DISAGREE_COMBINED_MATCHES_PHYS`: Physicochemical differs from ProtBERT; Combined model aligns with Physicochemical.
- `ALL_DISAGREE`: All three modalities predict different classes.

### D. Probability and Decision Margin Analysis
For probabilistic classifiers, the predictive margin is defined as:
$$\text{Margin} = P(\text{Top Class}) - P(\text{Second-Top Class})$$
Distinguishes high-confidence errors ($\text{Margin} \gg 0.30$) from low-margin, ambiguous predictions ($\text{Margin} \le 0.15$).

### E. Descriptive Physicochemical Error Analysis
For all 25 features, descriptive statistics (mean, SD, standardized mean difference Cohen's $d$) are computed comparing `Consistently_Correct` proteins against `Frequently/Consistently_Misclassified` proteins. No feature selection or causal claims are made.

### F. Cross-Phase Concordance (Phases 11 & 12)
Evaluates concordance between Phase 11 unsupervised clustering cohesion ($K=3$ subsampling co-clustering frequency) and Phase 12 supervised classification correctness:
- **Category A:** Correctly classified ($\ge 80\%$) and high clustering cohesion ($\ge 0.50$).
- **Category B:** Correctly classified ($\ge 80\%$) but cluster-ambiguous ($< 0.50$).
- **Category C:** Misclassified ($< 80\%$) but high clustering cohesion ($\ge 0.50$).
- **Category D:** Misclassified ($< 80\%$) and cluster-ambiguous ($< 0.50$).

---

## 5. Files Created & Verified

### Tabular Artifacts (`results/tables/`)
1. `phase13_protein_prediction_consistency.csv` — Protein-by-protein accuracy, class probabilities, margins, and consistency tiers across all primary models.
2. `phase13_class_confusion_analysis.csv` — Full $3 \times 3$ transition counts and class fractions.
3. `phase13_representation_disagreement.csv` — Modality concordance mappings across all 74 proteins.
4. `phase13_physicochemical_error_summary.csv` — 25-feature descriptive comparison (means, SDs, Cohen's $d$).
5. `phase13_error_profile.csv` — Model-level error summaries, class-specific error counts, and lists of consistently misclassified genes.
6. `phase13_outlier_analysis.csv` — Standardized geometric and k-NN embedding distances cross-referenced with classification accuracy.
7. `phase13_cross_phase_comparison.csv` — Cross-phase classification vs clustering concordance.

### Publication Figures (`results/figures/phase13/`)
1. `prediction_consistency_by_protein.png` — Per-protein out-of-fold accuracy barplot.
2. `representation_disagreement_matrix.png` — Modality agreement/disagreement frequency distribution.
3. `confusion_pattern_visualization.png` — Transition frequency profiles across primary models.
4. `pca_physicochemical_correctness.png` — Physicochemical PCA biplot dual-colored by true class and correctness.
5. `pca_protbert_correctness.png` — ProtBERT PCA biplot dual-colored by true class and correctness.
6. `pca_combined_correctness.png` — Combined Equal-Block PCA biplot colored by prediction consistency.
7. `predictive_margin_distribution.png` — Predictive probability margin boxplots across consistency tiers.

### Reports & Logs
- `results/logs/phase13_validation_report.txt`
- `docs/phases/PHASE_13_METHODOLOGY_AND_CHANGES.md`

---

## 6. Scientific Boundaries & Prohibited Interpretations
- **No Biological Relabeling:** Classification disagreements do not prove that an annotation is biologically incorrect.
- **No Causal Mechanism Claims:** Physicochemical differences between correct and incorrect predictions reflect computational classifier geometry, not molecular mechanisms of gene regulation.
- **Exploratory Status:** Disagreement proteins are designated as candidate targets for future experimental follow-up.

---

## 7. Final Phase 13 Audit Corrections

During the final audit of Phase 13, the following methodological checks and clarifications were performed and formalized:

1. **High-Confidence Incorrect Threshold Enforcement ($\Delta p > 0.35$):**
   - The threshold for "High-Confidence Incorrect" classification requires out-of-fold probability margin $\Delta p > 0.35$.
   - Protein `UL36` (Large tegument protein) has an out-of-fold margin of $\Delta p = 0.333 \le 0.35$. In the audit, `UL36` is strictly verified and categorized under *Moderate-Margin Incorrect* ($0.15 \le \Delta p \le 0.35$) rather than high-confidence error.
   - Proteins legitimately meeting $\Delta p > 0.35$ among misclassified proteins for the primary Combined Logistic Regression model are: `UL15` (0.849), `RL1` (0.771), `US12` (0.669), `UL11` (0.608), `UL8` (0.424), `UL24` (0.361), and `UL41` (0.352).

2. **Out-of-Fold Prediction Counting Terminology:**
   - In `phase13_class_confusion_analysis.csv`, the column `total_true_class_predictions` was renamed to `total_oof_prediction_instances` to avoid confusing out-of-fold prediction counts across 5 CV repeats ($25$ for IE, $75$ for Early, $270$ for Late) with unique biological protein counts ($5, 15, 54$).
   - Repeated out-of-fold validation predictions are documented as validation instances, not independent biological observations.

3. **Frozen 25-Feature Physicochemical Representation Verification:**
   - Verified that `phase13_physicochemical_error_summary.csv` exclusively evaluates the 25 frozen Phase 7 features (`sequence_length`, `molecular_weight`, `aromaticity`, `instability_index`, `isoelectric_point`, and 20 amino acid fractions `aa_A`..`aa_Y`).
   - Confirmed that `Gravy` is not part of the frozen Phase 7 feature matrix and is not present in Phase 13 tables.

4. **Descriptive Effect Sizes:**
   - Standardized mean differences (Cohen's $d$) are explicitly defined as descriptive summary metrics comparing computational representation geometry, not causal determinants or statistical significance tests.

5. **Cross-Phase & Outlier Interpretations:**
   - Unsupervised clustering structure (Phase 11) and supervised classification accuracy (Phase 12) are maintained as conceptually distinct analytical views.
   - Outlier proteins (`UL36`, `US12`, `UL41`) in representation space are treated descriptively without label alteration or sample pruning.

