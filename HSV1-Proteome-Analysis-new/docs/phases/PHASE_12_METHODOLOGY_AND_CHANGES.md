# Phase 12 Methodology: Leakage-Controlled Supervised Classification with Class-Balanced Learning

## 1. Scientific Objective
Phase 12 addresses the following computational biology question:
> *"How much information about the predefined HSV-1 temporal classes (Immediate-Early=5, Early=15, Late=54) can be recovered from the physicochemical, ProtBERT, and combined representations under leakage-controlled supervised learning?"*

This phase evaluates whether supervised decision boundaries trained on sequence-derived features can distinguish the viral transcriptional cascade while strictly controlling for severe class imbalance and data leakage.

---

## 2. Class Imbalance and Majority-Class Baseline

### The HSV-1 Temporal Class Distribution ($N=74$)
- **Immediate-Early ($\alpha$ / IE):** $N = 5$ ($6.76\%$)
- **Early ($\beta$ / E):** $N = 15$ ($20.27\%$)
- **Late ($\gamma$ / L):** $N = 54$ ($72.97\%$)

### Majority-Class Baseline
A naive classifier predicting the dominant `Late` class for all proteins achieves:
- **Raw Accuracy:** $\frac{54}{74} \approx 72.97\%$
- **Balanced Accuracy:** $\frac{0.0 + 0.0 + 1.0}{3} \approx 33.33\%$ (Chance level)
- **IE Recall:** $0.0\%$
- **Early Recall:** $0.0\%$
- **Macro-F1:** $\approx 0.281$

Therefore, standard classification accuracy is inherently misleading in this setting. **Balanced Accuracy** (the unweighted mean of class-specific recalls) and **Macro-F1** serve as the primary evaluation criteria.

---

## 3. Strict Leakage-Controlled Preprocessing Architecture

All learned feature transformations and statistics are computed **strictly inside each training fold**:

```
5-Fold Stratified Cross-Validation (Repeated across 5 seeds: 42, 123, 456, 789, 2026)
 │
 ├── Training Partition (N = 59: 4 IE, 12 Early, ~43 Late)
 │    ├── 1. Fit StandardScaler on Training Partition (Mean, Std)
 │    ├── 2. Compute Equal-Block Scaling Factors (1/sqrt(25), 1/sqrt(1024))
 │    ├── 3. Compute In-Fold Class Weights: w_c = N_train / (3 * n_c_train)
 │    └── 4. Train Supervised Classifier with Regularization
 │
 └── Validation Partition (N = 15: exactly 1 IE, 3 Early, ~11 Late)
      ├── 1. Apply Frozen Training Scaler & Block Weights
      └── 2. Generate Out-of-Fold Prediction & Probability Scores
```

### Evaluated Representation Spaces
1. **Physicochemical ($25$-D):** `StandardScaler()` fitted in-fold on the 25 sequence descriptors.
2. **ProtBERT ($1024$-D):** Native contextual transformer embeddings (evaluated without scaling as embeddings occupy unit-variance space).
3. **Combined Feature-Standardized ($1049$-D):** Concatenation of physicochemical and ProtBERT features, standardized in-fold with `StandardScaler()`.
4. **Combined Equal-Block-Weighted ($1049$-D):** In-fold independent standardization of each block, followed by scaling by $\frac{1}{\sqrt{25}}$ and $\frac{1}{\sqrt{1024}}$ to enforce equal $1.0$ aggregate variance contribution.

---

## 4. Class-Balanced Learning Protocol

### Inverse-Frequency Class Weighting
To give equal optimization priority to each biological temporal class without fabricating synthetic data:
$$w_c = \frac{N_{\text{train}}}{K \times n_{c, \text{train}}}$$
where $K = 3$. For a typical training fold ($N_{\text{train}} = 59$):
- $w_{\text{IE}} = \frac{59}{3 \times 4} \approx 4.9167$
- $w_{\text{Early}} = \frac{59}{3 \times 12} \approx 1.6389$
- $w_{\text{Late}} = \frac{59}{3 \times 43} \approx 0.4574$

### Experimental Conditions
- **Condition A (Unweighted Baseline):** Models fitted with `class_weight=None` (sample weight $= 1.0$) to establish the raw impact of class imbalance.
- **Condition B (Class-Balanced Primary Model):** Models fitted with in-fold inverse-frequency class weights (`class_weight='balanced'`).

---

## 5. Supervised Classifier Portfolio

Five regularized, robust baseline classification algorithms were evaluated across all configurations:
1. **Logistic Regression:** L2-regularized multinomial logistic regression ($C=1.0$, L-BFGS solver).
2. **Support Vector Machine (Linear):** Linear Support Vector Classifier ($C=1.0$, Platt scaling calibration).
3. **Support Vector Machine (RBF):** Radial Basis Function Support Vector Classifier ($C=1.0, \gamma=\text{'scale'}$) as a nonlinear sensitivity baseline.
4. **Random Forest:** Ensemble of 100 decision trees ($\text{max\_depth}=5$) with in-fold bootstrap reweighting.
5. **XGBoost:** Gradient boosted decision trees (50 estimators, $\text{max\_depth}=3$, $\eta=0.1$, subsample $= 0.8$, colsample $= 0.8$) fitted with training-fold sample weights.

---

## 6. Small-Sample Safeguards & Evaluation Metrics

Given $N=74$ and only $5$ Immediate-Early proteins:
- **No SMOTE / Synthetic Data:** No synthetic interpolation was used in primary modeling to avoid unrealistic manifold artifacts in high-dimensional embedding space ($D=1024$).
- **Uncertainty Reporting:** All metrics are reported as Mean $\pm$ Standard Deviation across all 25 cross-validation folds (5 folds $\times$ 5 random seeds).
- **Comprehensive Metric Suite:**
  - Accuracy, Balanced Accuracy, Macro Precision, Macro Recall, Macro F1, Weighted F1, Matthews Correlation Coefficient (MCC).
  - Class-specific Precision, Recall, and F1-score for Immediate-Early, Early, and Late.
  - Raw and row-normalized $3 \times 3$ confusion matrices in fixed class order `['Immediate-Early', 'Early', 'Late']`.

---

## 7. Files Created & Modified

### Files Created
1. `scripts/12_supervised_classification.py`
2. `tests/test_phase12_classification.py`
3. `results/tables/phase12_cv_results.csv` (Fold-level metric logs)
4. `results/tables/phase12_summary.csv` (Aggregated Mean $\pm$ SD summary table)
5. `results/tables/phase12_per_class_metrics.csv` (Class-specific precision/recall/F1 metrics)
6. `results/tables/phase12_confusion_matrices.csv` (Aggregated raw and row-normalized confusion matrices)
7. `results/tables/phase12_out_of_fold_predictions.csv` (Complete sample-level OOF predictions with class probabilities)
8. `results/figures/phase12/balanced_accuracy_comparison.png`
9. `results/figures/phase12/representation_balanced_accuracy_comparison.png`
10. `results/figures/phase12/minority_ie_f1_comparison.png`
11. `results/figures/phase12/confusion_matrix_heatmaps_logistic_regression.png`
12. `results/logs/phase12_validation_report.txt`
13. `docs/phases/PHASE_12_METHODOLOGY_AND_CHANGES.md`

### Files Modified
- *None.* All upstream datasets, annotations, embeddings, and Phase 1–11 results remain unmodified and frozen.

---

## 8. Scientific Limitations & Boundaries
- **Minority Sample Size ($N=5$):** Immediate-Early performance estimates exhibit wide variance across individual folds due to small sample size ($1$ test sample per fold).
- **Static Sequence Representations:** Supervised classification evaluates whether sequence composition and language model semantics correlate with temporal categories; it does not model dynamic viral transcription factors, epigenetic chromatin regulation, or host-pathway feedback loops.
- **Single Strain Scope:** Results are restricted to the reference genome HSV-1 strain 17 (`NC_001806.2`).
