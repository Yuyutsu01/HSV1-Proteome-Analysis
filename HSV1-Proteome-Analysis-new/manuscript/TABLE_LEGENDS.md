# Publication Table Legends and Summary

**Manuscript:** Computational Representation, Unsupervised Clustering, and Leakage-Controlled Classification of the Herpes Simplex Virus Type 1 Proteome  
**Dataset:** Human herpesvirus 1 strain 17 (RefSeq `NC_001806.2`, $N=74$ unique canonical proteins)  
**Date:** September 2026  

---

### Table 1. Dataset Composition and Temporal Class Distribution of the HSV-1 Proteome
Summary of the 74 unique canonical proteins of HSV-1 strain 17 extracted from NCBI RefSeq `NC_001806.2` following exact deduplication of terminal repeat copies (`RL2`, `RL1`, `RS1`). Columns list the temporal expression class, protein count, proteome percentage, representative genes, and primary biological functions.

| Temporal Class | Protein Count | Percentage (%) | Representative Genes | Primary Biological Function |
| :--- | :--- | :--- | :--- | :--- |
| **Immediate-Early ($\alpha$)** | 5 | 6.76% | `RL2` (ICP0), `RS1` (ICP4), `UL54` (ICP27), `US1` (ICP22), `US12` (ICP47) | Master transcriptional transactivation, post-transcriptional splicing, host immune evasion |
| **Early ($\beta$)** | 15 | 20.27% | `UL5`, `UL8`, `UL9`, `UL12`, `UL23` (TK), `UL29` (ICP8), `UL30` (Pol), `UL39` (ICP6), `UL42`, `UL52` | Viral DNA replication, nucleotide metabolism, proofreading |
| **Late ($\gamma$)** | 54 | 72.97% | `UL19` (VP5), `UL27` (gB), `UL22`/`UL1` (gH/gL), `US6` (gD), `UL36` (VP1/2), `UL48` (VP16) | Capsid structural shell, envelope glycoproteins, tegument matrix, DNA packaging |

---

### Table 2. Descriptive Statistics of the 25 Physicochemical Features
Summary statistics (Mean, Standard Deviation, Minimum, Maximum, and Units) for the exact 25 physicochemical features computed across all 74 proteins using `Bio.SeqUtils.ProtParam`.

- **Sequence Length:** Mean $= 522.9 \pm 438.8$ residues (range: 88 to 3,139 aa)
- **Molecular Weight:** Mean $= 56,842 \pm 48,154$ Da (range: 9,792 to 336,126 Da)
- **Aromaticity:** Mean $= 0.076 \pm 0.021$ (range: 0.027 to 0.136)
- **Instability Index:** Mean $= 48.7 \pm 10.9$ (range: 24.7 to 72.6; $81.1\%$ predicted unstable $>40$)
- **Isoelectric Point (pI):** Mean $= 7.82 \pm 1.84$ (range: 4.57 to 11.83)
- **Amino Acid Frequencies ($f_{\text{A}} \dots f_{\text{Y}}$):** Dominated by Alanine ($11.9\% \pm 2.9\%$), Proline ($9.3\% \pm 3.6\%$), and Leucine ($8.9\% \pm 2.2\%$), reflecting genomic GC richness.

---

### Table 3. Representation Spaces and Dimensionality Properties
Comparison of the four feature representation matrices evaluated in the study, including total feature dimensions, underlying data types, in-fold normalization methods, and number of principal components required to capture $95\%$ cumulative variance.

| Representation Name | Dimensions | Feature Type | In-Fold Normalization | PCs for 95% Variance |
| :--- | :--- | :--- | :--- | :--- |
| **Physicochemical** | 25 | Primary biophysical & compositional descriptors | StandardScaler ($Z_{\text{phys}}$) | 14 |
| **ProtBERT** | 1024 | Contextual transformer embeddings (`Rostlab/prot_bert`) | Native mean-pooled embeddings ($Z_{\text{pb}}$) | 18 |
| **Combined Feature-Standardized** | 1049 | Direct raw feature concatenation | Global StandardScaler | 28 |
| **Combined Equal-Block** | 1049 | Block-weighted multi-modal concatenation | Equal block scaling $[Z_{\text{phys}}/\sqrt{25} \;\|\; Z_{\text{pb}}/\sqrt{1024}]$ | 30 |

---

### Table 4. Unsupervised Clustering Summary across Representations and Algorithms
Summary of unsupervised clustering performance across representations and linkage methods for $K=2..10$. Outlines internal cluster validity (Silhouette, Inertia), stability across 100 bootstrap subsamples, and external temporal concordance (ARI, NMI) evaluated relative to 1,000 label permutations.

---

### Table 5. Supervised Classification Performance Benchmarks (Phase 12)
Cross-validation performance (Mean $\pm$ Standard Deviation across 25 folds) for 5 classifiers across the 4 representation spaces under leakage-controlled repeated Stratified Cross-Validation with inverse-frequency class-balancing.

- **Combined Equal-Block Logistic Regression (Primary Model):**
  - Accuracy: $76.51\% \pm 10.2\%$
  - Balanced Accuracy: $75.86\% \pm 12.3\%$
  - Macro-F1: $0.684 \pm 0.12$
  - Matthews Correlation Coefficient (MCC): $0.533 \pm 0.19$
  - Immediate-Early Recall: $80.0\%$
  - Early Recall: $69.3\%$
  - Late Recall: $78.3\%$
- **Linear SVM (Equal-Block):** Balanced Accuracy $= 73.23\% \pm 11.8\%$, Macro-F1 $= 0.671 \pm 0.11$.
- **RBF SVM (Equal-Block):** Balanced Accuracy $= 70.14\% \pm 13.0\%$, Macro-F1 $= 0.648 \pm 0.13$.

---

### Table 6. Per-Class Classification Metrics for Primary Models
Detailed breakdown of Precision, Recall, F1-score, and Support for Immediate-Early, Early, and Late temporal classes under the primary Combined Equal-Block Logistic Regression model.

| Class Label | Support (Instances) | Recall (%) | Precision (%) | F1-Score |
| :--- | :--- | :--- | :--- | :--- |
| **Immediate-Early ($\alpha$)** | 25 (5 folds $\times$ 5 repeats) | 80.0% | 43.5% | 0.561 |
| **Early ($\beta$)** | 75 (15 folds $\times$ 5 repeats) | 69.3% | 68.4% | 0.687 |
| **Late ($\gamma$)** | 270 (54 folds $\times$ 5 repeats) | 78.3% | 93.4% | 0.852 |
| **Macro Average** | 370 total instances | **75.86%** | **68.43%** | **0.684** |

---

### Table 7. Protein-Level Error and Decision Margin Taxonomy
Categorization of the 11 computationally difficult proteins from out-of-fold cross-validation predictions, reporting gene symbol, protein ID, sequence length, true class, frequent predicted class, mean decision margin ($\Delta p$), and assigned confidence tier.

- **High-Confidence Incorrect ($\Delta p > 0.35$):** `UL15` ($\Delta p = 0.849$), `RL1` ($\Delta p = 0.771$), `US12` ($\Delta p = 0.669$), `UL11` ($\Delta p = 0.608$), `UL8` ($\Delta p = 0.424$), `UL24` ($\Delta p = 0.361$), `UL41` ($\Delta p = 0.352$).
- **Moderate-Margin Incorrect ($0.15 \le \Delta p \le 0.35$):** `UL36` ($\Delta p = 0.333$), `UL49` ($\Delta p = 0.335$), `UL52` ($\Delta p = 0.274$), `UL13` ($\Delta p = 0.254$).

---

### Table 8. Methodological Robustness and Sensitivity Suite Summary
Summary of primary model performance across the 5 sensitivity dimensions evaluated in Phase 14:
1. **CV Seed Stability:** Original seeds ($0.758 \pm 0.020$) vs. Independent seeds ($0.754 \pm 0.017$) $\rightarrow$ Status: **ROBUST**.
2. **Class-Weight Sensitivity:** Unweighted (IE Recall $= 0.0\%$) vs. Inverse-Frequency (IE Recall $= 80.0\%$) $\rightarrow$ Status: **ROBUST**.
3. **Multi-Modal Scaling:** Equal-Block (Bal Acc $= 0.758$) vs. Raw Concat (Bal Acc $= 0.712$) $\rightarrow$ Status: **ROBUST**.
4. **In-Fold PCA Sensitivity:** PCA 95% variance (Bal Acc $= 0.749$) vs. Full 1049 dims (Bal Acc $= 0.758$) $\rightarrow$ Status: **ROBUST**.
5. **Outlier Sensitivity:** Full dataset (Bal Acc $= 0.758$) vs. Outlier-excluded (Bal Acc $= 0.782$) $\rightarrow$ Status: **ROBUST**.

---

### Table 9. Literature-Anchored Biological Contextualization of Difficult Proteins
Structured mapping linking the 11 computationally difficult proteins to documented biological functions, localization, virion packaging, and primary literature citations with PubMed IDs (PMIDs) and DOIs.

---

### Table 10. Evidence Hierarchy and Claim Traceability Matrix
Final synthesis matrix establishing evidence types (Documented Fact, Empirical Null Comparison, Leakage-Controlled CV, Descriptive Concordance, Literature Hypothesis), robustness status, biological interpretations, and manuscript safe wording for all major study findings.
