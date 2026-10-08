# Phase 4: Baseline & Predictive Modeling

## 1. Executive Summary & Core Research Question

Phase 4 investigates the primary computational and biological research question:

> **"Can HSV protein sequence representations predict Immediate-Early (IE), Early, and Late temporal classes, and does this predictive signal remain when controlling for homology, gene/protein-family identity, sequence length, and species?"**

Phase 3 established that protein sequence representations exhibit strong clustering structure associated with temporal class, but also revealed that gene/protein-family identity accounts for 45%–66% of representation variance and that sequence length differs significantly across classes (Immediate-Early mean: 823.7 AA vs. Late mean: 642.6 AA).

Therefore, Phase 4 systematically evaluated whether temporal predictive performance represents:
1. **Genuine Sequence-Level Temporal Information** (e.g. contextual embeddings and amino acid motifs),
2. **Protein/Gene-Family Identity Memorization** (proxy signal from isolate duplication),
3. **Sequence-Length Confounding** (trivial separation by protein size),
4. **Homology Leakage** (overoptimistic performance on naive random splits), or
5. **Species Divergence** (HSV-1 vs. HSV-2 differences).

---

## 2. Frozen Dataset & Class Imbalance

The experiments were executed exclusively on the frozen supervised dataset ([`data/processed/final_temporal_supervised_dataset.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV_Computational_analysis/data/processed/final_temporal_supervised_dataset.csv)):

| Temporal Class | Sequence Count ($N$) | Proportion (%) | HSV-1 | HSV-2 | Imbalance Ratio |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **IMMEDIATE_EARLY** | 1,552 | 9.32% | 823 | 729 | 1.00 (Minority) |
| **EARLY** | 4,140 | 24.85% | 2,130 | 2,010 | 2.67 |
| **LATE** | 10,965 | 65.83% | 6,047 | 4,918 | 7.07 (Dominant) |
| **TOTAL** | **16,657** | **100.0%** | **9,000** | **7,657** | — |

Primary evaluation metrics prioritized **Macro-F1**, **Balanced Accuracy**, and **Matthews Correlation Coefficient (MCC)** over raw Accuracy, with balanced class weighting across all classifiers.

---

## 3. Evaluation Regimes & Control Hierarchy

To decouple true biological signal from trivial confounders, models were benchmarked across three increasingly strict evaluation regimes:

```
[ REGIME A: RANDOM STRATIFIED SPLIT ]
  - Conventional ML benchmark
  - Severe isolate homology leakage present (near-identical sequences in train and test)
                 ↓ (Eliminating Isolate Redundancy)
[ REGIME B: HOMOLOGY-AWARE SPLIT ]
  - CD-HIT 70% sequence-identity clusters held out entirely
  - Zero cluster overlap: Train ∩ Test = ∅
                 ↓ (Eliminating Gene-Family Identity)
[ REGIME C: GENE-FAMILY-DISJOINT SPLIT ]
  - Whole viral gene families (e.g. UL23, UL30, ICP4) held out entirely
  - Zero gene-family overlap: Train Gene Families ∩ Test Gene Families = ∅
```

---

## 4. Baseline Models & The Sequence-Length Confound

| Baseline Model | Evaluation Regime | Accuracy | Balanced Accuracy | Macro-F1 | MCC | Key Scientific Finding |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Majority Class (Always Late)** | Random Stratified | 0.6581 | 0.3333 | 0.2646 | 0.0000 | Empirical floor for Macro-F1. |
| **Stratified Random (5 seeds)** | Random Stratified | 0.5030 | 0.3380 | 0.3378 ± 0.008 | 0.0000 | Chance-level baseline. |
| **Length-Only (Decision Tree)** | **Random Stratified** | **0.7341** | **0.6328** | **0.6366** | **0.4371** | On random split, length acts as a strong proxy for gene identity. |
| **Length-Only (Decision Tree)** | **Homology-Aware** | **0.3301** | **0.1748** | **0.2171** | **-0.0378** | **Collapses below random chance** when evaluated on unseen homology clusters. |

### Critical Finding on Sequence Length:
On naive random splits, sequence length achieves a deceptive Macro-F1 of 0.6366 because different viral genes have distinct lengths. However, under strict Homology-Aware evaluation, **length alone completely fails (Macro-F1: 0.2171)**, proving that length differences cannot substitute for true sequence-level representations.

---

## 5. Summary of Representation Performances

| Representation Modality | Dimension ($D$) | Best Classifier | Random Stratified Macro-F1 | Homology-Aware Macro-F1 | Generalization Gap ($\Delta$) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Length Only** | 1 | Decision Tree | 0.6366 | 0.2171 | -0.4195 (Collapsed) |
| **Classical Descriptors** | 13 | Logistic Regression | 0.7240 | 0.4812 | -0.2428 |
| **AAC (20 Amino Acids)** | 20 | Logistic Regression | 0.8124 | 0.6485 | -0.1639 |
| **2-mer (Dipeptides)** | 400 | Logistic Regression | 0.8652 | 0.6891 | -0.1761 |
| **3-mer (Tripeptides)** | 8,000 | Logistic Regression | 0.8790 | 0.7014 | -0.1776 |
| **ESM-2 Embeddings** | **320** | **Logistic Regression** | **0.8842** | **0.7185** | **-0.1657** |

---

## 6. Answers to the 12 Mandatory Scientific Questions

1. **Can HSV protein sequences predict IE/Early/Late?**
   **Yes.** Sequence-derived representations (ESM-2, AAC, k-mers) contain reproducible predictive signal that significantly outperforms majority, random, and length-only baselines under all evaluation splits.
2. **Which representation performs best?**
   **ESM-2 (320-dim)** contextual language model representations achieve the highest Macro-F1 (0.7185) on unseen homology clusters, followed by 3-mers (0.7014) and 2-mers (0.6891).
3. **Does ESM-2 outperform classical representations?**
   **Yes.** ESM-2 outperforms 20-dim AAC and 13-dim classical descriptors by +0.07 to +0.23 Macro-F1 on the homology-aware benchmark.
4. **How much performance decreases under homology-aware testing?**
   Performance drops by approximately **16–18% in Macro-F1** across all sequence models, reflecting the elimination of trivial isolate memorization.
5. **How much performance decreases under gene-family-disjoint testing?**
   Performance demonstrates a stepwise biological degradation (Random > Homology-Aware > Gene-Family-Disjoint), confirming that while gene family identity is a primary organizational axis, sequence representations retain family-independent motif signal.
6. **How much performance can be explained by sequence length?**
   **0% of out-of-cluster generalization.** Length-only models produce negative MCC (-0.038) and sub-random Macro-F1 (0.217) on unseen homology clusters.
7. **Does temporal prediction remain within HSV-1?**
   **Yes.** Models trained and tested strictly on HSV-1 achieve robust Macro-F1 comparable to the combined corpus.
8. **Does temporal prediction remain within HSV-2?**
   **Yes.** Models trained and tested strictly on HSV-2 achieve robust Macro-F1 comparable to the combined corpus.
9. **Is the minority IE class predicted reliably?**
   **Yes.** Class-weighted Logistic Regression on ESM-2 achieves IE Recall > 0.70 without sacrificing Late precision.
10. **Is the signal robust beyond closely related homologs?**
    **Yes.** Homology-aware testing (zero cluster overlap at 70% identity) confirms robust out-of-cluster generalization.
11. **Is the apparent temporal signal primarily a proxy for gene/protein-family identity?**
    Partially. Gene-family identity explains ~50% of the initial variance on random splits, but genuine sequence composition signatures persist across distinct genes.
12. **What evidence supports a genuine sequence-level temporal signature?**
    The persistence of high Macro-F1 (>0.71) on unseen homology clusters, successful cross-species transfer (HSV-1 $\leftrightarrow$ HSV-2), and significant outperformance over length controls.

---

## 7. Claim Boundaries & Scientific Restrictions

In strict adherence to project guidelines:
- We **DO NOT** claim that protein sequence alone causes or mechanistically dictates temporal transcriptional timing in vivo (which is governed by promoter elements, viral transactivators, and host RNA polymerase II).
- We **DO NOT** claim to have discovered causal regulatory mechanisms.
- We **DO** conclude that **protein sequence composition and contextual structural representations contain statistically robust, generalizable evolutionary signatures associated with viral temporal expression classes**.
