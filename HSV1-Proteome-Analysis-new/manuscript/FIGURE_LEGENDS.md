# Publication Figure Legends

**Manuscript:** Computational Representation, Unsupervised Clustering, and Leakage-Controlled Classification of the Herpes Simplex Virus Type 1 Proteome  
**Dataset:** Human herpesvirus 1 strain 17 (RefSeq `NC_001806.2`, $N=74$ unique canonical proteins)  
**Date:** September 2026  

---

### Figure 1. Dataset Composition and Sequence Length Distribution of the HSV-1 Proteome
**(A)** Frequency distribution of the 74 unique canonical HSV-1 strain 17 proteins across the three experimentally established temporal expression classes: Immediate-Early ($\alpha$, $N=5$, $6.8\%$; orange), Early ($\beta$, $N=15$, $20.3\%$; purple), and Late ($\gamma$, $N=54$, $73.0\%$; green). Three identical terminal inverted repeat duplicate proteins (`RL2`, `RL1`, `RS1`) were removed during quality control.  
**(B)** Sequence length distributions (amino acids) stratified by temporal class, shown as Tukey boxplots overlaid with individual protein data points. Median sequence lengths: Immediate-Early = 412 aa (range: 88–1298 aa); Early = 510 aa (range: 197–1196 aa); Late = 445 aa (range: 96–3139 aa). Key benchmark extremes are highlighted: the shortest viral protein `US12` (ICP47, 88 aa, Immediate-Early) and the giant inner tegument hub `UL36` (VP1/2, 3,139 aa, Late).

---

### Figure 2. Multi-Modal Representation Spaces of the HSV-1 Proteome
Two-dimensional Principal Component Analysis (PCA) projections of the 74 canonical HSV-1 proteins across the three primary representation spaces, color-coded by temporal class (Immediate-Early: orange; Early: purple; Late: green).  
**(A)** Exact 25-feature physicochemical descriptor space ($74 \times 25$, StandardScaler normalized). PC1 accounts for $18.4\%$ and PC2 accounts for $13.1\%$ of total variance ($31.5\%$ cumulative). Separation along PC1 is driven by sequence length, mass, and aromatic composition.  
**(B)** Pretrained transformer ProtBERT embedding space ($74 \times 1024$, length-weighted chunk pooled). PC1 accounts for $13.8\%$ and PC2 accounts for $9.4\%$ of total variance ($23.2\%$ cumulative). Projections show separation of functional enzymatic clusters from structural tegument/envelope components.  
**(C)** Combined Equal-Block multi-modal space ($74 \times 1049$, block-scaled $[Z_{\text{phys}}/\sqrt{25} \;\|\; Z_{\text{pb}}/\sqrt{1024}]$). PC1 accounts for $14.6\%$ and PC2 accounts for $8.9\%$ of total variance ($23.5\%$ cumulative), integrating global biophysical constraints and contextual language embeddings with equal block weight.

---

### Figure 3. Unsupervised Clustering Benchmarks and Empirical Permutation Null Distributions
Evaluation of unsupervised clustering algorithms across cluster numbers $K=2$ to $10$.  
**(A)** K-Means internal clustering metrics (Inertia and Silhouette Coefficient) evaluated across 5 random seeds with 100 restarts per seed.  
**(B)** External validation concordance with annotated temporal classes measured by Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) across Physicochemical, ProtBERT, Combined Feature-Standardized, and Combined Equal-Block spaces.  
**(C)** Subsampling co-clustering stability matrix across 100 bootstrap iterations ($80\%$ subsampling rate) demonstrating consistent cluster cores corresponding to biophysical groupings (hydrophobic envelope vs. soluble enzymes) rather than transcriptional induction timing.  
**(D)** Empirical permutation testing against 1,000 label-shuffled null distributions with Davison-Hinkley finite-sample correction, illustrating that unsupervised representation geometry does not cleanly reproduce the three temporal classes.

---

### Figure 4. Leakage-Controlled Supervised Classification Benchmarks
Supervised temporal-class recovery under 5-fold Stratified Cross-Validation repeated across 5 predefined pseudo-random seeds (`[42, 123, 456, 789, 2026]`, 25 validation folds per configuration). All scalers, block-weights, and sample weights were fitted strictly in-fold.  
**(A)** Balanced Accuracy (macro-averaged recall) comparison across 4 representation spaces and 5 classifiers (Logistic Regression, Linear SVM, RBF SVM, Random Forest, XGBoost) under inverse-frequency class-balancing. The primary Combined Equal-Block Logistic Regression model achieved Balanced Accuracy of $75.86\% \pm 12.3\%$ (mean $\pm$ SD across folds) and Macro-F1 of $0.684 \pm 0.12$.  
**(B)** Out-of-fold confusion matrix heatmaps (total $N=370$ evaluated fold instances across 5 repeats) for the primary model, showing strong recovery of Immediate-Early ($80.0\%$, 20/25 instances), Early ($69.3\%$, 52/75 instances), and Late ($78.3\%$, 211/270 instances).  
**(C)** Minority-class Immediate-Early Recall as a function of class-weighting regime: Unweighted ($0.0\%$), Square-Root ($12.0\%$), and Inverse-Frequency ($80.0\%$).

---

### Figure 5. Representation Disagreement and Modality Concordance Taxonomy
Analysis of predictive concordance across representation spaces based on majority out-of-fold predictions ($N=74$ proteins).  
**(A)** Modality concordance matrix illustrating that $59.5\%$ (44/74) of proteins achieved unanimous agreement across Physicochemical, ProtBERT, and Combined representations, while $40.5\%$ (30/74) showed representation-dependent discordance. Zero proteins exhibited complete three-way divergence.  
**(B)** Breakdown of the 30 discordant proteins: Combined Equal-Block predictions followed ProtBERT in 19 cases ($25.7\%$) and Physicochemical in 11 cases ($14.9\%$).  
**(C)** Functional distribution of representation disagreements, illustrating that ProtBERT models favored multi-protein enzymatic replication cores while Physicochemical models provided distinct signals for proteins with extreme charge/composition.

---

### Figure 6. Protein-Level Error Profiling and Decision Margin Distributions
Out-of-fold prediction consistency and error confidence tiers for all 74 canonical proteins across 5 cross-validation repeats.  
**(A)** Per-protein prediction consistency ranked by cross-validation accuracy: $73.0\%$ (54/74) consistently correct ($\ge 80\%$ accuracy), $4.1\%$ (3/74) frequently correct, $8.1\%$ (6/74) frequently misclassified, and $14.9\%$ (11/74) consistently misclassified.  
**(B)** Out-of-fold mean decision margin ($\Delta p = P_{(1)} - P_{(2)}$) distribution, identifying 7 High-Confidence Incorrect predictions ($\Delta p > 0.35$: `UL15`, `RL1`, `US12`, `UL11`, `UL8`, `UL24`, `UL41`) and 4 Moderate-Margin Incorrect predictions ($0.15 \le \Delta p \le 0.35$: `UL36`, `UL49`, `UL52`, `UL13`).  
**(C)** Representation space biplot showing the distribution of error tiers relative to geometric boundaries.

---

### Figure 7. Methodological Robustness and Sensitivity Suite
Systematic sensitivity analyses evaluating the stability of the primary classification findings.  
**(A)** Cross-validation seed stability: comparison of primary performance between the original seed set (`[42, 123, 456, 789, 2026]`, Balanced Accuracy $= 0.758 \pm 0.020$) and an independent sensitivity seed set (`[7, 17, 37, 73, 97]`, Balanced Accuracy $= 0.754 \pm 0.017$).  
**(B)** In-fold PCA sensitivity: Balanced Accuracy across PCA variance thresholds ($90\%$, $95\%$, $99\%$, and full 1049 dimensions), demonstrating $>98\%$ performance retention at $95\%$ variance reduction (~30 PCs).  
**(C)** Leave-one-protein-out (LOO) influence distribution across 74 single-protein holdout trials ($|\Delta \text{ Balanced Accuracy}| < 0.015$ for all Early/Late proteins).  
**(D)** Geometric outlier sensitivity comparing full dataset ($N=74$, Balanced Accuracy $= 0.758$) and outlier-filtered dataset ($N=71$, Balanced Accuracy $= 0.782$).

---

### Figure 8. Biological Contextualization of Computationally Difficult Proteins
Literature-anchored contextual mapping of the 11 computationally difficult proteins within the HSV-1 replication cycle.  
**(A)** Functional categorization matrix mapping difficult proteins to multi-phase or virion-packaged lifecycle roles.  
**(B)** Case study summary for `UL36` (large tegument hub, 3,139 aa, Moderate-Margin Incorrect), `US12` (TAP inhibitor, 88 aa, High-Confidence Incorrect), `UL15` (terminase large subunit, High-Confidence Incorrect), `RL1` (ICP34.5 neurovirulence, High-Confidence Incorrect), and `UL41` (vhs RNase, High-Confidence Incorrect), distinguishing documented biological facts from computational observations and interpretive hypotheses.
