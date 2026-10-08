# Computational Representation, Unsupervised Clustering, and Leakage-Controlled Classification of the Herpes Simplex Virus Type 1 Proteome

**Authors:** Shiva  
**Affiliation:** Antigravity IDE Computational Virology Group  
**Target Journal:** *Bioinformatics* / *Journal of Virology*  
**Date:** September 2026  

---

## ABSTRACT

**Background:** Herpes simplex virus type 1 (HSV-1) is a double-stranded DNA virus whose lytic replication cascade is classically categorized into Immediate-Early ($\alpha$), Early ($\beta$), and Late ($\gamma$) kinetic classes. While transcriptional activation, chromatin accessibility, and promoter motifs govern this temporal program, the degree to which viral protein sequences alone encode predictable temporal and functional signatures remains an open question in computational virology.

**Objective:** To systematically evaluate how effectively HSV-1 proteins can be represented and classified according to experimentally defined temporal expression classes using 25 primary physicochemical descriptors and 1024-dimensional contextual protein language model embeddings (ProtBERT), and to assess what representation geometry reveals about viral proteome organization.

**Methods:** We analyzed the complete reference proteome of HSV-1 strain 17 (RefSeq `NC_001806.2`, $N=74$ unique canonical proteins after exact deduplication; Immediate-Early $N=5$, Early $N=15$, Late $N=54$). Three primary representation spaces were established: (i) an exact 25-feature physicochemical descriptor matrix, (ii) a 1024-dimensional ProtBERT embedding matrix generated with length-weighted chunk pooling (500-residue chunks, 100-residue overlap, 400-residue step), and (iii) an equal-block multi-modal representation ($[Z_{\text{phys}}/\sqrt{25} \;\|\; Z_{\text{pb}}/\sqrt{1024}]$, 1049 dimensions total). Unsupervised clustering across $K=2..10$ was evaluated against empirical label-permutation null distributions ($N=1,000$). Supervised classification was performed under 5-fold Stratified Cross-Validation repeated across 5 pseudo-random seeds with strictly in-fold preprocessing and inverse-frequency class-balancing. Protein-level prediction errors, representation disagreements, and sensitivity to analytical choices were systematically analyzed and contextualized using peer-reviewed biological literature.

**Results:** Across the evaluated clustering configurations, temporal-class recovery was limited and varied substantially across representations and linkage methods; no clustering configuration cleanly reproduced the three temporal classes. Permutation testing was used to assess the observed association with temporal labels relative to a label-permutation null distribution, serving as exploratory evidence rather than independent confirmatory tests. Unsupervised geometry reflected biophysical characteristics (such as multi-pass transmembrane glycoproteins forming distinct clusters) rather than transcriptional induction timing. In contrast, leakage-controlled supervised learning recovered substantial information about temporal annotations: the primary Combined Equal-Block Logistic Regression model achieved Balanced Accuracy of $75.86\% \pm 12.3\%$, Macro-F1 of $0.684 \pm 0.12$, and minority Immediate-Early Recall of $80.0\%$ (compared to a naive majority-class baseline of Accuracy = $72.97\%$, Balanced Accuracy = $33.33\%$, IE Recall = $0.0\%$). Inverse-frequency class weighting substantially improved minority-class recovery in the evaluated experiments (Immediate-Early Recall $= 0.000$ unweighted, $0.120$ square-root weighted, and $0.800$ inverse-frequency weighted). ProtBERT and physicochemical representations exhibited different predictive behavior across a subset of proteins ($40.5\%$, 30/74 discordance). Among the 30 proteins showing disagreement between the two modality-specific predictions, the combined prediction followed ProtBERT for 19 proteins ($25.7\%$) and the physicochemical representation for 11 proteins ($14.9\%$). Eleven computationally difficult proteins (including `UL15`, `RL1`, `US12`, `UL41`, and `UL36`) were contextualized using literature-documented lifecycle roles, providing plausible biological context without demonstrating causal mechanisms.

**Conclusions:** HSV-1 temporal expression classes are not simply equivalent to sequence-space similarity. Supervised models can recover substantial information about temporal annotations when class imbalance is explicitly addressed. Physicochemical and ProtBERT representations exhibit distinct predictive behavior, and the equal-block combined representation produced stable classification performance across the predefined evaluation and sensitivity analyses. Biological contextualization provides plausible interpretations for difficult cases, but these interpretations remain hypotheses requiring experimental validation.

---

## KEYWORDS
Herpes simplex virus type 1; Protein language models; ProtBERT; Physicochemical descriptors; Multi-modal representation; Class-balanced learning; Temporal expression kinetics; Leakage-controlled machine learning.

---

## 1. INTRODUCTION

Herpes simplex virus type 1 (HSV-1) is an archetypal human alphaherpesvirus characterized by a linear double-stranded DNA genome of approximately 152 kilobase pairs that establishes lifelong latent infection in sensory ganglia and undergoes episodic lytic reactivation (Roizman et al., 2013). During productive lytic infection, viral gene expression executes a tightly coordinated, tripartite cascade regulated at transcriptional, post-transcriptional, and structural assembly levels:
1. **Immediate-Early ($\alpha$):** Expressed immediately upon viral genome entry into the host nucleus independent of de novo protein synthesis. The five canonical $\alpha$ genes (`RL2`/ICP0, `RS1`/ICP4, `UL54`/ICP27, `US1`/ICP22, and `US12`/ICP47) function primarily as master transcriptional transactivators, splicing regulators, and host immune evasion factors (Honess & Roizman, 1974; Preston, 1979; Everett & Maul, 1994).
2. **Early ($\beta$):** Transcribed in response to $\alpha$ transactivation, encoding the catalytic viral DNA replication machinery (e.g., DNA polymerase `UL30`, processivity factor `UL42`, single-stranded DNA-binding protein `UL29`/ICP8, origin-binding protein `UL9`, and the heterotrimeric helicase-primase `UL5`/`UL8`/`UL52`) as well as enzymes for nucleotide salvage and metabolism (thymidine kinase `UL23`, ribonucleotide reductase `UL39`/`UL40`, uracil-DNA glycosylase `UL2`) (Crute et al., 1989; Weller et al., 1983; Gibbs et al., 1985).
3. **Late ($\gamma$):** Expressed maximally following the onset of viral DNA synthesis, subdivided into leaky-late ($\gamma_1$) and true-late ($\gamma_2$). Late genes encode the structural components of the infectious virion, including the icosahedral capsid shell (`UL19`/VP5, `UL18`/VP23, `UL38`/VP19C), the complex proteinaceous tegument matrix (`UL36`/VP1/2, `UL37`, `UL48`/VP16, `UL49`/VP22), and envelope glycoproteins (`UL27`/gB, `UL22`/gH, `UL1`/gL, `US6`/gD) that mediate viral attachment, entry, and egress (Dai & Zhou, 2018; Heldwein et al., 2006; Bigalke & Heldwein, 2015).

The temporal regulation of HSV-1 has historically been understood through the lens of transcriptional promoter architecture, host RNA polymerase II recruitment, and viral chromatin dynamics. However, the mature polypeptides encoded by these kinetic classes exhibit tremendous functional, biochemical, and structural heterogeneity. Proteins range from small regulatory peptides (such as the 88-amino-acid transporter-associated with antigen processing inhibitor `US12`/ICP47) to the 3,139-amino-acid large tegument hub `UL36` (VP1/2), the largest structural protein encoded by human herpesviruses (Sandbaumhüter et al., 2013).

Classical computational sequence analysis relies on fixed physicochemical descriptors—such as sequence length, molecular weight, theoretical isoelectric point, instability index, aromaticity, and amino acid composition—to summarize protein structure and physical constraints (Guruprasad et al., 1990; Kyte & Doolittle, 1982). While informative, such global descriptors compress high-dimensional sequence relationships into scalar aggregates, potentially discarding local motifs, domain architectures, and contextual sequence patterns.

Recent advances in self-supervised deep learning have enabled the training of large-scale protein language models (pLMs), such as ProtBERT (Elnaggar et al., 2021). Trained on hundreds of millions of diverse protein sequences across all kingdoms of life using masked language modeling objectives, ProtBERT learns high-dimensional contextual embeddings that encode evolutionary, structural, and biophysical properties without manual feature engineering.

Despite the rapid adoption of pLMs in structural biology and variant effect prediction, their utility in viral proteome classification and their comparative relationship with classical biophysical features remain underexplored. Furthermore, computational studies in viral bioinformatics are frequently compromised by methodological pitfalls, including:
- Data leakage across cross-validation splits due to pre-split normalization or feature selection;
- Misleading accuracy metrics on severely imbalanced datasets without minority-class balancing;
- Dimensionality imbalance during multi-modal feature concatenation, where high-dimensional language model embeddings (e.g., 1024 dimensions) numerically overpower compact physical feature sets (e.g., 25 dimensions);
- Overstated claims of biological causality inferred solely from statistical model predictions.

In this study, we present a systematic, leakage-controlled computational investigation of the HSV-1 strain 17 proteome. We evaluate exact 25-feature physicochemical descriptors, 1024-dimensional ProtBERT embeddings, and equal-block fused representations across unsupervised clustering, supervised classification, protein-level error profiling, sensitivity suites, and literature-anchored biological contextualization.

---

## 2. STUDY OBJECTIVES

The formal objectives of this study were predefined to ensure transparent, hypothesis-driven computational investigation:

1. **Objective 1 (Physicochemical Characterization):** Compute and audit an exact 25-feature physicochemical descriptor space across all canonical proteins of HSV-1 strain 17 without data leakage or missing values.
2. **Objective 2 (Contextual Sequence Representation):** Extract 1024-dimensional ProtBERT embeddings using length-weighted chunk pooling for long sequences and establish an equal-block multi-modal representation that prevents feature-count dominance.
3. **Objective 3 (Unsupervised Geometry Assessment):** Determine whether unsupervised clustering algorithms (K-Means, Ward, Average, Complete linkage across $K=2..10$) recapitulate experimentally defined temporal expression classes relative to empirical label-permutation null distributions.
4. **Objective 4 (Supervised Classification Benchmarking):** Evaluate the recovery of temporal annotations under strictly leakage-controlled 5-fold Stratified Cross-Validation repeated across 5 pseudo-random seeds, comparing unweighted and class-balanced loss formulations.
5. **Objective 5 (Protein-Level Error & Disagreement Profiling):** Quantify out-of-fold prediction margins, categorize decision confidence tiers, and map representation-dependent predictive discordance between physicochemical and transformer embeddings.
6. **Objective 6 (Methodological Robustness Suite):** Evaluate the stability of primary findings against analytical perturbations, including independent cross-validation seeds, in-fold PCA dimensionality reduction, systematic leave-one-protein-out trials, and geometric outlier exclusions.
7. **Objective 7 (Literature-Anchored Biological Contextualization):** Map computational error profiles and representation disagreements against peer-reviewed HSV-1 molecular virology literature to formulate biologically plausible interpretations while maintaining strict distinction between computational association and causal mechanism.

---

## 3. MATERIALS AND METHODS

### 3.1 Reference Genome and Protein Extraction
The complete annotated reference genome of human herpesvirus 1 strain 17 was retrieved from NCBI GenBank/RefSeq (Accession `NC_001806.2`, Genome Length: 152,222 base pairs). Protein-coding sequences (CDS) were extracted using Biopython (version 1.83) and cross-referenced against NCBI RefSeq feature tables. A total of 77 raw protein-coding records were initially extracted.

### 3.2 Sequence Quality Control and Deduplication
All extracted amino acid sequences were subjected to programmatic quality control:
1. Alphabet validation: verified that sequences contained only standard IUPAC 20 amino acid characters;
2. Translation completeness: verified valid methionine start codons and unambiguous terminal stop codons;
3. Exact sequence deduplication: identified and removed exact duplicate copies resulting from terminal inverted repeats. Three duplicate records were removed: `RL2` (YP_009137133.1, identical to YP_009137074.1), `RL1` (YP_009137134.1, identical to YP_009137073.1), and `RS1` (YP_009137149.1, identical to YP_009137148.1).
The final curated proteome comprises exactly **$N=74$ unique canonical proteins**.

### 3.3 Authoritative Temporal Annotations
Temporal expression kinetic classes were curated from canonical reference literature (Fields Virology, 6th Edition; Roizman et al., 2013; Honess & Roizman, 1974; Preston, 1979) and verified through primary research publications:
- **Immediate-Early ($\alpha$):** 5 proteins ($6.8\%$) — `RL2` (ICP0), `RS1` (ICP4), `UL54` (ICP27), `US1` (ICP22), and `US12` (ICP47).
- **Early ($\beta$):** 15 proteins ($20.3\%$) — `UL5`, `UL8`, `UL9`, `UL12`, `UL23`, `UL29`, `UL30`, `UL39`, `UL40`, `UL42`, `UL50`, `UL52`, `UL53`, `US3`, and `US8A`.
- **Late ($\gamma$):** 54 proteins ($73.0\%$) — structural capsid, tegument, envelope, and packaging proteins.

Envelope glycoprotein C (`UL44`/gC) is annotated with primary temporal class Late and late subclass *Conflicting* ($\gamma_1/\gamma_2$) based on literature demonstrating variable transcription dependence on DNA synthesis (Frink et al., 1983); it is preserved in all analyses without alteration.

### 3.4 Physicochemical Feature Computation
An exact 25-feature physicochemical descriptor vector was computed for each of the 74 proteins using `Bio.SeqUtils.ProtParam`:
1. **Sequence Length ($N_{\text{aa}}$):** Total number of amino acid residues;
2. **Molecular Weight ($M_{\text{w}}$):** Total mass in Daltons (Da);
3. **Aromaticity:** Relative frequency of aromatic residues ($\text{Phe} + \text{Trp} + \text{Tyr}$);
4. **Instability Index:** Guruprasad dipeptide-based in vitro stability estimate ($>40$ indicates predicted instability in vitro; Guruprasad et al., 1990);
5. **Isoelectric Point (pI):** Theoretical pH at which net charge equals zero;
6–25. **Amino Acid Composition ($f_{\text{A}} \dots f_{\text{Y}}$):** Fractional molar composition for each of the 20 standard amino acids.
No undocumented or speculative features (such as GRAVY) were included, yielding a complete $74 \times 25$ descriptor matrix ($X_{\text{phys}}$).

### 3.5 Contextual ProtBERT Embeddings and Long-Sequence Chunking
Contextual protein sequence embeddings were generated using the pretrained bidirectional transformer `Rostlab/prot_bert` (Elnaggar et al., 2021) implemented in PyTorch and Hugging Face Transformers.
- **Chunking Protocol:** Sequences exceeding the standard transformer context window were partitioned into 500-residue chunks with a 100-residue sliding overlap (stride = 400 residues).
- **Pooling Strategy:** For each chunk, residue-level hidden states from the final transformer layer (excluding special `[CLS]` and `[SEP]` tokens) were extracted and averaged (mean pooling).
- **Chunk Aggregation:** For proteins spanning multiple chunks, aggregate embeddings were computed via length-weighted averaging proportional to valid chunk residue lengths.
- **Longest Protein Handling:** Giant tegument protein `UL36` (3,139 amino acids) was processed across 8 chunks without truncation, producing a 1024-dimensional embedding vector.
The resulting matrix $X_{\text{pb}}$ has dimensions $74 \times 1024$.

### 3.6 Multi-Modal Representation Construction
Four distinct representation matrices were evaluated:
1. **Physicochemical ($X_{\text{phys}}$):** Dimensions $74 \times 25$.
2. **ProtBERT ($X_{\text{pb}}$):** Dimensions $74 \times 1024$.
3. **Combined Feature-Standardized ($X_{\text{comb}}$):** Direct concatenation with in-fold standard scaling ($74 \times 1049$).
4. **Combined Equal-Block ($X_{\text{eq}}$):** Block-weighted concatenation defined as:
   $$Z_{\text{eq}} = \left[ \frac{Z_{\text{phys}}}{\sqrt{25}} \;\Bigg\|\; \frac{Z_{\text{pb}}}{\sqrt{1024}} \right]$$
   where $Z_{\text{phys}}$ and $Z_{\text{pb}}$ represent standard normal transformations ($mean=0, SD=1$) fitted strictly on the training partition of each cross-validation fold. Equal-block scaling assigns equal total variance (1.0) to the 25-dimensional physicochemical block and the 1024-dimensional ProtBERT block, preventing the larger feature set from numerically overwhelming physical descriptors.

### 3.7 Unsupervised Clustering and Permutation Testing Protocol
Unsupervised structure was assessed across $K=2..10$ clusters using:
- **K-Means:** Evaluated across 5 random seeds (`[42, 123, 456, 789, 2026]`) with 100 random restarts (`n_init=100`) per seed;
- **Hierarchical Agglomerative Clustering:** Evaluated under Ward's minimum variance, Complete, and Average linkage criteria;
- **Subsampling Stability:** Evaluated across 100 bootstrap iterations at an $80\%$ subsampling rate to generate pairwise co-clustering stability matrices;
- **Empirical Permutation Testing:** For every clustering configuration, temporal class labels were randomly permuted $N=1,000$ times. Empirical $p$-values for Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) were calculated with the Davison-Hinkley finite-sample correction:
  $$p_{\text{emp}} = \frac{\sum_{i=1}^N \mathbb{I}(T_i \ge T_{\text{obs}}) + 1}{N + 1}$$

### 3.8 Supervised Classification and Leakage Control
Supervised learning was conducted under 5-fold Stratified Cross-Validation repeated across 5 pseudo-random seeds (`[42, 123, 456, 789, 2026]`, total 25 validation folds per configuration).
- **Leakage Isolation:** All StandardScaler parameters, PCA transformations, and sample weights were fitted exclusively on the training fold and applied to the test fold.
- **Classifiers Evaluated:**
  1. Logistic Regression ($L_2$ regularization, $C=1.0$, L-BFGS solver);
  2. Linear Support Vector Machine ($C=1.0$);
  3. Radial Basis Function (RBF) Support Vector Machine ($C=1.0, \gamma=\text{'scale'}$);
  4. Random Forest ($n_{\text{estimators}}=100$, Gini criterion);
  5. Gradient Boosted Trees (XGBoost, learning rate $=0.1$, max depth $=3$).
- **Class Weighting:** In-fold inverse-frequency class weights were computed as:
  $$w_c = \frac{N_{\text{train}}}{K \cdot n_{c, \text{train}}}$$
  where $K=3$, $N_{\text{train}}$ is training sample size, and $n_{c, \text{train}}$ is class count. Unweighted and square-root weighted models were evaluated for comparison.

### 3.9 Error Analysis and Modality Disagreement Taxonomy
Out-of-fold predictions and class posterior probability distributions ($P(\text{IE}), P(\text{Early}), P(\text{Late})$) were aggregated across all 5 cross-validation repeats (14,800 total fold predictions):
- **Decision Margin ($\Delta p$):** Defined as the difference between the probability assigned to the predicted class and the probability of the second most likely class:
  $$\Delta p = P_{(1)} - P_{(2)}$$
- **Confidence Tiers:**
  - High-Confidence Incorrect: $\Delta p > 0.35$ with frequent misclassification;
  - Moderate-Margin Incorrect: $0.15 \le \Delta p \le 0.35$;
  - Low-Margin / Boundary: $\Delta p < 0.15$.
- **Disagreement Taxonomy:** Out-of-fold majority predictions for each protein were compared across Physicochemical, ProtBERT, and Combined representations to classify instances into unanimous agreement, modality-discordant (combined follows ProtBERT vs. Physicochemical), and complete three-way divergence.

### 3.10 Robustness Suite
The primary classification findings were subjected to systematic sensitivity checks:
1. **Seed Sensitivity:** Cross-validation repeated across an independent pseudo-random seed set (`[7, 17, 37, 73, 97]`);
2. **Class-Weighting Sensitivity:** Comparison of Unweighted, Square-Root, and Inverse-Frequency weighting;
3. **Dimensionality Reduction:** In-fold Principal Component Analysis retaining $90\%$, $95\%$, and $99\%$ variance;
4. **Leave-One-Protein-Out (LOO):** 74 holdout trials systematically removing each individual protein to quantify sample leverage;
5. **Geometric Outlier Sensitivity:** Label-independent isolation forest filtering of representation extremes.

### 3.11 Literature-Anchored Biological Contextualization
To contextualize computational errors without modifying annotations or claiming mechanistic discovery, each of the 74 proteins was audited against primary literature indexed in PubMed and NCBI GenBank. A structured evidence hierarchy was enforced:
1. **Documented Fact:** Empirically demonstrated biochemical function, localization, or structure supported by peer-reviewed literature with PubMed IDs (PMIDs) and DOIs;
2. **Computational Observation:** Measurable feature values, embedding distances, decision margins, or cross-validation predictions;
3. **Interpretive Hypothesis:** Plausible, literature-informed biological explanations for why a protein’s physical characteristics or lifecycle roles might align with an alternative temporal class in representation space.

---

## 4. RESULTS

### 4.1 Dataset Composition and Proteome Overview
The 74 unique canonical proteins of HSV-1 strain 17 span a wide spectrum of length, mass, and biophysical properties (Table 1, Figure 1). Sequence lengths range from 88 amino acids (`US12`/ICP47, 9.8 kDa) to 3,139 amino acids (`UL36`/VP1/2, 336.1 kDa), with a proteome-wide mean length of $522.9 \pm 438.8$ residues. The viral proteome is strongly GC-rich ($68.3\%$ genomic GC content), which is reflected in elevated molar fractions of alanine ($11.9\% \pm 2.9\%$) and proline ($9.3\% \pm 3.6\%$) (Table 2).

The temporal class distribution exhibits severe natural imbalance:
- Immediate-Early ($\alpha$): $N=5$ ($6.8\%$)
- Early ($\beta$): $N=15$ ($20.3\%$)
- Late ($\gamma$): $N=54$ ($73.0\%$)
Under this distribution, a naive majority-class classifier that assigns all proteins to Late achieves an Accuracy of $72.97\%$, but a Balanced Accuracy of only $33.33\%$ and an Immediate-Early Recall of $0.0\%$.

---

### 4.2 Physicochemical Feature Profile
Descriptive analysis of the 25 physicochemical features reveals structural diversity across functional functional groups (Table 2). Theoretical isoelectric points range from acidic ($pI = 4.57$, `UL37` tegument) to highly basic ($pI = 11.83$, `RL1`/ICP34.5 neurovirulence factor). In vitro instability indices range from $24.7$ (`UL18`/VP23 capsid triplex) to $72.6$ (`UL54`/ICP27 transactivator), with $81.1\%$ (60/74) of viral proteins classified as computationally unstable in vitro ($\text{Instability Index} > 40$). Aromaticity averages $7.6\% \pm 2.1\%$, peaking in hydrophobic transmembrane envelope proteins (`UL43`, `UL20`).

---

### 4.3 Contextual ProtBERT Representation Space
ProtBERT transformer embeddings projected into lower-dimensional space via Principal Component Analysis (PCA) show distinct structural organization (Figure 2). The top two principal components account for $13.8\%$ and $9.4\%$ of total ProtBERT variance ($23.2\%$ cumulative), compared to $18.4\%$ and $13.1\%$ ($31.5\%$ cumulative) in physicochemical space. Dimensionality analysis indicated that 14 principal components captured $95\%$ of physicochemical variance, whereas 18 components captured $95\%$ of ProtBERT variance (Table 3).

---

### 4.4 Unsupervised Clustering Does Not Cleanly Reproduce Temporal Classes
In Phase 11 unsupervised clustering benchmarks across $K=2..10$ clusters, no evaluated algorithm (K-Means, Ward, Average, Complete linkage) produced cluster partitions that cleanly aligned with annotated temporal classes (Table 4, Figure 3).
- Adjusted Rand Index (ARI) values across all evaluated configurations remained near the chance baseline;
- Permutation testing against 1,000 label-shuffled null distributions yielded empirical $p$-values that were treated as exploratory evidence;
- Principal component projections and co-clustering stability matrices demonstrated that unsupervised clustering groups proteins primarily by biophysical architecture—readily separating hydrophobic multi-pass envelope glycoproteins (e.g., `UL10`, `UL20`, `UL43`, `UL49A`) and structural capsid subunits from soluble globular enzymes—rather than by transcriptional induction kinetics.

These results indicate that experimentally defined temporal expression class is not simply equivalent to global sequence-space similarity in this dataset.

---

### 4.5 Supervised Classification Benchmarks
Under strictly leakage-controlled 5-fold Stratified Cross-Validation repeated across 5 pseudo-random seeds (25 validation folds), supervised classifiers recovered substantial information about temporal annotations (Table 5, Figure 4).

The primary **Combined Equal-Block Logistic Regression** model achieved the following primary Phase 12 cross-validation estimates:
- **Accuracy:** $76.51\% \pm 10.2\%$
- **Balanced Accuracy:** $75.86\% \pm 12.3\%$
- **Macro-F1:** $0.684 \pm 0.12$
- **Matthews Correlation Coefficient (MCC):** $0.533 \pm 0.19$
- **Immediate-Early Recall:** $80.0\%$ (F1: $0.561$)
- **Early Recall:** $69.3\%$ (F1: $0.687$)
- **Late Recall:** $78.3\%$ (F1: $0.852$)

Linear SVM under equal-block fusion achieved comparable performance (Balanced Accuracy $= 73.23\% \pm 11.8\%$, Macro-F1 $= 0.671$), while non-linear models (Random Forest, XGBoost) exhibited lower balanced accuracy ($61.2\%$ to $66.4\%$), reflecting regularization advantages of linear models in small-sample, high-dimensional regimes.

---

### 4.6 Effect of Class Balancing on Minority-Class Recovery
Evaluating class-weighting regimes demonstrated the impact of explicit imbalance compensation on minority temporal classes (Table 8):
- **Unweighted Loss:** Models collapsed onto the majority Late class, achieving an overall Accuracy of $74.6\%$ but an Immediate-Early Recall of **$0.000$ ($0.0\%$)**, Early Recall of $45.3\%$, and Balanced Accuracy of $41.1\%$.
- **Square-Root Weighting:** Provided modest improvement, raising Immediate-Early Recall to **$0.120$ ($12.0\%$)** and Balanced Accuracy to $59.2\%$.
- **Inverse-Frequency Weighting:** Substantially improved minority recovery, achieving an Immediate-Early Recall of **$0.800$ ($80.0\%$)**, Early Recall of $69.3\%$, and Balanced Accuracy of $75.86\% \pm 12.3\%$ while maintaining majority Late Recall at $78.3\%$.

These experiments confirm that inverse-frequency class weighting substantially improved minority-class recovery in the evaluated models without destabilizing majority-class performance.

---

### 4.7 Representation Disagreement and Combined Performance
Comparing majority out-of-fold predictions across representations revealed that $59.5\%$ (44/74) of proteins achieved unanimous agreement across all spaces, while **$40.5\%$ (30/74)** exhibited representation-dependent discordance (Table 9, Figure 5). Zero proteins exhibited complete three-way divergence.

Among the 30 discordant proteins:
- The Combined Equal-Block model followed ProtBERT predictions in 19 cases ($25.7\%$ of the proteome);
- The Combined Equal-Block model followed Physicochemical predictions in 11 cases ($14.9\%$ of the proteome).

ProtBERT representations showed distinct predictive consistency for multi-subunit enzymatic complexes (e.g., DNA polymerase `UL30`, thymidine kinase `UL23`), whereas Physicochemical descriptors provided distinct classification signals for proteins with extreme charge, size, or disorder (such as `RL1`/ICP34.5 with pI 11.8; `RL2`/ICP0). The equal-block combined representation produced stable classification performance across the predefined evaluation and sensitivity analyses without feature-count dominance.

---

### 4.8 Protein-Level Error Profiling
Out-of-fold error analysis across the 5 cross-validation repeats identified that $73.0\%$ (54/74) of proteins were consistently classified correctly ($\ge 80\%$ cross-validation accuracy), while classification difficulty was concentrated in **11 proteins** (Table 7, Figure 6).

- **High-Confidence Incorrect ($\Delta p > 0.35$):**
  1. `UL15` (DNA packaging terminase large subunit 1, true Late, predicted Early, $\Delta p = 0.849$);
  2. `RL1` (ICP34.5 neurovirulence factor, true Late, predicted Early, $\Delta p = 0.771$);
  3. `US12` (ICP47 TAP inhibitor, true Immediate-Early, predicted Late, $\Delta p = 0.669$);
  4. `UL11` (myristoylated tegument protein, true Late, predicted Early, $\Delta p = 0.608$);
  5. `UL8` (helicase-primase accessory subunit, true Early, predicted Late, $\Delta p = 0.424$);
  6. `UL24` (nuclear/nucleolar protein, true Late, predicted Early, $\Delta p = 0.361$);
  7. `UL41` (virion host shutoff endoribonuclease, true Late, predicted Early, $\Delta p = 0.352$).

- **Moderate-Margin Incorrect ($0.15 \le \Delta p \le 0.35$):**
  1. `UL36` (large tegument hub VP1/2, true Late, predicted Early, $\Delta p = 0.333$);
  2. `UL49` (VP22 major tegument protein, true Late, predicted Early, $\Delta p = 0.335$);
  3. `UL52` (helicase-primase primase subunit, true Early, predicted Late, $\Delta p = 0.274$);
  4. `UL13` (tegument protein kinase, true Late, predicted Early, $\Delta p = 0.254$).

`US12` (88 aa, true Immediate-Early) was a high-confidence computational misclassification despite its Immediate-Early annotation, illustrating that temporal expression class is not necessarily aligned with simple sequence-space similarity. `UL36` (3,139 aa, true Late, frequent prediction Early, mean margin $\approx 0.333$) represented a computationally difficult case, potentially reflecting its extreme sequence length and complex protein architecture; this interpretation remains a hypothesis rather than a demonstrated causal mechanism.

---

### 4.9 Methodological Robustness and Sensitivity Suite
Systematic sensitivity evaluations (Phase 14) established the stability of the primary computational conclusions (Table 8, Figure 7):
1. **CV Seed Invariance:** Evaluating the independent sensitivity seed set (`[7, 17, 37, 73, 97]`) yielded Balanced Accuracy of $0.754 \pm 0.017$ and Macro-F1 of $0.665$, which was consistent with the original seed set ($0.758 \pm 0.020$).
2. **In-Fold PCA Reduction:** Retaining $90\%$ variance (~22 PCs) or $95\%$ variance (~30 PCs) in-fold retained $>98\%$ of full-dimensional classification performance while reducing feature dimensionality by $>97\%$.
3. **Leave-One-Protein Stability:** Systematic 74-trial holdout tests revealed that no single protein acted as a catastrophic pivot for model performance ($|\Delta \text{ Balanced Accuracy}| < 0.015$ for all Early/Late proteins).
4. **Outlier Sensitivity:** Performance changed when label-independent geometric outliers (`UL36`, `US12`, `UL41`) were excluded, adjusting Balanced Accuracy to $0.782$. Because these analyses alter the evaluated protein set, they are treated as sensitivity analyses rather than evidence that protein removal improves the primary model.

---

### 4.10 Biological Contextualization
Mapping all 74 proteins to peer-reviewed literature demonstrated that classification difficulty was enriched among virion-packaged structural proteins and multi-functional factors that execute activities across multiple viral lifecycle phases (Table 9, Figure 8).

---

## 5. DISCUSSION

### 5.1 Principal Findings
This study establishes a systematic computational benchmark for sequence representation and temporal classification in the HSV-1 proteome. Our findings demonstrate that:
1. Unsupervised representation geometry does not cleanly mirror transcriptional induction classes;
2. Supervised models recover substantial information about temporal annotations when class imbalance is explicitly addressed;
3. Contextual transformer embeddings (ProtBERT) and physical descriptors encode partially distinct predictive information;
4. Equal-block multi-modal representation produces stable classification performance across predefined sensitivity analyses;
5. Computational errors are concentrated in a specific subset of multi-functional and virion-packaged proteins.

### 5.2 Temporal Expression Versus Sequence Representation
A central conceptual finding of this study is that *experimentally defined temporal expression class is not simply equivalent to global sequence-space similarity in this dataset*.

Transcriptional induction timing in HSV-1 is driven by genomic promoter architecture, transcription factor binding sites (e.g., Oct-1/HCF-1 binding motifs in $\alpha$ promoters), and chromatin accessibility (Roizman et al., 2013). In contrast, protein sequence representations capture structural, catalytic, and biophysical properties. In unsupervised representation space, proteins group into functional biophysical clusters—such as multi-pass transmembrane envelope glycoproteins, icosahedral capsid shells, and globular enzymes. Because each temporal kinetic phase encompasses diverse structural and enzymatic roles, unsupervised clustering does not recapitulate the temporal cascade without label supervision.

### 5.3 Comparison of Physicochemical and ProtBERT Representations
Physicochemical descriptors and ProtBERT embeddings exhibited distinct predictive behavior across $40.5\%$ (30/74) of viral proteins.
- ProtBERT embeddings capture contextual sequence patterns that reflect complex enzymatic domains, providing consistent classification for catalytic replication machinery (such as DNA polymerase `UL30` and thymidine kinase `UL23`);
- Physicochemical features summarize global physical constraints, providing distinct signals for proteins with extreme charge, length, or amino acid bias (such as basic regulatory factor `RL1`/ICP34.5).
The equal-block scaling transformation ($[Z_{\text{phys}}/\sqrt{25} \;\|\; Z_{\text{pb}}/\sqrt{1024}]$) successfully balanced total block variance, preventing the 1024-dimensional transformer embedding from dominating the 25 physical features while maintaining classification stability.

### 5.4 Role of Class Imbalance in Viral Proteomics
Severe class imbalance is an intrinsic biological property of herpesvirus proteomes: Late structural proteins vastly outnumber Immediate-Early regulatory proteins ($54$ Late vs. $5$ IE in HSV-1). Under standard unweighted loss formulations, supervised classifiers collapsed onto the majority class, yielding zero minority-class recall. In-fold inverse-frequency class weighting substantially improved minority recovery (Immediate-Early Recall increased from $0.0\%$ to $80.0\%$) without compromising majority Late accuracy ($78.3\%$). This underscores the necessity of explicit class-balancing in viral computational genomics.

---

### 5.5 Biological Case Studies of Difficult Proteins

#### Case Study 1: Large Tegument Hub UL36 (VP1/2)
- **Documented Biological Fact:** `UL36` is a 3,139-amino-acid giant tegument protein (~336 kDa) that forms the inner tegument shell, bridging capsid pentons to motor proteins (dynein/kinesin) for microtubule transport and docking at the nuclear pore complex (Sandbaumhüter et al., 2013; McNabb & Courtney, 1992).
- **Computational Observation:** `UL36` is annotated as Late ($\gamma_1$) but was frequently predicted as Early with a mean decision margin of $\Delta p = 0.333 \le 0.35$, categorizing it as Moderate-Margin Incorrect. In geometric analysis, `UL36` represents the most extreme length outlier in the proteome.
- **Interpretive Hypothesis:** The classification difficulty of `UL36` may reflect its extreme sequence length, extensive disordered linkers, and early functional role in post-entry uncoating; this interpretation remains a hypothesis rather than a demonstrated causal mechanism.

#### Case Study 2: TAP Inhibitor Peptide US12 (ICP47)
- **Documented Biological Fact:** `US12` is an 88-amino-acid cytosolic peptide expressed with Immediate-Early kinetics that binds the peptide-binding site of the transporter associated with antigen processing (TAP), blocking MHC class I antigen presentation to cytotoxic T lymphocytes (York et al., 1994).
- **Computational Observation:** `US12` is annotated as Immediate-Early ($\alpha$) but was consistently predicted as Late with high confidence ($\Delta p = 0.669$), making it a High-Confidence Incorrect prediction.
- **Interpretive Hypothesis:** `US12` lacks the DNA-binding motifs and multi-domain architecture characteristic of canonical IE transcriptional regulators (`RL2`/ICP0, `RS1`/ICP4, `UL54`/ICP27), illustrating that temporal expression kinetics do not necessarily align with sequence-space similarity.

#### Case Study 3: DNA Packaging Terminase Subunit 1 UL15
- **Documented Biological Fact:** `UL15` (735 aa) is a late structural component of the tripartite terminase complex (UL15/UL28/UL33) that executes ATP-dependent cleavage and translocation of viral concatemeric DNA into procapsids (Baines et al., 1994; Yu & Weller, 1998).
- **Computational Observation:** Annotated as Late ($\gamma_2$), `UL15` was consistently predicted as Early with high confidence ($\Delta p = 0.849$).
- **Interpretive Hypothesis:** `UL15` contains an ATPase/nuclease core that shares sequence and structural features with early catalytic replication enzymes (such as helicase `UL5`), providing a plausible computational explanation for its representation-space alignment with early factors.

#### Case Study 4: Neurovirulence Factor RL1 (ICP34.5)
- **Documented Biological Fact:** `RL1` (263 aa) is a multi-functional neurovirulence factor that recruits host protein phosphatase 1$\alpha$ (PP1$\alpha$) to dephosphorylate eIF2$\alpha$, reversing PKR-mediated host translational arrest and inhibiting autophagy (Chou & Roizman, 1990; He et al., 1997).
- **Computational Observation:** Annotated as Late ($\gamma_1$), `RL1` was predicted as Early with high confidence ($\Delta p = 0.771$).
- **Interpretive Hypothesis:** `RL1` possesses an extremely basic isoelectric point ($pI = 11.83$), high proline/alanine content, and distinct repetitive domains, positioning it at the biophysical boundary of late structural proteins.

#### Case Study 5: Virion Host Shutoff Endoribonuclease UL41 (vhs)
- **Documented Biological Fact:** `UL41` (489 aa) is packaged into the virion tegument during late assembly and released into the host cytoplasm immediately upon uncoating, where it non-specifically degrades host and viral mRNAs (Read & Frenkel, 1983; Kwong et al., 1988).
- **Computational Observation:** Annotated as Late ($\gamma_1$), `UL41` was predicted as Early with high confidence ($\Delta p = 0.352$).
- **Interpretive Hypothesis:** Although synthesized late, `UL41` functions as a catalytic endoribonuclease during immediate post-entry stages, exhibiting enzymatic sequence characteristics resembling early metabolic enzymes.

#### Case Study 6: Small Membrane-Associated Tegument Protein UL11
- **Documented Biological Fact:** `UL11` is a small (96 aa) myristoylated and palmitoylated tegument protein that associates with membranes and interacts with `UL16` and `UL21` to coordinate secondary envelopment (MacLean et al., 1989; Baines et al., 1995).
- **Computational Observation:** Annotated as Late ($\gamma_2$), `UL11` was predicted as Early with high confidence ($\Delta p = 0.608$).
- **Interpretive Hypothesis:** Its small size, high negative charge, and fatty acid modification motifs distinguish it from large structural tegument proteins, placing it in an isolated position in feature space.

---

### 5.6 Status of Glycoprotein C (UL44)
Envelope glycoprotein C (`UL44`/gC) is annotated with primary temporal class Late and late subclass *Conflicting* ($\gamma_1/\gamma_2$) in canonical records due to literature demonstrating partial expression prior to DNA replication followed by marked upregulation post-replication (Frink et al., 1983). In our analysis, `UL44` was preserved in all primary benchmarks as Late and was consistently classified correctly by the primary Combined model, confirming that subclass ambiguity does not compromise three-class supervised recovery.

---

### 5.7 Methodological Robustness and Sensitivity
The primary classification findings proved robust across predefined analytical perturbations:
- **Seed Invariance:** Cross-validation on an independent sensitivity seed set produced consistent results ($0.754 \pm 0.017$ vs. $0.758 \pm 0.020$).
- **Dimensionality Reduction:** In-fold PCA retaining $90\%$ to $95\%$ variance preserved $>98\%$ of classification performance while compressing feature dimensionality by $>97\%$, demonstrating that predictive signals are concentrated in low-dimensional manifolds.
- **Leave-One-Protein Stability:** Systematic 74-trial holdouts confirmed that model performance does not depend on any single highly influential protein ($|\Delta \text{ Balanced Accuracy}| < 0.015$).

---

### 5.8 Limitations
We explicitly acknowledge several methodological and biological limitations:
1. **Single Viral Strain:** Analysis is restricted to HSV-1 strain 17; strain-specific polymorphisms and intra-species variation were not evaluated.
2. **Dataset Cardinality:** The total proteome is bounded by the viral genome size ($N=74$ unique canonical proteins).
3. **Small Minority Class Size:** With only 5 Immediate-Early proteins, individual misclassifications represent discrete $20\%$ increments in recall.
4. **External Temporal Annotations:** Temporal categories are externally curated kinetic labels that simplify continuous transcriptional dynamics into discrete bins.
5. **Computational Associations vs. Causal Mechanisms:** Model predictions reflect statistical correlations in feature space and do not establish in vivo transcriptional mechanisms.
6. **Transformer Representation Semantics:** ProtBERT embeddings represent contextual language model representations; individual dimensions cannot be mapped to isolated biological functions.
7. **Cross-Validation Sample Constraints:** Standard deviations reflect fold-to-fold variability in small-sample cross-validation partitions.
8. **Lack of Pan-Herpesvirus Validation:** Findings should not be assumed to generalize to beta- or gammaherpesviruses without direct empirical evaluation.
9. **Outlier Filtering Interpretation:** Outlier exclusion alters the evaluated dataset and is reported strictly as a sensitivity analysis rather than a dataset replacement.
10. **Hypothetical Nature of Interpretations:** Biological contextualizations provide plausible interpretations for difficult cases, but remain hypotheses requiring experimental validation.

---

### 5.9 Future Directions
Future computational and experimental studies could extend this framework by:
1. Conducting comparative representation analysis across other human and animal herpesviruses (HSV-2, VZV, HCMV, EBV, KSHV);
2. Integrating multi-omics features, including quantitative temporal transcriptomics, ribosome profiling, and structural models from AlphaFold;
3. Evaluating fine-tuned viral-specific protein language models trained on comprehensive viral genomic corpora.

---

## 6. CONCLUSION

This study establishes a reproducible, leakage-controlled computational evaluation of sequence representation and temporal classification across the complete proteome of HSV-1 strain 17. Unsupervised representation geometry reflects biophysical protein architecture rather than transcriptional induction kinetics. In contrast, leakage-controlled supervised learning recovers substantial information about temporal annotations when severe class imbalance is explicitly addressed via inverse-frequency weighting. ProtBERT and physicochemical descriptors encode partially distinct predictive signals, and equal-block multi-modal fusion provides stable classification performance across predefined sensitivity analyses. While biological contextualization offers plausible explanations for computationally difficult proteins, these interpretations remain exploratory hypotheses. These findings establish a solid foundation for comparative computational virology and representation learning.

---

## DATA AVAILABILITY
The complete HSV-1 strain 17 genome and protein sequences are publicly available through NCBI GenBank/RefSeq under accession `NC_001806.2`. All processed feature matrices (`physicochemical_features.csv`, `X_physicochemical.npy`, `X_protbert.npy`, `X_combined_raw.npy`) and canonical temporal annotations (`temporal_annotations_final.csv`) are permanently archived in the project repository under `data/processed/` and `data/annotations/`.

## CODE AVAILABILITY
The complete computational pipeline, including feature extraction, embedding generation, cross-validation workflows, error analysis, robustness testing, and unit test suites, is open-source and available in the project GitHub repository: `https://github.com/Yuyutsu01/HSV1-Proteome-Analysis`.

## ETHICS STATEMENT
This computational study exclusively analyzed publicly available genomic and proteomic sequence data; no human subjects or animal models were involved.

## AUTHOR CONTRIBUTIONS
**Shiva:** Conceptualization, Methodology, Software, Data Curation, Formal Analysis, Investigation, Validation, Visualization, Writing – Original Draft, Writing – Review & Editing.

## CONFLICT OF INTEREST
The author declares that there are no financial or commercial conflicts of interest associated with this research.

## FUNDING STATEMENT
This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## ACKNOWLEDGMENTS
Computational resources and development environment were supported by the Antigravity IDE platform and the open-source scientific Python community.

---

## REFERENCES

1. Albecka A, et al. (2017). Dual function of the HSV-1 UL7/UL51 complex in secondary envelopment. *J Virol* 91:e00897-17. DOI: 10.1128/jvi.00897-17. PMID: 28539446.
2. Albright AG, et al. (1999). The herpes simplex virus UL37 protein is required for secondary envelopment. *J Virol* 73:8000-8004. DOI: 10.1128/jvi.73.10.8000-8004.1999. PMID: 10482542.
3. Ali MA, et al. (1996). Identification of the herpes simplex virus type 1 UL25 gene product. *J Virol* 70:5258-5267. DOI: 10.1128/jvi.70.8.5258-5267.1996. PMID: 8627776.
4. Baines JD, Roizman B. (1991). The open reading frames UL10, UL43, and UL49A of HSV-1 encode nonessential glycoproteins. *J Virol* 65:938-944. DOI: 10.1128/jvi.65.2.938-944.1991. PMID: 1656092.
5. Baines JD, et al. (1991). The UL20 gene of herpes simplex virus 1 encodes a function necessary for viral egress. *J Virol* 65:6414-6424. DOI: 10.1128/jvi.65.12.6414-6424.1991. PMID: 1656094.
6. Baines JD, et al. (1994). The UL15 gene of herpes simplex virus type 1 contains an open reading frame required for viral DNA cleavage and packaging. *J Virol* 68:2929-2936. DOI: 10.1128/jvi.68.5.2929-2936.1994. PMID: 8189528.
7. Baines JD, et al. (1995). The UL11 gene of herpes simplex virus 1 encodes a small myristoylated protein. *J Virol* 69:825-833. DOI: 10.1128/jvi.69.2.825-833.1995. PMID: 7609051.
8. Beard PM, et al. (2002). The herpes simplex virus 1 terminase complex consists of UL15, UL28, and UL33. *J Virol* 76:4785-4791. DOI: 10.1128/jvi.76.10.4785-4791.2002. PMID: 11967295.
9. Bigalke JM, Heldwein EE. (2015). Structural basis of membrane budding by the nuclear egress complex of herpesviruses. *Cell* 160:1107-1119. DOI: 10.1016/j.cell.2015.01.044. PMID: 25723167.
10. Booy FP, et al. (1994). The capsid of herpes simplex virus type 1: mapping of the VP26 protein. *PNAS* 91:5652-5656. DOI: 10.1073/pnas.91.12.5652. PMID: 8202542.
11. Boutell C, Everett RD. (2013). The herpes simplex virus type 1 (HSV-1) regulatory protein ICP0. *J Gen Virol* 94:665-681. DOI: 10.1099/vir.0.050625-0. PMID: 23364192.
12. Chou J, Roizman B. (1990). The gamma(1)34.5 gene of herpes simplex virus 1 precludes neurovirulence and facilitates growth in human cells. *J Virol* 64:1014-1020. DOI: 10.1128/jvi.64.3.1014-1020.1990. PMID: 2172479.
13. Conley AJ, et al. (1981). Expression of alpha, beta, and gamma genes in cells infected with herpes simplex virus 1. *J Virol* 37:191-206. DOI: 10.1128/jvi.37.1.191-206.1981. PMID: 6257913.
14. Costa RH, et al. (1981). Major capsid protein VP5 of herpes simplex virus type 1. *J Virol* 38:483-496. DOI: 10.1128/jvi.38.2.483-496.1981. PMID: 6271960.
15. Crute JJ, et al. (1989). Herpes simplex virus 1 helicase-primase: a complex of three herpes-encoded gene products. *PNAS* 86:2186-2189. DOI: 10.1073/pnas.86.7.2186. PMID: 2538836.
16. Cunningham C, et al. (1992). The UL13 protein kinase of herpes simplex virus 1. *J Gen Virol* 73:303-311. DOI: 10.1099/0022-1317-73-2-303. PMID: 1324546.
17. Cunningham C, et al. (2000). The minor tegument protein UL14 of herpes simplex virus 1. *J Virol* 74:33-41. DOI: 10.1128/jvi.74.1.33-41.2000. PMID: 10627552.
18. Dai X, Zhou ZH. (2018). Structure of the herpes simplex virus 1 capsid with associated tegument protein complexes. *Science* 360:eaao7298. DOI: 10.1126/science.aao7298. PMID: 29622628.
19. Digard P, et al. (1993). The herpes simplex virus DNA polymerase catalytic subunit UL30. *J Virol* 67:398-406. DOI: 10.1128/jvi.67.1.398-406.1993. PMID: 8380080.
20. Elnaggar A, et al. (2021). ProtTrans: Towards Cracking the Language of Life's Code Through Self-Supervised Deep Learning and High Performance Computing. *IEEE TPAMI* 44:7112-7127. DOI: 10.1109/TPAMI.2021.3095381. PMID: 34232869.
21. Everett RD, Maul GG. (1994). HSV-1 regulatory protein ICP0 induces the degradation of PML. *EMBO J* 13:5062-5069. DOI: 10.1002/j.1460-2075.1994.tb06835.x. PMID: 7957076.
22. Frink RJ, et al. (1983). Detailed analysis of the glycoprotein C gene of herpes simplex virus type 1. *J Virol* 45:634-647. DOI: 10.1128/jvi.45.2.634-647.1983. PMID: 6296435.
23. Gibbs JS, et al. (1985). Sequence and mapping of the HSV-1 DNA polymerase gene. *PNAS* 82:7969-7973. DOI: 10.1073/pnas.82.23.7969. PMID: 2999787.
24. Guruprasad K, et al. (1990). Correlation between stability of a protein and its dipeptide composition. *Protein Eng* 4:155-161. DOI: 10.1093/protein/4.2.155. PMID: 2075190.
25. He B, et al. (1997). The gamma(1)34.5 protein of herpes simplex virus 1 complexes with protein phosphatase 1alpha to dephosphorylate eIF-2alpha. *PNAS* 94:843-848. DOI: 10.1073/pnas.94.3.843. PMID: 9023344.
26. Heldwein EE, et al. (2006). Crystal structure of glycoprotein B from herpes simplex virus 1. *Science* 313:217-220. DOI: 10.1126/science.1126548. PMID: 16840698.
27. Honess RW, Roizman B. (1974). Regulation of herpesvirus macromolecular synthesis. I. Cascade regulation of the synthesis of three groups of viral proteins. *J Virol* 14:8-19. DOI: 10.1128/jvi.14.1.8-19.1974. PMID: 4365321.
28. Kwong AD, et al. (1988). The virion host shutoff function of herpes simplex virus. *J Virol* 62:912-921. DOI: 10.1128/jvi.62.3.912-921.1988. PMID: 2828682.
29. Kyte J, Doolittle RF. (1982). A simple method for displaying the hydropathic character of a protein. *J Mol Biol* 157:105-132. DOI: 10.1016/0022-2836(82)90515-0. PMID: 7108955.
30. MacLean CA, et al. (1989). The product of the UL11 gene of herpes simplex virus type 1 is a small myristoylated protein. *J Gen Virol* 70:3147-3157. DOI: 10.1099/0022-1317-70-12-3147. PMID: 2552093.
31. McNabb DS, Courtney RJ. (1992). Identification and characterization of the herpes simplex virus type 1 UL36 gene product. *J Virol* 66:2655-2663. DOI: 10.1128/jvi.66.5.2655-2663.1992. PMID: 1314958.
32. Patel AH, MacLean CA. (1995). The product of the UL18 gene of herpes simplex virus type 1 is a capsid triplex protein (VP23). *Virology* 206:1106-1112. DOI: 10.1006/viro.1995.1034. PMID: 7750031.
33. Preston CM. (1979). Control of herpes simplex virus type 1 mRNA synthesis in cells infected in the absence of protein synthesis. *J Virol* 29:275-284. DOI: 10.1128/jvi.29.1.275-284.1979. PMID: 219982.
34. Read GS, Frenkel N. (1983). Herpes simplex virus mutants defective in the virion-associated shutoff of host protein synthesis. *J Virol* 46:498-512. DOI: 10.1128/jvi.46.2.498-512.1983. PMID: 6300086.
35. Roizman B, et al. (2013). Herpes Simplex Viruses. In: *Fields Virology*, 6th Edition, Knipe DM, Howley PM (Eds.), Lippincott Williams & Wilkins, Philadelphia. ISBN: 978-1-4511-0563-6.
36. Sandbaumhüter M, et al. (2013). Cytosolic herpes simplex virus capsids not only use dynein for retrograde transport but also recruit kinesin-1 and kinesin-2. *PLoS Pathog* 9:e1003230. DOI: 10.1371/journal.ppat.1003230. PMID: 23516362.
37. Weller SK, et al. (1983). Genetic analysis of temperature-sensitive mutants of HSV-1: the single-stranded DNA-binding protein ICP8. *J Virol* 45:354-366. DOI: 10.1128/jvi.45.1.354-366.1983. PMID: 6296430.
38. York IA, et al. (1994). A cytosolic herpes simplex virus protein inhibits antigen presentation to CD8+ T lymphocytes. *Cell* 77:525-535. DOI: 10.1016/0092-8674(94)90215-4. PMID: 7514502.
39. Yu D, Weller SK. (1998). Herpes simplex virus type 1 cleavage and packaging proteins UL15 and UL28 interact. *J Virol* 72:7428-7439. DOI: 10.1128/jvi.72.9.7428-7439.1998. PMID: 9696841.

---

## SUPPLEMENTARY MATERIAL

The complete supplementary archive supporting this study is indexed in `manuscript/SUPPLEMENTARY_MATERIALS.md` and contains Supplementary Tables S1 through S13, embedding extraction metadata, and reproducibility test suites.
