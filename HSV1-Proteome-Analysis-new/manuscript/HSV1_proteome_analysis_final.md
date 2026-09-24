# Computational Representation, Unsupervised Clustering, and Leakage-Controlled Classification of the Herpes Simplex Virus Type 1 Proteome

**Authors:** Shiva  
**Affiliation:** Antigravity IDE Computational Virology Group  
**Target Submission:** *Journal of Virology* / *Bioinformatics*  
**Date:** September 2026  

---

## ABSTRACT

**Background:** Herpes simplex virus type 1 (HSV-1) executes a tightly coordinated temporal gene expression cascade traditionally categorized into Immediate-Early ($\alpha$), Early ($\beta$), and Late ($\gamma$) classes. While temporal kinetics are driven by viral transactivators, promoter architectures, and host interactions, the degree to which viral protein sequences alone encode predictable temporal and functional signatures remains a fundamental computational question.

**Objective:** To systematically evaluate how effectively HSV-1 proteins can be represented and classified according to experimentally defined temporal expression classes using physicochemical descriptors and contextual protein language model embeddings (ProtBERT), and to determine what representation geometry reveals about viral protein organization.

**Methods:** We analyzed the complete proteome of HSV-1 strain 17 (RefSeq `NC_001806.2`, $N=74$ unique canonical proteins; IE=5, Early=15, Late=54). Three primary representation spaces were constructed: (i) an exact 25-feature physicochemical descriptor matrix, (ii) a 1024-dimensional ProtBERT embedding matrix with length-weighted chunk pooling (500 aa chunk size, 100 aa overlap, 400 aa step), and (iii) an equal-block multi-modal representation ($Z_{\text{phys}}/\sqrt{25} \,\|\, Z_{\text{pb}}/\sqrt{1024}$, 1049 dimensions total). Unsupervised clustering ($K=2..10$) was evaluated against label-permutation null distributions ($N=1,000$). Supervised classification was performed under 5-fold Stratified Cross-Validation repeated across 5 pseudo-random seeds with strictly in-fold preprocessing and inverse-frequency class-balancing. Protein-level prediction errors, representation disagreements, and sensitivity to analytical choices were systematically analyzed and contextualized using peer-reviewed biological literature.

**Results:** Across the evaluated clustering configurations, temporal-class recovery was limited and varied substantially across representations and linkage methods; no clustering configuration cleanly reproduced the three temporal classes. Permutation testing was used to assess the observed association with the temporal labels relative to a label-permutation null distribution. These permutation results were treated as exploratory evidence rather than as independent confirmatory tests. Unsupervised geometry reflected biophysical characteristics (e.g. membrane glycoproteins vs globular enzymes) rather than induction kinetics. In contrast, leakage-controlled supervised learning recovered substantial information about temporal annotations: the primary Combined Equal-Block Logistic Regression model achieved Balanced Accuracy of $75.86\% \pm 12.3\%$, Macro-F1 of $0.684 \pm 0.12$, and minority Immediate-Early Recall of $80.0\%$ (compared to a majority-class baseline of Accuracy = $72.97\%$, Balanced Accuracy = $33.33\%$, IE Recall = $0.0\%$). Inverse-frequency class weighting substantially improved minority-class recovery in the evaluated experiments (IE Recall = $0.000$ unweighted, $0.120$ square-root weighted, $0.800$ inverse-frequency weighted). ProtBERT and physicochemical representations exhibited different predictive behavior across a subset of proteins ($40.5\%$, 30/74 discordance). Among the 30 proteins showing disagreement between the two modality-specific predictions, the combined prediction followed ProtBERT for 19 proteins and the physicochemical representation for 11 proteins. Eleven computationally difficult proteins (e.g., `UL15`, `RL1`, `US12`, `UL41`, `UL36`) were contextualized using literature-documented lifecycle roles, providing plausible biological context without demonstrating causal mechanisms.

**Conclusions:** HSV-1 temporal expression classes are not simply equivalent to sequence-space similarity. Supervised models can recover substantial information about temporal annotations when class imbalance is explicitly addressed. Physicochemical and ProtBERT representations exhibit distinct predictive behavior, and the equal-block combined representation produced stable classification performance across the predefined evaluation and sensitivity analyses. Biological contextualization provides plausible interpretations for difficult cases, but these interpretations remain hypotheses requiring experimental validation.

---

## 1. INTRODUCTION

Herpes simplex virus type 1 (HSV-1) is a ubiquitous human alphaherpesvirus with a double-stranded DNA genome of approximately 152 kilobase pairs encoding over 70 distinct proteins (Roizman et al., 2013). During lytic infection, viral gene expression proceeds through a tightly regulated, sequential cascade:
1. **Immediate-Early ($\alpha$):** Five viral genes (`RL2`/ICP0, `RS1`/ICP4, `UL54`/ICP27, `US1`/ICP22, `US12`/ICP47) expressed immediately upon host cell entry independent of de novo viral protein synthesis, primarily functioning as master transcriptional transactivators and innate immune modulators.
2. **Early ($\beta$):** Genes transcribed in response to immediate-early transactivation, encoding the core catalytic replication machinery (DNA polymerase, helicase-primase, origin-binding protein) and nucleotide metabolic enzymes.
3. **Late ($\gamma$):** Genes whose maximal expression requires viral DNA synthesis, subdivided into leaky-late ($\gamma_1$) and true-late ($\gamma_2$), encoding the structural building blocks of the virion (capsid shell, tegument matrix, and envelope glycoproteins).

While the transcriptional promoter motifs, chromatin states, and transcription factor cascades governing this temporal sequence are extensively characterized, a fundamental computational biology question remains unresolved: *To what extent is temporal classification predictable directly from mature protein sequences, and how do biophysical descriptors compare to modern deep protein language models?*

Prior computational studies in viral proteomics have frequently suffered from methodological limitations, including improper duplicate handling of identical terminal inverted repeats, unweighted learning on severely imbalanced datasets resulting in minority class suppression, feature-count dominance in multi-modal concatenation, and data leakage across cross-validation partitions.

In this study, we present a leakage-controlled computational analysis of the HSV-1 proteome. Using the authoritative strain 17 RefSeq genome (`NC_001806.2`), we evaluate exact 25-feature physicochemical descriptors, 1024-dimensional ProtBERT embeddings, and equal-block multi-modal representations across unsupervised clustering, supervised classification, protein-level error analysis, systematic robustness testing, and literature-anchored biological contextualization.

---

## 2. MATERIALS AND METHODS

### 2.1 Authoritative Dataset and Quality Control
The complete reference genome of HSV-1 strain 17 was retrieved from NCBI GenBank/RefSeq (Accession `NC_001806.2`). Protein-coding sequences were extracted and audited. Exact terminal repeat duplicates (`RL2` YP_009137133.1, `RL1` YP_009137134.1, and `RS1` YP_009137149.1) were removed, yielding an authoritative, non-redundant dataset of **$N=74$ unique canonical proteins** (from 77 raw extracted records). All 74 sequences passed strict quality control (valid IUPAC amino acid alphabet, complete translation, unambiguous start/stop codons).

### 2.2 Temporal Class Annotations
Authoritative temporal annotations were assigned based on peer-reviewed experimental literature and verified against canonical reference texts (Fields Virology; Roizman et al., 2013). The 74 proteins comprise:
- **Immediate-Early ($\alpha$):** 5 proteins ($6.8\%$)
- **Early ($\beta$):** 15 proteins ($20.3\%$)
- **Late ($\gamma$):** 54 proteins ($73.0\%$)

Envelope glycoprotein C (`UL44`/gC) is annotated with primary temporal class Late and late subclass *Conflicting* ($\gamma_1/\gamma_2$) in accordance with documented literature evidence, and is preserved in all analyses.

### 2.3 Physicochemical Feature Representation
An exact 25-feature physicochemical descriptor vector was computed for each protein:
1. Sequence Length ($N_{\text{aa}}$)
2. Molecular Weight (Da)
3. Aromaticity (frequency of F, W, Y)
4. Instability Index (Guruprasad metric)
5. Isoelectric Point (pI)
6–25. 20 standard amino acid molar fractions ($f_{\text{A}}, f_{\text{C}}, \dots, f_{\text{Y}}$).

### 2.4 ProtBERT Embeddings
Contextual sequence representations were generated using the pretrained bidirectional transformer `Rostlab/prot_bert` (Elnaggar et al., 2021). Sequences were processed in 500-residue chunks with a 100-residue overlap (stride = 400). Mean pooling across residue-level hidden states (excluding `[CLS]` and `[SEP]`) was performed for each chunk, followed by length-weighted chunk aggregation. The longest protein, `UL36` (3,139 amino acids), was processed across 8 chunks. Each protein is represented by a 1024-dimensional embedding vector.

### 2.5 Multi-Modal Representation Fusion
To combine the 25 physicochemical features and 1024 ProtBERT dimensions without allowing the massive ProtBERT dimensionality to overwhelm physical descriptors, we implemented an **Equal-Block Scaling Transformation**:
$$Z_{\text{combined}} = \left[ \frac{Z_{\text{phys}}}{\sqrt{25}} \;\Bigg\|\; \frac{Z_{\text{pb}}}{\sqrt{1024}} \right]$$
where $Z_{\text{phys}}$ and $Z_{\text{pb}}$ represent standard normal transformations ($mean=0, SD=1$) fitted strictly inside each training fold, yielding a 1049-dimensional fused representation.

### 2.6 Unsupervised Clustering Protocol
Unsupervised structure was evaluated across $K=2..10$ using K-Means (100 random restarts per seed across 5 predefined seeds), Hierarchical Agglomerative Clustering (Ward, Average, Complete linkage), and subsampling co-clustering stability (100 subsamples at $80\%$ sampling rate). Empirical statistical evaluation against temporal class annotations was conducted using $N=1,000$ label permutations with Davison-Hinkley $(r+1)/(N+1)$ finite-sample correction.

### 2.7 Leakage-Controlled Supervised Classification
Supervised classification was evaluated under 5-fold Stratified Cross-Validation repeated across 5 pseudo-random seeds (`[42, 123, 456, 789, 2026]`, 25 validation folds total). Preprocessing scalers, block-weight transformations, and sample weights were fitted exclusively on the training partition of each fold. Five primary class-balanced classifiers were evaluated: Logistic Regression ($C=1.0$), Linear SVM ($C=1.0$), RBF SVM ($C=1.0$), Random Forest ($n=100$), and XGBoost. Loss functions utilized in-fold inverse-frequency class weights:
$$w_c = \frac{N_{\text{train}}}{K \cdot n_{c, \text{train}}}$$

### 2.8 Statistical Reporting and Robustness
Metrics reported include Accuracy, Balanced Accuracy (macro-averaged recall), Macro-F1, Matthews Correlation Coefficient (MCC), and per-class Recall and F1. Sensitivity analyses tested independent seed sets (`[7, 17, 37, 73, 97]`), alternative weighting regimes (Unweighted, Square-Root, Inverse-Frequency), in-fold PCA dimensionality reduction ($90\%, 95\%, 99\%$ variance), 74 systematic leave-one-protein ablation trials, and label-independent geometric outlier exclusions.

---

## 3. RESULTS

### 3.1 Dataset Characteristics and Feature Space
The 74 HSV-1 proteins exhibit substantial structural diversity, ranging from the 88-amino-acid regulatory peptide `US12` (ICP47, 9.8 kDa) to the 3,139-amino-acid giant tegument protein `UL36` (336.1 kDa) (Table 1, Figure 1). Mean proteome length is $522.9 \pm 438.8$ residues. The proteome is characteristically GC-rich, reflected in high alanine ($11.9\% \pm 2.9\%$) and proline ($9.3\% \pm 3.6\%$) fractions.

### 3.2 Unsupervised Clustering Evaluation and Temporal Structure
Across the evaluated clustering configurations (K-Means, Ward, Average, Complete across $K=2..10$), temporal-class recovery was limited and varied substantially across representations and linkage methods; no clustering configuration cleanly reproduced the three temporal classes (Table 4, Figure 3). Permutation testing was used to assess the observed association with the temporal labels relative to a label-permutation null distribution. These permutation results were treated as exploratory evidence rather than as independent confirmatory tests. Principal component projections reveal that unsupervised geometry is organized by biophysical characteristics (e.g. multi-pass transmembrane glycoproteins forming distinct clusters) rather than transcriptional induction timing (Figure 2). The results indicate that experimentally defined temporal expression class is not simply equivalent to global sequence-space similarity in this dataset.

### 3.3 Supervised Classification Performance
Under leakage-controlled repeated stratified cross-validation, supervised models recovered substantial information about temporal annotations (Table 5, Figure 4). The primary **Combined Equal-Block Logistic Regression** model achieved the following primary Phase 12 cross-validation estimates:
- **Accuracy:** $76.51\% \pm 10.2\%$
- **Balanced Accuracy:** $75.86\% \pm 12.3\%$
- **Macro-F1:** $0.684 \pm 0.12$
- **Matthews Correlation Coefficient (MCC):** $0.533 \pm 0.19$
- **Immediate-Early Recall:** $80.0\%$ (F1: $0.561$)
- **Early Recall:** $69.3\%$ (F1: $0.687$)
- **Late Recall:** $78.3\%$ (F1: $0.852$)

Linear SVM under equal-block fusion achieved comparable performance (Balanced Accuracy $= 73.23\% \pm 11.8\%$, Macro-F1 $= 0.671$).

### 3.4 Effect of Class Balancing on Minority-Class Recovery
Because Late proteins constitute $73.0\%$ of the dataset ($N=54$, with Early $N=15$ and IE $N=5$), a naive majority baseline achieves an Accuracy of $72.97\%$ but a Balanced Accuracy of only $33.33\%$ and an Immediate-Early Recall of $0.0\%$. Under unweighted training, supervised models collapsed onto the majority class, achieving an overall accuracy of $74.6\%$ but an Immediate-Early Recall of **$0.000$ ($0.0\%$)** and Balanced Accuracy of $41.1\%$ (Table 8). Square-root weighting provided modest improvement (**IE Recall $= 0.120$ ($12.0\%$)**). Inverse-frequency class weighting substantially improved minority-class recovery in the evaluated experiments (**IE Recall $= 0.800$ ($80.0\%$)**, Early Recall $= 69.3\%$, Late Recall $= 78.3\%$, Balanced Accuracy $= 75.86\% \pm 12.3\%$).

### 3.5 Representation Disagreement and Combined Performance
Comparing majority out-of-fold predictions across representations revealed that $59.5\%$ (44/74) of proteins achieved unanimous agreement, while **$40.5\%$ (30/74)** exhibited representation-dependent discordance, with zero proteins showing complete three-way divergence (Table 9, Figure 5). Among the 30 proteins showing disagreement between the two modality-specific predictions, the combined prediction followed ProtBERT for 19 proteins ($25.7\%$) and the physicochemical representation for 11 proteins ($14.9\%$). These disagreements indicate that the representation types provide partially distinct predictive information. Physicochemical and ProtBERT representations exhibited distinct predictive behavior for proteins with different sequence characteristics: ProtBERT models predicted enzymatic and replication-core proteins with higher consistency, while Physicochemical models exhibited distinct decision patterns for proteins with extreme charge or composition (such as `RL1`/ICP34.5 with pI 11.8). The equal-block combined representation produced stable classification performance across the predefined evaluation and sensitivity analyses without feature-count dominance.

### 3.6 Protein-Level Error Profiling
Out-of-fold error analysis across the 5 cross-validation repeats identified that $73.0\%$ (54/74) of proteins were consistently classified correctly ($\ge 80\%$ accuracy across folds), while classification difficulty was concentrated in **11 proteins** (Table 7, Figure 6).
- **High-Confidence Incorrect ($\Delta p > 0.35$):** `UL15` (terminase large subunit, $\Delta p = 0.849$), `RL1` (ICP34.5, $\Delta p = 0.771$), `US12` (ICP47, $\Delta p = 0.669$), `UL11` (myristoylated tegument, $\Delta p = 0.608$), `UL8` (helicase-primase accessory, $\Delta p = 0.424$), `UL24` (nucleolar regulator, $\Delta p = 0.361$), and `UL41` (vhs RNase, $\Delta p = 0.352$).
- **Moderate-Margin Incorrect ($0.15 \le \Delta p \le 0.35$):** `UL36` (large tegument hub, $\Delta p = 0.333$), `UL49` (VP22 tegument, $\Delta p = 0.335$), `UL52` (primase, $\Delta p = 0.274$), and `UL13` (protein kinase, $\Delta p = 0.254$).

`US12` (88 aa, true Immediate-Early) was a high-confidence computational misclassification despite its Immediate-Early annotation, illustrating that temporal expression class is not necessarily aligned with simple sequence-space similarity. `UL36` (3,139 aa, true Late, frequent prediction Early, mean margin $\approx 0.333$, tier = Moderate-Margin Incorrect) represented a computationally difficult case, potentially reflecting its extreme sequence length and complex protein architecture; this interpretation remains a hypothesis rather than a demonstrated causal mechanism.

### 3.7 Methodological Robustness and Sensitivity
1. **CV Seed Invariance:** The primary combined-model classification result was stable across the predefined independent seed set (`[7, 17, 37, 73, 97]`), yielding Balanced Accuracy of $0.754 \pm 0.017$ and Macro-F1 of $0.665$, compared to the original seed set Balanced Accuracy of $0.758 \pm 0.020$ (Table 8, Figure 7), while PCA and leave-one-protein sensitivity analyses produced broadly consistent conclusions.
2. **In-Fold PCA Reduction:** Retaining $90\%$ variance (~22 PCs) or $95\%$ variance (~30 PCs) in-fold retained $>98\%$ of full-dimensional classification performance while reducing feature dimensionality by $>97\%$.
3. **Leave-One-Protein Stability:** Systematic 74-trial holdout tests revealed that no single protein acted as a catastrophic pivot for model performance ($|\Delta \text{ Balanced Accuracy}| < 0.015$ for all Early/Late proteins).
4. **Outlier Sensitivity:** Performance changed when label-independent geometric outliers (`UL36`, `US12`, `UL41`) were excluded, indicating that representation-space extremes influence classification behavior (Balanced Accuracy adjusted to $0.782$). Because these analyses alter the evaluated protein set, they are treated as sensitivity analyses rather than evidence that protein removal improves the primary model.

---

## 4. DISCUSSION

### 4.1 Sequence Space vs. Temporal Kinetics
Our findings indicate that experimentally defined temporal expression class is not simply equivalent to global sequence-space similarity in this dataset. In unsupervised space, clustering algorithms group proteins according to biophysical characteristics—separating hydrophobic envelope glycoproteins, highly ordered icosahedral capsid subunits, and globular metabolic enzymes. Because each temporal expression phase in HSV-1 encompasses multiple functional and structural classes, unsupervised clustering does not reconstruct the temporal cascade.

### 4.2 Difficult Proteins and Biological Context
Biological contextualization provides plausible explanations for several computationally difficult proteins whose lifecycle activities or physical properties span conventional boundaries:
- **`UL41` (Virion Host Shutoff):** Annotated as Late ($\gamma_1$), `UL41` is packaged into the virion tegument and functions as an endoribonuclease immediately upon viral uncoating. Its globular ribonuclease fold computationally resembles early enzymatic proteins, which may contribute to its classification margin ($\Delta p = 0.352$).
- **`UL15` (Terminase Subunit 1):** Annotated as Late ($\gamma_2$) for viral DNA cleavage and packaging, `UL15` possesses an ATPase and nuclease core structurally homologous to early helicases, providing a plausible computational context for its misclassification ($\Delta p = 0.849$).
- **`US12` (ICP47):** Annotated as Immediate-Early, `US12` functions as a short cytosolic inhibitor (88 aa) of TAP-mediated peptide loading rather than a transcriptional transactivator. Lacking the regulatory domains of canonical IE transactivators (`ICP0`, `ICP4`, `ICP27`), `US12` represents a high-confidence misclassification ($\Delta p = 0.669$).
- **`UL36` (Large Tegument Protein VP1/2):** At 3,139 amino acids, `UL36` is an extreme sequence outlier that functions as a structural nexus in the tegument, coordinating capsid transport and nuclear docking. Its unique sequence properties and modular linkers may contribute to its moderate-margin classification difficulty ($\Delta p = 0.333$).

These biological contextualizations offer plausible insights into computational behavior but remain exploratory interpretations rather than demonstrated causal mechanisms.

### 4.3 Methodological Implications
The equal-block multi-modal representation produced stable classification performance across the predefined evaluation and sensitivity analyses, preventing feature-count dominance while maintaining descriptive diversity. Inverse-frequency class weighting substantially improved minority-class recovery across all evaluated models, highlighting the importance of explicit class balancing in viral proteomics datasets characterized by severe majority class dominance.

### 4.4 Limitations
1. **Single Viral Strain:** Analysis is restricted to HSV-1 strain 17; findings should not be assumed to generalize across all herpesviruses or other clinical strains without direct evaluation.
2. **Minority Class Size:** With only 5 IE proteins in the HSV-1 genome, individual prediction changes cause discrete $20\%$ shifts in recall.
3. **Computational Scope:** Model predictions represent statistical associations within computed feature spaces and do not constitute direct experimental validation of transcriptional mechanisms, clinical utility, or therapeutic efficacy.

---

## 5. CONCLUSION

This study provides a rigorous, leakage-controlled computational evaluation of sequence representation and temporal classification in the HSV-1 proteome. The results demonstrate that temporal expression classes are not simply equivalent to sequence-space similarity, yet supervised models recover substantial temporal information when class imbalance is explicitly addressed. Physicochemical and ProtBERT representations exhibit distinct predictive behavior, and the equal-block combined representation provides stable performance across sensitivity analyses. While biological contextualization offers plausible explanations for difficult cases, computational associations should not be conflated with demonstrated biological mechanisms. These findings establish a reproducible computational baseline for viral proteome analysis.

---

## REFERENCES

1. Albecka A, et al. (2017). Dual function of the HSV-1 UL7/UL51 complex in secondary envelopment. *J Virol* 91:e00897-17.
2. Baines JD, Roizman B. (1991). The open reading frames UL10, UL43, and UL49A of HSV-1 encode nonessential glycoproteins. *J Virol* 65:938-944.
3. Bigalke JM, Heldwein EE. (2015). Structural basis of membrane budding by the nuclear egress complex of herpesviruses. *Cell* 160:1107-1119.
4. Chou J, Roizman B. (1990). The gamma(1)34.5 gene of herpes simplex virus 1 precludes neurovirulence and facilitates growth in human cells. *J Virol* 64:1014-1020.
5. Crute JJ, et al. (1989). Herpes simplex virus 1 helicase-primase: a complex of three herpes-encoded gene products. *PNAS* 86:2186-2189.
6. Dai X, Zhou ZH. (2018). Structure of the herpes simplex virus 1 capsid with associated tegument protein complexes. *Science* 360:eaao7298.
7. Elnaggar A, et al. (2021). ProtTrans: Towards Cracking the Language of Life's Code Through Self-Supervised Deep Learning and High Performance Computing. *IEEE TPAMI* 44:7112-7127.
8. Everett RD, Maul GG. (1994). HSV-1 regulatory protein ICP0 induces the degradation of PML. *EMBO J* 13:5062-5069.
9. Heldwein EE, et al. (2006). Crystal structure of glycoprotein B from herpes simplex virus 1. *Science* 313:217-220.
10. Honess RW, Roizman B. (1974). Regulation of herpesvirus macromolecular synthesis. I. Cascade regulation of the synthesis of three groups of viral proteins. *J Virol* 14:8-19.
11. Prekeris R, et al. (1997). The UL11 protein of herpes simplex virus 1 is packaged into the tegument. *J Virol* 71:4394-4402.
12. Preston CM. (1979). Control of herpes simplex virus type 1 mRNA synthesis in cells infected in the absence of protein synthesis. *J Virol* 29:275-284.
13. Read GS, Frenkel N. (1983). Herpes simplex virus mutants defective in the virion-associated shutoff of host protein synthesis. *J Virol* 46:498-512.
14. Roizman B, et al. (2013). Herpes Simplex Viruses. In: *Fields Virology*, 6th Edition, Knipe DM, Howley PM (Eds.), Lippincott Williams & Wilkins, Philadelphia.
15. Sandbaumhüter M, et al. (2013). Cytosolic herpes simplex virus capsids not only use dynein for retrograde transport but also recruit kinesin-1 and kinesin-2. *PLoS Pathog* 9:e1003230.
16. York IA, et al. (1994). A cytosolic herpes simplex virus protein inhibits antigen presentation to CD8+ T lymphocytes. *Cell* 77:525-535.

---

## SUPPLEMENTARY MATERIAL

Supplementary Tables S1–S13, complete feature matrices, embedding extraction metadata, and reproducibility test suites are available in the project repository and indexed in `supplementary/README.md`.
