# Supplementary Materials Index and Technical Details

**Manuscript:** Computational Representation, Unsupervised Clustering, and Leakage-Controlled Classification of the Herpes Simplex Virus Type 1 Proteome  
**Dataset:** Human herpesvirus 1 strain 17 (RefSeq `NC_001806.2`, $N=74$ unique canonical proteins)  
**Repository:** `https://github.com/Yuyutsu01/HSV1-Proteome-Analysis`  
**Date:** September 2026  

---

## 1. Supplementary Tables Index

| Table ID | File Name | Canonical Path | Description |
| :--- | :--- | :--- | :--- |
| **Table S1** | `temporal_annotations_final.csv` | `data/annotations/temporal_annotations_final.csv` | Complete 74-protein biological annotation table with gene names, accession numbers, length, temporal classes, late subclasses, and literature citations |
| **Table S2** | `physicochemical_features.csv` | `data/processed/physicochemical_features.csv` | Complete 25-feature physicochemical descriptor matrix for all 74 proteins without missing values |
| **Table S3** | `protbert_embedding_metadata.csv` | `data/processed/protbert_embedding_metadata.csv` | ProtBERT embedding extraction tracking, chunking counts, residue coverage, and sequence boundaries |
| **Table S4** | `phase11_clustering_metrics.csv` | `results/tables/phase11_internal_metrics.csv` | Full unsupervised clustering evaluation metrics across $K=2..10$ for KMeans, Ward, Complete, and Average linkage |
| **Table S5** | `phase11_permutation_tests.csv` | `results/tables/phase11_permutation_tests.csv` | Empirical permutation null distribution p-values for clustering ARI and NMI based on $N=1,000$ label shuffles |
| **Table S6** | `phase12_cv_results.csv` | `results/tables/phase12_cv_results.csv` | Fold-by-fold cross-validation performance across 4 representations, 5 classifiers, and 2 weighting regimes (25 folds each) |
| **Table S7** | `phase12_out_of_fold_predictions.csv` | `results/tables/phase12_out_of_fold_predictions.csv` | Complete out-of-fold validation predictions and class probability distributions (14,800 records) |
| **Table S8** | `phase13_protein_prediction_consistency.csv` | `results/tables/phase13_protein_prediction_consistency.csv` | Per-protein prediction consistency and error confidence tiers across cross-validation repeats |
| **Table S9** | `phase13_representation_disagreement.csv` | `results/tables/phase13_representation_disagreement.csv` | Modality concordance and disagreement taxonomy for all 74 proteins |
| **Table S10** | `phase14_cv_seed_sensitivity.csv` | `results/tables/phase14_cv_seed_sensitivity.csv` | Cross-validation seed stability comparison (Original seeds vs. Independent sensitivity seeds) |
| **Table S11** | `phase14_leave_one_protein_sensitivity.csv` | `results/tables/phase14_leave_one_protein_sensitivity.csv` | 74 systematic leave-one-protein ablation trials evaluating individual protein leverage on macro metrics |
| **Table S12** | `phase15_protein_biological_context.csv` | `results/tables/phase15_protein_biological_context.csv` | Comprehensive biological context table with subcellular localizations, molecular functions, and lifecycle roles |
| **Table S13** | `phase15_source_audit.csv` | `results/tables/phase15_source_audit.csv` | Traceable source literature audit linking all biological claims to peer-reviewed primary publications, PMIDs, and DOIs |

---

## 2. Supplementary Technical Notes

### Note S1: Exact Sequence Deduplication Protocol
During genome curation, three pairs of protein coding sequences located within the terminal inverted repeats ($TR_L/IR_L$ and $TR_S/IR_S$) were found to have $100\%$ identical amino acid sequences:
- `RL2` (ICP0): YP_009137074.1 and YP_009137133.1 (775 aa)
- `RL1` (ICP34.5): YP_009137073.1 and YP_009137134.1 (263 aa)
- `RS1` (ICP4): YP_009137148.1 and YP_009137149.1 (1298 aa)
To prevent cross-validation data leakage (identical sequences partitioned across training and validation folds), only the canonical first copy of each duplicate was retained in the authoritative 74-protein dataset.

### Note S2: ProtBERT Long-Sequence Chunking and Aggregation
Standard transformer self-attention mechanisms exhibit quadratic time and memory complexity $\mathcal{O}(L^2)$ with respect to sequence length $L$. For proteins exceeding the 512-token context window of `Rostlab/prot_bert`, sequences were partitioned into 500-amino-acid chunks with a 100-amino-acid sliding overlap (step $= 400$ amino acids). Let $C_k$ denote the $k$-th chunk of length $L_k$, and $E_k \in \mathbb{R}^{1024}$ denote the mean-pooled embedding of chunk $C_k$. The final aggregated embedding $E_{\text{protein}}$ is given by:
$$E_{\text{protein}} = \sum_{k=1}^M \left( \frac{L_k}{\sum_{j=1}^M L_j} \right) E_k$$
This length-weighted aggregation ensures that longer chunks contribute proportionally to the aggregate representation without token truncation.

### Note S3: Leakage-Controlled In-Fold Preprocessing Pipeline
To guarantee complete isolation between training and validation data, all feature transformations were encapsulated inside cross-validation fold loops:
1. Training data $X_{\text{train}}$ and validation data $X_{\text{val}}$ are partitioned according to StratifiedKFold;
2. Feature scalers (e.g. `StandardScaler`) compute mean $\mu_{\text{train}}$ and standard deviation $\sigma_{\text{train}}$ strictly from $X_{\text{train}}$;
3. Transformations are applied to $X_{\text{val}}$ using the precomputed parameters: $X_{\text{val, scaled}} = (X_{\text{val}} - \mu_{\text{train}}) / \sigma_{\text{train}}$;
4. Equal-block weights ($1/\sqrt{25}$ and $1/\sqrt{1024}$) are applied after standardized scaling;
5. In-fold inverse-frequency sample weights are calculated from training label frequencies $n_{c, \text{train}}$ and supplied to the classifier loss function.
Zero validation data statistics are accessible to the model during training.
