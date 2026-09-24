# Supplementary Material: Computational Representation and Classification of the HSV-1 Proteome

This directory and referenced canonical data archives contain the complete, transparent supplementary records for the study.

## Supplementary Tables Index

| Table ID | File Name | Canonical Location | Description |
| :--- | :--- | :--- | :--- |
| **Table S1** | `temporal_annotations_final.csv` | `data/annotations/` | Complete 74-protein biological annotation table with literature citations and PMIDs |
| **Table S2** | `physicochemical_features.csv` | `data/processed/` | Complete 25-feature physicochemical descriptor matrix for all 74 proteins |
| **Table S3** | `protbert_embedding_metadata.csv` | `data/processed/` | ProtBERT embedding extraction metadata, chunking tracking, and residue counts |
| **Table S4** | `phase11_clustering_metrics.csv` | `results/tables/` | Full unsupervised clustering evaluation metrics across K=2..10 (KMeans, Ward, Avg, Complete) |
| **Table S5** | `phase11_permutation_tests.csv` | `results/tables/` | Empirical permutation null distribution p-values for clustering ARI and NMI (N=1,000) |
| **Table S6** | `phase12_cv_results.csv` | `results/tables/` | Fold-by-fold cross-validation performance across 4 representations, 5 models, and 2 weighting regimes |
| **Table S7** | `phase12_out_of_fold_predictions.csv` | `results/tables/` | Complete out-of-fold validation predictions and class probabilities (14,800 records) |
| **Table S8** | `phase13_protein_prediction_consistency.csv` | `results/tables/` | Per-protein prediction consistency and error confidence tiers across CV repeats |
| **Table S9** | `phase13_representation_disagreement.csv` | `results/tables/` | Modality concordance and disagreement taxonomy for all 74 proteins |
| **Table S10** | `phase14_cv_seed_sensitivity.csv` | `results/tables/` | Cross-validation seed stability comparison (Original vs Independent sensitivity seeds) |
| **Table S11** | `phase14_leave_one_protein_sensitivity.csv` | `results/tables/` | 74 systematic leave-one-protein ablation trials evaluating individual protein leverage |
| **Table S12** | `phase15_protein_biological_context.csv` | `results/tables/` | Comprehensive biological context table with localizations, functions, and lifecycle roles |
| **Table S13** | `phase15_source_audit.csv` | `results/tables/` | Traceable source literature audit linking all biological claims to DOIs and PMIDs |
