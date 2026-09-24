# Computational Representation and Classification Analysis of the HSV-1 Proteome

An end-to-end, leakage-controlled computational study evaluating physicochemical descriptors, deep protein language model embeddings (ProtBERT), and equal-block multi-modal representations across the 74 unique canonical proteins of Herpes Simplex Virus Type 1 (strain 17, RefSeq `NC_001806.2`).

## Key Findings
1. **Unsupervised Representation Geometry (Phase 11):** Across the evaluated clustering configurations (K-Means, Ward, Complete, Average across $K=2..10$), temporal-class recovery was limited and varied substantially across representations and linkage methods; no clustering configuration cleanly reproduced the three temporal classes. Representation space geometry clusters proteins primarily by biophysical characteristics rather than strict transcriptional timing.
2. **Supervised Temporal Recovery (Phase 12):** Under strictly leakage-controlled 5-fold Stratified Cross-Validation repeated across 5 seeds, the primary Combined Equal-Block Logistic Regression model achieved **Accuracy = 76.51% +- 10.2%**, **Balanced Accuracy = 75.86% +- 12.3%**, **Macro-F1 = 0.684 +- 0.12**, and **IE Recall = 80.0%** (compared to a majority-class baseline of Accuracy = 72.97%, Balanced Accuracy = 33.33%, IE Recall = 0.0%).
3. **Effect of Class Balancing on Minority-Class Recovery (Phases 12 & 14):** Inverse-frequency class weighting substantially improved minority-class recovery in the evaluated experiments (IE Recall = 80.0% with inverse-frequency weighting vs. 12.0% with square-root weighting and 0.0% unweighted).
4. **Representation Disagreement (Phase 13):** 40.5% (30/74) of proteins show representation-dependent predictive discordance between Physicochemical and ProtBERT models. The equal-block combined representation produced stable classification performance across the predefined evaluation and sensitivity analyses.
5. **Difficult Proteins and Biological Context (Phases 13 & 15):** Classification difficulty was concentrated in 11 proteins (e.g. `UL15` terminase, `RL1` ICP34.5, `US12` ICP47, `UL41` vhs, `UL36` large tegument hub), whose biological contexts span multi-phase or virion-packaged lifecycle roles.

## Dataset Summary
- **Total Unique Canonical Proteins:** 74 (77 raw records minus 3 exact duplicates)
- **Immediate-Early (alpha):** 5 proteins (`RL2`, `RS1`, `UL54`, `US1`, `US12`)
- **Early (beta):** 15 proteins
- **Late (gamma):** 54 proteins

## Directory Structure
```
├── data/
│   ├── raw/                  # Downloaded NCBI GenBank records
│   ├── processed/            # 25 features, ProtBERT embeddings, representation matrices
│   └── annotations/          # Authoritative temporal annotations (temporal_annotations_final.csv)
├── scripts/                  # Reproducible pipeline scripts (01 through 15 & final_synthesis.py)
├── tests/                    # Pytest test suite covering all phases
├── results/
│   ├── tables/               # Canonical manuscript tables and CSV artifacts
│   ├── figures/              # Publication figures across all phases
│   └── logs/                 # Detailed execution and audit logs
├── docs/                     # Methodological records and phase documentation
├── manuscript/               # Final scientific manuscript (HSV1_proteome_analysis_final.md) & figures
└── supplementary/            # Supplementary material index and documentation
```

## Validation & Testing
To verify full computational reproducibility:
```bash
pytest tests/ -v
```
