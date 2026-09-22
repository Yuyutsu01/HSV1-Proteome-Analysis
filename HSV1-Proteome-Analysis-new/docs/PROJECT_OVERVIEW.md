# Project Overview

## Scientific Purpose and Framework

### 1. Scientific Objective
The primary objective of this project is to conduct a systematic, reproducible, and mathematically rigorous computational characterization and representation comparison of the Human alphaherpesvirus 1 (HSV-1) proteome.

The HSV-1 proteome is expressed in a tightly coordinated temporal cascade upon viral entry:
1. **Immediate-Early ($\alpha$ / IE):** Regulatory proteins transcribed independently of *de novo* viral protein synthesis, responsible for initiating viral transcription and subverting host innate immunity.
2. **Early ($\beta$ / E):** Enzymes and factors required for viral DNA replication and nucleotide metabolism.
3. **Late ($\gamma$ / L):** Structural virion proteins (capsid, tegument, and envelope glycoproteins) and assembly factors expressed predominantly after viral genome replication onset.

This study investigates whether these temporally distinct functional classes exhibit characteristic sequence-level signatures, structural properties, or semantic patterns that can be detected through computational representations.

---

### 2. Evaluated Modalities
The project evaluates three primary representations across the entire HSV-1 unique proteome ($N=74$):

1. **Physicochemical Sequence Descriptors ($25$-D):**
   - 5 global sequence properties: Sequence Length, Molecular Weight, Isoelectric Point (pI), Instability Index, and Aromaticity.
   - 20 canonical amino acid composition fractions.
2. **Contextual Protein Language Model Representations ($1024$-D):**
   - Dense, context-aware embeddings extracted from the pretrained `Rostlab/prot_bert` transformer architecture.
   - Sequence-length-weighted pooling of sliding-window chunks ($C=500, O=100$) to prevent sequence truncation for long viral polypeptides (e.g., UL36 at 3,139 aa).
3. **Combined Multi-Modal Representations ($1049$-D):**
   - Concatenation of the 25 physicochemical features and 1024 ProtBERT dimensions.
   - Preserved in raw unscaled format to enable controlled downstream weighting analyses.

---

### 3. Core Scientific Principles and Boundaries

> [!IMPORTANT]
> **No Presumed Superiority:**
> This study does **NOT** assume in advance that ProtBERT contextual embeddings are superior to classical physicochemical features or amino acid composition. The purpose is strictly empirical and comparative evaluation.

- **Zero Data Fabrication / Simulation:** Every result, metric, table, and figure is computed directly from biological sequences downloaded from authoritative public databases (NCBI RefSeq).
- **Strict Isolation of Representation from Classification:** Preprocessing, representation generation, dimensionality reduction, and unsupervised clustering are performed in an entirely unsupervised manner. Temporal class annotations are frozen and used exclusively for post-hoc external validation.
- **Small-Sample Safety ($N=74$):** Because the HSV-1 proteome consists of 74 unique proteins, all transformations, feature scaling, and dimensionality reduction for supervised tasks must strictly occur inside cross-validation folds to prevent data leakage.
- **Neutral Reporting:** Visual separation in low-dimensional projections (PCA, UMAP, t-SNE) is treated as exploratory hypothesis generation rather than proof of biological partitioning. Formal statistical evaluation is reserved for clustering and supervised modeling.
