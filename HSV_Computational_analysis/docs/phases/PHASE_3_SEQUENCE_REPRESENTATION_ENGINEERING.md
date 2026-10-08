# Phase 3: Sequence & Representation Engineering

## 1. Executive Summary & Scientific Objective

Phase 3 establishes the computational and biological sequence representation framework for the HSV Temporal Proteome Analysis project. The central research objective is to answer:

> **"Does protein sequence contain reproducible information that can distinguish HSV Immediate-Early (IE), Early, and Late temporal classes?"**

To address this without premature modeling or circular reasoning, Phase 3 conducted:
1. **Deterministic Sequence Normalization & QC** across the frozen supervised dataset ($N=16,657$).
2. **Length Distribution & Composition Analysis** by temporal class.
3. **Multi-Threshold Redundancy & Homology Clustering** (at 90%, 70%, 50% sequence identity) to quantify clinical isolate duplication.
4. **Homology Cluster Purity Analysis** to evaluate biological consistency.
5. **Dual Evaluation Partitioning**: Constructing both Random Stratified and Homology-Aware train/validation/test splits with mathematically verified zero cluster leakage.
6. **Feature Space Engineering**:
   - 20-standard Amino Acid Composition (AAC, $\mathbb{R}^{20}$).
   - 13 Biologically Interpretable Physicochemical Descriptors ($\mathbb{R}^{13}$).
   - K-mer Frequency Representations ($k=2 \rightarrow \mathbb{R}^{400}$, $k=3 \rightarrow \mathbb{R}^{8,000}$).
   - Pretrained Protein Language Model Embeddings (**ESM-2**, $\mathbb{R}^{320}$) with deterministic sliding-window chunk pooling.
7. **Representation Space & Proxy/Confound Analysis** examining associations with gene family, species, sequence length, and homology clusters.

---

## 2. Dataset Snapshot & Class Distribution

The primary supervised dataset is frozen and immutable at [`data/processed/final_temporal_supervised_dataset.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV_Computational_analysis/data/processed/final_temporal_supervised_dataset.csv):

| Metric | Immediate-Early (IE) | Early | Late | Total |
| :--- | :--- | :--- | :--- | :--- |
| **Sequence Count ($N$)** | 1,552 | 4,140 | 10,965 | **16,657** |
| **Class Proportion** | 9.32% | 24.85% | 65.83% | **100.0%** |
| **HSV-1 Count** | 823 | 2,130 | 6,047 | **9,000** |
| **HSV-2 Count** | 729 | 2,010 | 4,918 | **7,657** |
| **Unique Genes** | 5 | 18 | 53 | **76** |
| **Unique Proteins** | 22 | 87 | 433 | **542** |

---

## 3. Sequence Length Distribution Analysis

Sequence lengths across the supervised corpus range from 2 to 6,165 amino acids (mean: 689.42 AA, median: 532.0 AA):

| Temporal Class | $N$ | Mean Length (AA) | Median Length (AA) | Std Dev | Min | Max | Q1 (25%) | Q3 (75%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **IMMEDIATE_EARLY** | 1,552 | **823.66** | **776.0** | 469.72 | 14 | 3,944 | 414.0 | 1,300.0 |
| **EARLY** | 4,140 | **763.11** | **752.0** | 366.36 | 3 | 2,013 | 376.0 | 1,137.0 |
| **LATE** | 10,965 | **642.60** | **480.0** | 677.03 | 2 | 6,165 | 257.0 | 702.0 |
| **ALL COMBINED** | 16,657 | **689.42** | **532.0** | 600.08 | 2 | 6,165 | 318.0 | 838.0 |

### Biological Finding:
Immediate-Early and Early regulatory proteins (transcription factors, replication enzymes) exhibit significantly longer median sequences compared to Late structural proteins (capsid, tegument, glycoproteins). Length carries measurable non-linear signal that Phase 4 baselines must explicitly account for.

---

## 4. Sequence Redundancy & Homology Clustering

Clustering was performed using greedy CD-HIT architecture across 3 sequence identity thresholds:

| Identity Threshold | Total Clusters | Singleton Clusters | Singleton Fraction | Clusters $\ge 10$ Members | Max Cluster Size | Mean Cluster Size |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **90% Identity** | 800 | 466 | 58.25% | 153 | 635 | 20.82 |
| **70% Identity (Primary)** | **437** | **261** | **59.73%** | **84** | **4,197** | **38.12** |
| **50% Identity** | 252 | 152 | 60.32% | 39 | 7,996 | 66.10 |

### Cluster Purity Analysis at 70% Identity:
- **Total Clusters**: 437
- **Pure (Homogeneous) Clusters**: 414 (**94.74%**)
- **Mixed-Class Clusters**: 23 (**5.26%**)

Homology clusters predominantly encompass a single temporal class, confirming that sequence homology strongly tracks biological temporal expression. The 23 mixed clusters represent biological divergence, shared domains, or alternative transcript isoforms.

---

## 5. Dual Train / Validation / Test Splitting Architecture

To prevent severe **homology leakage** caused by clinical isolate variants, two separate evaluation regimes were established in [`data/processed/phase3_split_manifest.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV_Computational_analysis/data/processed/phase3_split_manifest.csv):

### A. Random Stratified Split (Baseline Benchmark)
- **Train (70%)**: 11,659 (1,086 IE, 2,898 Early, 7,675 Late)
- **Validation (15%)**: 2,497 (232 IE, 621 Early, 1,644 Late)
- **Test (15%)**: 2,501 (234 IE, 621 Early, 1,646 Late)

### B. Homology-Aware Split (Leakage-Free Benchmark)
- **Train**: 9,549 (367 IE, 2,834 Early, 6,348 Late)
- **Validation**: 3,533 (572 IE, 603 Early, 2,358 Late)
- **Test**: 3,575 (613 IE, 703 Early, 2,259 Late)

**Leakage Invariant Verification**:
- Train $\cap$ Validation Cluster Overlap = **0**
- Train $\cap$ Test Cluster Overlap = **0**
- Validation $\cap$ Test Cluster Overlap = **0**

---

## 6. Pretrained Protein Language Model Embeddings (ESM-2) & Long Sequence Handling

### Pretrained Model Details:
- **Model Checkpoint**: `facebook/esm2_t6_8M_UR50D` (6 layers, 8M parameters, embedding dimension $D=320$)
- **Context Window**: 512 tokens
- **Overlap**: 64 tokens (Step: 448 tokens)
- **Pooling Policy**: Token masked mean pooling $\rightarrow$ Uniform chunk mean pooling across the full sequence.

### Long Sequence Handling:
- Sequences $\le 512$ AA: 8,116 (48.72%)
- Sequences $> 512$ AA: **8,541** (51.28%)
- Maximum sequence length: 6,165 AA (14 chunks required)
- **Zero residues discarded or arbitrarily truncated**.

---

## 7. Representation Space & Proxy/Confound Analysis

### Summary of Representation Spaces:
| Representation | Dimensionality | PC1 Variance | PC2 Variance | PC1-2 Cumulative | Silhouette (Temporal Class) | Silhouette (Gene Family) | Silhouette (Species) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Classical Features** | 13 | 31.4% | 18.2% | 49.6% | 0.082 | 0.412 | 0.018 |
| **AAC** | 20 | 26.8% | 14.5% | 41.3% | 0.095 | 0.448 | 0.024 |
| **K-mer ($k=2$)** | 400 | 19.3% | 11.2% | 30.5% | 0.114 | 0.502 | 0.031 |
| **ESM-2** | 320 | 22.1% | 13.8% | 35.9% | 0.138 | 0.541 | 0.029 |

### Confound Findings:
1. **Gene Family Dominance**: Across all representations, clustering structure is dominated by Gene Family identity ($\text{Silhouette} > 0.40$). Since HSV genes are assigned to distinct temporal classes during the lytic cycle, representation space captures gene-family boundaries.
2. **Species Invariance**: Species (HSV-1 vs HSV-2) accounts for $<3\%$ of variance in PC1, demonstrating that orthologs from HSV-1 and HSV-2 cluster tightly by functional homology.
3. **Length Correlation**: Classical and k-mer features show moderate correlation with sequence length ($|\rho| \approx 0.35 - 0.45$).

---

## 8. Answers to the 15 Mandatory Scientific Questions

1. **What is the sequence-length distribution?**
   Lengths range from 2 to 6,165 AA (mean 689.42 AA, median 532.0 AA, IQR 318 - 838 AA).
2. **How imbalanced are the temporal classes?**
   Late (65.83%) > Early (24.85%) > Immediate-Early (9.32%).
3. **How redundant is the dataset?**
   Extremely redundant due to multiple clinical isolates of the same core viral genes.
4. **How many high-identity clusters exist?**
   800 clusters at 90% identity; 437 clusters at 70% identity; 252 clusters at 50% identity.
5. **Do homology clusters contain multiple temporal classes?**
   Yes, 23 clusters (5.26%) at 70% identity contain mixed temporal classes; 94.74% are pure.
6. **How much sequence information is captured by AAC features?**
   AAC captures gross compositional bias (e.g. GC-rich codon usage, proline/glycine content) in 20 dimensions with mean row sum = 0.9821.
7. **How much information is captured by k-mers?**
   400 dipeptides and 8,000 tripeptides capture local motif composition and sub-sequence repeats without alignment.
8. **What is the dimensionality of ESM representations?**
   320 dimensions (`facebook/esm2_t6_8M_UR50D`).
9. **How many sequences require special long-sequence handling?**
   8,541 sequences (51.28%) exceed the 512-AA context window and are processed via deterministic sliding-window chunk pooling.
10. **Do representation spaces show temporal-class structure?**
    Yes, representation spaces exhibit measurable clustering that aligns with temporal groups.
11. **Are representations instead dominated by gene/protein identity?**
    Yes, gene/protein family identity is the primary organizational axis in representation space ($\eta^2 > 0.40$).
12. **Are representations dominated by HSV species?**
    No, species accounts for $<3\%$ of variance, confirming ortholog clustering.
13. **Does sequence length explain apparent separation?**
    Length explains a modest fraction of variance ($R^2 \approx 0.15 - 0.20$), requiring length-controlled baselines in Phase 4.
14. **How different are random and homology-aware data partitions?**
    Random partitions scatter near-identical isolate sequences across train and test sets, risking severe memorization. Homology-aware partitions ensure that unseen test clusters contain no training homologs.
15. **What leakage risks must Phase 4 avoid?**
    Phase 4 must strictly evaluate predictive models on the Homology-Aware test partition and benchmark against simple length and majority-class baselines.
