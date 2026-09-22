# Research Methodology

## Computational Biology Pipeline Workflow

The research pipeline is structured as a sequential, multi-phase computational biology framework. Each phase enforces strict input/output boundaries, automated unit tests, and zero data leakage.

```
Phase 01: NCBI Data Acquisition (RefSeq NC_001806.2 GenBank)
   │
Phase 02: CDS Translation & Protein Extraction (77 records)
   │
Phase 03: Quality Control (77/77 Validated)
   │
Phase 04: Exact Sequence Deduplication (3 Inverted Repeat Copies Removed -> 74 Unique)
   │
Phase 05: Annotation Template Construction
   │
Phase 06: Biological Annotation Audit & Finalization (IE=5, E=15, L=54)
   │
   ├──────────────────────────────────────────┬──────────────────────────────────────────┐
   ▼                                          ▼                                          ▼
Phase 07: Physicochemical Features      Phase 08: ProtBERT Embeddings              [Temporal Annotations]
  (25 Features: 5 Global + 20 AA)         (1024-D Mean-Pooled Transformer)           (Isolated Ground Truth)
   │                                          │                                          │
   └──────────────────────────────────────────┴──────────────────────────────────────────┘
                                              │
Phase 09: Multi-Modal Representation Matrix Construction
  (X_physicochemical: 74x25, X_protbert: 74x1024, X_combined_raw: 74x1049)
                                              │
Phase 10: Unsupervised Exploratory Representation Analysis
  (PCA Scree/SVD, Multi-Parameter UMAP & t-SNE Manifold Exploration, Zero-Leakage Projections)
                                              │
Phase 11: Unsupervised Clustering & Validation [PLANNED]
  (K-Means K=2..10, Hierarchical Linkages, Stability Protocols, External Validation ARI/NMI)
```

---

### Core Methodological Rules

#### 1. Temporal Annotation Independence
No machine learning, representation construction, feature selection, or dimensionality reduction was performed before temporal-expression annotations were frozen and validated in Phase 6.

#### 2. Strict Unsupervised Representation Generation
- Feature extraction (Biopython `ProtParam`) and language model inferences (`Rostlab/prot_bert`) operate strictly on amino acid sequences without access to class labels.
- ProtBERT long-sequence handling uses a deterministic sliding window ($C=500, O=100, \text{step}=400$) with sequence-length-weighted mean pooling across chunk embeddings, excluding special tokens (`[CLS]`, `[SEP]`, `[PAD]`).

#### 3. Small-Sample Integrity ($N=74$)
- **Raw representations are unscaled:** `X_combined_raw.npy` preserves raw feature scales ($74 \times 1049$).
- **No global scaling before cross-validation:** In supervised phases, any StandardScaler, MinMax, or PCA transformation must be fitted strictly on training folds to prevent optimistic bias.
- **Exploratory analysis isolation:** Projections in Phase 10 (PCA, UMAP, t-SNE) are computed strictly on feature matrices. Class labels are overlaid post-hoc exclusively for neutral visual inspection.

#### 4. External Validation vs Internal Clustering
In unsupervised clustering (Phase 11), internal cluster quality metrics (Silhouette, Calinski-Harabasz, Davies-Bouldin) and stability metrics guide cluster structure assessment. External biological metrics (Adjusted Rand Index, Normalized Mutual Information, Purity) measure concordance with temporal biology without driving cluster optimization.
