# Scientific Decision Record

## Methodological Rationale and Foundational Assumptions

This document maintains the chronological record of core scientific and design decisions established across Phases 1 through 10.

---

### Chronological Decision Log

| ID | Decision Title | Phase | Rationale & Scientific Justification |
| :---: | :--- | :---: | :--- |
| **D01** | **Reference Genome Selection** | Phase 1 | Human alphaherpesvirus 1 strain 17 (`NC_001806.2`) was selected as the reference proteome because it represents the most extensively curated and experimentally verified HSV-1 laboratory strain. Multi-strain clinical variants were excluded to establish a single-strain reference baseline. |
| **D02** | **Exact Deduplication Strategy** | Phase 4 | The HSV-1 inverted repeat regions produce identical diploid translations (`RL1`, `RL2`, `RS1`). Only exact sequence hash matches ($100\%$ identity) were deduplicated. Divergent gene family members were preserved to prevent information loss. |
| **D03** | **Authoritative 74-Protein Dataset** | Phase 4 | 74 single-copy unique polypeptides constitute the canonical proteome dataset, eliminating artificial diploid gene weighting in downstream machine learning. |
| **D04** | **Three Primary Temporal Classes** | Phase 6 | The viral replication cascade is partitioned into three discrete primary stages: Immediate-Early (5), Early (15), and Late (54), reflecting standard virological kinetics. |
| **D05** | **Strict 25 Physicochemical Features** | Phase 7 | While original literature text mentioned 26 features, only 25 were explicitly specified in the formal feature table. The pipeline strictly implements the 25 documented descriptors without inventing an undocumented 26th feature. |
| **D06** | **Annotation Freezing Prior to ML** | Phase 6 | Biological temporal annotations were finalized and frozen prior to representation extraction and learning to ensure complete separation between ground truth curation and computational experimentation. |
| **D07** | **ProtBERT 1024-D Representations** | Phase 8 | Pretrained `Rostlab/prot_bert` was selected to evaluate whether context-aware protein language representations encode temporal or functional viral signatures beyond classical amino acid composition. |
| **D08** | **Sliding-Window Chunking for Long Sequences** | Phase 8 | Large viral structural polypeptides exceed standard transformer contexts. A chunking scheme ($C=500, O=100, S=400$) with length-weighted mean pooling was implemented to achieve 100% sequence coverage (0% truncation). |
| **D09** | **Preservation of Raw Combined Matrices** | Phase 9 | The multi-modal matrix $X_{\text{combined\_raw}}$ was stored unscaled ($74 \times 1049$) to prevent premature scaling bias and allow cross-validation-contained transformations in supervised stages. |
| **D10** | **Exploratory Status of Manifold Projections** | Phase 10 | Dimensionality reduction techniques (PCA, UMAP, t-SNE) are treated strictly as exploratory visual tools rather than statistical proof of biological separability. |
| **D11** | **Comparative Neutrality Principle** | General | No representation is declared "superior" or "optimal" *a priori*. All representations are evaluated under identical benchmark conditions. |
