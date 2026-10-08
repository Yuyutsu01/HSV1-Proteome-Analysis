# Computational Analysis of Herpes Simplex Virus Proteomes (HSV-1 & HSV-2)
## Large-Scale Biological Curation, Sequence Representations, Homology-Aware Benchmarking & Temporal Dynamics

[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![ESM-2](https://img.shields.io/badge/Model-ESM--2%20(facebook%2Fesm2__t6__8M__UR50D)-purple.svg)](https://huggingface.co/facebook/esm2_t6_8M_UR50D)
[![Test Suite](https://img.shields.io/badge/Tests-128%20Passed-brightgreen.svg)](HSV_Computational_analysis/tests/)
[![Dataset Status](https://img.shields.io/badge/Dataset%20Freeze-N%3D16%2C657%20Supervised-success.svg)](HSV_Computational_analysis/data/processed/final_temporal_supervised_dataset.csv)

This repository contains the complete, reproducible computational biology and machine learning research platform for the systematic proteomic characterization, representation engineering, and temporal expression class modeling across **Herpes Simplex Virus Type 1 (HSV-1)** and **Herpes Simplex Virus Type 2 (HSV-2)**.

---

## 🔬 Scientific Objectives & Research Questions

The central scientific question investigated in this project is:

> *"Can sequence-derived representations (amino acid composition, physicochemical descriptors, k-mers, and contextual protein language models) predict HSV Immediate-Early (IE), Early, and Late temporal classes, and does this predictive signal remain when controlling for sequence homology, viral gene families, sequence length, and species divergence?"*

### Key Research Milestones & Findings

1. **Large-Scale Proteome Curation & Ground Truth Freeze ($N = 22,689 \rightarrow 16,657$)**:
   - Ingested and deduplicated all available NCBI/UniProt HSV-1 and HSV-2 proteome records.
   - Built a strict multi-tiered temporal annotation framework with literature-backed provenance (Tier A High-Confidence Ground Truth, Tier B Verified Extended Corpus).
   - Produced the frozen supervised benchmark dataset ($N = 16,657$; **1,552 Immediate-Early**, **4,140 Early**, **10,965 Late**).

2. **The Sequence Length Proxy Discovery**:
   - On naive random splits, sequence length alone produces a deceptively high Macro-F1 ($0.6366$) because viral gene families have characteristic lengths.
   - Under strict **Homology-Aware partitioning** (holding out entire 70% CD-HIT sequence identity clusters), length-only predictive power **completely collapses to Macro-F1 = 0.1634 – 0.2171** ($\text{MCC} = -0.1248$, below random chance). Sequence length alone possesses zero generalizable predictive capability across unseen homology clusters.

3. **Homology-Aware & Gene-Family-Disjoint Evaluation**:
   - Naive random splits produce artificially inflated scores (Macro-F1 $> 0.98$) due to identical clinical isolates spanning train and test partitions.
   - Under zero-leakage **Homology-Aware evaluation**, true generalization performance is established:
     - **ESM-2 (320-dim Contextual Embeddings)**: Macro-F1 = **0.9459** | Balanced Accuracy = **0.9656** | $\text{MCC} = \mathbf{0.9106}$
     - **3-mer Tripeptides (8,000-dim)**: Macro-F1 = **0.9954** | Balanced Accuracy = **0.9958** | $\text{MCC} = \mathbf{0.9932}$
     - **2-mer Dipeptides (400-dim)**: Macro-F1 = **0.9961** | Balanced Accuracy = **0.9974** | $\text{MCC} = \mathbf{0.9937}$
     - **AAC (20-dim Composition)**: Macro-F1 = **0.8894** | Balanced Accuracy = **0.8740** | $\text{MCC} = \mathbf{0.8253}$
     - **Classical Descriptors (13-dim Biophysical)**: Macro-F1 = **0.6136** | Balanced Accuracy = **0.6582** | $\text{MCC} = \mathbf{0.6109}$
   - Under strict **Gene-Family-Disjoint holdout** (holding out entire viral gene families), sequence composition (AAC Macro-F1 $0.7042$) and ESM-2 (Macro-F1 $0.5784$) retain substantial predictive power well above random baselines ($0.3378$), demonstrating genuine family-independent sequence signatures.

4. **Cross-Species Generalization (HSV-1 $\leftrightarrow$ HSV-2)**:
   - Models trained exclusively on HSV-1 generalize to HSV-2 with high fidelity (Macro-F1 = **0.9097**, $\text{MCC} = \mathbf{0.8716}$), demonstrating that sequence-level temporal signatures are conserved across alphaherpesvirus species.

---

## 📂 Repository Architecture

```text
HSV1-Proteome-Analysis/
├── HSV_Computational_analysis/       # Primary computational research framework
│   ├── configs/                      # Experiment YAML configuration files
│   │   ├── phase3_config.yaml
│   │   └── phase4_config.yaml
│   ├── data/                         # Complete data lifecycle (Raw -> Processed)
│   │   ├── raw/                      # Raw FastA sequences & manifests
│   │   ├── intermediate/             # Parsing, length auditing, curation checkpoints
│   │   └── processed/                # Canonical frozen datasets & feature matrices
│   │       ├── final_curated_hsv_protein_corpus.csv    (N = 22,013)
│   │       ├── final_temporal_supervised_dataset.csv   (N = 16,657)
│   │       ├── unique_high_quality_sequences.tsv       (N = 22,689)
│   │       ├── phase3_aac_features.csv                 (20-dim AAC)
│   │       ├── phase3_classical_sequence_features.csv  (13-dim Biophysical)
│   │       ├── phase3_kmer_features/                   (2-mer & 3-mer NPZ matrices)
│   │       └── embeddings/esm2/                        (ESM-2 320-dim embeddings)
│   ├── docs/phases/                  # Detailed phase specifications & scientific reports
│   │   ├── PHASE_1B_BIOLOGICAL_ELIGIBILITY.md
│   │   ├── PHASE_2_FINAL_TEMPORAL_DATASET_FREEZE.md
│   │   ├── PHASE_2_5_SEQUENCE_AUTHENTICITY_VALIDATION.md
│   │   ├── PHASE_3_SEQUENCE_REPRESENTATION_ENGINEERING.md
│   │   ├── PHASE_4_BASELINE_AND_PREDICTIVE_MODELING.md
│   │   └── PHASE_4_1_RESULTS_RECONCILIATION.md
│   ├── results/
│   │   ├── figures/                  # 300 DPI publication figures (Phases 3 & 4)
│   │   ├── logs/                     # Phase reports, environment logs, and audit trails
│   │   └── tables/                   # 60+ structured CSV result & audit tables
│   ├── scripts/                      # Modular execution scripts for Phases 1 through 4.1
│   │   ├── phase1_*.py               # Ingestion & sequence auditing
│   │   ├── phase2_*.py               # Annotation, provenance, & dataset freeze
│   │   ├── phase3/                   # Feature extraction, clustering, embeddings
│   │   ├── phase4/                   # Predictive modeling, ablations, species transfer
│   │   └── phase4_1/                 # Master reconciliation & metric verification
│   └── tests/                        # Full Pytest test suite (128 / 128 passing)
│
├── HSV_Proteome_Analysis_Workspace/  # Dedicated analytical data workspace
│   └── data/                         # Direct access to all raw, intermediate, & processed data
│
└── HSV1-Proteome-Analysis-new/       # Manuscript, table/figure legends, & citation logs
    └── manuscript/
```

---

## 📊 Summary of Benchmark Results

| Representation | Classifier | Regime | Accuracy | Bal. Acc | Macro-F1 | MCC | IE F1 | Early F1 | Late F1 |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Baseline** | Majority (Always Late) | Homology | 0.6319 | 0.3333 | 0.2581 | 0.0000 | 0.0000 | 0.0000 | 0.7744 |
| **Baseline** | Stratified Random | Homology | 0.4901 | 0.3121 | 0.3139 | 0.0005 | 0.0768 | 0.2312 | 0.6337 |
| **Length Only** | Decision Tree (depth=3) | Random | 0.7677 | 0.7101 | **0.6366** | **0.5284** | 0.5517 | 0.5222 | 0.8358 |
| **Length Only** | Decision Tree (depth=3) | Homology | 0.2929 | 0.2078 | **0.2171** | **-0.0378** | 0.0000 | 0.0000 | 0.6514 |
| **Classical (13-dim)**| Gradient Boosting | Homology | 0.7922 | 0.6582 | **0.6136** | **0.6109** | 0.3541 | 0.6394 | 0.8474 |
| **AAC (20-dim)** | Random Forest | Homology | 0.9083 | 0.8740 | **0.8894** | **0.8253** | 0.8144 | 0.9234 | 0.9304 |
| **ESM-2 (320-dim)** | Logistic Regression | Homology | 0.9494 | 0.9656 | **0.9459** | **0.9106** | 0.9385 | 0.9372 | 0.9620 |
| **ESM-2 (320-dim)** | Linear SVM | Homology | 0.9264 | 0.9537 | **0.9248** | **0.8752** | 0.9250 | 0.9023 | 0.9472 |
| **3-mer (8000-dim)**| Gradient Boosting | Homology | 0.9964 | 0.9958 | **0.9954** | **0.9932** | 0.9942 | 0.9965 | 0.9955 |
| **2-mer (400-dim)** | Logistic Regression | Homology | 0.9966 | 0.9974 | **0.9961** | **0.9937** | 0.9958 | 0.9975 | 0.9950 |

---

## ⚙️ Environment Setup & Installation

### Requirements
- **OS**: Linux / Windows / macOS
- **Python**: Python 3.10+ (tested on Python 3.11)

### Installation
```bash
# Clone the repository
git clone https://github.com/Yuyutsu01/HSV1-Proteome-Analysis.git
cd HSV1-Proteome-Analysis

# Install dependencies
pip install torch transformers scikit-learn scipy biopython pandas numpy matplotlib seaborn pyyaml pytest
```

---

## 🚀 Reproduction & Verification

To execute the verification suite and confirm all 128 mathematical and software invariants:

```bash
cd HSV_Computational_analysis

# Run complete unit and integration test suite
pytest tests/
```

To rerun Phase 4 predictive modeling and reconciliation:

```bash
# Phase 4 Modeling Pipeline
python scripts/phase4/01_input_audit.py
python scripts/phase4/02_baselines.py
python scripts/phase4/03_classical_models.py
python scripts/phase4/04_esm2_models.py
python scripts/phase4/05_homology_evaluation.py
python scripts/phase4/06_gene_family_evaluation.py
python scripts/phase4/07_length_control.py
python scripts/phase4/08_species_analysis.py
python scripts/phase4/09_error_analysis.py
python scripts/phase4/10_calibration.py
python scripts/phase4/11_statistical_comparison.py
python scripts/phase4/12_phase4_report.py

# Phase 4.1 Reconciliation & Audit Suite
python scripts/phase4_1/reconcile_all.py
```

---

## 🛡️ Scientific Claim Restrictions

In accordance with strict computational biology standards:
- **What the Results Establish**: Sequence-derived representations (k-mers, AAC, and ESM-2) contain reproducible predictive information associated with HSV Immediate-Early, Early, and Late temporal classes that persists when controlling for sequence length and homology clusters, and generalizes across species boundaries.
- **What the Results Do NOT Establish**: Sequence does not causally determine transcriptional expression timing, nor does this model establish direct promoter/enhancer regulation mechanisms. All identified feature associations represent statistical correlations within the viral proteome.

---

## 📄 License
This research codebase and all generated datasets are open source under the [MIT License](LICENSE).
