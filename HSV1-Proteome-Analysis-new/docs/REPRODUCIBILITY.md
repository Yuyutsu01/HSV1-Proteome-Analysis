# Reproducibility Framework

## Deterministic Execution and Integrity Protocol

To ensure full reproducibility across computational platforms, the pipeline implements deterministic random seeding, explicit package dependencies, cryptographic hash checks, and automated validation suites.

---

### 1. Cryptographic Hashes and Raw Data Verification

| File | Type | Expected SHA-256 Checksum |
| :--- | :--- | :--- |
| `data/raw/ncbi/NC_001806.2.gbk` | GenBank Record | `aaee9fc50cfaf8f5bf9f2b8b939f8adba061ca81d11ff982f6dbd9fef038b3a0` |
| `data/processed/unique_proteins.fasta` | Protein Sequences | Verified 74 entries, exactly matching NCBI translations |
| `data/annotations/temporal_annotations_final.csv` | Final Annotations | Frozen (74 rows, 10/10 audit checks passed) |

---

### 2. Software Environment and Dependencies

- **Python Version:** 3.11+
- **Core Computational Packages:**
  - `biopython`: Sequence parsing, CDS translation, and `ProtParam` descriptors.
  - `torch` & `transformers`: `Rostlab/prot_bert` model inference.
  - `numpy` & `pandas`: Data structures and matrix manipulation.
  - `scikit-learn`: Standard scaling, PCA, and metric evaluation.
  - `umap-learn`: Uniform Manifold Approximation and Projection.
  - `matplotlib` & `seaborn`: Visualization.

*See [`requirements.txt`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/requirements.txt) for exact environment specifications.*

---

### 3. Random Seeds and Determinism

| Phase | Operation | Random Seed | Determinism Mechanism |
| :--- | :--- | :---: | :--- |
| **Phase 07** | Physicochemical extraction | N/A | Deterministic analytical formulas |
| **Phase 08** | ProtBERT inference | `42` | `torch.no_grad()`, `model.eval()`, CPU/GPU float32 deterministic execution |
| **Phase 09** | Matrix concatenation | N/A | Direct NumPy array slicing and concatenation |
| **Phase 10** | PCA | N/A | Deterministic SVD (`svd_solver='full'`) |
| **Phase 10** | UMAP | `42` | Fixed seed (`random_state=42`) |
| **Phase 10** | t-SNE | `42` | Fixed seed (`random_state=42`, exact method) |

---

### 4. Automated Verification Test Suites

Every pipeline phase includes an automated test module in `tests/` with 100% pass verification:
- `tests/test_phase7_features.py`: 8/8 tests passed (finite values, valid ranges, AA sum = 1.0, length match).
- `tests/test_phase8_protbert.py`: 6/6 tests passed (74x1024 shape, zero NaNs, exact protein ID alignment).
- `tests/test_phase9_representations.py`: 6/6 tests passed (bijective 1:1 mapping, exact dimension preservation).
- `tests/test_phase10_exploration.py`: 6/6 tests passed (PCA variance properties, coordinate stability, figure outputs).
