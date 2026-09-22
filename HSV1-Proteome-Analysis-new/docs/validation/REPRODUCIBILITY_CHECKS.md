# Reproducibility Checks and Re-run Logs

## Verification of Deterministic Execution Across Pipeline Phases

This document records the exact scripts, test suites, and deterministic verification runs confirming computational stability across all phases.

---

### 1. Phase-by-Phase Determinism Verification

| Phase | Script Responsible | Test Suite | Determinism Metric / Check | Outcome |
| :---: | :--- | :--- | :--- | :---: |
| **Phase 01** | `scripts/01_download_ncbi.py` | `download_manifest.json` | SHA-256 hash match on raw `.gbk` | **VERIFIED** |
| **Phase 02** | `scripts/02_extract_proteins.py` | Ingestion check | Deterministic CDS sequence extraction | **VERIFIED** |
| **Phase 03** | `scripts/03_quality_control.py` | `qc_summary.json` | 77/77 proteins pass syntax rules | **VERIFIED** |
| **Phase 04** | `scripts/04_exact_deduplication.py` | `duplicates.csv` | Exact SHA-256 sequence hash matching | **VERIFIED** |
| **Phase 05** | `scripts/05_annotation_template.py` | Schema validator | Template generation consistency | **VERIFIED** |
| **Phase 06** | `scripts/06_validate_annotations.py` | 10-point audit | 10/10 automated audit checks pass | **VERIFIED** |
| **Phase 07** | `scripts/07_physicochemical_features.py` | `tests/test_phase7_features.py` | Deterministic recomputation: $\max \|\Delta\| = 0.0$ | **VERIFIED (8/8 Pass)** |
| **Phase 08** | `scripts/08_protbert_embeddings.py` | `tests/test_phase8_protbert.py` | Model deterministic mode: $\max \|\Delta\| = 0.0$ | **VERIFIED (6/6 Pass)** |
| **Phase 09** | `scripts/09_construct_representations.py` | `tests/test_phase9_representations.py` | Exact array slice equality: $\max \|\Delta\| = 0.0$ | **VERIFIED (6/6 Pass)** |
| **Phase 10** | `scripts/10_exploratory_representation_analysis.py` | `tests/test_phase10_exploration.py` | Fixed seeds (`random_state=42`) | **VERIFIED (6/6 Pass)** |

---

### 2. Test Execution Summary

All test suites were executed in the local environment and passed unconditionally:

```bash
pytest tests/test_phase7_features.py
# Result: 8 passed in 0.42s

pytest tests/test_phase8_protbert.py
# Result: 6 passed in 0.85s

pytest tests/test_phase9_representations.py
# Result: 6 passed in 0.38s

pytest tests/test_phase10_exploration.py
# Result: 6 passed in 2.14s
```

**Total Automated Tests:** **26 / 26 Passed (100.0%)**
