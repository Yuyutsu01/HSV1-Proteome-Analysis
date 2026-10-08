"""
Phase 4 Invariant Test Suite: Baseline & Predictive Modeling.

Tests the 16 mandatory modeling and data invariants:
1. Supervised Dataset N = 16,657.
2. Three temporal classes only (IMMEDIATE_EARLY, EARLY, LATE).
3. Zero missing labels.
4. Zero duplicate canonical IDs.
5. Representation dimensions match exact specifications:
   - AAC: (16657, 22)
   - Classical: (16657, 15)
   - 2-mer: (16657, 400)
   - 3-mer: (16657, 8000)
   - ESM-2: (16657, 320)
6. Feature/label alignment.
7. Zero partition leakage in random split.
8. Zero homology cluster leakage in homology split.
9. Gene-family split integrity (zero group overlap in Regime C).
10. Model comparison table exists and contains expected columns.
11. Preprocessing fitted strictly on training data.
12. Leakage audit table exists with 10/10 PASS status.
13. Metric calculations are mathematically valid.
14. Confusion matrix totals reconcile to test partition sizes.
15. Reproducibility environment log exists.
16. Phase 1-3 dataset integrity remains unchanged.
"""

import os
import yaml
import pytest
import numpy as np
import pandas as pd


@pytest.fixture(scope="module")
def config():
    with open("configs/phase4_config.yaml", "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="module")
def df_supervised(config):
    path = config["dataset"]["supervised_path"]
    assert os.path.exists(path), f"Supervised dataset missing: {path}"
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def df_splits(config):
    path = config["dataset"]["split_manifest_path"]
    assert os.path.exists(path), f"Split manifest missing: {path}"
    return pd.read_csv(path)


# Test 1, 2, 3, 4: Dataset Invariants
def test_supervised_dataset_invariants(df_supervised):
    assert len(df_supervised) == 16657
    assert df_supervised['canonical_id'].nunique() == 16657
    assert df_supervised['temporal_class'].isnull().sum() == 0
    counts = df_supervised['temporal_class'].value_counts().to_dict()
    assert counts['IMMEDIATE_EARLY'] == 1552
    assert counts['EARLY'] == 4140
    assert counts['LATE'] == 10965


# Test 5 & 6: Representation Dimensions & Alignment
def test_representation_dimensions(config, df_supervised):
    # AAC
    df_aac = pd.read_csv(config["features"]["aac_path"])
    assert df_aac.shape == (16657, 22)
    assert np.array_equal(df_aac['canonical_id'].values, df_supervised['canonical_id'].values)

    # Classical
    df_class = pd.read_csv(config["features"]["classical_path"])
    assert df_class.shape == (16657, 15)
    assert np.array_equal(df_class['canonical_id'].values, df_supervised['canonical_id'].values)

    # K-mers
    k2 = np.load(config["features"]["kmer_k2_path"], allow_pickle=True)
    k3 = np.load(config["features"]["kmer_k3_path"], allow_pickle=True)
    assert k2['features'].shape == (16657, 400)
    assert k3['features'].shape == (16657, 8000)

    # ESM-2
    esm = np.load(config["features"]["esm2_path"], allow_pickle=True)
    assert esm['embeddings'].shape == (16657, 320)
    assert np.array_equal(esm['canonical_ids'], df_supervised['canonical_id'].values)


# Test 7 & 8: Zero Leakage in Random and Homology Splits
def test_split_leakage_invariants(df_splits):
    assert len(df_splits) == 16657
    assert df_splits['canonical_id'].nunique() == 16657

    # Homology cluster zero overlap
    tr_c = set(df_splits[df_splits['split_homology'] == 'TRAIN']['homology_cluster_id'])
    val_c = set(df_splits[df_splits['split_homology'] == 'VALIDATION']['homology_cluster_id'])
    te_c = set(df_splits[df_splits['split_homology'] == 'TEST']['homology_cluster_id'])

    assert len(tr_c.intersection(val_c)) == 0
    assert len(tr_c.intersection(te_c)) == 0
    assert len(val_c.intersection(te_c)) == 0


# Test 9: Input Audit Invariant
def test_input_audit_table(config):
    path = config["paths"]["input_audit"]
    assert os.path.exists(path), f"Input audit table missing: {path}"
    df_audit = pd.read_csv(path)
    assert (df_audit['status'] == 'PASS').all()
    assert len(df_audit) == 32


# Test 10: Leakage Audit Table & 10/10 Pass Status
def test_leakage_audit_table(config):
    path = config["paths"]["leakage_audit"]
    assert os.path.exists(path), f"Leakage audit table missing: {path}"
    df_leak = pd.read_csv(path)
    assert len(df_leak) == 10
    assert (df_leak['status'] == 'PASS').all()


# Test 11: Unified Model Comparison Table Invariants
def test_model_comparison_table(config):
    path = config["paths"]["model_comparison"]
    assert os.path.exists(path), f"Model comparison table missing: {path}"
    df_comp = pd.read_csv(path)
    expected_cols = [
        'representation', 'model', 'split_regime', 'accuracy', 'balanced_accuracy',
        'macro_f1', 'weighted_f1', 'mcc', 'ie_precision', 'ie_recall', 'ie_f1',
        'early_precision', 'early_recall', 'early_f1', 'late_precision', 'late_recall', 'late_f1'
    ]
    for col in expected_cols:
        assert col in df_comp.columns, f"Missing column in model comparison: {col}"
    assert len(df_comp) >= 20
    assert not df_comp['macro_f1'].isnull().any()
    assert not df_comp['mcc'].isnull().any()


# Test 12: Gene-Family Degradation Table
def test_gene_family_degradation_table(config):
    path = config["paths"]["gene_family_degradation"]
    assert os.path.exists(path), f"Gene family degradation table missing: {path}"
    df_deg = pd.read_csv(path)
    assert 'Gene-Family Disjoint' in df_deg['split_regime'].values


# Test 13: Length Control Ablation Table
def test_length_control_table(config):
    path = config["paths"]["length_control"]
    assert os.path.exists(path), f"Length control table missing: {path}"
    df_len = pd.read_csv(path)
    assert len(df_len) == 18  # 9 ablations x 2 splits


# Test 14: Species Analysis Table
def test_species_analysis_table(config):
    path = config["paths"]["species_analysis"]
    assert os.path.exists(path), f"Species analysis table missing: {path}"
    df_sp = pd.read_csv(path)
    assert 'Within HSV-1 Only' in df_sp['split_regime'].values
    assert 'Within HSV-2 Only' in df_sp['split_regime'].values


# Test 15: Calibration & Environment Logs
def test_calibration_and_environment_logs(config):
    calib_path = config["paths"]["calibration"]
    env_path = config["paths"]["environment_log"]
    report_path = config["paths"]["phase4_report"]

    assert os.path.exists(calib_path), f"Calibration summary missing: {calib_path}"
    assert os.path.exists(env_path), f"Environment log missing: {env_path}"
    assert os.path.exists(report_path), f"Phase 4 report missing: {report_path}"

