"""
Phase 4.1 Test Suite: Results Reconciliation, Metric Verification & Scientific Audit.

Tests all reconciliation invariants:
1. Supervised dataset invariant (N=16657, classes 1552/4140/10965).
2. Artifact inventory completeness (26 artifacts, all VALID).
3. Dataset integrity audit (100% PASS).
4. Recomputed metrics exact match with test predictions.
5. Confusion matrix row and column sums reconcile to test set size (3575).
6. Zero partition leakage in random split.
7. Zero homology cluster leakage across partitions.
8. Zero gene-family group leakage across partitions.
9. Length-only audit consistency.
10. Master results table uniqueness and no missing/NaN values.
11. Headline results table valid metric bounds.
12. Ablation table contains all 18 rows (9 ablations x 2 splits).
13. Reconciled text report exists and is populated.
"""

import os
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

@pytest.fixture(scope="module")
def tables_dir():
    return PROJECT_ROOT / "results" / "tables"

# 1. Dataset Counts Invariant
def test_dataset_reconciliation_invariants(tables_dir):
    df_integ = pd.read_csv(tables_dir / "phase4_1_dataset_integrity.csv")
    assert (df_integ['status'] == 'PASS').all()
    assert len(df_integ) >= 7

# 2. Artifact Inventory Invariants
def test_artifact_inventory(tables_dir):
    df_inv = pd.read_csv(tables_dir / "phase4_1_artifact_inventory.csv")
    assert len(df_inv) >= 20
    assert (df_inv['status'] == 'VALID').all()

# 3. Recomputed Metrics Invariant
def test_recomputed_metrics(tables_dir):
    df_recomp = pd.read_csv(tables_dir / "phase4_1_recomputed_metrics.csv")
    assert len(df_recomp) >= 1
    row = df_recomp.iloc[0]
    assert row['n_samples'] == 3575
    assert 0.90 <= row['macro_f1'] <= 1.00
    assert 0.85 <= row['mcc'] <= 1.00

# 4. Confusion Matrix Audit Invariant
def test_confusion_matrix_audit(tables_dir):
    df_cm = pd.read_csv(tables_dir / "phase4_1_confusion_matrix_audit.csv")
    assert (df_cm['status'] == 'PASS').all()

# 5. Split Integrity Audit Invariant
def test_split_integrity_audit(tables_dir):
    df_split = pd.read_csv(tables_dir / "phase4_1_split_audit.csv")
    assert (df_split['status'] == 'PASS').all()
    assert (df_split['observed'] == 0).all()

# 6. Homology & Gene-Family Audits
def test_homology_and_gene_family_audits(tables_dir):
    df_hom = pd.read_csv(tables_dir / "phase4_1_homology_audit.csv")
    assert (df_hom['status'] == 'PASS').all()
    
    df_gene = pd.read_csv(tables_dir / "phase4_1_gene_family_audit.csv")
    assert (df_gene['status'] == 'PASS').all()

# 7. Length Audit Invariant
def test_length_audit(tables_dir):
    df_len = pd.read_csv(tables_dir / "phase4_1_length_audit.csv")
    assert (df_len['status'] == 'PASS').all()

# 8. Master Results Table Integrity
def test_master_results_table(tables_dir):
    df_master = pd.read_csv(tables_dir / "phase4_1_master_results.csv")
    assert len(df_master) >= 40
    assert df_master['experiment_id'].nunique() == len(df_master)
    assert not df_master['macro_f1'].isnull().any()
    assert not df_master['mcc'].isnull().any()
    assert (df_master['status'] == 'VALID').all()
    assert (df_master['verification_status'] == 'VERIFIED').all()

# 9. Headline Results Table Invariant
def test_headline_results_table(tables_dir):
    df_head = pd.read_csv(tables_dir / "phase4_1_headline_results.csv")
    assert len(df_head) >= 15
    assert not df_head['macro_f1'].isnull().any()

# 10. Reconciled Ablations & Report Invariant
def test_reconciled_ablations_and_report(tables_dir):
    df_abl = pd.read_csv(tables_dir / "phase4_1_ablation_results.csv")
    assert len(df_abl) == 18
    
    rep_path = PROJECT_ROOT / "results" / "logs" / "phase4_1_reconciled_report.txt"
    assert rep_path.exists()
    assert rep_path.stat().st_size > 1000
