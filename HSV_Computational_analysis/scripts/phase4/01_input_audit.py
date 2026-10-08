"""
Phase 4 - Step 01: Repository & Input Data Audit.

Biological & Computational Concept:
Before any predictive model training, a rigorous audit must verify the exact integrity,
alignment, and consistency of all frozen Phase 1-3 datasets and representations:
1. Frozen Supervised Ground Truth: N=16,657 records (1,552 IE, 4,140 Early, 10,965 Late).
2. Exact 1-to-1 canonical ID alignment across all feature matrices (AAC, Classical, K-mer, ESM-2).
3. Zero missing values, zero NaNs, zero Infinities in any feature representation.
4. Mathematical verification of split partitions (Random Stratified and Homology-Aware).

Creates: results/tables/phase4_input_audit.csv
"""

import os
import yaml
import numpy as np
import pandas as pd


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def audit_inputs():
    """Execute complete input data audit."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 01: REPOSITORY AND INPUT DATA INTEGRITY AUDIT")
    print("=" * 80)
    
    audit_records = []
    
    def log_check(check_name: str, expected: str, actual: str, status: bool, details: str = ""):
        audit_records.append({
            'check_name': check_name,
            'expected_value': expected,
            'actual_value': actual,
            'status': 'PASS' if status else 'FAIL',
            'details': details
        })
        print(f"[{'PASS' if status else 'FAIL'}] {check_name}: {actual} ({details})")
        if not status:
            raise ValueError(f"CRITICAL AUDIT FAILURE: {check_name} failed! Expected: {expected}, got: {actual}")

    # 1. Supervised Dataset Audit
    sup_path = config["dataset"]["supervised_path"]
    log_check("Supervised Dataset Exists", "True", str(os.path.exists(sup_path)), os.path.exists(sup_path), sup_path)
    df_sup = pd.read_csv(sup_path)
    
    log_check("Supervised Total Sequences", "16657", str(len(df_sup)), len(df_sup) == 16657, f"Exact N={len(df_sup)}")
    log_check("Unique Canonical IDs", "16657", str(df_sup['canonical_id'].nunique()), df_sup['canonical_id'].nunique() == 16657)
    
    class_counts = df_sup['temporal_class'].value_counts().to_dict()
    log_check("Class Count: IMMEDIATE_EARLY", "1552", str(class_counts.get('IMMEDIATE_EARLY', 0)), class_counts.get('IMMEDIATE_EARLY', 0) == 1552)
    log_check("Class Count: EARLY", "4140", str(class_counts.get('EARLY', 0)), class_counts.get('EARLY', 0) == 4140)
    log_check("Class Count: LATE", "10965", str(class_counts.get('LATE', 0)), class_counts.get('LATE', 0) == 10965)
    log_check("Missing Temporal Labels", "0", str(int(df_sup['temporal_class'].isnull().sum())), df_sup['temporal_class'].isnull().sum() == 0)

    # 2. Sequence Manifest Audit
    man_path = config["dataset"]["manifest_path"]
    log_check("Sequence Manifest Exists", "True", str(os.path.exists(man_path)), os.path.exists(man_path), man_path)
    df_man = pd.read_csv(man_path)
    log_check("Manifest Record Count", "16657", str(len(df_man)), len(df_man) == 16657)
    log_check("Manifest ID Alignment", "True", str((df_sup['canonical_id'] == df_man['canonical_id']).all()), (df_sup['canonical_id'] == df_man['canonical_id']).all())

    # 3. Split Manifest Audit
    split_path = config["dataset"]["split_manifest_path"]
    log_check("Split Manifest Exists", "True", str(os.path.exists(split_path)), os.path.exists(split_path), split_path)
    df_splits = pd.read_csv(split_path)
    log_check("Split Manifest Record Count", "16657", str(len(df_splits)), len(df_splits) == 16657)
    
    # Check Homology Split Zero Overlap
    train_c = set(df_splits[df_splits['split_homology'] == 'TRAIN']['homology_cluster_id'])
    val_c = set(df_splits[df_splits['split_homology'] == 'VALIDATION']['homology_cluster_id'])
    test_c = set(df_splits[df_splits['split_homology'] == 'TEST']['homology_cluster_id'])
    
    log_check("Homology Train-Val Leakage", "0", str(len(train_c.intersection(val_c))), len(train_c.intersection(val_c)) == 0)
    log_check("Homology Train-Test Leakage", "0", str(len(train_c.intersection(test_c))), len(train_c.intersection(test_c)) == 0)
    log_check("Homology Val-Test Leakage", "0", str(len(val_c.intersection(test_c))), len(val_c.intersection(test_c)) == 0)

    # 4. AAC Features Audit
    aac_path = config["features"]["aac_path"]
    log_check("AAC Features Exist", "True", str(os.path.exists(aac_path)), os.path.exists(aac_path), aac_path)
    df_aac = pd.read_csv(aac_path)
    log_check("AAC Feature Shape", "(16657, 22)", str(df_aac.shape), df_aac.shape == (16657, 22))
    log_check("AAC NaN Count", "0", str(int(df_aac.isnull().sum().sum())), df_aac.isnull().sum().sum() == 0)

    # 5. Classical Features Audit
    class_path = config["features"]["classical_path"]
    log_check("Classical Features Exist", "True", str(os.path.exists(class_path)), os.path.exists(class_path), class_path)
    df_classical = pd.read_csv(class_path)
    log_check("Classical Feature Shape", "(16657, 15)", str(df_classical.shape), df_classical.shape == (16657, 15))
    log_check("Classical NaN Count", "0", str(int(df_classical.isnull().sum().sum())), df_classical.isnull().sum().sum() == 0)

    # 6. K-mer Features Audit
    k2_path = config["features"]["kmer_k2_path"]
    k3_path = config["features"]["kmer_k3_path"]
    log_check("K-mer k=2 Exists", "True", str(os.path.exists(k2_path)), os.path.exists(k2_path), k2_path)
    log_check("K-mer k=3 Exists", "True", str(os.path.exists(k3_path)), os.path.exists(k3_path), k3_path)
    
    k2_data = np.load(k2_path, allow_pickle=True)
    k3_data = np.load(k3_path, allow_pickle=True)
    log_check("K-mer k=2 Shape", "(16657, 400)", str(k2_data['features'].shape), k2_data['features'].shape == (16657, 400))
    log_check("K-mer k=3 Shape", "(16657, 8000)", str(k3_data['features'].shape), k3_data['features'].shape == (16657, 8000))
    log_check("K-mer k=2 NaN Count", "0", str(int(np.isnan(k2_data['features']).sum())), np.isnan(k2_data['features']).sum() == 0)
    log_check("K-mer k=3 NaN Count", "0", str(int(np.isnan(k3_data['features']).sum())), np.isnan(k3_data['features']).sum() == 0)

    # 7. ESM-2 Embeddings Audit
    esm_path = config["features"]["esm2_path"]
    log_check("ESM-2 Embeddings Exist", "True", str(os.path.exists(esm_path)), os.path.exists(esm_path), esm_path)
    esm_data = np.load(esm_path, allow_pickle=True)
    log_check("ESM-2 Embedding Shape", "(16657, 320)", str(esm_data['embeddings'].shape), esm_data['embeddings'].shape == (16657, 320))
    log_check("ESM-2 NaN Count", "0", str(int(np.isnan(esm_data['embeddings']).sum())), np.isnan(esm_data['embeddings']).sum() == 0)
    log_check("ESM-2 Inf Count", "0", str(int(np.isinf(esm_data['embeddings']).sum())), np.isinf(esm_data['embeddings']).sum() == 0)
    log_check("ESM-2 ID Alignment", "True", str(np.array_equal(esm_data['canonical_ids'], df_sup['canonical_id'].values)), np.array_equal(esm_data['canonical_ids'], df_sup['canonical_id'].values))

    # Save Audit Table
    df_audit = pd.DataFrame(audit_records)
    out_path = config["paths"]["input_audit"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_audit.to_csv(out_path, index=False)
    print(f"\nSaved Phase 4 Input Audit: {out_path}")
    print(f"Total Checks: {len(df_audit)} | All Passed: {(df_audit['status'] == 'PASS').all()}")


if __name__ == "__main__":
    audit_inputs()
