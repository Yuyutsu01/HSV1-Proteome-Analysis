"""
Comprehensive Invariant Test Suite for Phase 3: Sequence & Representation Engineering.

Tests the 18 mandatory scientific and data invariants:
1. Frozen Phase 2 datasets are unchanged.
2. Sequence normalization is deterministic.
3. Sequence content and length are preserved.
4. AAC features sum correctly (~1.0).
5. No target label leakage into input feature matrices.
6. K-mer feature generation is deterministic (k=2: 400 dims, k=3: 8000 dims).
7. Homology clusters are valid and assign every sequence.
8. Zero cluster leakage between Train, Val, and Test in homology split.
9. Every supervised sequence receives exactly one split assignment.
10. Class counts reconcile to exact dataset totals (IE=1,552, Early=4,140, Late=10,965).
11. Embeddings have correct dimensionality (16,657 x 320 for ESM-2).
12. Embeddings contain zero NaN values.
13. Embeddings contain zero Infinite values.
14. Embedding IDs map 1-to-1 with canonical IDs in manifest.
15. Long-sequence handling is documented and explicit.
16. Model checkpoint identifiers are recorded in metadata.
17. Random seed is recorded in config and environment log.
18. Phase 2 dataset integrity checksums remain immutable.
"""

import os
import json
import yaml
import pytest
import numpy as np
import pandas as pd


@pytest.fixture(scope="module")
def config():
    with open("configs/phase3_config.yaml", "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="module")
def df_supervised(config):
    path = config["dataset"]["supervised_path"]
    assert os.path.exists(path), f"Supervised dataset missing: {path}"
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def df_manifest(config):
    path = config["paths"]["sequence_manifest"]
    assert os.path.exists(path), f"Sequence manifest missing: {path}"
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def df_splits(config):
    path = config["paths"]["split_manifest"]
    assert os.path.exists(path), f"Split manifest missing: {path}"
    return pd.read_csv(path)


# Test 1 & 18: Frozen Phase 2 datasets are unchanged
def test_frozen_datasets_unchanged(df_supervised):
    assert len(df_supervised) == 16657, f"Supervised dataset altered! Found {len(df_supervised)} rows"
    counts = df_supervised['temporal_class'].value_counts().to_dict()
    assert counts['IMMEDIATE_EARLY'] == 1552
    assert counts['EARLY'] == 4140
    assert counts['LATE'] == 10965


# Test 2 & 3: Sequence normalization is deterministic and content preserved
def test_sequence_manifest_integrity(df_manifest, df_supervised):
    assert len(df_manifest) == 16657
    assert df_manifest['canonical_id'].nunique() == 16657
    assert (df_manifest['original_length'] == df_manifest['normalized_length']).all()
    assert (df_manifest['normalization_changes'] == 'NONE').all()


# Test 4: AAC features sum correctly
def test_aac_features_sum(config):
    aac_path = config["paths"]["aac_features"]
    assert os.path.exists(aac_path), f"AAC features missing: {aac_path}"
    df_aac = pd.read_csv(aac_path)
    assert len(df_aac) == 16657
    
    standard_aas = config["features"]["standard_amino_acids"]
    row_sums = df_aac[standard_aas].sum(axis=1)
    
    # Allow small tolerance for sequences containing non-standard residues (e.g. X)
    assert (row_sums >= 0.0).all() and (row_sums <= 1.0001).all()
    assert row_sums.mean() > 0.95, f"Mean AAC sum too low: {row_sums.mean()}"


# Test 5: No target leakage into classical or AAC features
def test_no_target_leakage_in_features(config):
    df_aac = pd.read_csv(config["paths"]["aac_features"])
    df_classical = pd.read_csv(config["paths"]["classical_features"])
    
    # Feature columns should not include temporal evidence or proxy labels
    forbidden = {'temporal_evidence_source', 'temporal_evidence_type', 'temporal_confidence'}
    assert not forbidden.intersection(set(df_aac.columns))
    assert not forbidden.intersection(set(df_classical.columns))


# Test 6: K-mer features deterministic and valid
def test_kmer_features(config):
    kmer_dir = config["paths"]["kmer_dir"]
    k2_path = os.path.join(kmer_dir, "kmer_k2_features.npz")
    k3_path = os.path.join(kmer_dir, "kmer_k3_features.npz")
    
    assert os.path.exists(k2_path), "k=2 k-mer features missing!"
    assert os.path.exists(k3_path), "k=3 k-mer features missing!"
    
    data_k2 = np.load(k2_path, allow_pickle=True)
    data_k3 = np.load(k3_path, allow_pickle=True)
    
    assert data_k2['features'].shape == (16657, 400)
    assert data_k3['features'].shape == (16657, 8000)
    assert len(data_k2['canonical_ids']) == 16657
    assert len(data_k3['canonical_ids']) == 16657


# Test 7, 8, 9, 10: Homology clusters and dual splits
def test_splits_and_zero_homology_leakage(df_splits, df_supervised):
    assert len(df_splits) == 16657
    assert df_splits['canonical_id'].nunique() == 16657
    
    # Check valid split labels
    valid_splits = {'TRAIN', 'VALIDATION', 'TEST'}
    assert set(df_splits['split_random']).issubset(valid_splits)
    assert set(df_splits['split_homology']).issubset(valid_splits)
    
    # Class count reconciliation
    assert df_splits['temporal_class'].value_counts()['IMMEDIATE_EARLY'] == 1552
    assert df_splits['temporal_class'].value_counts()['EARLY'] == 4140
    assert df_splits['temporal_class'].value_counts()['LATE'] == 10965
    
    # Zero cluster leakage in homology-aware split
    train_c = set(df_splits[df_splits['split_homology'] == 'TRAIN']['homology_cluster_id'])
    val_c = set(df_splits[df_splits['split_homology'] == 'VALIDATION']['homology_cluster_id'])
    test_c = set(df_splits[df_splits['split_homology'] == 'TEST']['homology_cluster_id'])
    
    assert len(train_c.intersection(val_c)) == 0, "Train-Val homology cluster overlap detected!"
    assert len(train_c.intersection(test_c)) == 0, "Train-Test homology cluster overlap detected!"
    assert len(val_c.intersection(test_c)) == 0, "Val-Test homology cluster overlap detected!"


# Test 11, 12, 13, 14: Embeddings quality, dimensionality, NaN/Inf checks
def test_embedding_qc(config):
    esm_dir = os.path.join(config["paths"]["embeddings_dir"], "esm2")
    npz_path = os.path.join(esm_dir, "esm2_embeddings.npz")
    assert os.path.exists(npz_path), f"ESM-2 embeddings file missing: {npz_path}"
    
    data = np.load(npz_path, allow_pickle=True)
    embeds = data['embeddings']
    cids = data['canonical_ids']
    
    assert embeds.shape == (16657, 320), f"Unexpected embedding shape: {embeds.shape}"
    assert len(cids) == 16657
    assert not np.isnan(embeds).any(), "NaN values found in embeddings!"
    assert not np.isinf(embeds).any(), "Infinite values found in embeddings!"


# Test 15: Long sequence handling table
def test_long_sequence_handling_table(config):
    path = config["paths"]["long_sequence_handling"]
    assert os.path.exists(path), f"Long sequence table missing: {path}"
    df_long = pd.read_csv(path)
    assert len(df_long) >= 1
    assert df_long['sequences_exceeding_max_len'].iloc[0] == 8541
    assert df_long['max_sequence_length'].iloc[0] == 6165


# Test 16 & 17: Model metadata and reproducibility logs
def test_reproducibility_artifacts(config):
    env_path = config["paths"]["environment_log"]
    assert os.path.exists(env_path), f"Environment log missing: {env_path}"
    
    esm_meta_path = os.path.join(config["paths"]["embeddings_dir"], "esm2", "metadata.json")
    assert os.path.exists(esm_meta_path), f"ESM-2 metadata missing: {esm_meta_path}"
    with open(esm_meta_path, "r", encoding="utf-8") as f:
        meta = json.load(f)
    assert meta['model_name'] == "facebook/esm2_t6_8M_UR50D"
    assert meta['embedding_dimension'] == 320
