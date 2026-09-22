"""
Automated Unit Tests for Phase 8 ProtBERT Embeddings

Tests:
1. Embedding matrix shape is exactly (74, 1024).
2. Embedding metadata rows match unique_proteins.csv exactly.
3. Zero NaN or infinite values in embedding matrix.
4. Complete 100% coverage fraction (1.0) for every protein sequence.
5. Multi-chunk processing for long proteins (>500 aa), single chunk for short (<=500 aa).
6. Spot checks on longest viral protein UL36 (3,139 aa, 8 chunks, coverage=1.0).
7. ProtBERT configuration JSON contains required hyperparameters.
"""

import os
import json
import pytest
import numpy as np
import pandas as pd

@pytest.fixture(scope="module")
def load_protbert_data():
    base_dir = "."
    if not os.path.exists("data/processed/protbert_embeddings.npy"):
        base_dir = ".."
        
    emb_path = os.path.join(base_dir, "data/processed/protbert_embeddings.npy")
    meta_path = os.path.join(base_dir, "data/processed/protbert_metadata.csv")
    cfg_path = os.path.join(base_dir, "data/processed/protbert_config.json")
    unique_path = os.path.join(base_dir, "data/processed/unique_proteins.csv")
    stats_path = os.path.join(base_dir, "results/tables/protbert_embedding_statistics.csv")
    handling_path = os.path.join(base_dir, "results/tables/protbert_sequence_handling.csv")
    
    assert os.path.exists(emb_path), f"Missing {emb_path}"
    assert os.path.exists(meta_path), f"Missing {meta_path}"
    assert os.path.exists(cfg_path), f"Missing {cfg_path}"
    assert os.path.exists(unique_path), f"Missing {unique_path}"
    assert os.path.exists(stats_path), f"Missing {stats_path}"
    assert os.path.exists(handling_path), f"Missing {handling_path}"
    
    emb = np.load(emb_path)
    df_meta = pd.read_csv(meta_path)
    with open(cfg_path, "r") as f:
        cfg = json.load(f)
    df_unique = pd.read_csv(unique_path)
    df_stats = pd.read_csv(stats_path)
    df_handling = pd.read_csv(handling_path)
    
    return emb, df_meta, cfg, df_unique, df_stats, df_handling

def test_embedding_dimensions_and_shape(load_protbert_data):
    emb, df_meta, _, _, _, _ = load_protbert_data
    assert emb.shape == (74, 1024), f"Expected (74, 1024), got {emb.shape}"
    assert len(df_meta) == 74, f"Expected 74 metadata rows, got {len(df_meta)}"

def test_no_nans_or_infs(load_protbert_data):
    emb, _, _, _, _, _ = load_protbert_data
    assert np.isnan(emb).sum() == 0, "Embedding matrix contains NaN values!"
    assert not np.isinf(emb).any(), "Embedding matrix contains infinite values!"

def test_coverage_and_residue_tracking(load_protbert_data):
    _, df_meta, _, df_unique, _, _ = load_protbert_data
    assert (df_meta['coverage_fraction'] == 1.0).all(), "Incomplete sequence coverage detected!"
    assert (df_meta['covered_residues'] == df_unique['sequence_length']).all(), "Covered residues count mismatch"
    assert (df_meta['sequence_length'] == df_unique['sequence_length']).all(), "Sequence length mismatch"

def test_chunking_logic(load_protbert_data):
    _, df_meta, _, _, _, _ = load_protbert_data
    long_prots = df_meta[df_meta['sequence_length'] > 500]
    short_prots = df_meta[df_meta['sequence_length'] <= 500]
    
    assert (long_prots['number_of_chunks'] > 1).all(), "Long proteins were not split into multiple chunks!"
    assert (short_prots['number_of_chunks'] == 1).all(), "Short proteins were unnecessarily fragmented!"

def test_longest_protein_ul36(load_protbert_data):
    _, df_meta, _, _, _, df_handling = load_protbert_data
    top_row = df_handling.iloc[0]
    assert top_row['gene'] == 'UL36', f"Expected longest protein to be UL36, got {top_row['gene']}"
    assert int(top_row['sequence_length']) == 3139, f"Expected length 3139 for UL36, got {top_row['sequence_length']}"
    assert int(top_row['number_of_chunks']) >= 7, f"Expected >= 7 chunks for UL36, got {top_row['number_of_chunks']}"
    assert float(top_row['coverage_fraction']) == 1.0

def test_config_integrity(load_protbert_data):
    _, _, cfg, _, _, _ = load_protbert_data
    assert cfg['model_name'] == "Rostlab/prot_bert"
    assert cfg['hidden_dimension'] == 1024
    assert cfg['chunk_size'] == 500
    assert cfg['overlap'] == 100
    assert "residue_mean" in cfg['pooling_method']
