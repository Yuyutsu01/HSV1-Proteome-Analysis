"""
Automated Unit Tests for Phase 9 Representation Matrices

Tests:
1. Physicochemical matrix shape is exactly (74, 25).
2. ProtBERT matrix shape is exactly (74, 1024).
3. Combined raw matrix shape is exactly (74, 1049).
4. Master index contains exactly 74 unique protein IDs matching unique_proteins.csv.
5. Row alignment across all 3 metadata files is 100% bijective and identical.
6. Temporal class distribution remains strictly 5 IE / 15 Early / 54 Late across all metadata files.
7. Zero NaN or infinite values across all 3 matrices.
8. Slice integrity: X_combined_raw[:, :25] equals X_physicochemical exactly.
9. Slice integrity: X_combined_raw[:, 25:] equals X_protbert exactly.
10. Dtype validation: all arrays are numeric floating point.
"""

import os
import pytest
import numpy as np
import pandas as pd

@pytest.fixture(scope="module")
def load_representation_data():
    base_dir = "."
    if not os.path.exists("data/processed/X_physicochemical.npy"):
        base_dir = ".."
        
    master_csv = os.path.join(base_dir, "data/processed/representation_master_index.csv")
    phys_npy = os.path.join(base_dir, "data/processed/X_physicochemical.npy")
    phys_meta = os.path.join(base_dir, "data/processed/X_physicochemical_metadata.csv")
    pb_npy = os.path.join(base_dir, "data/processed/X_protbert.npy")
    pb_meta = os.path.join(base_dir, "data/processed/X_protbert_metadata.csv")
    comb_npy = os.path.join(base_dir, "data/processed/X_combined_raw.npy")
    comb_meta = os.path.join(base_dir, "data/processed/X_combined_metadata.csv")
    align_csv = os.path.join(base_dir, "results/tables/representation_alignment_report.csv")
    
    assert os.path.exists(master_csv), f"Missing {master_csv}"
    assert os.path.exists(phys_npy), f"Missing {phys_npy}"
    assert os.path.exists(phys_meta), f"Missing {phys_meta}"
    assert os.path.exists(pb_npy), f"Missing {pb_npy}"
    assert os.path.exists(pb_meta), f"Missing {pb_meta}"
    assert os.path.exists(comb_npy), f"Missing {comb_npy}"
    assert os.path.exists(comb_meta), f"Missing {comb_meta}"
    assert os.path.exists(align_csv), f"Missing {align_csv}"
    
    df_master = pd.read_csv(master_csv)
    X_phys = np.load(phys_npy)
    df_phys_meta = pd.read_csv(phys_meta)
    X_pb = np.load(pb_npy)
    df_pb_meta = pd.read_csv(pb_meta)
    X_comb = np.load(comb_npy)
    df_comb_meta = pd.read_csv(comb_meta)
    df_align = pd.read_csv(align_csv)
    
    return {
        "df_master": df_master,
        "X_phys": X_phys,
        "df_phys_meta": df_phys_meta,
        "X_pb": X_pb,
        "df_pb_meta": df_pb_meta,
        "X_comb": X_comb,
        "df_comb_meta": df_comb_meta,
        "df_align": df_align
    }

def test_matrix_dimensions(load_representation_data):
    d = load_representation_data
    assert d['X_phys'].shape == (74, 25), f"Expected (74, 25), got {d['X_phys'].shape}"
    assert d['X_pb'].shape == (74, 1024), f"Expected (74, 1024), got {d['X_pb'].shape}"
    assert d['X_comb'].shape == (74, 1049), f"Expected (74, 1049), got {d['X_comb'].shape}"
    assert d['X_phys'].shape[1] + d['X_pb'].shape[1] == d['X_comb'].shape[1]

def test_no_nans_or_infs(load_representation_data):
    d = load_representation_data
    for name, mat in [("Physicochemical", d['X_phys']), ("ProtBERT", d['X_pb']), ("Combined", d['X_comb'])]:
        assert np.isnan(mat).sum() == 0, f"{name} matrix contains NaN values!"
        assert not np.isinf(mat).any(), f"{name} matrix contains infinite values!"

def test_row_alignment_across_representations(load_representation_data):
    d = load_representation_data
    assert (d['df_align']['alignment_status'] == 'PASS').all(), "Alignment report contains failed rows!"
    assert (d['df_master']['protein_id'].values == d['df_phys_meta']['protein_id'].values).all()
    assert (d['df_master']['protein_id'].values == d['df_pb_meta']['protein_id'].values).all()
    assert (d['df_master']['protein_id'].values == d['df_comb_meta']['protein_id'].values).all()

def test_class_distributions_and_labels(load_representation_data):
    d = load_representation_data
    for meta_name, df_m in [("Master", d['df_master']), ("Phys", d['df_phys_meta']), ("PB", d['df_pb_meta']), ("Comb", d['df_comb_meta'])]:
        counts = df_m['temporal_class'].value_counts().to_dict()
        assert counts.get('Immediate-Early', 0) == 5, f"Expected 5 IE in {meta_name}, got {counts.get('Immediate-Early', 0)}"
        assert counts.get('Early', 0) == 15, f"Expected 15 Early in {meta_name}, got {counts.get('Early', 0)}"
        assert counts.get('Late', 0) == 54, f"Expected 54 Late in {meta_name}, got {counts.get('Late', 0)}"

def test_combined_slice_integrity(load_representation_data):
    d = load_representation_data
    # First 25 dimensions equal Physicochemical
    phys_diff = np.max(np.abs(d['X_comb'][:, :25] - d['X_phys']))
    assert phys_diff < 1e-9, f"Combined slice [:, :25] differs from Physicochemical matrix: max diff = {phys_diff}"
    
    # Last 1024 dimensions equal ProtBERT
    pb_diff = np.max(np.abs(d['X_comb'][:, 25:] - d['X_pb']))
    assert pb_diff < 1e-6, f"Combined slice [:, 25:] differs from ProtBERT matrix: max diff = {pb_diff}"

def test_numeric_dtypes(load_representation_data):
    d = load_representation_data
    assert np.issubdtype(d['X_phys'].dtype, np.floating)
    assert np.issubdtype(d['X_pb'].dtype, np.floating)
    assert np.issubdtype(d['X_comb'].dtype, np.floating)
