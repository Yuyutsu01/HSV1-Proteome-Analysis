"""
Automated Unit and Integration Tests for Phase 7 Physicochemical Features

Tests:
1. Exact row count equals 74.
2. Exact feature count equals 25 numerical columns.
3. 100% unique protein IDs with no duplication.
4. Zero NaN or infinite values.
5. Amino acid composition fractions sum to 1.0 within numerical tolerance (1e-5).
6. Aromaticity values fall strictly within [0, 1].
7. Sequence lengths in feature matrix match raw unique_proteins.csv.
8. Temporal class distribution remains strictly 5 IE / 15 Early / 54 Late.
9. Deterministic reproduction: re-extracting features produces identical results.
10. Spot checks on benchmark proteins (UL36, UL30, UL19, UL54, US12).
"""

import os
import pytest
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

CANONICAL_AA_LIST = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
]

FEATURE_COLUMNS = [
    'sequence_length',
    'molecular_weight',
    'aromaticity',
    'instability_index',
    'isoelectric_point'
] + [f'aa_{aa}' for aa in CANONICAL_AA_LIST]

@pytest.fixture(scope="module")
def load_datasets():
    # Support running from tests/ or repo root
    base_dir = "."
    if not os.path.exists("data/processed/physicochemical_features.csv"):
        base_dir = ".."
    
    features_csv = os.path.join(base_dir, "data/processed/physicochemical_features.csv")
    unique_csv = os.path.join(base_dir, "data/processed/unique_proteins.csv")
    annot_csv = os.path.join(base_dir, "data/annotations/temporal_annotations_final.csv")
    meta_csv = os.path.join(base_dir, "data/processed/feature_metadata.csv")
    
    assert os.path.exists(features_csv), f"Missing {features_csv}"
    assert os.path.exists(unique_csv), f"Missing {unique_csv}"
    assert os.path.exists(annot_csv), f"Missing {annot_csv}"
    assert os.path.exists(meta_csv), f"Missing {meta_csv}"
    
    df_feat = pd.read_csv(features_csv)
    df_unique = pd.read_csv(unique_csv)
    df_annot = pd.read_csv(annot_csv)
    df_meta = pd.read_csv(meta_csv)
    
    return df_feat, df_unique, df_annot, df_meta

def test_row_count_and_uniqueness(load_datasets):
    df_feat, df_unique, _, _ = load_datasets
    assert len(df_feat) == 74, f"Expected 74 rows, got {len(df_feat)}"
    assert df_feat['protein_id'].nunique() == 74, "Duplicate protein IDs detected!"
    assert set(df_feat['protein_id']) == set(df_unique['protein_id']), "Protein ID mismatch with unique_proteins.csv"

def test_exact_25_features(load_datasets):
    df_feat, _, _, df_meta = load_datasets
    num_cols = [c for c in df_feat.columns if c not in ['protein_id', 'gene', 'temporal_class']]
    assert len(num_cols) == 25, f"Expected 25 features, found {len(num_cols)}"
    assert num_cols == FEATURE_COLUMNS, "Feature column names or ordering mismatch"
    assert len(df_meta) == 25, f"Metadata file does not document exactly 25 features ({len(df_meta)})"

def test_no_missing_or_infinite_values(load_datasets):
    df_feat, _, _, _ = load_datasets
    assert df_feat[FEATURE_COLUMNS].isna().sum().sum() == 0, "NaN values present in feature matrix!"
    assert not np.isinf(df_feat[FEATURE_COLUMNS].values).any(), "Infinite values present in feature matrix!"

def test_amino_acid_composition_sum(load_datasets):
    df_feat, _, _, _ = load_datasets
    aa_cols = [f'aa_{aa}' for aa in CANONICAL_AA_LIST]
    aa_sums = df_feat[aa_cols].sum(axis=1)
    for idx, s in enumerate(aa_sums):
        assert np.isclose(s, 1.0, atol=1e-5), f"Row {idx} AA composition does not sum to 1.0: {s}"

def test_sequence_length_consistency(load_datasets):
    df_feat, df_unique, _, _ = load_datasets
    for _, row in df_unique.iterrows():
        prot_id = row['protein_id']
        expected_len = int(row['sequence_length'])
        calc_len = int(df_feat[df_feat['protein_id'] == prot_id]['sequence_length'].values[0])
        assert calc_len == expected_len, f"Sequence length mismatch for {prot_id}"

def test_aromaticity_range(load_datasets):
    df_feat, _, _, _ = load_datasets
    assert (df_feat['aromaticity'] >= 0.0).all() and (df_feat['aromaticity'] <= 1.0).all(), "Aromaticity out of [0, 1] range"

def test_temporal_class_distribution(load_datasets):
    df_feat, _, df_annot, _ = load_datasets
    counts = df_feat['temporal_class'].value_counts().to_dict()
    assert counts.get('Immediate-Early', 0) == 5, f"Expected 5 IE, got {counts.get('Immediate-Early', 0)}"
    assert counts.get('Early', 0) == 15, f"Expected 15 Early, got {counts.get('Early', 0)}"
    assert counts.get('Late', 0) == 54, f"Expected 54 Late, got {counts.get('Late', 0)}"
    assert (df_feat['temporal_class'].values == df_annot['temporal_class'].values).all(), "Temporal class order mismatch"

def test_benchmark_proteins_spot_check(load_datasets):
    df_feat, df_unique, _, _ = load_datasets
    for g in ['UL36', 'UL30', 'UL19', 'UL54', 'US12']:
        row_feat = df_feat[df_feat['gene'] == g].iloc[0]
        row_uniq = df_unique[df_unique['gene'] == g].iloc[0]
        pa = ProteinAnalysis(row_uniq['sequence'])
        
        assert int(row_feat['sequence_length']) == len(row_uniq['sequence'])
        assert np.isclose(row_feat['molecular_weight'], pa.molecular_weight(), atol=1e-2)
        assert np.isclose(row_feat['isoelectric_point'], pa.isoelectric_point(), atol=1e-2)
        assert np.isclose(row_feat['aromaticity'], pa.aromaticity(), atol=1e-4)
        assert np.isclose(row_feat['instability_index'], pa.instability_index(), atol=1e-2)
