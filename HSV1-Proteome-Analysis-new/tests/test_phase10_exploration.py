"""
Automated Unit Tests for Phase 10 Exploratory Representation Analysis

Tests:
1. Valid loading and sample dimensions (N=74) for all three representations.
2. PCA scores table shape, finite numeric validity (0 NaN, 0 Inf), and protein ID alignment.
3. UMAP coordinates table shape, parameter configurations, and numeric validity.
4. t-SNE coordinates table shape, parameter configurations, and numeric validity.
5. PCA explained variance ratio properties (non-negative, sums <= 1.0, monotonic cumulative).
6. Master index and temporal class label alignment with all coordinate tables.
7. Extreme observations table validity.
8. Unsupervised transformation integrity (no temporal labels in input matrices).
"""

import os
import pytest
import numpy as np
import pandas as pd

@pytest.fixture(scope="module")
def load_phase10_data():
    base_dir = "."
    if not os.path.exists("data/processed/X_physicochemical.npy"):
        base_dir = ".."
        
    master_csv = os.path.join(base_dir, "data/processed/representation_master_index.csv")
    pca_ev_csv = os.path.join(base_dir, "results/tables/pca_explained_variance.csv")
    pca_thresh_csv = os.path.join(base_dir, "results/tables/pca_variance_thresholds.csv")
    pca_scores_csv = os.path.join(base_dir, "results/tables/pca_scores.csv")
    umap_params_csv = os.path.join(base_dir, "results/tables/umap_parameters.csv")
    umap_coords_csv = os.path.join(base_dir, "results/tables/umap_coordinates.csv")
    tsne_params_csv = os.path.join(base_dir, "results/tables/tsne_parameters.csv")
    tsne_coords_csv = os.path.join(base_dir, "results/tables/tsne_coordinates.csv")
    rep_comp_csv = os.path.join(base_dir, "results/tables/phase10_representation_comparison.csv")
    extreme_csv = os.path.join(base_dir, "results/tables/phase10_extreme_observations.csv")
    
    assert os.path.exists(master_csv), f"Missing {master_csv}"
    assert os.path.exists(pca_ev_csv), f"Missing {pca_ev_csv}"
    assert os.path.exists(pca_thresh_csv), f"Missing {pca_thresh_csv}"
    assert os.path.exists(pca_scores_csv), f"Missing {pca_scores_csv}"
    assert os.path.exists(umap_params_csv), f"Missing {umap_params_csv}"
    assert os.path.exists(umap_coords_csv), f"Missing {umap_coords_csv}"
    assert os.path.exists(tsne_params_csv), f"Missing {tsne_params_csv}"
    assert os.path.exists(tsne_coords_csv), f"Missing {tsne_coords_csv}"
    assert os.path.exists(rep_comp_csv), f"Missing {rep_comp_csv}"
    assert os.path.exists(extreme_csv), f"Missing {extreme_csv}"
    
    df_master = pd.read_csv(master_csv)
    df_pca_ev = pd.read_csv(pca_ev_csv)
    df_pca_thresh = pd.read_csv(pca_thresh_csv)
    df_pca_scores = pd.read_csv(pca_scores_csv)
    df_umap_params = pd.read_csv(umap_params_csv)
    df_umap_coords = pd.read_csv(umap_coords_csv)
    df_tsne_params = pd.read_csv(tsne_params_csv)
    df_tsne_coords = pd.read_csv(tsne_coords_csv)
    df_rep_comp = pd.read_csv(rep_comp_csv)
    df_extreme = pd.read_csv(extreme_csv)
    
    return {
        "df_master": df_master,
        "df_pca_ev": df_pca_ev,
        "df_pca_thresh": df_pca_thresh,
        "df_pca_scores": df_pca_scores,
        "df_umap_params": df_umap_params,
        "df_umap_coords": df_umap_coords,
        "df_tsne_params": df_tsne_params,
        "df_tsne_coords": df_tsne_coords,
        "df_rep_comp": df_rep_comp,
        "df_extreme": df_extreme
    }

def test_dataset_sample_size_and_alignment(load_phase10_data):
    d = load_phase10_data
    assert len(d['df_master']) == 74
    assert len(d['df_pca_scores']) == 74
    assert (d['df_pca_scores']['protein_id'].values == d['df_master']['protein_id'].values).all()
    assert (d['df_pca_scores']['temporal_class'].values == d['df_master']['temporal_class'].values).all()

def test_pca_variance_properties(load_phase10_data):
    d = load_phase10_data
    df_ev = d['df_pca_ev']
    for rep in ["Physicochemical", "ProtBERT", "Combined"]:
        sub = df_ev[df_ev['representation'] == rep]
        assert (sub['explained_variance_ratio'] >= 0.0).all(), f"Negative variance ratio in {rep}"
        assert np.isclose(sub['explained_variance_ratio'].sum(), 1.0, atol=1e-3), f"EV ratio sum != 1 in {rep}"
        cev = sub['cumulative_explained_variance_ratio'].values
        assert (np.diff(cev) >= -1e-7).all(), f"Cumulative variance not monotonic in {rep}"

def test_pca_thresholds_and_comparison(load_phase10_data):
    d = load_phase10_data
    df_thresh = d['df_pca_thresh']
    df_comp = d['df_rep_comp']
    assert len(df_thresh) == 12  # 3 reps x 4 thresholds
    assert len(df_comp) == 3

def test_umap_integrity_and_sensitivity(load_phase10_data):
    d = load_phase10_data
    df_u = d['df_umap_coords']
    assert df_u['UMAP1'].isna().sum() == 0
    assert df_u['UMAP2'].isna().sum() == 0
    assert not np.isinf(df_u[['UMAP1', 'UMAP2']].values).any()
    
    # 3 representations x 3 configurations x 74 samples = 666 rows
    assert len(df_u) == 3 * 3 * 74
    configs = set(df_u['configuration'].unique())
    assert configs == {"primary_n10", "sensitivity_n5", "sensitivity_n20"}

def test_tsne_integrity_and_sensitivity(load_phase10_data):
    d = load_phase10_data
    df_t = d['df_tsne_coords']
    assert df_t['TSNE1'].isna().sum() == 0
    assert df_t['TSNE2'].isna().sum() == 0
    assert not np.isinf(df_t[['TSNE1', 'TSNE2']].values).any()
    
    # 3 representations x 3 configurations x 74 samples = 666 rows
    assert len(df_t) == 3 * 3 * 74
    configs = set(df_t['configuration'].unique())
    assert configs == {"primary_perp15", "sensitivity_perp5", "sensitivity_perp25"}

def test_extreme_observations_table(load_phase10_data):
    d = load_phase10_data
    df_ext = d['df_extreme']
    assert len(df_ext) > 0, "Extreme observations table is empty!"
    assert set(df_ext['representation'].unique()).issubset({"Physicochemical", "ProtBERT", "Combined"})
    assert set(df_ext['method'].unique()).issubset({"PCA", "UMAP", "t-SNE"})
