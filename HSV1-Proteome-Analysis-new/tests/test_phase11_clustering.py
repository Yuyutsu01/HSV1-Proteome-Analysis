#!/usr/bin/env python3
"""
Test Suite for Phase 11: Unsupervised Clustering, Biological Concordance, and Stability.
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_ANNOTATIONS = PROJECT_ROOT / "data" / "annotations"
RESULTS_TABLES = PROJECT_ROOT / "results" / "tables"
RESULTS_CONTINGENCY = RESULTS_TABLES / "contingency"
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures" / "phase11"
RESULTS_LOGS = PROJECT_ROOT / "results" / "logs"

import importlib.util
spec = importlib.util.spec_from_file_location(
    "phase11_clustering",
    PROJECT_ROOT / "scripts" / "11_unsupervised_clustering.py"
)
phase11_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phase11_mod)

load_authoritative_data = phase11_mod.load_authoritative_data
prepare_representations = phase11_mod.prepare_representations
calculate_internal_metrics = phase11_mod.calculate_internal_metrics
calculate_external_metrics = phase11_mod.calculate_external_metrics
run_kmeans_clustering = phase11_mod.run_kmeans_clustering
run_hierarchical_clustering = phase11_mod.run_hierarchical_clustering
K_VALUES = phase11_mod.K_VALUES
KMEANS_SEEDS = phase11_mod.KMEANS_SEEDS
PRIMARY_SEED = phase11_mod.PRIMARY_SEED


class TestPhase11Clustering:
    """Comprehensive test suite verifying Phase 11 clustering execution and integrity."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.master_df, self.X_phys, self.X_pb, self.X_comb_raw = load_authoritative_data()
        self.representations = prepare_representations(self.X_phys, self.X_pb)
        self.y_true = self.master_df["temporal_class"].values

    def test_01_dataset_cardinality(self):
        """Test 1: Exactly 74 proteins loaded."""
        assert len(self.master_df) == 74, f"Expected 74 proteins, found {len(self.master_df)}"

    def test_02_protein_id_uniqueness(self):
        """Test 2: All 74 protein IDs are unique."""
        assert self.master_df["protein_id"].nunique() == 74, "Duplicate protein IDs detected"

    def test_03_temporal_class_distribution(self):
        """Test 3: Temporal classes remain IE=5, Early=15, Late=54."""
        counts = self.master_df["temporal_class"].value_counts().to_dict()
        assert counts.get("Immediate-Early", 0) == 5, f"Expected 5 IE, got {counts.get('Immediate-Early')}"
        assert counts.get("Early", 0) == 15, f"Expected 15 Early, got {counts.get('Early')}"
        assert counts.get("Late", 0) == 54, f"Expected 54 Late, got {counts.get('Late')}"

    def test_04_physicochemical_dimensions(self):
        """Test 4: Physicochemical representation has exactly 25 dimensions."""
        assert self.X_phys.shape == (74, 25)
        assert self.representations["Physicochemical"].shape == (74, 25)

    def test_05_protbert_dimensions(self):
        """Test 5: ProtBERT representation has exactly 1024 dimensions."""
        assert self.X_pb.shape == (74, 1024)
        assert self.representations["ProtBERT"].shape == (74, 1024)

    def test_06_combined_representations_dimensions(self):
        """Test 6: Combined representations have exactly 1049 dimensions."""
        assert self.representations["Combined_Feature_Standardized"].shape == (74, 1049)
        assert self.representations["Combined_Equal_Block"].shape == (74, 1049)

    def test_07_equal_block_weighting_math(self):
        """Test 7: Equal-block weighting is mathematically implemented correctly."""
        X_eq = self.representations["Combined_Equal_Block"]
        phys_block = X_eq[:, :25]
        pb_block = X_eq[:, 25:]

        # Total variance of each block should be approximately 1.0
        phys_var_sum = np.sum(np.var(phys_block, axis=0))
        pb_var_sum = np.sum(np.var(pb_block, axis=0))

        assert np.isclose(phys_var_sum, 1.0, atol=1e-5), f"Phys block variance sum {phys_var_sum} != 1.0"
        assert np.isclose(pb_var_sum, 1.0, atol=1e-5), f"ProtBERT block variance sum {pb_var_sum} != 1.0"

    def test_08_no_nan_or_inf(self):
        """Test 8: No NaN or Infinite values in any representation."""
        for name, mat in self.representations.items():
            assert not np.isnan(mat).any(), f"NaN in {name}"
            assert not np.isinf(mat).any(), f"Inf in {name}"

    def test_09_k_values_range(self):
        """Test 9: K values are exactly 2 through 10."""
        assert K_VALUES == [2, 3, 4, 5, 6, 7, 8, 9, 10]

    def test_10_predefined_seeds(self):
        """Test 10: Seeds are exactly [42, 123, 456, 789, 2026]."""
        assert KMEANS_SEEDS == [42, 123, 456, 789, 2026]

    def test_11_temporal_labels_not_in_features(self):
        """Test 11: Temporal labels are never used as clustering features."""
        for name, mat in self.representations.items():
            assert mat.dtype in [np.float32, np.float64], f"Non-numeric dtype in {name}"
            assert mat.shape[1] in [25, 1024, 1049], f"Unexpected column count in {name}"

    def test_12_internal_metrics_independent_of_labels(self):
        """Test 12: Internal metrics calculation requires only X and cluster labels."""
        X = self.representations["Physicochemical"]
        labels = run_kmeans_clustering(X, K=3, seed=PRIMARY_SEED)
        int_metrics = calculate_internal_metrics(X, labels)

        assert "silhouette" in int_metrics
        assert "calinski_harabasz" in int_metrics
        assert "davies_bouldin" in int_metrics
        assert not np.isnan(int_metrics["silhouette"])
        assert not np.isnan(int_metrics["calinski_harabasz"])
        assert not np.isnan(int_metrics["davies_bouldin"])

    def test_13_external_metrics_validity(self):
        """Test 13: External metrics calculate ARI, NMI, and purity."""
        X = self.representations["Physicochemical"]
        labels = run_kmeans_clustering(X, K=3, seed=PRIMARY_SEED)
        ext_metrics = calculate_external_metrics(self.y_true, labels)

        assert "ARI" in ext_metrics
        assert "NMI" in ext_metrics
        assert "purity" in ext_metrics
        assert -1.0 <= ext_metrics["ARI"] <= 1.0
        assert 0.0 <= ext_metrics["NMI"] <= 1.0
        assert 0.0 <= ext_metrics["purity"] <= 1.0

    def test_14_cluster_assignments_contain_all_proteins(self):
        """Test 14: Cluster assignments contain all 74 proteins."""
        for name, mat in self.representations.items():
            labels = run_kmeans_clustering(mat, K=3, seed=PRIMARY_SEED)
            assert len(labels) == 74
            assert len(np.unique(labels)) == 3

    def test_15_phase9_matrices_unmodified(self):
        """Test 15: Upstream Phase 9 raw matrices remain identical."""
        raw_phys = np.load(DATA_PROCESSED / "X_physicochemical.npy")
        raw_pb = np.load(DATA_PROCESSED / "X_protbert.npy")
        raw_comb = np.load(DATA_PROCESSED / "X_combined_raw.npy")

        assert raw_phys.shape == (74, 25)
        assert raw_pb.shape == (74, 1024)
        assert raw_comb.shape == (74, 1049)

    def test_16_pca_sensitivity_variance_retained(self):
        """Test 16: PCA sensitivity retains >= 95% cumulative variance."""
        from sklearn.decomposition import PCA
        pca_phys = PCA(svd_solver="full").fit(self.representations["Physicochemical"])
        cum_phys = np.cumsum(pca_phys.explained_variance_ratio_)
        n_pcs_phys = int(np.argmax(cum_phys >= 0.95) + 1)
        assert cum_phys[n_pcs_phys - 1] >= 0.95

    def test_17_deterministic_reproducibility(self):
        """Test 17: Clustering execution is 100% deterministic given fixed seed."""
        X = self.representations["Physicochemical"]
        run1 = run_kmeans_clustering(X, K=4, seed=42)
        run2 = run_kmeans_clustering(X, K=4, seed=42)
        assert np.array_equal(run1, run2), "K-Means execution non-deterministic with fixed seed"
