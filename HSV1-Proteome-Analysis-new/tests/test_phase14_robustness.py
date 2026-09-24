#!/usr/bin/env python3
"""
Test Suite for Phase 14: Robustness and Sensitivity Analysis.
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
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures" / "phase14"
RESULTS_LOGS = PROJECT_ROOT / "results" / "logs"

import importlib.util
spec = importlib.util.spec_from_file_location(
    "phase14_robustness",
    PROJECT_ROOT / "scripts" / "14_robustness_analysis.py"
)
phase14_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phase14_mod)

load_authoritative_data = phase14_mod.load_authoritative_data
compute_sample_weights = phase14_mod.compute_sample_weights
ORIGINAL_SEEDS = phase14_mod.ORIGINAL_SEEDS
INDEPENDENT_SEEDS = phase14_mod.INDEPENDENT_SEEDS


class TestPhase14Robustness:
    """Comprehensive test suite verifying Phase 14 robustness analysis execution and integrity."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.annot_df, self.X_phys, self.X_pb, self.X_comb_raw, self.feat_df = load_authoritative_data()

    def test_01_authoritative_data_cardinality(self):
        """Test 1: Authoritative dataset remains exactly 74 proteins with IE=5, Early=15, Late=54."""
        assert len(self.annot_df) == 74
        assert self.annot_df["protein_id"].nunique() == 74
        counts = self.annot_df["temporal_class"].value_counts().to_dict()
        assert counts["Immediate-Early"] == 5
        assert counts["Early"] == 15
        assert counts["Late"] == 54

    def test_02_upstream_matrices_unmodified(self):
        """Test 2: Upstream representation arrays from Phase 9 remain unmodified."""
        raw_phys = np.load(DATA_PROCESSED / "X_physicochemical.npy")
        raw_pb = np.load(DATA_PROCESSED / "X_protbert.npy")
        raw_comb = np.load(DATA_PROCESSED / "X_combined_raw.npy")
        assert raw_phys.shape == (74, 25)
        assert raw_pb.shape == (74, 1024)
        assert raw_comb.shape == (74, 1049)

    def test_03_cv_seed_sensitivity_table(self):
        """Test 3: CV seed sensitivity table has 10 rows (2 seed sets x 5 models) with no NaNs."""
        p = RESULTS_TABLES / "phase14_cv_seed_sensitivity.csv"
        assert p.exists(), "phase14_cv_seed_sensitivity.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 10
        assert set(df["seed_set"]) == {"Original_Seeds", "Independent_Sensitivity_Seeds"}
        assert not df.isnull().any().any()
        assert df["balanced_accuracy_mean"].between(0.0, 1.0).all()

    def test_04_class_weight_sensitivity_table(self):
        """Test 4: Class weight sensitivity table covers 3 weighting schemes across 5 models (15 rows)."""
        p = RESULTS_TABLES / "phase14_class_weight_sensitivity.csv"
        assert p.exists(), "phase14_class_weight_sensitivity.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 15
        assert set(df["weighting_scheme"]) == {"Unweighted", "Inverse_Frequency", "Square_Root"}
        assert not df.isnull().any().any()

    def test_05_combined_representation_sensitivity_table(self):
        """Test 5: Combined representation table compares Equal-Block vs Feature-Standardized (4 rows)."""
        p = RESULTS_TABLES / "phase14_combined_representation_sensitivity.csv"
        assert p.exists(), "phase14_combined_representation_sensitivity.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 4
        assert set(df["representation_strategy"]) == {"Combined_Feature_Standardized", "Combined_Equal_Block"}
        assert not df.isnull().any().any()

    def test_06_pca_supervised_sensitivity_table(self):
        """Test 6: In-fold PCA sensitivity table has 8 rows with valid retained dimension counts."""
        p = RESULTS_TABLES / "phase14_pca_supervised_sensitivity.csv"
        assert p.exists(), "phase14_pca_supervised_sensitivity.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 8
        assert set(df["pca_variance_threshold"]) == {"Full_Dimensional", "PCA_90pct", "PCA_95pct", "PCA_99pct"}
        assert not df.isnull().any().any()

    def test_07_leave_one_protein_sensitivity_table(self):
        """Test 7: Leave-one-protein-out table contains exactly 74 rows (one per protein)."""
        p = RESULTS_TABLES / "phase14_leave_one_protein_sensitivity.csv"
        assert p.exists(), "phase14_leave_one_protein_sensitivity.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 74
        assert df["excluded_protein_id"].nunique() == 74
        assert not df.isnull().any().any()
        assert df["remaining_sample_size"].unique() == [73]

    def test_08_outlier_sensitivity_table(self):
        """Test 8: Outlier sensitivity table evaluates all 5 predefined cohorts."""
        p = RESULTS_TABLES / "phase14_outlier_sensitivity.csv"
        assert p.exists(), "phase14_outlier_sensitivity.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 5
        assert not df.isnull().any().any()

    def test_09_robustness_summary_synthesis(self):
        """Test 9: Robustness summary classifies all 5 observations A-E."""
        p = RESULTS_TABLES / "phase14_robustness_summary.csv"
        assert p.exists(), "phase14_robustness_summary.csv missing"
        df = pd.read_csv(p)
        assert len(df) == 5
        expected_obs = {"Observation_A", "Observation_B", "Observation_C", "Observation_D", "Observation_E"}
        assert set(df["observation_id"]) == expected_obs
        assert df["status"].isin(["ROBUST", "SENSITIVITY-DEPENDENT", "INCONCLUSIVE"]).all()

    def test_10_all_figures_generated(self):
        """Test 10: All 6 publication figures exist with non-zero file sizes."""
        fig_names = [
            "cv_seed_sensitivity.png",
            "class_weight_sensitivity.png",
            "combined_representation_sensitivity.png",
            "pca_supervised_sensitivity.png",
            "leave_one_protein_influence.png",
            "overall_robustness_summary.png"
        ]
        for name in fig_names:
            fig_p = RESULTS_FIGURES / name
            assert fig_p.exists(), f"Missing figure {name}"
            assert fig_p.stat().st_size > 1000, f"Figure {name} is empty or corrupted"

    def test_11_sample_weights_normalization(self):
        """Test 11: Square-root and inverse-frequency weights strictly preserve mean training weight."""
        dummy_y = np.array([0, 0, 1, 1, 1, 2, 2, 2, 2, 2]) # 2 IE, 3 Early, 5 Late
        w_sqrt = compute_sample_weights(dummy_y, "Square_Root")
        assert np.isclose(np.mean(w_sqrt), 1.0)
        w_inv = compute_sample_weights(dummy_y, "Inverse_Frequency")
        # In Inverse_Frequency, mean weight is also 1.0 (since sum(w_i) = n_classes * (N/n_classes) = N)
        assert np.isclose(np.mean(w_inv), 1.0)
