#!/usr/bin/env python3
"""
Test Suite for Phase 12: Leakage-Controlled Supervised Classification.
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
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures" / "phase12"
RESULTS_LOGS = PROJECT_ROOT / "results" / "logs"

import importlib.util
spec = importlib.util.spec_from_file_location(
    "phase12_classification",
    PROJECT_ROOT / "scripts" / "12_supervised_classification.py"
)
phase12_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phase12_mod)

load_authoritative_data = phase12_mod.load_authoritative_data
transform_train_val = phase12_mod.transform_train_val
build_model = phase12_mod.build_model
compute_sample_weights = phase12_mod.compute_sample_weights
evaluate_predictions = phase12_mod.evaluate_predictions
CLASS_ORDER = phase12_mod.CLASS_ORDER
SEEDS = phase12_mod.SEEDS
REPRESENTATIONS = phase12_mod.REPRESENTATIONS
WEIGHTING_CONDITIONS = phase12_mod.WEIGHTING_CONDITIONS
MODEL_NAMES = phase12_mod.MODEL_NAMES


class TestPhase12Classification:
    """Comprehensive test suite verifying Phase 12 supervised classification integrity."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.master_df, self.X_phys, self.X_pb, self.X_comb_raw = load_authoritative_data()
        self.y_str = self.master_df["temporal_class"].values
        self.y_int = np.array([phase12_mod.LABEL_TO_INT[s] for s in self.y_str])

    def test_01_dataset_cardinality(self):
        """Test 1: Exactly 74 proteins loaded."""
        assert len(self.master_df) == 74

    def test_02_protein_id_uniqueness(self):
        """Test 2: No duplicate protein IDs."""
        assert self.master_df["protein_id"].nunique() == 74

    def test_03_class_distribution(self):
        """Test 3: Temporal classes are exactly IE=5, Early=15, Late=54."""
        counts = self.master_df["temporal_class"].value_counts().to_dict()
        assert counts.get("Immediate-Early", 0) == 5
        assert counts.get("Early", 0) == 15
        assert counts.get("Late", 0) == 54

    def test_04_representation_dimensions(self):
        """Test 4: Input representation dimensions are correct."""
        assert self.X_phys.shape == (74, 25)
        assert self.X_pb.shape == (74, 1024)
        assert self.X_comb_raw.shape == (74, 1049)

    def test_05_no_nan_or_inf(self):
        """Test 5: No NaN or Infinite values in input matrices."""
        assert not np.isnan(self.X_phys).any()
        assert not np.isinf(self.X_phys).any()
        assert not np.isnan(self.X_pb).any()
        assert not np.isinf(self.X_pb).any()
        assert not np.isnan(self.X_comb_raw).any()
        assert not np.isinf(self.X_comb_raw).any()

    def test_06_representation_label_alignment(self):
        """Test 6: Master index matches row order across matrices."""
        master_ids = self.master_df["protein_id"].tolist()
        phys_meta = pd.read_csv(DATA_PROCESSED / "X_physicochemical_metadata.csv")
        assert phys_meta["protein_id"].tolist() == master_ids

    def test_07_stratified_split_minority_class_presence(self):
        """Test 7: Every 5-fold split contains exactly 1 IE sample in validation fold."""
        from sklearn.model_selection import StratifiedKFold
        for seed in SEEDS:
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
            for fold, (train_idx, val_idx) in enumerate(skf.split(self.X_phys, self.y_int)):
                val_ie = np.sum(self.y_int[val_idx] == 0)
                train_ie = np.sum(self.y_int[train_idx] == 0)
                assert val_ie == 1, f"Seed {seed} Fold {fold} val IE count {val_ie} != 1"
                assert train_ie == 4, f"Seed {seed} Fold {fold} train IE count {train_ie} != 4"

    def test_08_in_fold_preprocessing_leakage_isolation(self):
        """Test 8: Changing validation data does NOT alter training transformation."""
        train_idx = np.arange(59)
        val_idx1 = np.arange(59, 74)
        val_idx2 = np.arange(59, 74)

        # Alter val_idx2 values
        X_phys_val_perturbed = self.X_phys[val_idx2] * 100.0

        tr1, _ = transform_train_val(
            "Physicochemical",
            self.X_phys[train_idx], self.X_pb[train_idx],
            self.X_phys[val_idx1], self.X_pb[val_idx1]
        )
        tr2, _ = transform_train_val(
            "Physicochemical",
            self.X_phys[train_idx], self.X_pb[train_idx],
            X_phys_val_perturbed, self.X_pb[val_idx2]
        )
        assert np.array_equal(tr1, tr2), "Training transformation leaked validation data"

    def test_09_in_fold_class_weights_calculation(self):
        """Test 9: In-fold class weights correctly implement inverse frequency."""
        y_train_mock = np.array([0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2])
        # N=12, K=3. Class 0: n=4 -> 12/(3*4)=1.0. Class 1: n=3 -> 12/(3*3)=1.333. Class 2: n=5 -> 12/(3*5)=0.8.
        weights = compute_sample_weights(y_train_mock)
        assert np.isclose(weights[0], 1.0)
        assert np.isclose(weights[4], 12.0 / 9.0)
        assert np.isclose(weights[7], 12.0 / 15.0)

    def test_10_equal_block_transformation_balance(self):
        """Test 10: In-fold equal block transformation scales each block variance to 1.0."""
        train_idx = np.arange(60)
        val_idx = np.arange(60, 74)
        X_tr, _ = transform_train_val(
            "Combined_Equal_Block",
            self.X_phys[train_idx], self.X_pb[train_idx],
            self.X_phys[val_idx], self.X_pb[val_idx]
        )
        var_phys = np.sum(np.var(X_tr[:, :25], axis=0))
        var_pb = np.sum(np.var(X_tr[:, 25:], axis=0))
        assert np.isclose(var_phys, 1.0, atol=1e-5)
        assert np.isclose(var_pb, 1.0, atol=1e-5)

    def test_11_metrics_calculation_validity(self):
        """Test 11: Metric calculation function computes balanced accuracy and MCC."""
        y_true = np.array([0, 0, 1, 1, 2, 2])
        y_pred = np.array([0, 1, 1, 1, 2, 2])
        metrics = evaluate_predictions(y_true, y_pred)
        assert "balanced_accuracy" in metrics
        assert "mcc" in metrics
        assert 0.0 <= metrics["accuracy"] <= 1.0
        assert 0.0 <= metrics["balanced_accuracy"] <= 1.0

    def test_12_class_order_consistency(self):
        """Test 12: Fixed class order is strictly ['Immediate-Early', 'Early', 'Late']."""
        assert CLASS_ORDER == ["Immediate-Early", "Early", "Late"]

    def test_13_no_clustering_labels_in_features(self):
        """Test 13: Supervised features contain only biochemical and embedding dimensions."""
        assert self.X_phys.shape[1] == 25
        assert self.X_pb.shape[1] == 1024
        assert self.X_comb_raw.shape[1] == 1049

    def test_14_reproducibility(self):
        """Test 14: Model instantiation with fixed seed is deterministic."""
        m1 = build_model("Random_Forest", "Class_Balanced", seed=42)
        m2 = build_model("Random_Forest", "Class_Balanced", seed=42)
        assert m1.random_state == m2.random_state == 42
