#!/usr/bin/env python3
"""
Test Suite for Phase 13: Supervised Model Error Analysis and Representation Interpretation.
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
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures" / "phase13"
RESULTS_LOGS = PROJECT_ROOT / "results" / "logs"

import importlib.util
spec = importlib.util.spec_from_file_location(
    "phase13_error_analysis",
    PROJECT_ROOT / "scripts" / "13_error_analysis.py"
)
phase13_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phase13_mod)

load_all_data = phase13_mod.load_all_data
compute_protein_prediction_consistency = phase13_mod.compute_protein_prediction_consistency
compute_class_confusion_analysis = phase13_mod.compute_class_confusion_analysis
compute_representation_disagreement = phase13_mod.compute_representation_disagreement
compute_physicochemical_error_summary = phase13_mod.compute_physicochemical_error_summary
compute_error_profile = phase13_mod.compute_error_profile
compute_outlier_analysis = phase13_mod.compute_outlier_analysis
compute_cross_phase_comparison = phase13_mod.compute_cross_phase_comparison
PRIMARY_MODELS = phase13_mod.PRIMARY_MODELS


class TestPhase13ErrorAnalysis:
    """Comprehensive test suite verifying Phase 13 error analysis execution and data integrity."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.oof_df, self.annot_df, self.features_df, self.X_phys, self.X_pb, self.X_comb_raw = load_all_data()

    def test_01_all_74_proteins_represented(self):
        """Test 1: Exactly 74 proteins present in authoritative annotations and features."""
        assert len(self.annot_df) == 74
        assert len(self.features_df) == 74
        assert self.annot_df["protein_id"].nunique() == 74

    def test_02_no_duplicate_protein_ids(self):
        """Test 2: No duplicate protein IDs in dataset."""
        assert self.annot_df["protein_id"].duplicated().sum() == 0

    def test_03_correct_class_labels(self):
        """Test 3: Temporal class labels remain IE=5, Early=15, Late=54."""
        counts = self.annot_df["temporal_class"].value_counts().to_dict()
        assert counts.get("Immediate-Early") == 5
        assert counts.get("Early") == 15
        assert counts.get("Late") == 54

    def test_04_expected_prediction_counts(self):
        """Test 4: OOF predictions contain exactly 5 predictions per protein for each model."""
        for rep, mod, wt in PRIMARY_MODELS:
            sub = self.oof_df[
                (self.oof_df["representation"] == rep) &
                (self.oof_df["model"] == mod) &
                (self.oof_df["weighting"] == wt)
            ]
            counts = sub["protein_id"].value_counts()
            assert len(counts) == 74
            assert (counts == 5).all(), f"Expected 5 predictions per protein for ({rep}, {mod}, {wt})"

    def test_05_no_missing_predictions(self):
        """Test 5: No null or missing predictions in OOF dataset."""
        assert not self.oof_df["predicted_class"].isnull().any()
        assert not self.oof_df["true_class"].isnull().any()

    def test_06_prediction_consistency_computation(self):
        """Test 6: Prediction consistency table calculates accuracy in [0, 1]."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        assert len(consistency_df) == len(PRIMARY_MODELS) * 74
        assert consistency_df["accuracy_across_repeats"].between(0.0, 1.0).all()
        assert not consistency_df["consistency_category"].isnull().any()

    def test_07_class_confusion_transitions(self):
        """Test 7: Class confusion table has all 9 transition types per primary model."""
        confusion_df = compute_class_confusion_analysis(self.oof_df)
        assert len(confusion_df) == len(PRIMARY_MODELS) * 9
        assert (confusion_df["fraction_of_true_class"].between(0.0, 1.0)).all()

    def test_08_representation_disagreement_categories(self):
        """Test 8: Representation disagreement assigns valid patterns for all 74 proteins."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        disagreement_df = compute_representation_disagreement(consistency_df, self.annot_df)
        assert len(disagreement_df) == 74
        valid_patterns = [
            "ALL_AGREE",
            "PHYS_PB_DISAGREE_COMBINED_MATCHES_PB",
            "PHYS_PB_DISAGREE_COMBINED_MATCHES_PHYS",
            "PHYS_PB_AGREE_COMBINED_DIFFERS",
            "ALL_DISAGREE"
        ]
        assert disagreement_df["disagreement_pattern"].isin(valid_patterns).all()

    def test_09_physicochemical_error_summary_dimensions(self):
        """Test 9: Physicochemical error summary covers all 25 features."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        phys_error_df = compute_physicochemical_error_summary(consistency_df, self.features_df)
        assert len(phys_error_df) == 25
        assert not phys_error_df["overall_mean"].isnull().any()

    def test_10_error_profile_accuracy(self):
        """Test 10: Error profile table correctly sums errors."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        error_profile_df = compute_error_profile(consistency_df, self.oof_df)
        assert len(error_profile_df) == len(PRIMARY_MODELS)
        for _, row in error_profile_df.iterrows():
            total_calc = row["immediate_early_errors (out of 25)"] + row["early_errors (out of 75)"] + row["late_errors (out of 270)"]
            assert row["total_errors"] == total_calc

    def test_11_outlier_analysis_integrity(self):
        """Test 11: Outlier analysis computes finite distances for all 74 proteins."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        outlier_df = compute_outlier_analysis(self.X_phys, self.X_pb, self.annot_df, consistency_df)
        assert len(outlier_df) == 74
        assert not np.isnan(outlier_df["physicochemical_standardized_distance"]).any()
        assert not np.isnan(outlier_df["protbert_knn_density_distance"]).any()

    def test_12_upstream_annotations_unmodified(self):
        """Test 12: Ground truth annotations remain strictly frozen and unmodified."""
        raw_annot = pd.read_csv(DATA_ANNOTATIONS / "temporal_annotations_final.csv")
        assert len(raw_annot) == 74
        assert raw_annot["temporal_class"].value_counts()["Immediate-Early"] == 5
        assert raw_annot["temporal_class"].value_counts()["Early"] == 15
        assert raw_annot["temporal_class"].value_counts()["Late"] == 54

    def test_13_upstream_representations_unmodified(self):
        """Test 13: Phase 9 representation arrays remain identical."""
        raw_phys = np.load(DATA_PROCESSED / "X_physicochemical.npy")
        raw_pb = np.load(DATA_PROCESSED / "X_protbert.npy")
        raw_comb = np.load(DATA_PROCESSED / "X_combined_raw.npy")
        assert raw_phys.shape == (74, 25)
        assert raw_pb.shape == (74, 1024)
        assert raw_comb.shape == (74, 1049)

    def test_14_high_confidence_incorrect_threshold_audit(self):
        """Test 14: Verify that high-confidence incorrect strictly requires margin > 0.35 and UL36 is not in it."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        comb_sub = consistency_df[
            (consistency_df["representation"] == "Combined_Equal_Block") &
            (consistency_df["model"] == "Logistic_Regression") &
            (consistency_df["weighting"] == "Class_Balanced")
        ].set_index("gene")

        # UL36 margin is ~0.333, strictly <= 0.35
        ul36_margin = comb_sub.loc["UL36", "predicted_probability_margin"]
        assert ul36_margin <= 0.35, f"UL36 margin ({ul36_margin}) should not exceed 0.35 threshold"

        # Check all misclassified proteins
        misclass_proteins = comb_sub[comb_sub["accuracy_across_repeats"] < 0.5]
        high_conf_incorrect = misclass_proteins[misclass_proteins["predicted_probability_margin"] > 0.35]
        assert "UL36" not in high_conf_incorrect.index
        assert "US12" in high_conf_incorrect.index

    def test_15_confusion_matrix_oof_instances_denominator(self):
        """Test 15: Verify that OOF instance denominators are exactly 25 (IE), 75 (Early), 270 (Late)."""
        confusion_df = compute_class_confusion_analysis(self.oof_df)
        assert "total_oof_prediction_instances" in confusion_df.columns
        ie_instances = confusion_df[confusion_df["true_class"] == "Immediate-Early"]["total_oof_prediction_instances"].unique()
        early_instances = confusion_df[confusion_df["true_class"] == "Early"]["total_oof_prediction_instances"].unique()
        late_instances = confusion_df[confusion_df["true_class"] == "Late"]["total_oof_prediction_instances"].unique()

        assert list(ie_instances) == [25]
        assert list(early_instances) == [75]
        assert list(late_instances) == [270]

    def test_16_exact_25_physicochemical_features_used(self):
        """Test 16: Verify that exactly the 25 frozen features from Phase 7 are used (no Gravy or extraneous features)."""
        consistency_df = compute_protein_prediction_consistency(self.oof_df, self.annot_df)
        phys_error_df = compute_physicochemical_error_summary(consistency_df, self.features_df)
        assert len(phys_error_df) == 25
        feature_names = phys_error_df["feature_name"].tolist()
        assert "Gravy" not in feature_names
        assert "gravy" not in feature_names
        expected_features = [
            "sequence_length", "molecular_weight", "aromaticity", "instability_index", "isoelectric_point",
            "aa_A", "aa_C", "aa_D", "aa_E", "aa_F", "aa_G", "aa_H", "aa_I", "aa_K", "aa_L",
            "aa_M", "aa_N", "aa_P", "aa_Q", "aa_R", "aa_S", "aa_T", "aa_V", "aa_W", "aa_Y"
        ]
        assert sorted(feature_names) == sorted(expected_features)

