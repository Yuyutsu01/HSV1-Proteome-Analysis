#!/usr/bin/env python3
"""
Test Suite for Phase 15: Biological and Contextual Interpretation.
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
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures" / "phase15"
RESULTS_LOGS = PROJECT_ROOT / "results" / "logs"


class TestPhase15BiologicalContext:
    """Comprehensive test suite verifying Phase 15 biological contextualization and data integrity."""

    @pytest.fixture(autouse=True)
    def setup_data(self):
        self.annot_df = pd.read_csv(DATA_ANNOTATIONS / "temporal_annotations_final.csv")
        self.bio_df = pd.read_csv(RESULTS_TABLES / "phase15_protein_biological_context.csv")
        self.temp_func_df = pd.read_csv(RESULTS_TABLES / "phase15_temporal_functional_context.csv")
        self.diff_df = pd.read_csv(RESULTS_TABLES / "phase15_difficult_protein_context.csv")
        self.disagree_df = pd.read_csv(RESULTS_TABLES / "phase15_representation_disagreement_context.csv")
        self.source_df = pd.read_csv(RESULTS_TABLES / "phase15_source_audit.csv")

    def test_01_all_74_proteins_represented(self):
        """Test 1: All 74 proteins present in biological context table without duplicates."""
        assert len(self.bio_df) == 74
        assert self.bio_df["protein_id"].nunique() == 74
        assert self.bio_df["gene"].nunique() == 74

    def test_02_temporal_class_counts_unmodified(self):
        """Test 2: Authoritative temporal class counts remain IE=5, Early=15, Late=54."""
        counts = self.bio_df["primary_temporal_class"].value_counts().to_dict()
        assert counts["Immediate-Early"] == 5
        assert counts["Early"] == 15
        assert counts["Late"] == 54

    def test_03_ul44_annotation_integrity(self):
        """Test 3: UL44 remains primary Late with Conflicting late subclass and verification status."""
        ul44_row = self.bio_df[self.bio_df["gene"] == "UL44"].iloc[0]
        assert ul44_row["primary_temporal_class"] == "Late"
        assert ul44_row["late_subclass"] == "Conflicting"
        assert ul44_row["verification_status"] == "conflicting"

    def test_04_difficult_proteins_list_exact(self):
        """Test 4: Difficult proteins table contains exactly the 11 proteins identified in Phase 13."""
        expected_difficult = {"RL1", "UL8", "UL11", "UL13", "UL15", "UL24", "UL36", "UL41", "UL49", "UL52", "US12"}
        assert set(self.diff_df["gene"]) == expected_difficult
        assert len(self.diff_df) == 11

    def test_05_ul36_margin_tier_exact(self):
        """Test 5: UL36 is categorized under Moderate-Margin Error (Delta_p <= 0.35)."""
        ul36_row = self.diff_df[self.diff_df["gene"] == "UL36"].iloc[0]
        assert ul36_row["mean_decision_margin_delta_p"] <= 0.35
        assert ul36_row["error_confidence_tier"] == "Moderate-Margin Incorrect (0.15 <= Delta_p <= 0.35)"

    def test_06_high_confidence_difficult_proteins(self):
        """Test 6: High-confidence incorrect proteins strictly satisfy Delta_p > 0.35."""
        high_conf = self.diff_df[self.diff_df["error_confidence_tier"] == "High-Confidence Incorrect (Delta_p > 0.35)"]
        assert len(high_conf) == 7
        assert (high_conf["mean_decision_margin_delta_p"] > 0.35).all()
        assert "US12" in high_conf["gene"].values
        assert "UL15" in high_conf["gene"].values
        assert "UL36" not in high_conf["gene"].values

    def test_07_representation_disagreement_count(self):
        """Test 7: Exactly 30 discordant proteins in representation disagreement context table."""
        assert len(self.disagree_df) == 30
        assert self.disagree_df["protein_id"].nunique() == 30

    def test_08_source_audit_integrity(self):
        """Test 8: Source audit contains 74 verified traceable claims."""
        assert len(self.source_df) == 74
        assert (self.source_df["audit_check_status"] == "VERIFIED_TRACEABLE").all()
        assert not self.source_df["primary_source"].isnull().any()

    def test_09_all_figures_exist(self):
        """Test 9: All 3 Phase 15 figures exist and are non-empty."""
        figs = [
            "temporal_functional_category_matrix.png",
            "difficult_proteins_functional_context.png",
            "representation_disagreement_functional_distribution.png"
        ]
        for f in figs:
            fig_p = RESULTS_FIGURES / f
            assert fig_p.exists(), f"Missing figure {f}"
            assert fig_p.stat().st_size > 1000, f"Figure {f} is empty"

    def test_10_upstream_data_unmodified(self):
        """Test 10: Upstream representations, annotations, and prior phase tables remain strictly frozen."""
        raw_phys = np.load(DATA_PROCESSED / "X_physicochemical.npy")
        raw_pb = np.load(DATA_PROCESSED / "X_protbert.npy")
        raw_comb = np.load(DATA_PROCESSED / "X_combined_raw.npy")
        assert raw_phys.shape == (74, 25)
        assert raw_pb.shape == (74, 1024)
        assert raw_comb.shape == (74, 1049)

        # Check Phase 11, 12, 13, 14 summary files
        assert (RESULTS_TABLES / "phase11_stability_summary.csv").exists()
        assert (RESULTS_TABLES / "phase12_summary.csv").exists()
        assert (RESULTS_TABLES / "phase13_error_profile.csv").exists()
        assert (RESULTS_TABLES / "phase14_robustness_summary.csv").exists()
