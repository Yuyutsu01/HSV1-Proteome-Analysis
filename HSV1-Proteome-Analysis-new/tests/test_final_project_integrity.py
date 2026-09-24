#!/usr/bin/env python3
"""
Test Suite for Final Project Integrity and Synthesis.
Verifies all frozen datasets, representations, model outputs, manuscript tables,
figures, reproducibility guide, and claims across all phases (1–15 and Final Synthesis).
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_ANNOTATIONS = PROJECT_ROOT / "data" / "annotations"
RESULTS_TABLES = PROJECT_ROOT / "results" / "tables"
RESULTS_FIGURES = PROJECT_ROOT / "results" / "figures"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
MANUSCRIPT_FIGS = MANUSCRIPT_DIR / "figures"
SUPPLEMENTARY_DIR = PROJECT_ROOT / "supplementary"


class TestFinalProjectIntegrity:
    """Master validation suite for overall project integrity and synthesis."""

    def test_01_dataset_cardinality_and_classes(self):
        """Verify N=74 unique proteins and exact class distribution: IE=5, Early=15, Late=54."""
        annot_df = pd.read_csv(DATA_ANNOTATIONS / "temporal_annotations_final.csv")
        assert len(annot_df) == 74
        assert annot_df["protein_id"].nunique() == 74
        counts = annot_df["temporal_class"].value_counts().to_dict()
        assert counts["Immediate-Early"] == 5
        assert counts["Early"] == 15
        assert counts["Late"] == 54

    def test_02_physicochemical_exact_25_features(self):
        """Verify exact 25 physicochemical features with zero NaNs or Infs."""
        feat_df = pd.read_csv(DATA_PROCESSED / "physicochemical_features.csv")
        assert len(feat_df) == 74
        feature_cols = [c for c in feat_df.columns if c not in ["protein_id", "gene", "temporal_class"]]
        assert len(feature_cols) == 25
        assert "Gravy" not in feature_cols
        assert "gravy" not in feature_cols
        assert not feat_df[feature_cols].isnull().any().any()
        assert not np.isinf(feat_df[feature_cols].values).any()

    def test_03_representation_matrices_dimensions(self):
        """Verify dimensions of representation matrices: Phys (74x25), PB (74x1024), Comb (74x1049)."""
        X_phys = np.load(DATA_PROCESSED / "X_physicochemical.npy")
        X_pb = np.load(DATA_PROCESSED / "X_protbert.npy")
        X_comb = np.load(DATA_PROCESSED / "X_combined_raw.npy")

        assert X_phys.shape == (74, 25)
        assert X_pb.shape == (74, 1024)
        assert X_comb.shape == (74, 1049)

        assert not np.isnan(X_phys).any()
        assert not np.isnan(X_pb).any()
        assert not np.isnan(X_comb).any()

    def test_04_ul44_annotation_preservation(self):
        """Verify UL44 is Late with Conflicting late subclass."""
        annot_df = pd.read_csv(DATA_ANNOTATIONS / "temporal_annotations_final.csv")
        ul44 = annot_df[annot_df["gene"] == "UL44"].iloc[0]
        assert ul44["temporal_class"] == "Late"
        assert ul44["late_subclass"] == "Conflicting"
        assert ul44["verification_status"] == "conflicting"

    def test_05_difficult_proteins_exact_eleven(self):
        """Verify difficult proteins list matches the exact 11 proteins from Phase 13 and 15."""
        diff_df = pd.read_csv(RESULTS_TABLES / "phase15_difficult_protein_context.csv")
        expected_difficult = {"RL1", "UL8", "UL11", "UL13", "UL15", "UL24", "UL36", "UL41", "UL49", "UL52", "US12"}
        assert set(diff_df["gene"]) == expected_difficult
        assert len(diff_df) == 11

    def test_06_ul36_moderate_margin_tier(self):
        """Verify UL36 is categorized under Moderate-Margin Error (Delta_p <= 0.35)."""
        diff_df = pd.read_csv(RESULTS_TABLES / "phase15_difficult_protein_context.csv")
        ul36 = diff_df[diff_df["gene"] == "UL36"].iloc[0]
        assert ul36["mean_decision_margin_delta_p"] <= 0.35
        assert "Moderate-Margin" in ul36["error_confidence_tier"]

    def test_07_robustness_summary_all_robust(self):
        """Verify Phase 14 robustness summary classifies all 5 observations A-E as ROBUST."""
        rob_df = pd.read_csv(RESULTS_TABLES / "phase14_robustness_summary.csv")
        assert len(rob_df) == 5
        assert (rob_df["status"] == "ROBUST").all()

    def test_08_canonical_tables_exist(self):
        """Verify all 10 canonical manuscript tables exist and are non-empty."""
        for i in range(1, 11):
            if i == 1:
                p = RESULTS_TABLES / "table1_dataset_composition.csv"
            elif i == 2:
                p = RESULTS_TABLES / "table2_physicochemical_features.csv"
            elif i == 3:
                p = RESULTS_TABLES / "table3_representation_dimensions.csv"
            elif i == 4:
                p = RESULTS_TABLES / "table4_unsupervised_clustering_summary.csv"
            elif i == 5:
                p = RESULTS_TABLES / "table5_supervised_classification_summary.csv"
            elif i == 6:
                p = RESULTS_TABLES / "table6_per_class_classification_metrics.csv"
            elif i == 7:
                p = RESULTS_TABLES / "table7_protein_level_error_analysis.csv"
            elif i == 8:
                p = RESULTS_TABLES / "table8_robustness_analysis.csv"
            elif i == 9:
                p = RESULTS_TABLES / "table9_biological_contextualization.csv"
            elif i == 10:
                p = RESULTS_TABLES / "table10_final_evidence_hierarchy.csv"

            assert p.exists(), f"Missing canonical table: {p.name}"
            df = pd.read_csv(p)
            assert len(df) > 0, f"Table {p.name} is empty"

    def test_09_all_manuscript_figures_exist(self):
        """Verify all 8 canonical manuscript figures exist and have valid file sizes."""
        for i in range(1, 9):
            pattern = f"Figure_{i}_*.png"
            matches = list(MANUSCRIPT_FIGS.glob(pattern))
            assert len(matches) > 0, f"Missing manuscript figure for index {i}"
            assert matches[0].stat().st_size > 1000, f"Figure {matches[0].name} is empty"

    def test_10_manuscript_document_integrity(self):
        """Verify manuscript document exists and contains complete required sections."""
        ms_path = MANUSCRIPT_DIR / "HSV1_proteome_analysis_final.md"
        assert ms_path.exists(), "Missing final manuscript file"
        content = ms_path.read_text(encoding="utf-8")
        assert "## ABSTRACT" in content
        assert "## 1. INTRODUCTION" in content
        assert "## 2. MATERIALS AND METHODS" in content
        assert "## 3. RESULTS" in content
        assert "## 4. DISCUSSION" in content
        assert "## 5. CONCLUSION" in content
        assert "## REFERENCES" in content
        assert len(content) > 10000

    def test_11_reproducibility_and_readme_exist(self):
        """Verify REPRODUCIBILITY.md and README.md exist and contain required documentation."""
        repro = PROJECT_ROOT / "REPRODUCIBILITY.md"
        readme = PROJECT_ROOT / "README.md"
        assert repro.exists()
        assert readme.exists()
        assert len(repro.read_text(encoding="utf-8")) > 1000
        assert len(readme.read_text(encoding="utf-8")) > 1000
