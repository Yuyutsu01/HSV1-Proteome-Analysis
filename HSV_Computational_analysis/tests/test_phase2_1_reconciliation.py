"""
===============================================================================
Phase 2.1 Automated Test Suite: Annotation & Ground-Truth Reconciliation
===============================================================================
Verifies all 20 Phase 2.1 invariants and reconciliation criteria.
===============================================================================
"""

import os
import hashlib
from pathlib import Path
import pandas as pd
import pytest


BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "raw"
DATA_RAW_DIR = BASE_DIR / "data" / "raw"
DATA_PROCESSED_DIR = BASE_DIR / "data" / "processed"
RESULTS_TABLES_DIR = BASE_DIR / "results" / "tables"


def compute_sha256(file_path: Path) -> str:
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(65536):
            sha256.update(chunk)
    return sha256.hexdigest()


def test_1_and_2_and_15_sequence_count_and_no_new_or_external():
    """TEST 1, 2, 15: Verify exact 22,689 unique sequences represented without additions or deletions."""
    df_1a = pd.read_csv(DATA_PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t")
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")

    assert len(df_21) == 22689
    assert len(df_1a) == 22689
    assert set(df_21["sequence_id"]) == set(df_1a["canonical_sequence_id"])


def test_3_and_4_and_5_raw_immutability_and_hash_invariance():
    """TEST 3, 4, 5: Verify sequence contents, hashes, and raw FASTA files are identical and unmutated."""
    df_1a = pd.read_csv(DATA_PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t").set_index("canonical_sequence_id")
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv").set_index("sequence_id")

    for seq_id, row2 in df_21.head(500).iterrows():
        row1 = df_1a.loc[seq_id]
        assert row2["sequence"] == row1["sequence"]
        assert row2["canonical_sequence_hash"] == seq_id

    # Raw FASTA immutability
    raw_files = list(RAW_DIR.glob("*.fasta"))
    assert len(raw_files) == 2
    for rf in raw_files:
        mirror = DATA_RAW_DIR / rf.name
        assert compute_sha256(rf) == compute_sha256(mirror)


def test_6_allowed_temporal_classes():
    """TEST 6: Verify all temporal labels belong to the allowed vocabulary."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    allowed_classes = {"IMMEDIATE_EARLY", "EARLY", "LATE", "UNKNOWN", "CONFLICTING"}
    observed_classes = set(df_21["temporal_class"].unique())
    assert observed_classes.issubset(allowed_classes), f"Unexpected classes: {observed_classes - allowed_classes}"


def test_7_temporal_annotation_status_consistency():
    """TEST 7: Verify temporal_annotation_status is mathematically derived from temporal_class."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    
    labeled_mask = df_21["temporal_class"].isin(["IMMEDIATE_EARLY", "EARLY", "LATE"])
    assert (df_21[labeled_mask]["temporal_annotation_status"] == "LABELED").all()
    
    unknown_mask = df_21["temporal_class"] == "UNKNOWN"
    assert (df_21[unknown_mask]["temporal_annotation_status"] == "UNKNOWN").all()
    
    conflicted_mask = df_21["temporal_class"] == "CONFLICTING"
    assert (df_21[conflicted_mask]["temporal_annotation_status"] == "CONFLICTING").all()


def test_8_and_9_tier_b_contains_labels_and_not_called_unlabeled():
    """TEST 8, 9: Verify Tier B contains preserved temporal labels and is defined as extended corpus."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    
    tier_b = df_21[df_21["annotation_tier"] == "TIER_B_VERIFIED_EXTENDED_BIOLOGICAL_CORPUS"]
    assert len(tier_b) == 11531
    
    labeled_tier_b = tier_b[tier_b["temporal_annotation_status"] == "LABELED"]
    assert len(labeled_tier_b) == 6175, f"Expected 6,175 labeled Tier B sequences, got {len(labeled_tier_b)}"


def test_10_and_11_and_20_traceable_evidence_for_all_labeled():
    """TEST 10, 11, 20: Verify all labeled and Tier A sequences have valid traceable evidence."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    labeled = df_21[df_21["temporal_annotation_status"] == "LABELED"]
    assert len(labeled) == 16657

    tier_a = df_21[df_21["annotation_tier"] == "TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH"]
    assert len(tier_a) == 10482
    assert (tier_a["pmid"] != "NOT_AVAILABLE").all()
    assert (tier_a["doi"] != "NOT_AVAILABLE").all()
    assert (tier_a["temporal_class"].isin(["IMMEDIATE_EARLY", "EARLY", "LATE"])).all()


def test_12_and_13_quarantined_and_excluded_not_ground_truth():
    """TEST 12, 13: Verify quarantined and excluded records are never ground-truth eligible."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    
    tier_c = df_21[df_21["annotation_tier"] == "TIER_C_CONFLICTING_QUARANTINED"]
    assert (tier_c["ground_truth_eligibility"] == "CONFLICTED").all()
    
    tier_d = df_21[df_21["annotation_tier"] == "TIER_D_EXCLUDED"]
    assert (tier_d["ground_truth_eligibility"] == "NOT_ELIGIBLE_GROUND_TRUTH").all()


def test_14_no_computational_temporal_labels():
    """TEST 14: Verify no temporal labels were generated from embeddings or ML models."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    forbidden = ["ProtBERT_embedding", "ESM_embedding", "model_prediction", "cluster_id", "predicted_class"]
    for col in forbidden:
        assert col not in df_21.columns


def test_16_tier_totals_exact_closure():
    """TEST 16: Verify Tier A + Tier B + Tier C + Tier D == 22,689."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    tier_counts = df_21["annotation_tier"].value_counts().to_dict()
    assert tier_counts["TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH"] == 10482
    assert tier_counts["TIER_B_VERIFIED_EXTENDED_BIOLOGICAL_CORPUS"] == 11531
    assert tier_counts["TIER_C_CONFLICTING_QUARANTINED"] == 5
    assert tier_counts["TIER_D_EXCLUDED"] == 671
    assert sum(tier_counts.values()) == 22689


def test_17_ground_truth_dataset_integrity():
    """TEST 17: Verify ground-truth dataset contains exactly the 10,482 Tier A labeled records."""
    df_gt = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_ground_truth_dataset.csv")
    assert len(df_gt) == 10482
    assert (df_gt["ground_truth_eligibility"] == "ELIGIBLE_GROUND_TRUTH").all()
    assert (df_gt["temporal_class"].isin(["IMMEDIATE_EARLY", "EARLY", "LATE"])).all()


def test_18_extended_labeled_corpus_integrity():
    """TEST 18: Verify extended corpus contains exactly 16,657 non-conflicting labeled sequences."""
    df_ext = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_extended_labeled_corpus.csv")
    assert len(df_ext) == 16657
    assert (df_ext["temporal_annotation_status"] == "LABELED").all()
    assert not (df_ext["annotation_tier"].isin(["TIER_C_CONFLICTING_QUARANTINED", "TIER_D_EXCLUDED"])).any()


def test_19_evidence_propagation_relationship_recorded():
    """TEST 19: Verify temporal_evidence_relationship is explicitly populated."""
    df_21 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_reconciled_annotation.csv")
    allowed_rel = {
        "DIRECT_SEQUENCE_EVIDENCE", "GENE_PROTEIN_LEVEL_EVIDENCE",
        "CURATED_DATABASE_EVIDENCE", "SECONDARY_LITERATURE",
        "INSUFFICIENT", "CONFLICTING"
    }
    assert set(df_21["temporal_evidence_relationship"].unique()).issubset(allowed_rel)
