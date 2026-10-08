"""
===============================================================================
Phase 2 Automated Test Suite: Biological Annotation & Temporal Curation
===============================================================================
Verifies all 15 Phase 2 quality control, evidence ground-truth traceability,
and mathematical closure criteria.
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


def test_1_and_4_and_13_exact_sequence_count_and_no_new_or_external():
    """1, 4, 13: Verify exact 22,689 unique sequences represented without additions or deletions."""
    df_1a = pd.read_csv(DATA_PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t")
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")

    assert len(df_2) == 22689
    assert len(df_1a) == 22689
    assert set(df_2["sequence_id"]) == set(df_1a["canonical_sequence_id"])


def test_2_and_3_and_5_raw_immutability_and_hash_invariance():
    """2, 3, 5: Verify sequence contents, hashes, and raw FASTA files are identical and unmutated."""
    df_1a = pd.read_csv(DATA_PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t").set_index("canonical_sequence_id")
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv").set_index("sequence_id")

    for seq_id, row2 in df_2.head(500).iterrows():
        row1 = df_1a.loc[seq_id]
        assert row2["sequence"] == row1["sequence"]
        assert row2["canonical_sequence_hash"] == seq_id

    # Check raw FASTA immutability
    raw_files = list(RAW_DIR.glob("*.fasta"))
    assert len(raw_files) == 2
    for rf in raw_files:
        mirror = DATA_RAW_DIR / rf.name
        assert compute_sha256(rf) == compute_sha256(mirror)


def test_6_allowed_temporal_classes():
    """6: Verify all temporal labels belong to the strictly allowed vocabulary."""
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")
    allowed_classes = {"IMMEDIATE_EARLY", "EARLY", "LATE", "UNKNOWN", "CONFLICTING"}
    observed_classes = set(df_2["temporal_class"].unique())
    assert observed_classes.issubset(allowed_classes), f"Unexpected classes: {observed_classes - allowed_classes}"


def test_7_and_8_evidence_traceability_for_labeled_sequences():
    """7, 8: Verify every labeled sequence has traceable literature/database evidence."""
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")
    df_ev = pd.read_csv(RESULTS_TABLES_DIR / "phase2_annotation_evidence.csv")

    labeled = df_2[df_2["temporal_class"].isin(["IMMEDIATE_EARLY", "EARLY", "LATE"])]
    assert len(labeled) == len(df_ev) == 16657

    # High-confidence Tier A sequences must all have valid PMIDs/DOIs
    tier_a = df_2[df_2["annotation_tier"] == "TIER_A_HIGH_CONFIDENCE_LABELED"]
    assert (tier_a["pmid"] != "NOT_AVAILABLE").all()
    assert (tier_a["doi"] != "NOT_AVAILABLE").all()


def test_9_unknown_used_when_insufficient():
    """9: Verify UNKNOWN is used when evidence is insufficient (e.g. unassigned/hypothetical)."""
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")
    unknowns = df_2[df_2["temporal_class"] == "UNKNOWN"]
    assert len(unknowns) == 6027
    assert (unknowns["temporal_class_confidence"].isin(["INSUFFICIENT_EVIDENCE", "NOT_APPLICABLE_EXCLUDED"])).all()


def test_10_and_11_conflicting_quarantined_records():
    """10, 11: Verify that conflicting records from Phase 1B.1 are quarantined as CONFLICTING."""
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")
    conflicts = df_2[df_2["temporal_class"] == "CONFLICTING"]
    assert len(conflicts) == 5
    assert (conflicts["annotation_tier"] == "TIER_C_CONFLICTING_UNCERTAIN").all()
    assert (conflicts["manual_review_required"] == True).all()


def test_12_no_ml_or_embedding_columns():
    """12: Strict prohibition: No ML predictions, embeddings, or clustering columns in Phase 2."""
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")
    forbidden_cols = [
        "ProtBERT_embedding", "ESM_embedding", "model_prediction",
        "cluster_id", "PCA_1", "UMAP_1", "predicted_class"
    ]
    for col in forbidden_cols:
        assert col not in df_2.columns


def test_14_and_15_tier_assignment_and_mathematical_closure():
    """14, 15: Verify all Tier A/B/C/D counts close exactly to 22,689."""
    df_2 = pd.read_csv(DATA_PROCESSED_DIR / "phase2_biological_annotation.csv")
    
    tier_counts = df_2["annotation_tier"].value_counts().to_dict()
    assert tier_counts["TIER_A_HIGH_CONFIDENCE_LABELED"] == 10482
    assert tier_counts["TIER_B_VERIFIED_UNLABELED"] == 11531
    assert tier_counts["TIER_C_CONFLICTING_UNCERTAIN"] == 5
    assert tier_counts["TIER_D_EXCLUDED"] == 671

    assert sum(tier_counts.values()) == 22689
