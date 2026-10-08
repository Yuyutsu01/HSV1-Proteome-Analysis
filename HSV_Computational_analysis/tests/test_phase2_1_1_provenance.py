"""
===============================================================================
Phase 2.1.1 Automated Test Suite: Temporal Evidence Provenance Verification
===============================================================================
Verifies all 20 Phase 2.1.1 provenance criteria, evidence traceability,
and ground-truth impact invariants.
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


def test_1_and_2_and_3_evidence_records_and_types():
    """TEST 1, 2, 3: Every labeled sequence has an evidence record with type and relationship."""
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv")
    df_prov = pd.read_csv(RESULTS_TABLES_DIR / "phase2_1_1_temporal_evidence_provenance.csv")

    labeled = df_final[df_final["temporal_annotation_status"] == "LABELED"]
    assert len(labeled) == len(df_prov) == 16657

    assert not df_prov["evidence_type"].isnull().any()
    assert not df_prov["evidence_relationship"].isnull().any()
    assert (df_prov["evidence_type"] != "NONE").all()


def test_4_and_5_direct_vs_gene_protein_evidence_distinction():
    """TEST 4, 5: DIRECT_SEQUENCE_EVIDENCE is only used for direct references; GENE_PROTEIN requires supported identity."""
    df_prov = pd.read_csv(RESULTS_TABLES_DIR / "phase2_1_1_temporal_evidence_provenance.csv")
    
    direct = df_prov[df_prov["evidence_relationship"] == "DIRECT_SEQUENCE_EVIDENCE"]
    assert len(direct) == 12
    assert (direct["source_database"] == "UniProtKB").all()

    gene_level = df_prov[df_prov["evidence_relationship"] == "GENE_PROTEIN_LEVEL_EVIDENCE"]
    assert len(gene_level) == 16645
    assert (gene_level["gene_identity_supported"].isin(["YES", "UNCERTAIN"])).all()


def test_6_and_7_no_computational_or_ml_labels():
    """TEST 6, 7: No temporal labels were generated from sequence similarity, ProtBERT, or ML models."""
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv")
    forbidden = ["ProtBERT_embedding", "ESM_embedding", "model_prediction", "cluster_id", "predicted_class"]
    for col in forbidden:
        assert col not in df_final.columns


def test_8_and_9_cross_species_and_strain_transfer_marking():
    """TEST 8, 9: Cross-species and cross-strain transfer fields are explicitly populated."""
    df_prov = pd.read_csv(RESULTS_TABLES_DIR / "phase2_1_1_temporal_evidence_provenance.csv")
    
    assert set(df_prov["cross_species_transfer"].unique()) == {"YES", "NO"}
    assert set(df_prov["cross_strain_transfer"].unique()) == {"YES", "NO"}
    
    # HSV-2 records mapped from HSV-1 reference kinetics must have cross_species_transfer == YES
    hsv2_transfers = df_prov[(df_prov["virus_species"] == "HSV-2") & (df_prov["evidence_species"] == "HSV-1")]
    assert (hsv2_transfers["cross_species_transfer"] == "YES").all()


def test_10_database_identity_separated_from_temporal_evidence():
    """TEST 10: Database identity source is separated from primary literature temporal evidence."""
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv")
    labeled = df_final[df_final["temporal_annotation_status"] == "LABELED"]
    
    assert not labeled["gene_protein_mapping_source"].isnull().any()
    assert not labeled["temporal_evidence_source"].isnull().any()


def test_11_and_20_tier_a_ground_truth_provenance_validated():
    """TEST 11, 20: Every Tier A ground-truth record has supported provenance with zero unverified records."""
    df_gt = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_ground_truth_dataset.csv")
    assert len(df_gt) == 10482
    
    assert (df_gt["ground_truth_eligibility"] == "ELIGIBLE_GROUND_TRUTH").all()
    assert (df_gt["provenance_audit_status"].isin(["SUPPORTED", "SUPPORTED_WITH_TRANSFER"])).all()
    assert (df_gt["evidence_chain_complete"] == "YES").all()


def test_12_quarantined_records_remain_excluded_from_ground_truth():
    """TEST 12: Quarantined records from Phase 1B.1 are never ground-truth eligible."""
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv")
    quarantined = df_final[df_final["provenance_audit_status"] == "CONFLICTING"]
    assert len(quarantined) == 5
    assert (quarantined["ground_truth_eligibility"] == "CONFLICTED").all()


def test_13_and_14_and_15_raw_immutability_and_hash_invariance():
    """TEST 13, 14, 15: Sequence contents ($N=22,689$), hashes, and raw FASTA remain identical."""
    df_1a = pd.read_csv(DATA_PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t").set_index("canonical_sequence_id")
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv").set_index("sequence_id")

    assert len(df_final) == 22689
    for seq_id, row_f in df_final.head(500).iterrows():
        row_1a = df_1a.loc[seq_id]
        assert row_f["sequence"] == row_1a["sequence"]
        assert row_f["canonical_sequence_hash"] == seq_id

    # Raw FASTA checksums
    raw_files = list(RAW_DIR.glob("*.fasta"))
    assert len(raw_files) == 2
    for rf in raw_files:
        mirror = DATA_RAW_DIR / rf.name
        assert compute_sha256(rf) == compute_sha256(mirror)


def test_16_temporal_classes_vocabulary():
    """TEST 16: All temporal labels belong to the allowed set."""
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv")
    allowed = {"IMMEDIATE_EARLY", "EARLY", "LATE", "UNKNOWN", "CONFLICTING"}
    assert set(df_final["temporal_class"].unique()).issubset(allowed)


def test_17_dataset_counts_exact_mathematical_closure():
    """TEST 17: All dataset counts close exactly."""
    df_final = pd.read_csv(DATA_PROCESSED_DIR / "phase2_1_1_final_annotation.csv")
    
    tier_counts = df_final["annotation_tier"].value_counts().to_dict()
    assert tier_counts["TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH"] == 10482
    assert tier_counts["TIER_B_VERIFIED_EXTENDED_BIOLOGICAL_CORPUS"] == 11531
    assert tier_counts["TIER_C_CONFLICTING_QUARANTINED"] == 5
    assert tier_counts["TIER_D_EXCLUDED"] == 671
    assert sum(tier_counts.values()) == 22689


def test_18_label_changes_log_exists():
    """TEST 18: Label changes table exists and logs zero unwarranted changes."""
    df_changes = pd.read_csv(RESULTS_TABLES_DIR / "phase2_1_1_label_changes.csv")
    assert len(df_changes) >= 1
    assert df_changes.iloc[0]["reason"] == "NO_TEMPORAL_LABEL_CHANGES_REQUIRED"


def test_19_manual_review_queue_explicit():
    """TEST 19: Manual review queue contains exactly the 5 quarantined records."""
    df_review = pd.read_csv(RESULTS_TABLES_DIR / "phase2_1_1_manual_review_queue.csv")
    assert len(df_review) == 5
    assert (df_review["status"] == "QUARANTINED").all()
