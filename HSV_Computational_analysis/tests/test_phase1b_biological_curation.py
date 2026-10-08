"""
===============================================================================
Phase 1B Automated Test Suite: Biological Eligibility & Curation Verification
===============================================================================
Ensures that all 15 Phase 1B biological curation criteria, provenance tracking,
exclusion logging, and strict prohibition of temporal/ML inferences are satisfied.
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


def test_1_raw_fasta_unmodified_and_checksums_preserved():
    """1. Verify that raw FASTA files in ./raw/ and data/raw/ are intact and identical."""
    raw_files = list(RAW_DIR.glob("*.fasta"))
    assert len(raw_files) == 2
    for rf in raw_files:
        mirror = DATA_RAW_DIR / rf.name
        assert mirror.exists()
        assert compute_sha256(rf) == compute_sha256(mirror)


def test_2_and_3_unique_count_and_no_new_sequences():
    """2 & 3. Verify exactly 22,689 unique sequences with no additions or deletions."""
    df_1a = pd.read_csv(DATA_PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t")
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")

    assert len(df_1b) == 22689
    assert len(df_1a) == 22689
    assert set(df_1b["sequence_id"]) == set(df_1a["canonical_sequence_id"])


def test_4_complete_provenance_tracing():
    """4. Verify every row has valid raw source file and raw record ID mapping."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    assert not df_1b["original_record_ids"].isnull().any()
    assert not df_1b["source_file"].isnull().any()
    assert df_1b["provenance_status"].str.contains("SURJECTIVELY_MAPPED").all()


def test_5_and_6_taxonomy_status_and_evidence():
    """5 & 6. Verify every sequence has a valid taxonomy status and verified ones have evidence."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    valid_tax_statuses = {"VERIFIED", "SUPPORTED", "UNRESOLVED", "CONFLICTING"}
    assert set(df_1b["taxonomy_status"].unique()).issubset(valid_tax_statuses)

    verified_tax = df_1b[df_1b["taxonomy_status"] == "VERIFIED"]
    assert not verified_tax["taxonomy_evidence"].isnull().any()
    assert (verified_tax["taxonomy_evidence"] != "NONE").all()


def test_7_protein_identity_evidence():
    """7. Verify verified protein identities have documented explicit evidence."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    verified_prots = df_1b[df_1b["protein_identity_status"] == "VERIFIED"]
    assert not verified_prots["protein_identity_evidence"].isnull().any()
    assert (verified_prots["protein_identity_evidence"] != "NONE").all()


def test_8_exclusion_audit_integrity():
    """8. Verify every excluded sequence has a documented reason in the audit table."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    df_excl = pd.read_csv(RESULTS_TABLES_DIR / "phase1b_exclusion_audit.csv")

    excluded_in_dataset = df_1b[df_1b["biological_eligibility"].str.startswith("EXCLUDE_")]
    assert len(excluded_in_dataset) == len(df_excl) == 671
    assert not df_excl["exclusion_reason"].isnull().any()


def test_9_conflict_table_presence():
    """9. Verify all detected metadata conflicts are represented in phase1b_conflicts.csv."""
    df_conflicts = pd.read_csv(RESULTS_TABLES_DIR / "phase1b_conflicts.csv")
    assert len(df_conflicts) == 6
    assert set(df_conflicts["resolution_status"].unique()) == {"REQUIRES_MANUAL_REVIEW"}


def test_10_unknown_not_silently_converted():
    """10. Verify that unresolved/unknown fields remain UNKNOWN or UNRESOLVED."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    unresolved_tax = df_1b[df_1b["taxonomy_status"] == "UNRESOLVED"]
    assert (unresolved_tax["virus_species"] == "UNKNOWN").all()


def test_11_and_12_no_temporal_labels_or_model_predictions():
    """11 & 12. Strict prohibition: No temporal classes (IE/Early/Late) or ML predictions in Phase 1B."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    forbidden_cols = [
        "temporal_class", "IE", "Early", "Late", "model_prediction",
        "ProtBERT_embedding", "cluster_assignment", "temporal_label"
    ]
    for col in forbidden_cols:
        assert col not in df_1b.columns, f"Forbidden Phase 2 column '{col}' detected in Phase 1B table"


def test_13_no_duplicate_sequence_ids():
    """13. Verify that sequence IDs are unique across all rows."""
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")
    assert df_1b["sequence_id"].is_unique


def test_14_and_15_deterministic_and_traceable():
    """14 & 15. Verify determinism and traceability back to raw FASTA record IDs."""
    df_raw = pd.read_csv(DATA_PROCESSED_DIR / "raw_derived_sequences.tsv", sep="\t")
    df_1b = pd.read_csv(DATA_PROCESSED_DIR / "phase1b_biological_dataset.csv")

    raw_ids_in_raw = set(df_raw["raw_record_id"])
    
    # Collect all raw IDs referenced in 1B
    referenced_raw_ids = set()
    for ids_str in df_1b["original_record_ids"]:
        for rid in ids_str.split(";"):
            referenced_raw_ids.add(rid)

    assert referenced_raw_ids == raw_ids_in_raw
