"""
===============================================================================
Phase 1 Automated Test Suite
===============================================================================
Verifies that all Phase 1 data integrity, provenance, determinism, QC terminology,
and immutability criteria are strictly met.
===============================================================================
"""

import os
import hashlib
from pathlib import Path
import pandas as pd
import pytest


# Define base project directories
BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "raw"
DATA_RAW_DIR = BASE_DIR / "data" / "raw"
INTERMEDIATE_DIR = BASE_DIR / "data" / "intermediate"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
RESULTS_DIR = BASE_DIR / "results"


def compute_sha256(file_path: Path) -> str:
    """Helper to compute SHA-256 digest of a file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(65536):
            sha256.update(chunk)
    return sha256.hexdigest()


def test_raw_files_exist_and_unmutated():
    """Ensure raw files exist in ./raw/ and match their copies in data/raw/."""
    assert RAW_DIR.exists(), "Raw source directory does not exist"
    raw_files = list(RAW_DIR.glob("*.fasta"))
    assert len(raw_files) == 2, f"Expected 2 raw FASTA files, found {len(raw_files)}"

    for rf in raw_files:
        dest_copy = DATA_RAW_DIR / rf.name
        assert dest_copy.exists(), f"Raw file mirror {dest_copy} missing"
        # Check that checksums are identical
        assert compute_sha256(rf) == compute_sha256(dest_copy), f"SHA-256 mismatch for {rf.name}"


def test_manifests_consistency():
    """Verify manifests exist, are non-empty, and report exact record counts."""
    intermediate_manifest = INTERMEDIATE_DIR / "raw_manifest.csv"
    data_raw_manifest = DATA_RAW_DIR / "RAW_DATA_MANIFEST.csv"

    assert intermediate_manifest.exists()
    assert data_raw_manifest.exists()

    df_manifest = pd.read_csv(intermediate_manifest)
    assert len(df_manifest) == 2
    assert "sha256" in df_manifest.columns
    assert "record_count" in df_manifest.columns
    assert df_manifest["record_count"].sum() == 22710


def test_raw_derived_sequences_provenance():
    """Verify that every raw sequence is recorded with valid provenance."""
    tsv_path = PROCESSED_DIR / "raw_derived_sequences.tsv"
    assert tsv_path.exists()

    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) == 22710, f"Expected 22,710 raw records, found {len(df)}"

    # Required provenance columns
    required_cols = [
        "source_file", "record_index", "raw_record_id", "raw_header",
        "accession", "protein_name", "organism", "sequence",
        "sequence_length", "sequence_sha256", "qc_status"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing required column {col}"

    # Verify no nulls in critical identifiers
    assert not df["raw_record_id"].isnull().any()
    assert not df["sequence"].isnull().any()
    assert not df["sequence_sha256"].isnull().any()


def test_qc_categories_strict_terminology():
    """Verify all sequences are assigned to strict Phase 1 QC terminology."""
    tsv_path = PROCESSED_DIR / "raw_derived_sequences.tsv"
    df = pd.read_csv(tsv_path, sep="\t")
    
    valid_qc_categories = {
        "CANONICAL_20AA", "AMBIGUOUS_RESIDUE", "NON_STANDARD_RESIDUE",
        "STOP_SYMBOL", "INVALID_CHARACTER", "EXCLUDED_EMPTY"
    }
    observed_categories = set(df["qc_status"].unique())
    assert observed_categories.issubset(valid_qc_categories), f"Unexpected QC categories: {observed_categories - valid_qc_categories}"
    assert len(df[df["qc_status"] == "CANONICAL_20AA"]) == 15229
    assert len(df[df["qc_status"] == "AMBIGUOUS_RESIDUE"]) == 7481


def test_surjective_provenance_mapping():
    """
    Verify provenance mapping topology:
    1. Every raw FASTA record maps to exactly one canonical sequence (or documented exclusion).
    2. Every canonical sequence maps to >= 1 raw source records.
    """
    df_raw = pd.read_csv(PROCESSED_DIR / "raw_derived_sequences.tsv", sep="\t")
    df_unique = pd.read_csv(PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t")
    df_prov = pd.read_csv(INTERMEDIATE_DIR / "phase1_provenance_audit.csv")

    # 1. Total raw records mapped
    assert len(df_prov) == len(df_raw) == 22710
    assert (df_prov["mapped_canonical_id"] != "EXCLUDED_OR_UNMAPPED").sum() == 22710

    # 2. Every canonical sequence is represented
    mapped_canonical_ids = set(df_prov["mapped_canonical_id"].unique())
    unique_canonical_ids = set(df_unique["canonical_sequence_id"].unique())
    assert mapped_canonical_ids == unique_canonical_ids

    # 3. Every canonical sequence maps to >= 1 raw source records
    for _, row in df_unique.iterrows():
        assert row["occurrence_count"] >= 1
        raw_srcs = row["source_record_ids"].split(";")
        assert len(raw_srcs) == row["occurrence_count"]


def test_hash_determinism_and_invariance():
    """Verify hashing algorithm determinism: SHA-256(uppercase sequence)[:12]."""
    df_unique = pd.read_csv(PROCESSED_DIR / "unique_high_quality_sequences.tsv", sep="\t")
    
    for _, row in df_unique.head(500).iterrows():
        expected_hash = hashlib.sha256(row["sequence"].strip().upper().encode("utf-8")).hexdigest()
        assert row["sequence_sha256"] == expected_hash
        assert row["canonical_sequence_id"] == f"CANONICAL_{expected_hash[:12]}"


def test_unique_high_quality_sequences_deduplication():
    """Verify deduplication integrity: unique count + duplicates == valid count."""
    unique_path = PROCESSED_DIR / "unique_high_quality_sequences.tsv"
    assert unique_path.exists()

    df_unique = pd.read_csv(unique_path, sep="\t")
    assert len(df_unique) == 22689, f"Expected 22,689 unique sequences, found {len(df_unique)}"

    # Verify occurrence counts sum to valid record count
    total_occurrences = df_unique["occurrence_count"].sum()
    assert total_occurrences == 22710, f"Occurrences sum ({total_occurrences}) must equal valid count (22710)"

    # Verify all sequence hashes in unique table are distinct
    assert df_unique["sequence_sha256"].nunique() == len(df_unique)


def test_fasta_export_integrity():
    """Verify generated unique FASTA exists and matches unique sequence count."""
    fasta_path = PROCESSED_DIR / "unique_high_quality_sequences.fasta"
    assert fasta_path.exists()

    header_count = 0
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                header_count += 1

    assert header_count == 22689, f"FASTA headers ({header_count}) must match unique count (22689)"
