"""
===============================================================================
Phase 1B.1 Automated Test Suite: Conflict Resolution & Quarantine Verification
===============================================================================
"""

import os
from pathlib import Path
import pandas as pd
import pytest


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"


def test_conflict_resolution_table_exists():
    """Verify phase1b_final_conflict_resolution.csv exists and contains exactly 7 records."""
    table_path = RESULTS_DIR / "tables" / "phase1b_final_conflict_resolution.csv"
    assert table_path.exists()
    df = pd.read_csv(table_path)
    assert len(df) == 7

    # Check required columns
    expected_cols = [
        "sequence_id", "accession_1", "accession_2", "source_1", "source_2",
        "verified_species", "verified_gene", "verified_protein", "verified_strain",
        "evidence_source", "resolution_status", "resolution_reason"
    ]
    for col in expected_cols:
        assert col in df.columns


def test_specific_resolutions_accuracy():
    """Verify specific resolution outcomes for the 7 audited records."""
    table_path = RESULTS_DIR / "tables" / "phase1b_final_conflict_resolution.csv"
    df = pd.read_csv(table_path).set_index("sequence_id")

    # 1. CANONICAL_40dd30c0e1cc: UL54/ICP27 HSV-2 HG52
    assert df.loc["CANONICAL_40dd30c0e1cc", "resolution_status"] == "RESOLVED"
    assert df.loc["CANONICAL_40dd30c0e1cc", "verified_species"] == "HSV-2"

    # 2. CANONICAL_08f271887ce9: 1-aa fragment
    assert df.loc["CANONICAL_08f271887ce9", "resolution_status"] == "EXCLUDED — DOCUMENTED_REASON"

    # 3. 5 Quarantined records
    quarantine_ids = [
        "CANONICAL_01eb55d7fe80", "CANONICAL_5082249f592d", "CANONICAL_a9996922e85c",
        "CANONICAL_e9476297d514", "CANONICAL_ec8d790b87fd"
    ]
    for qid in quarantine_ids:
        assert df.loc[qid, "resolution_status"] == "UNRESOLVED — QUARANTINED"
        assert "AUTHORITATIVE RESOLUTION UNAVAILABLE" in df.loc[qid, "resolution_reason"]


def test_dataset_reconciliation_exact_totals():
    """Verify that after Phase 1B.1, all categories sum to exactly 22,689."""
    df = pd.read_csv(DATA_DIR / "processed" / "phase1b_biological_dataset.csv")
    assert len(df) == 22689

    counts = df["biological_eligibility"].value_counts().to_dict()
    assert counts.get("ELIGIBLE_PRIMARY") == 13321
    assert counts.get("ELIGIBLE_SEQUENCE_ANALYSIS_ONLY") == 8692
    assert counts.get("REVIEW_REQUIRED") == 5
    assert counts.get("EXCLUDE_NON_HSV") == 411
    assert counts.get("EXCLUDE_TECHNICAL") == 147
    assert counts.get("EXCLUDE_ENGINEERED") == 50
    assert counts.get("EXCLUDE_SYNTHETIC") == 49
    assert counts.get("EXCLUDE_UNRESOLVED_IDENTITY") == 14

    assert sum(counts.values()) == 22689
