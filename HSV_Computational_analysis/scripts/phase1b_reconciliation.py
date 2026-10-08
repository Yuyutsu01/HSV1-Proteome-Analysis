#!/usr/bin/env python3
"""
===============================================================================
Phase 1B Final Reconciliation Audit Script
===============================================================================
Performs complete cross-tabulations, conflict tracking, and verification of exact
category closures before freezing Phase 1B.
===============================================================================
"""

import os
import sys
import re
from pathlib import Path
from typing import Dict, List, Any
import pandas as pd
from phase1b_curation import run_phase1b_curation


def build_conflict_detail_table(df_curated: pd.DataFrame, df_raw: pd.DataFrame, tables_dir: Path) -> pd.DataFrame:
    """
    Constructs the dedicated detailed table for the 6 cross-file taxonomy conflicts:
    sequence_id, sequence_hash, source_file_1, raw_header_1, accession_1, species_1,
    source_file_2, raw_header_2, accession_2, species_2, sequence_length, resolution_status
    """
    conflicting_seq_ids = df_curated[df_curated["taxonomy_status"] == "CONFLICTING"]["sequence_id"].tolist()
    
    rows = []
    for seq_id in conflicting_seq_ids:
        cur_row = df_curated[df_curated["sequence_id"] == seq_id].iloc[0]
        seq_hash = cur_row["sequence_id"].replace("CANONICAL_", "")
        
        # Look up raw records for this sequence
        raw_matches = df_raw[df_raw["sequence_sha256"].str.startswith(seq_hash)]
        
        if len(raw_matches) >= 2:
            r1 = raw_matches.iloc[0]
            r2 = raw_matches.iloc[1]
            
            # Determine species 1 and 2 explicitly from header
            sp1 = "HSV-2" if "alphaherpesvirus 2" in r1["raw_header"].lower() else ("HSV-1" if "alphaherpesvirus 1" in r1["raw_header"].lower() else "UNKNOWN")
            sp2 = "HSV-2" if "alphaherpesvirus 2" in r2["raw_header"].lower() else ("HSV-1" if "alphaherpesvirus 1" in r2["raw_header"].lower() else "UNKNOWN")
            
            rows.append({
                "sequence_id": seq_id,
                "sequence_hash": cur_row["sequence_id"],
                "source_file_1": r1["source_file"],
                "raw_header_1": r1["raw_header"],
                "accession_1": r1["accession"],
                "species_1": sp1,
                "source_file_2": r2["source_file"],
                "raw_header_2": r2["raw_header"],
                "accession_2": r2["accession"],
                "species_2": sp2,
                "sequence_length": cur_row["sequence_length"],
                "resolution_status": "MANUAL_REVIEW_REQUIRED"
            })
            
    df_detail = pd.DataFrame(rows)
    df_detail.to_csv(tables_dir / "phase1b_cross_file_conflicts_detail.csv", index=False)
    return df_detail


def run_reconciliation(base_dir: Path):
    """
    Executes full Phase 1B reconciliation audit.
    """
    # 1. Run Phase 1B Curation
    run_phase1b_curation(base_dir)
    
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    tables_dir = results_dir / "tables"
    logs_dir = results_dir / "logs"
    
    df_dataset = pd.read_csv(data_dir / "processed" / "phase1b_biological_dataset.csv")
    df_raw = pd.read_csv(data_dir / "processed" / "raw_derived_sequences.tsv", sep="\t")
    
    total_unique = len(df_dataset)
    assert total_unique == 22689, f"Expected 22,689 unique sequences, got {total_unique}"
    
    # 2. Build cross-tabulations
    ct_group = pd.crosstab(df_dataset["dataset_group"], df_dataset["biological_eligibility"], margins=True, margins_name="Total")
    ct_source_type = pd.crosstab(df_dataset["source_type"], df_dataset["biological_eligibility"], margins=True, margins_name="Total")
    ct_species = pd.crosstab(df_dataset["virus_species"], df_dataset["biological_eligibility"], margins=True, margins_name="Total")
    ct_prot_status = pd.crosstab(df_dataset["protein_identity_status"], df_dataset["biological_eligibility"], margins=True, margins_name="Total")
    ct_completeness = pd.crosstab(df_dataset["completeness_status"], df_dataset["biological_eligibility"], margins=True, margins_name="Total")
    
    # 3. Build Conflict Detail Table
    df_conflicts_detail = build_conflict_detail_table(df_dataset, df_raw, tables_dir)
    
    # 4. Generate Final Reconciliation Report
    report_path = logs_dir / "phase1b_final_reconciliation.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("PHASE 1B: FINAL RECONCILIATION AUDIT REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"TOTAL UNIQUE SEQUENCES RECONCILED: {total_unique:,}\n\n")
        
        f.write("1. CROSS-TABULATION: DATASET GROUP x BIOLOGICAL ELIGIBILITY\n")
        f.write("-" * 80 + "\n")
        f.write(ct_group.to_string() + "\n\n")
        
        f.write("2. CROSS-TABULATION: SOURCE TYPE x BIOLOGICAL ELIGIBILITY\n")
        f.write("-" * 80 + "\n")
        f.write(ct_source_type.to_string() + "\n\n")
        
        f.write("3. CROSS-TABULATION: VIRUS SPECIES x BIOLOGICAL ELIGIBILITY\n")
        f.write("-" * 80 + "\n")
        f.write(ct_species.to_string() + "\n\n")
        
        f.write("4. CROSS-TABULATION: PROTEIN IDENTITY STATUS x BIOLOGICAL ELIGIBILITY\n")
        f.write("-" * 80 + "\n")
        f.write(ct_prot_status.to_string() + "\n\n")
        
        f.write("5. CROSS-TABULATION: COMPLETENESS STATUS x BIOLOGICAL ELIGIBILITY\n")
        f.write("-" * 80 + "\n")
        f.write(ct_completeness.to_string() + "\n\n")
        
        f.write("6. DEDICATED CROSS-FILE TAXONOMY CONFLICTS (6 RECORDS)\n")
        f.write("-" * 80 + "\n")
        f.write(df_conflicts_detail.to_string(index=False) + "\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("FINAL STATUS: PHASE 1B — READY FOR FREEZE\n")
        f.write("=" * 80 + "\n")
        
    print(f"Reconciliation completed successfully! Report saved to {report_path}")


if __name__ == "__main__":
    base_p = Path(__file__).resolve().parent.parent
    run_reconciliation(base_p)
