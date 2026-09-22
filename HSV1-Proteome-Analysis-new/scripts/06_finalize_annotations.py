"""
Script 06 Finalize: Final Annotation Dataset Creation and Integrity Validation

Scientific Objective & Concept:
--------------------------------
Constructs the authoritative biological ground truth dataset (`data/annotations/temporal_annotations_final.csv`)
for all downstream supervised learning, unsupervised cluster validation, and representation comparisons.

Handling Biological Ambiguity (UL44 / glycoprotein C):
- In herpesvirus virology, viral glycoproteins are frequently grouped broadly as leaky late (Gamma-1).
- However, classical experimental inhibitor assays (PAA-block experiments; Conley et al. 1981, Homa et al. 1986)
  demonstrate stringent DNA-replication dependence for UL44 (gC), characterizing it as true late (Gamma-2).
- Rather than forcing an ad-hoc or arbitrary subclass decision, the primary class remains unambiguously 'Late',
  while late_subclass is explicitly labeled 'Conflicting' and verification_status is set to 'conflicting'.
- This preserves biological ground truth fidelity without compromising the primary 3-class target
  (Immediate-Early vs Early vs Late).

Validation Checks:
1. Exactly 74 non-redundant protein records.
2. 100% bijective match with `data/processed/unique_proteins.csv`.
3. No duplicate rows or missing values in primary class.
4. Class counts: IE=5, Early=15, Late=54 (Sum=74).
5. Explicit verification of UL44 status.
6. Immutability check on historical raw/audited annotation files.

Author: Computational Biology Pipeline
"""

import os
import sys
import yaml
import hashlib
import pandas as pd

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def compute_file_hash(filepath):
    """Computes SHA-256 hash for file immutability tracking."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()

def finalize_and_validate_annotations(config_path="config.yaml"):
    cfg = load_config(config_path)
    proc_dir = cfg['paths']['data_processed']
    annot_dir = cfg['paths']['data_annotations']
    tables_dir = cfg['paths']['results_tables']
    logs_dir = cfg['paths']['results_logs']
    
    os.makedirs(annot_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    
    unique_csv = os.path.join(proc_dir, "unique_proteins.csv")
    audited_csv = os.path.join(annot_dir, "temporal_annotations_audited.csv")
    orig_csv = os.path.join(annot_dir, "temporal_annotations.csv")
    final_csv = os.path.join(annot_dir, "temporal_annotations_final.csv")
    
    assert os.path.exists(unique_csv), f"Missing {unique_csv}"
    assert os.path.exists(audited_csv), f"Missing {audited_csv}"
    assert os.path.exists(orig_csv), f"Missing {orig_csv}"
    
    # Record baseline hashes for immutability check
    orig_hash_before = compute_file_hash(orig_csv)
    audited_hash_before = compute_file_hash(audited_csv)
    
    df_unique = pd.read_csv(unique_csv)
    df_audited = pd.read_csv(audited_csv)
    
    print("=" * 60)
    print("PHASE 6: FINALIZING BIOLOGICAL ANNOTATION DATASET")
    print("=" * 60)
    
    # Build final dataset by applying explicit UL44 rule
    df_final = df_audited.copy()
    
    # Update UL44 row specifically
    ul44_mask = df_final['protein_id'] == 'YP_009137119.1'
    assert ul44_mask.sum() == 1, "Error: UL44 (YP_009137119.1) not found in audited dataset!"
    
    df_final.loc[ul44_mask, 'temporal_class'] = 'Late'
    df_final.loc[ul44_mask, 'late_subclass'] = 'Conflicting'
    df_final.loc[ul44_mask, 'verification_status'] = 'conflicting'
    
    # Save final authoritative annotation CSV
    df_final.to_csv(final_csv, index=False)
    print(f"[Saved] Authoritative Dataset: {final_csv}")
    
    # Verify file immutability of predecessors
    orig_hash_after = compute_file_hash(orig_csv)
    audited_hash_after = compute_file_hash(audited_csv)
    assert orig_hash_before == orig_hash_after, "CRITICAL ERROR: Original temporal_annotations.csv was modified!"
    assert audited_hash_before == audited_hash_after, "CRITICAL ERROR: Audited temporal_annotations_audited.csv was modified!"
    
    # ==========================================
    # AUTOMATED INTEGRITY VALIDATION CHECKS
    # ==========================================
    checks_passed = True
    validation_log = []
    
    validation_log.append("=" * 60)
    validation_log.append("HSV-1 FINAL BIOLOGICAL ANNOTATION VALIDATION REPORT")
    validation_log.append("=" * 60)
    
    # Check 1: Exactly 74 unique protein IDs
    num_proteins = len(df_final['protein_id'].unique())
    c1 = (num_proteins == 74) and (len(df_final) == 74)
    checks_passed = checks_passed and c1
    validation_log.append(f"Check 1 - Exactly 74 Unique Protein IDs: {'PASS' if c1 else 'FAIL'} (Count = {num_proteins})")
    
    # Check 2: Every protein ID exists in unique_proteins.csv
    set_unique = set(df_unique['protein_id'])
    set_final = set(df_final['protein_id'])
    c2 = (set_unique == set_final)
    checks_passed = checks_passed and c2
    validation_log.append(f"Check 2 - 100% Match with unique_proteins.csv: {'PASS' if c2 else 'FAIL'}")
    
    # Check 3: No duplicate annotation rows
    c3 = (df_final['protein_id'].duplicated().sum() == 0)
    checks_passed = checks_passed and c3
    validation_log.append(f"Check 3 - Zero Duplicate Annotation Rows: {'PASS' if c3 else 'FAIL'}")
    
    # Check 4: No missing primary temporal class
    c4 = (df_final['temporal_class'].isna().sum() == 0) and ((df_final['temporal_class'] == '').sum() == 0)
    checks_passed = checks_passed and c4
    validation_log.append(f"Check 4 - No Missing Primary Temporal Class Values: {'PASS' if c4 else 'FAIL'}")
    
    # Check 5: Only three primary classes exist
    allowed_classes = {'Immediate-Early', 'Early', 'Late'}
    actual_classes = set(df_final['temporal_class'].unique())
    c5 = (actual_classes == allowed_classes)
    checks_passed = checks_passed and c5
    validation_log.append(f"Check 5 - Only Canonical 3 Primary Classes Exist: {'PASS' if c5 else 'FAIL'} ({actual_classes})")
    
    # Check 6: Class counts are IE=5, Early=15, Late=54
    ie_c = (df_final['temporal_class'] == 'Immediate-Early').sum()
    e_c = (df_final['temporal_class'] == 'Early').sum()
    l_c = (df_final['temporal_class'] == 'Late').sum()
    c6 = (ie_c == 5) and (e_c == 15) and (l_c == 54) and (ie_c + e_c + l_c == 74)
    checks_passed = checks_passed and c6
    validation_log.append(f"Check 6 - Exact Class Distribution (IE=5, Early=15, Late=54; Sum=74): {'PASS' if c6 else 'FAIL'} (IE={ie_c}, Early={e_c}, Late={l_c})")
    
    # Check 7: UL44 remains in dataset
    ul44_row = df_final[df_final['gene'] == 'UL44']
    c7 = len(ul44_row) == 1
    checks_passed = checks_passed and c7
    validation_log.append(f"Check 7 - UL44 (YP_009137119.1) Present in Dataset: {'PASS' if c7 else 'FAIL'}")
    
    # Check 8: UL44 primary class = Late
    c8 = (ul44_row['temporal_class'].values[0] == 'Late')
    checks_passed = checks_passed and c8
    validation_log.append(f"Check 8 - UL44 Primary Class is 'Late': {'PASS' if c8 else 'FAIL'}")
    
    # Check 9: UL44 subclass = Conflicting & status = conflicting
    c9 = (ul44_row['late_subclass'].values[0] == 'Conflicting') and (ul44_row['verification_status'].values[0] == 'conflicting')
    checks_passed = checks_passed and c9
    validation_log.append(f"Check 9 - UL44 Subclass is 'Conflicting' & Status 'conflicting': {'PASS' if c9 else 'FAIL'}")
    
    # Check 10: Immutability of predecessor files
    c10 = (orig_hash_before == orig_hash_after) and (audited_hash_before == audited_hash_after)
    checks_passed = checks_passed and c10
    validation_log.append(f"Check 10 - Predecessor Files Immutability Maintained: {'PASS' if c10 else 'FAIL'}")
    
    # Subclass breakdown summary for documentation
    gamma1_c = (df_final['late_subclass'] == 'Gamma-1').sum()
    gamma2_c = (df_final['late_subclass'] == 'Gamma-2').sum()
    confl_c = (df_final['late_subclass'] == 'Conflicting').sum()
    none_c = (df_final['late_subclass'] == 'None').sum()
    
    validation_log.append("\n--------------------------------------------------")
    validation_log.append("DETAILED BREAKDOWN OF AUTHORITATIVE DATASET")
    validation_log.append("--------------------------------------------------")
    validation_log.append(f"Immediate-Early (alpha):  {ie_c:2d} ({ie_c/74*100:5.2f}%)")
    validation_log.append(f"Early (beta):             {e_c:2d} ({e_c/74*100:5.2f}%)")
    validation_log.append(f"Late (gamma):              {l_c:2d} ({l_c/74*100:5.2f}%)")
    validation_log.append(f"  - Gamma-1 (Leaky Late): {gamma1_c:2d}")
    validation_log.append(f"  - Gamma-2 (True Late):  {gamma2_c:2d}")
    validation_log.append(f"  - Conflicting Subclass: {confl_c:2d} (UL44 / gC)")
    validation_log.append(f"  - Non-Late Subclass:    {none_c:2d} (IE & Early genes)")
    validation_log.append("--------------------------------------------------")
    validation_log.append(f"OVERALL VALIDATION STATUS: {'ALL CHECKS PASSED' if checks_passed else 'VALIDATION FAILED'}")
    validation_log.append("==================================================")
    
    # Generate Results Table: final_annotation_summary.csv
    summary_data = [
        {"metric": "Total Unique Proteins", "value": str(num_proteins), "expected": "74", "status": "PASS" if c1 else "FAIL"},
        {"metric": "ID Match with Unique Sequences", "value": "100%", "expected": "100%", "status": "PASS" if c2 else "FAIL"},
        {"metric": "Duplicate Rows", "value": "0", "expected": "0", "status": "PASS" if c3 else "FAIL"},
        {"metric": "Immediate-Early Count", "value": str(ie_c), "expected": "5", "status": "PASS" if ie_c == 5 else "FAIL"},
        {"metric": "Early Count", "value": str(e_c), "expected": "15", "status": "PASS" if e_c == 15 else "FAIL"},
        {"metric": "Late Count", "value": str(l_c), "expected": "54", "status": "PASS" if l_c == 54 else "FAIL"},
        {"metric": "Total Annotated Sum", "value": str(ie_c + e_c + l_c), "expected": "74", "status": "PASS" if (ie_c + e_c + l_c == 74) else "FAIL"},
        {"metric": "UL44 Primary Class", "value": str(ul44_row['temporal_class'].values[0]), "expected": "Late", "status": "PASS" if c8 else "FAIL"},
        {"metric": "UL44 Late Subclass", "value": str(ul44_row['late_subclass'].values[0]), "expected": "Conflicting", "status": "PASS" if c9 else "FAIL"},
        {"metric": "UL44 Verification Status", "value": str(ul44_row['verification_status'].values[0]), "expected": "conflicting", "status": "PASS" if c9 else "FAIL"},
        {"metric": "Historical Files Immutability", "value": "Preserved", "expected": "Preserved", "status": "PASS" if c10 else "FAIL"},
        {"metric": "Overall Pipeline Readiness", "value": "Verified", "expected": "Verified", "status": "PASS" if checks_passed else "FAIL"}
    ]
    df_summary = pd.DataFrame(summary_data)
    summary_path = os.path.join(tables_dir, "final_annotation_summary.csv")
    df_summary.to_csv(summary_path, index=False)
    print(f"[Saved] Summary Table: {summary_path}")
    
    # Save Validation Text Log
    log_text = "\n".join(validation_log)
    log_path = os.path.join(logs_dir, "final_annotation_validation.txt")
    with open(log_path, "w") as f:
        f.write(log_text)
    print(f"[Saved] Validation Log: {log_path}")
    
    print("\n" + log_text)
    
    if not checks_passed:
        print("\nCRITICAL WARNING: Integrity validation failed!")
        sys.exit(1)
        
    print("\n[STOP CONDITION REACHED] Phase 6 finalized and validated. Standing by for next phase instructions.")

if __name__ == "__main__":
    finalize_and_validate_annotations()
