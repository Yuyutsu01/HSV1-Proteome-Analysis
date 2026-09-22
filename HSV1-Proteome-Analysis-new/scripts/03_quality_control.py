"""
Script 03: Automated Quality Control & Sequence Integrity Verification

Scientific Objective:
Performs rigorous sequence validation on all raw extracted CDS records:
1. Validates sequence length against minimum threshold (default >= 20 aa).
2. Verifies presence of exclusively canonical amino acid characters (ACDEFGHIKLMNPQRSTVWY).
3. Detects ambiguous/non-standard residues (X, B, Z, J, U, O).
4. Detects internal stop codons (*) or malformed translations.
5. Identifies duplicate accession identifiers and exact sequence duplicates.
6. Generates comprehensive audit trails in JSON and CSV format without silent data deletion.

Author: Computational Biology Pipeline
"""

import os
import json
import yaml
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def run_quality_control(config_path="config.yaml"):
    """
    Execute comprehensive QC inspection on raw extracted protein records.
    """
    cfg = load_config(config_path)
    inter_dir = cfg['paths']['data_intermediate']
    log_dir = cfg['paths']['results_logs']
    
    os.makedirs(inter_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    raw_csv_path = os.path.join(inter_dir, "raw_proteins.csv")
    assert os.path.exists(raw_csv_path), f"Error: raw_proteins.csv not found at {raw_csv_path}"
    
    df_raw = pd.read_csv(raw_csv_path)
    print(f"[Phase 3] Running Quality Control on {len(df_raw)} extracted protein records...")
    
    allowed_aa = set(cfg['qc'].get('allowed_amino_acids', "ACDEFGHIKLMNPQRSTVWY"))
    min_len = cfg['qc'].get('min_sequence_length', 20)
    
    qc_records = []
    seen_ids = set()
    duplicate_ids = []
    
    seen_seqs = {}
    
    passed_count = 0
    flagged_count = 0
    
    for idx, row in df_raw.iterrows():
        prot_id = str(row['protein_id'])
        seq = str(row['sequence']).upper()
        gene = str(row['gene'])
        seq_len = len(seq)
        
        flags = []
        is_passed = True
        
        # Check 1: Empty or null sequence
        if not seq or seq == "NONE" or seq == "NAN":
            flags.append("EMPTY_SEQUENCE")
            is_passed = False
            
        # Check 2: Minimum sequence length
        if seq_len < min_len:
            flags.append(f"SHORT_SEQUENCE (<{min_len}aa)")
            is_passed = False
            
        # Check 3: Invalid or ambiguous amino acids
        invalid_chars = sorted(list(set([aa for aa in seq if aa not in allowed_aa])))
        if invalid_chars:
            flags.append(f"NON_CANONICAL_RESIDUES({','.join(invalid_chars)})")
            # If only ambiguous (e.g. X), log flag
            
        # Check 4: Internal stop codons
        if "*" in seq[:-1]:
            flags.append("INTERNAL_STOP_CODON")
            is_passed = False
            
        # Check 5: Duplicate Protein ID
        if prot_id in seen_ids:
            flags.append("DUPLICATE_PROTEIN_ID")
            duplicate_ids.append(prot_id)
        seen_ids.add(prot_id)
        
        # Check 6: Exact Sequence Duplication tracking
        is_exact_duplicate = False
        duplicate_of = ""
        if seq in seen_seqs:
            is_exact_duplicate = True
            duplicate_of = seen_seqs[seq]
            flags.append(f"EXACT_SEQUENCE_DUPLICATE_OF_{duplicate_of}")
        else:
            seen_seqs[seq] = prot_id
            
        # Check 7: Missing descriptive metadata
        if pd.isna(row['gene']) or str(row['gene']).startswith("UNKNOWN"):
            flags.append("MISSING_GENE_NAME")
            
        if is_passed:
            passed_count += 1
        else:
            flagged_count += 1
            
        qc_records.append({
            "protein_id": prot_id,
            "gene": gene,
            "product": row['product'],
            "sequence_length": seq_len,
            "qc_status": "PASS" if is_passed else "FAIL",
            "is_exact_duplicate": is_exact_duplicate,
            "duplicate_of": duplicate_of,
            "non_canonical_residues": len(invalid_chars),
            "invalid_char_list": ",".join(invalid_chars),
            "qc_flags": "; ".join(flags) if flags else "NONE"
        })
        
    df_qc = pd.DataFrame(qc_records)
    qc_report_path = os.path.join(inter_dir, "qc_report.csv")
    df_qc.to_csv(qc_report_path, index=False)
    
    # Save QC passed records
    df_passed = df_raw[df_raw['protein_id'].isin(df_qc[df_qc['qc_status'] == 'PASS']['protein_id'])].copy()
    passed_csv_path = os.path.join(inter_dir, "qc_passed_proteins.csv")
    df_passed.to_csv(passed_csv_path, index=False)
    
    # Summary dictionary
    qc_summary = {
        "total_records_evaluated": len(df_raw),
        "qc_passed_records": passed_count,
        "qc_failed_records": flagged_count,
        "exact_sequence_duplicates_detected": sum(df_qc['is_exact_duplicate']),
        "unique_amino_acid_sequences": len(seen_seqs),
        "min_sequence_length_observed": int(df_raw['sequence_length'].min()),
        "max_sequence_length_observed": int(df_raw['sequence_length'].max()),
        "mean_sequence_length_observed": float(df_raw['sequence_length'].mean()),
        "ambiguous_or_invalid_residue_count": int(df_qc['non_canonical_residues'].sum()),
        "qc_report_path": qc_report_path,
        "qc_passed_path": passed_csv_path
    }
    
    summary_json_path = os.path.join(log_dir, "qc_summary.json")
    with open(summary_json_path, "w") as f:
        json.dump(qc_summary, f, indent=2)
        
    print(f"[Phase 3] QC Complete:")
    print(f"  - Total Records Evaluated: {len(df_raw)}")
    print(f"  - Passed Structural QC: {passed_count}")
    print(f"  - Exact Duplicate Copies Detected: {sum(df_qc['is_exact_duplicate'])}")
    print(f"  - Unique Sequence Count: {len(seen_seqs)}")
    print(f"  - Summary log saved to: {summary_json_path}")
    
    return df_qc, qc_summary

if __name__ == "__main__":
    run_quality_control()
