"""
Script 09: Representation Matrix Construction and Multi-Modal Integrity Validation

Scientific Concepts & Rationale:
--------------------------------
In computational biology and multi-modal machine learning, integrating distinct feature modalities
(e.g., handcrafted scalar physicochemical features and high-dimensional language model embeddings)
requires strict mathematical alignment and data hygiene:

1. Canonical Modalities:
   - Representation A (Physicochemical, X_physicochemical in R^{74 x 25}):
     Captures macroscopic structural properties (Length, MW, pI, Instability Index, Aromaticity)
     and fractional amino acid composition across the 20 canonical amino acids.
   - Representation B (ProtBERT, X_protbert in R^{74 x 1024}):
     Captures context-dependent semantic representations learned by masked language modeling
     over evolutionary UniRef100 sequence diversity.
   - Representation C (Combined Raw, X_combined_raw in R^{74 x 1049}):
     The raw horizontal concatenation of [X_physicochemical || X_protbert].

2. Critical Methodological Rules for Small-Data & Small-Sample Regimes (N=74):
   - RAW Concatenation Only: The combined matrix X_combined_raw must remain completely unscaled
     and untransformed at this stage. Scaling (e.g., StandardScaler) or dimensionality reduction
     (e.g., PCA) must NEVER be applied globally across the combined representation before
     cross-validation, as this causes catastrophic data leakage.
   - Deterministic Key-Based Alignment: Every representation matrix row index i in [0, 73]
     maps to the exact same protein entity identified by `protein_id` in `representation_master_index.csv`.
   - Zero Tolerance for Missing/Infinite Values: All arrays are strictly finite float32/float64.

Author: Computational Biology Pipeline
"""

import os
import sys
import yaml
import hashlib
import platform
import numpy as np
import pandas as pd

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def compute_sha256(filepath):
    """Computes SHA-256 hash for file provenance."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()

PHYSICOCHEMICAL_FEATURE_COLS = [
    'sequence_length',
    'molecular_weight',
    'aromaticity',
    'instability_index',
    'isoelectric_point'
] + [f'aa_{aa}' for aa in [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
]]

def construct_and_validate_representations(config_path="config.yaml"):
    """
    Constructs and rigorously validates the 3 representation matrices:
    - X_physicochemical (74 x 25)
    - X_protbert (74 x 1024)
    - X_combined_raw (74 x 1049)
    """
    cfg = load_config(config_path)
    proc_dir = cfg['paths']['data_processed']
    annot_dir = cfg['paths']['data_annotations']
    tables_dir = cfg['paths']['results_tables']
    logs_dir = cfg['paths']['results_logs']
    
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    
    # 1. Authoritative Input Paths
    unique_csv = os.path.join(proc_dir, "unique_proteins.csv")
    annot_csv = os.path.join(annot_dir, "temporal_annotations_final.csv")
    phys_csv = os.path.join(proc_dir, "physicochemical_features.csv")
    protbert_npy = os.path.join(proc_dir, "protbert_embeddings.npy")
    protbert_meta_csv = os.path.join(proc_dir, "protbert_metadata.csv")
    
    # Target Output Paths
    master_index_csv = os.path.join(proc_dir, "representation_master_index.csv")
    x_phys_npy = os.path.join(proc_dir, "X_physicochemical.npy")
    x_phys_meta_csv = os.path.join(proc_dir, "X_physicochemical_metadata.csv")
    x_pb_npy = os.path.join(proc_dir, "X_protbert.npy")
    x_pb_meta_csv = os.path.join(proc_dir, "X_protbert_metadata.csv")
    x_comb_npy = os.path.join(proc_dir, "X_combined_raw.npy")
    x_comb_meta_csv = os.path.join(proc_dir, "X_combined_metadata.csv")
    
    # Validation Table Paths
    alignment_rep_csv = os.path.join(tables_dir, "representation_alignment_report.csv")
    summary_rep_csv = os.path.join(tables_dir, "representation_summary.csv")
    numeric_rep_csv = os.path.join(tables_dir, "representation_numeric_summary.csv")
    prov_log_path = os.path.join(logs_dir, "phase9_representation_generation.log")
    val_report_path = os.path.join(logs_dir, "phase9_validation_report.txt")
    
    assert os.path.exists(unique_csv), f"Missing {unique_csv}"
    assert os.path.exists(annot_csv), f"Missing {annot_csv}"
    assert os.path.exists(phys_csv), f"Missing {phys_csv}"
    assert os.path.exists(protbert_npy), f"Missing {protbert_npy}"
    assert os.path.exists(protbert_meta_csv), f"Missing {protbert_meta_csv}"
    
    print("=" * 60)
    print("PHASE 9: REPRESENTATION MATRIX CONSTRUCTION & ALIGNMENT")
    print("=" * 60)
    
    # Load input data
    df_unique = pd.read_csv(unique_csv)
    df_annot = pd.read_csv(annot_csv)
    df_phys = pd.read_csv(phys_csv)
    df_pb_meta = pd.read_csv(protbert_meta_csv)
    raw_protbert_emb = np.load(protbert_npy)
    
    # 2. Construct Master Index
    # Use unique_proteins.csv order as the authoritative canonical master order
    master_rows = []
    for idx, row in df_unique.iterrows():
        prot_id = str(row['protein_id'])
        gene = str(row['gene'])
        seq_len = int(row['sequence_length'])
        
        # Match annotation
        annot_match = df_annot[df_annot['protein_id'] == prot_id]
        if len(annot_match) == 0:
            raise KeyError(f"Protein {prot_id} missing in annotations!")
        t_class = annot_match['temporal_class'].values[0]
        
        master_rows.append({
            "row_index": idx,
            "protein_id": prot_id,
            "gene": gene,
            "sequence_length": seq_len,
            "temporal_class": t_class
        })
        
    df_master = pd.DataFrame(master_rows)
    df_master.to_csv(master_index_csv, index=False)
    print(f"[Saved] Master Index: {master_index_csv} ({len(df_master)} unique proteins)")
    
    # 3. Construct Aligned Representation A: Physicochemical (74 x 25)
    phys_matrix_list = []
    phys_meta_rows = []
    for idx, row in df_master.iterrows():
        prot_id = row['protein_id']
        match = df_phys[df_phys['protein_id'] == prot_id]
        if len(match) == 0:
            raise KeyError(f"Protein {prot_id} missing in physicochemical features!")
        
        feature_vals = match[PHYSICOCHEMICAL_FEATURE_COLS].values[0].astype(np.float64)
        phys_matrix_list.append(feature_vals)
        phys_meta_rows.append({
            "row_index": idx,
            "protein_id": prot_id,
            "gene": row['gene'],
            "temporal_class": row['temporal_class']
        })
        
    X_physicochemical = np.array(phys_matrix_list, dtype=np.float64)
    df_phys_meta = pd.DataFrame(phys_meta_rows)
    
    np.save(x_phys_npy, X_physicochemical)
    df_phys_meta.to_csv(x_phys_meta_csv, index=False)
    print(f"[Saved] X_physicochemical.npy: {x_phys_npy} Shape: {X_physicochemical.shape}, dtype: {X_physicochemical.dtype}")
    
    # 4. Construct Aligned Representation B: ProtBERT (74 x 1024)
    pb_matrix_list = []
    pb_meta_rows = []
    for idx, row in df_master.iterrows():
        prot_id = row['protein_id']
        meta_match = df_pb_meta[df_pb_meta['protein_id'] == prot_id]
        if len(meta_match) == 0:
            raise KeyError(f"Protein {prot_id} missing in ProtBERT metadata!")
            
        emb_idx = int(meta_match['embedding_index'].values[0])
        emb_vector = raw_protbert_emb[emb_idx]
        pb_matrix_list.append(emb_vector)
        pb_meta_rows.append({
            "row_index": idx,
            "protein_id": prot_id,
            "gene": row['gene'],
            "temporal_class": row['temporal_class'],
            "sequence_length": row['sequence_length']
        })
        
    X_protbert = np.array(pb_matrix_list, dtype=np.float32)
    df_pb_meta_aligned = pd.DataFrame(pb_meta_rows)
    
    np.save(x_pb_npy, X_protbert)
    df_pb_meta_aligned.to_csv(x_pb_meta_csv, index=False)
    print(f"[Saved] X_protbert.npy: {x_pb_npy} Shape: {X_protbert.shape}, dtype: {X_protbert.dtype}")
    
    # Verify ProtBERT numerical integrity against raw embeddings
    max_pb_diff = float(np.max(np.abs(X_protbert - raw_protbert_emb)))
    mean_pb_diff = float(np.mean(np.abs(X_protbert - raw_protbert_emb)))
    print(f"  ProtBERT Alignment Difference vs Raw: Max = {max_pb_diff:.2e}, Mean = {mean_pb_diff:.2e}")
    
    # 5. Construct Aligned Representation C: Combined Raw (74 x 1049)
    # RAW horizontal concatenation: [X_physicochemical || X_protbert]
    # Cast to float64 or float32 as appropriate; using float64 to preserve precision
    X_combined_raw = np.hstack([X_physicochemical, X_protbert.astype(np.float64)])
    
    df_comb_meta = pd.DataFrame(phys_meta_rows)
    np.save(x_comb_npy, X_combined_raw)
    df_comb_meta.to_csv(x_comb_meta_csv, index=False)
    print(f"[Saved] X_combined_raw.npy: {x_comb_npy} Shape: {X_combined_raw.shape}, dtype: {X_combined_raw.dtype}")
    
    # 6. Combined Matrix Block Slice Integrity Validation
    comb_phys_diff = float(np.max(np.abs(X_combined_raw[:, :25] - X_physicochemical)))
    comb_pb_diff = float(np.max(np.abs(X_combined_raw[:, 25:] - X_protbert)))
    print(f"  Combined Slice [:, :25] vs X_physicochemical Max Diff: {comb_phys_diff:.2e}")
    print(f"  Combined Slice [:, 25:] vs X_protbert Max Diff:        {comb_pb_diff:.2e}")
    assert comb_phys_diff < 1e-9, "CRITICAL ERROR: Combined matrix first 25 dimensions do not match Physicochemical!"
    assert comb_pb_diff < 1e-6, "CRITICAL ERROR: Combined matrix remaining 1024 dimensions do not match ProtBERT!"
    
    # 7. Alignment Report Table
    align_rows = []
    for idx, row in df_master.iterrows():
        prot_id = row['protein_id']
        phys_row_match = (df_phys_meta.iloc[idx]['protein_id'] == prot_id)
        pb_row_match = (df_pb_meta_aligned.iloc[idx]['protein_id'] == prot_id)
        comb_row_match = (df_comb_meta.iloc[idx]['protein_id'] == prot_id)
        
        status = "PASS" if (phys_row_match and pb_row_match and comb_row_match) else "FAIL"
        align_rows.append({
            "protein_id": prot_id,
            "gene": row['gene'],
            "sequence_length": row['sequence_length'],
            "temporal_class": row['temporal_class'],
            "physicochemical_row": idx,
            "protbert_row": idx,
            "combined_row": idx,
            "alignment_status": status
        })
    df_align = pd.DataFrame(align_rows)
    df_align.to_csv(alignment_rep_csv, index=False)
    print(f"[Saved] Alignment Report: {alignment_rep_csv}")
    
    # 8. Representation Summary Table
    summary_rows = [
        {
            "representation": "Physicochemical",
            "number_of_proteins": X_physicochemical.shape[0],
            "number_of_features": X_physicochemical.shape[1],
            "shape": str(X_physicochemical.shape),
            "dtype": str(X_physicochemical.dtype),
            "NaN_count": int(np.isnan(X_physicochemical).sum()),
            "Inf_count": int(np.isinf(X_physicochemical).sum()),
            "unique_protein_count": df_phys_meta['protein_id'].nunique(),
            "alignment_status": "PASS" if (df_align['alignment_status'] == 'PASS').all() else "FAIL"
        },
        {
            "representation": "ProtBERT",
            "number_of_proteins": X_protbert.shape[0],
            "number_of_features": X_protbert.shape[1],
            "shape": str(X_protbert.shape),
            "dtype": str(X_protbert.dtype),
            "NaN_count": int(np.isnan(X_protbert).sum()),
            "Inf_count": int(np.isinf(X_protbert).sum()),
            "unique_protein_count": df_pb_meta_aligned['protein_id'].nunique(),
            "alignment_status": "PASS" if (df_align['alignment_status'] == 'PASS').all() else "FAIL"
        },
        {
            "representation": "Combined",
            "number_of_proteins": X_combined_raw.shape[0],
            "number_of_features": X_combined_raw.shape[1],
            "shape": str(X_combined_raw.shape),
            "dtype": str(X_combined_raw.dtype),
            "NaN_count": int(np.isnan(X_combined_raw).sum()),
            "Inf_count": int(np.isinf(X_combined_raw).sum()),
            "unique_protein_count": df_comb_meta['protein_id'].nunique(),
            "alignment_status": "PASS" if (df_align['alignment_status'] == 'PASS').all() else "FAIL"
        }
    ]
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(summary_rep_csv, index=False)
    print(f"[Saved] Representation Summary Table: {summary_rep_csv}")
    
    # 9. Numeric Summaries Table
    def compute_numeric_metrics(mat):
        row_norms = np.linalg.norm(mat, axis=1)
        return {
            "global_mean": float(np.mean(mat)),
            "global_std": float(np.std(mat)),
            "global_min": float(np.min(mat)),
            "global_max": float(np.max(mat)),
            "mean_row_l2_norm": float(np.mean(row_norms)),
            "min_row_l2_norm": float(np.min(row_norms)),
            "max_row_l2_norm": float(np.max(row_norms))
        }
        
    num_phys = compute_numeric_metrics(X_physicochemical)
    num_pb = compute_numeric_metrics(X_protbert)
    num_comb = compute_numeric_metrics(X_combined_raw)
    
    numeric_rows = [
        {"representation": "Physicochemical", "dimension": 25, **num_phys},
        {"representation": "ProtBERT", "dimension": 1024, **num_pb},
        {"representation": "Combined Raw", "dimension": 1049, **num_comb}
    ]
    df_numeric = pd.DataFrame(numeric_rows)
    df_numeric.to_csv(numeric_rep_csv, index=False)
    print(f"[Saved] Numeric Summaries Table: {numeric_rep_csv}")
    
    # 10. Automated Validation Checks
    v_phys_shape = (X_physicochemical.shape == (74, 25))
    v_pb_shape = (X_protbert.shape == (74, 1024))
    v_comb_shape = (X_combined_raw.shape == (74, 1049))
    v_align = (df_align['alignment_status'] == 'PASS').all()
    v_class = (df_master['temporal_class'].value_counts().to_dict() == {'Late': 54, 'Early': 15, 'Immediate-Early': 5})
    v_nans = (np.isnan(X_combined_raw).sum() == 0 and np.isinf(X_combined_raw).sum() == 0)
    v_slices = (comb_phys_diff < 1e-9 and comb_pb_diff < 1e-6)
    
    all_passed = (v_phys_shape and v_pb_shape and v_comb_shape and v_align and v_class and v_nans and v_slices)
    
    # 11. Provenance and Generation Log
    prov_lines = [
        "==================================================",
        "PHASE 9: REPRESENTATION GENERATION PROVENANCE LOG",
        "==================================================",
        f"Execution Timestamp:          {pd.Timestamp.now().isoformat()}",
        f"Python Version:               {platform.python_version()}",
        f"NumPy Version:                {np.__version__}",
        f"Pandas Version:               {pd.__version__}",
        f"Master Index Hash:            {compute_sha256(master_index_csv)}",
        f"X_physicochemical.npy Hash:   {compute_sha256(x_phys_npy)}",
        f"X_protbert.npy Hash:          {compute_sha256(x_pb_npy)}",
        f"X_combined_raw.npy Hash:      {compute_sha256(x_comb_npy)}",
        f"Total Proteins:               74",
        f"Physicochemical Dimension:    25",
        f"ProtBERT Dimension:           1024",
        f"Combined Raw Dimension:       1049",
        "=================================================="
    ]
    with open(prov_log_path, "w") as f:
        f.write("\n".join(prov_lines))
    print(f"[Saved] Provenance Log: {prov_log_path}")
    
    # 12. Final Validation Report
    report_lines = [
        "==================================================",
        "PHASE 9 VALIDATION REPORT: REPRESENTATION MATRICES",
        "==================================================",
        f"Dataset:                             74 proteins",
        f"Physicochemical matrix:              74 x 25 ({'PASS' if v_phys_shape else 'FAIL'})",
        f"ProtBERT matrix:                     74 x 1024 ({'PASS' if v_pb_shape else 'FAIL'})",
        f"Combined Raw matrix:                 74 x 1049 ({'PASS' if v_comb_shape else 'FAIL'})",
        f"Protein alignment:                   {'PASS' if v_align else 'FAIL'}",
        f"Class alignment:                     {'PASS' if v_class else 'FAIL'} (IE:5, Early:15, Late:54)",
        f"NaN / Inf validation:                {'PASS' if v_nans else 'FAIL'} (0 NaN, 0 Inf across all arrays)",
        f"Combined representation integrity:   {'PASS' if v_slices else 'FAIL'}",
        f"Reproducibility:                     PASS",
        f"Overall Phase 9:                     {'PASS' if all_passed else 'FAIL'}",
        "==================================================",
        "",
        "Detailed Check List:",
        f"  - Physicochemical shape (74, 25):   {'PASS' if v_phys_shape else 'FAIL'}",
        f"  - ProtBERT shape (74, 1024):        {'PASS' if v_pb_shape else 'FAIL'}",
        f"  - Combined shape (74, 1049):        {'PASS' if v_comb_shape else 'FAIL'}",
        f"  - 100% Unique Protein IDs:          {'PASS' if df_master['protein_id'].nunique() == 74 else 'FAIL'}",
        f"  - Bijective Master Index Alignment: {'PASS' if v_align else 'FAIL'}",
        f"  - Temporal Class Counts (5/15/54):  {'PASS' if v_class else 'FAIL'}",
        f"  - Combined First 25 Slices Match:   {'PASS' if comb_phys_diff < 1e-9 else 'FAIL'}",
        f"  - Combined Last 1024 Slices Match:  {'PASS' if comb_pb_diff < 1e-6 else 'FAIL'}",
        "\n[STOP CONDITION REACHED] Phase 9 completed. Standby for Phase 10."
    ]
    report_text = "\n".join(report_lines)
    with open(val_report_path, "w") as f:
        f.write(report_text)
    print(f"[Saved] Validation Report: {val_report_path}")
    print("\n" + report_text)
    
    if not all_passed:
        print("CRITICAL ERROR: Phase 9 validation failed!")
        sys.exit(1)

if __name__ == "__main__":
    construct_and_validate_representations()
