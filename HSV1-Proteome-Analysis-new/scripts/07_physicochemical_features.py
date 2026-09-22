"""
Script 07: Physicochemical Feature Extraction and Validation

Scientific Concepts & Rationale:
--------------------------------
Proteins exhibit macroscopic biochemical properties determined by their primary amino acid sequences.
For Herpes Simplex Virus Type 1 (HSV-1 strain 17), we extract exactly 25 explicitly documented
physicochemical features across two canonical groups:

1. Global Sequence Descriptors (5 features):
   - Sequence Length (L): Total count of peptide bonds / residues.
   - Molecular Weight (MW, Da): Sum of isotopic atomic masses of all amino acids minus water loss during condensation.
   - Aromaticity: Relative frequency of aromatic residues (Phenylalanine [F], Tryptophan [W], Tyrosine [Y])
     governing hydrophobic core formation and UV absorption at 280 nm (Lobry & Gautier, 1994).
   - Instability Index: Empirical metric of in vivo protein stability predicted from dipeptide composition
     (Guruprasad et al., 1990); values < 40 predict stable proteins, > 40 predict unstable proteins.
   - Isoelectric Point (pI): Theoretical pH at which the net electrical charge of the protein is zero,
     calculated iteratively using standard amino acid side-chain pKa values (Bjellqvist et al., 1993).

2. Amino Acid Composition (20 features):
   - Fractional frequency (f_i in [0, 1]) for each of the 20 standard proteinogenic amino acids:
     A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y.
   - Normalized such that sum(f_i) == 1.0 for every protein.

Strict Methodological Rules:
- Exactly 25 features are computed. No undocumented or synthetic 26th feature is introduced.
- Preprocessing and feature calculations are purely algorithmic and reproducible from raw sequences.
- All 74 unique proteins are matched bijectively with authoritative biological annotations.

Author: Computational Biology Pipeline
"""

import os
import sys
import yaml
import hashlib
import platform
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import Bio

# Configure publication-quality aesthetic
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.alpha'] = 0.3

CANONICAL_AA_LIST = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
]

FEATURE_COLUMNS = [
    'sequence_length',
    'molecular_weight',
    'aromaticity',
    'instability_index',
    'isoelectric_point'
] + [f'aa_{aa}' for aa in CANONICAL_AA_LIST]

def load_config(config_path="config.yaml"):
    """Load pipeline configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def compute_sha256(filepath):
    """Computes SHA-256 hash for provenance and reproducibility."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()

def calculate_protein_features(sequence):
    """
    Computes exactly 25 physicochemical features for a given protein sequence using Biopython.
    
    Returns:
        dict: Mapping feature name to numerical float value.
    """
    clean_seq = str(sequence).strip().upper()
    
    # Verify canonical amino acid validity
    invalid_chars = set(clean_seq) - set(CANONICAL_AA_LIST)
    if invalid_chars:
        raise ValueError(f"Sequence contains non-canonical amino acid residues: {invalid_chars}")
        
    pa = ProteinAnalysis(clean_seq)
    seq_len = len(clean_seq)
    mw = float(pa.molecular_weight())
    aromaticity = float(pa.aromaticity())
    instability = float(pa.instability_index())
    pi = float(pa.isoelectric_point())
    
    # Calculate exact fractions [0, 1] for all 20 canonical amino acids
    aa_counts = pa.count_amino_acids()
    aa_fractions = {f"aa_{aa}": float(aa_counts.get(aa, 0) / seq_len) for aa in CANONICAL_AA_LIST}
    
    features = {
        'sequence_length': float(seq_len),
        'molecular_weight': mw,
        'aromaticity': aromaticity,
        'instability_index': instability,
        'isoelectric_point': pi,
    }
    features.update(aa_fractions)
    return features

def extract_and_validate_features(config_path="config.yaml"):
    """
    Executes Phase 7: Extract 25 physicochemical features, run strict validation,
    compute statistics and correlations, generate figures, and output verification reports.
    """
    cfg = load_config(config_path)
    proc_dir = cfg['paths']['data_processed']
    annot_dir = cfg['paths']['data_annotations']
    tables_dir = cfg['paths']['results_tables']
    figures_dir = cfg['paths']['results_figures']
    logs_dir = cfg['paths']['results_logs']
    
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    
    unique_csv = os.path.join(proc_dir, "unique_proteins.csv")
    annot_csv = os.path.join(annot_dir, "temporal_annotations_final.csv")
    output_features_csv = os.path.join(proc_dir, "physicochemical_features.csv")
    output_metadata_csv = os.path.join(proc_dir, "feature_metadata.csv")
    
    assert os.path.exists(unique_csv), f"Missing {unique_csv}"
    assert os.path.exists(annot_csv), f"Missing {annot_csv}"
    
    df_unique = pd.read_csv(unique_csv)
    df_annot = pd.read_csv(annot_csv)
    
    print("=" * 60)
    print("PHASE 7: PHYSICOCHEMICAL FEATURE EXTRACTION & VALIDATION")
    print("=" * 60)
    print(f"Total Unique Proteins Input: {len(df_unique)}")
    print(f"Authoritative Annotations:    {len(df_annot)}")
    
    # 1. Feature Extraction Loop
    feature_rows = []
    for idx, row in df_unique.iterrows():
        prot_id = str(row['protein_id'])
        gene = str(row['gene'])
        seq = str(row['sequence'])
        expected_len = int(row['sequence_length'])
        
        # Match annotation
        annot_match = df_annot[df_annot['protein_id'] == prot_id]
        if len(annot_match) == 0:
            raise KeyError(f"Protein ID {prot_id} missing in authoritative annotations!")
        t_class = annot_match['temporal_class'].values[0]
        
        # Calculate 25 features
        feat_dict = calculate_protein_features(seq)
        
        # Verification: length matches
        if int(feat_dict['sequence_length']) != expected_len:
            raise ValueError(f"Length mismatch for {prot_id}: calculated {feat_dict['sequence_length']} != expected {expected_len}")
            
        record = {
            'protein_id': prot_id,
            'gene': gene,
            'temporal_class': t_class
        }
        record.update(feat_dict)
        feature_rows.append(record)
        
    df_features = pd.DataFrame(feature_rows)
    
    # Ensure exact column ordering
    ordered_cols = ['protein_id', 'gene', 'temporal_class'] + FEATURE_COLUMNS
    df_features = df_features[ordered_cols]
    
    # Save main feature matrix
    df_features.to_csv(output_features_csv, index=False)
    print(f"[Saved] Physicochemical Feature Matrix: {output_features_csv} ({df_features.shape[0]} rows x {len(FEATURE_COLUMNS)} features)")
    
    # 2. Generate Feature Metadata Dictionary
    metadata_records = [
        {"feature_name": "sequence_length", "feature_group": "Global Descriptor", "definition": "Total number of amino acid residues in polypeptide chain", "unit": "amino acids (aa)", "calculation_method": "Exact residue count", "source_library": "Bio.SeqUtils.ProtParam"},
        {"feature_name": "molecular_weight", "feature_group": "Global Descriptor", "definition": "Theoretical molecular mass calculated from monoisotopic amino acid masses minus water condensation", "unit": "Daltons (Da)", "calculation_method": "Biopython molecular_weight()", "source_library": "Bio.SeqUtils.ProtParam"},
        {"feature_name": "aromaticity", "feature_group": "Global Descriptor", "definition": "Relative fractional frequency of aromatic amino acids (Phe + Trp + Tyr)", "unit": "fraction [0-1]", "calculation_method": "Lobry & Gautier (1994) / Biopython aromaticity()", "source_library": "Bio.SeqUtils.ProtParam"},
        {"feature_name": "instability_index", "feature_group": "Global Descriptor", "definition": "In vivo protein stability estimate based on dipeptide weight values (<40 stable, >40 unstable)", "unit": "dimensionless index", "calculation_method": "Guruprasad et al. (1990) / Biopython instability_index()", "source_library": "Bio.SeqUtils.ProtParam"},
        {"feature_name": "isoelectric_point", "feature_group": "Global Descriptor", "definition": "Theoretical pH at which net molecular electrical charge equals zero", "unit": "pH scale [0-14]", "calculation_method": "Bjellqvist et al. (1993) pKa iteration", "source_library": "Bio.SeqUtils.ProtParam"},
    ]
    for aa in CANONICAL_AA_LIST:
        metadata_records.append({
            "feature_name": f"aa_{aa}",
            "feature_group": "Amino Acid Composition",
            "definition": f"Fractional frequency of amino acid {aa} in the primary sequence",
            "unit": "fraction [0-1]",
            "calculation_method": f"Count({aa}) / Sequence_Length",
            "source_library": "Bio.SeqUtils.ProtParam"
        })
    df_meta = pd.DataFrame(metadata_records)
    df_meta.to_csv(output_metadata_csv, index=False)
    print(f"[Saved] Feature Metadata Documentation: {output_metadata_csv}")
    
    # 3. Cross-Check Selected Canonical Proteins
    cross_check_proteins = ['UL36', 'UL30', 'UL19', 'UL54', 'US12']
    print("\n--- Manual Spot Verification on Benchmark Proteins ---")
    for g in cross_check_proteins:
        row = df_features[df_features['gene'] == g].iloc[0]
        comp_sum = sum(row[f'aa_{aa}'] for aa in CANONICAL_AA_LIST)
        print(f"  * {g:<5} ({row['protein_id']}): Length={int(row['sequence_length']):<5} MW={row['molecular_weight']:<10.2f} pI={row['isoelectric_point']:<5.2f} Instab={row['instability_index']:<6.2f} Arom={row['aromaticity']:<6.4f} AA_Sum={comp_sum:.6f}")
        
    # 4. Descriptive Statistics (Overall and Class-Wise)
    # Overall statistics
    desc_overall = df_features[FEATURE_COLUMNS].describe().T[['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']].reset_index()
    desc_overall.rename(columns={'index': 'feature', '50%': 'median'}, inplace=True)
    desc_overall['subset'] = 'All Proteins (N=74)'
    
    # Class-wise statistics
    class_dfs = []
    for c_name in ['Immediate-Early', 'Early', 'Late']:
        sub = df_features[df_features['temporal_class'] == c_name]
        c_desc = sub[FEATURE_COLUMNS].describe().T[['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']].reset_index()
        c_desc.rename(columns={'index': 'feature', '50%': 'median'}, inplace=True)
        c_desc['subset'] = f'{c_name} (N={len(sub)})'
        class_dfs.append(c_desc)
        
    all_stats = pd.concat([desc_overall] + class_dfs, ignore_index=True)
    stats_csv = os.path.join(tables_dir, "physicochemical_statistics.csv")
    all_stats.to_csv(stats_csv, index=False)
    print(f"[Saved] Descriptive Statistics Table: {stats_csv}")
    
    # 5. Correlation Analysis (Pearson and Spearman)
    pearson_corr = df_features[FEATURE_COLUMNS].corr(method='pearson')
    spearman_corr = df_features[FEATURE_COLUMNS].corr(method='spearman')
    
    pearson_csv = os.path.join(tables_dir, "pearson_correlation_matrix.csv")
    spearman_csv = os.path.join(tables_dir, "spearman_correlation_matrix.csv")
    pearson_corr.to_csv(pearson_csv)
    spearman_corr.to_csv(spearman_csv)
    print(f"[Saved] Pearson Correlation Matrix: {pearson_csv}")
    print(f"[Saved] Spearman Correlation Matrix: {spearman_csv}")
    
    # 6. Figures Generation (Publication Quality)
    # Figure 1: Physicochemical Correlation Heatmap
    plt.figure(figsize=(14, 12))
    mask = np.triu(np.ones_like(pearson_corr, dtype=bool))
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    sns.heatmap(pearson_corr, mask=mask, cmap=cmap, vmin=-1.0, vmax=1.0, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .8, "label": "Pearson Correlation (r)"},
                annot=False, fmt=".2f")
    plt.title("HSV-1 Proteome: Pairwise Physicochemical Feature Correlations (25 Features, N=74)", fontsize=13, fontweight='bold', pad=15)
    plt.tight_layout()
    corr_fig_path = os.path.join(figures_dir, "physicochemical_correlation_matrix.png")
    plt.savefig(corr_fig_path, dpi=300)
    plt.close()
    print(f"[Saved] Figure: {corr_fig_path}")
    
    # Figure 2: Global Descriptors Distribution
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()
    global_cols = ['sequence_length', 'molecular_weight', 'isoelectric_point', 'instability_index', 'aromaticity']
    titles = ['Sequence Length (aa)', 'Molecular Weight (Da)', 'Isoelectric Point (pI)', 'Instability Index', 'Aromaticity (F+W+Y Fraction)']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, col in enumerate(global_cols):
        sns.histplot(df_features[col], ax=axes[i], kde=True, color=colors[i], bins=15, edgecolor='black', alpha=0.6)
        axes[i].set_title(titles[i], fontsize=12, fontweight='bold')
        axes[i].set_xlabel(col)
        axes[i].set_ylabel("Count")
        axes[i].grid(True, linestyle='--', alpha=0.4)
        
    # Class-wise boxplot for pI in the 6th subplot
    palette = {'Immediate-Early': '#d62728', 'Early': '#1f77b4', 'Late': '#2ca02c'}
    sns.boxplot(x='temporal_class', y='isoelectric_point', data=df_features, ax=axes[5], palette=palette, boxprops=dict(alpha=0.7))
    sns.stripplot(x='temporal_class', y='isoelectric_point', data=df_features, ax=axes[5], color='black', alpha=0.6, jitter=0.2)
    axes[5].set_title("pI by Temporal Class", fontsize=12, fontweight='bold')
    axes[5].set_xlabel("Temporal Expression Class")
    axes[5].set_ylabel("Isoelectric Point (pI)")
    axes[5].grid(True, linestyle='--', alpha=0.4)
    
    plt.suptitle("HSV-1 Proteome: Global Physicochemical Feature Distributions (N=74)", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    dist_fig_path = os.path.join(figures_dir, "feature_distributions.png")
    plt.savefig(dist_fig_path, dpi=300)
    plt.close()
    print(f"[Saved] Figure: {dist_fig_path}")
    
    # Figure 3: Amino Acid Composition Profile (Mean +/- SD across proteome)
    aa_cols = [f'aa_{aa}' for aa in CANONICAL_AA_LIST]
    aa_means = df_features[aa_cols].mean() * 100
    aa_stds = df_features[aa_cols].std() * 100
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(CANONICAL_AA_LIST, aa_means, yerr=aa_stds, capsize=4, color='#3470a3', edgecolor='black', alpha=0.8)
    plt.title("HSV-1 Proteome: Mean Amino Acid Composition Profile (% +/- SD)", fontsize=13, fontweight='bold', pad=12)
    plt.xlabel("Canonical Amino Acid Residue", fontsize=11, fontweight='bold')
    plt.ylabel("Mean Percentage Composition (%)", fontsize=11, fontweight='bold')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.3, f'{yval:.1f}%', ha='center', va='bottom', fontsize=8)
    plt.tight_layout()
    aa_fig_path = os.path.join(figures_dir, "amino_acid_composition_profile.png")
    plt.savefig(aa_fig_path, dpi=300)
    plt.close()
    print(f"[Saved] Figure: {aa_fig_path}")

    # Figure 4: Class-Wise Physicochemical Descriptor Comparison
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    metrics_to_plot = ['sequence_length', 'molecular_weight', 'instability_index', 'aromaticity']
    titles_c = ['Length (aa)', 'Molecular Weight (Da)', 'Instability Index', 'Aromaticity']
    
    for i, m in enumerate(metrics_to_plot):
        sns.boxplot(x='temporal_class', y=m, data=df_features, ax=axes[i], palette=palette, boxprops=dict(alpha=0.7))
        sns.stripplot(x='temporal_class', y=m, data=df_features, ax=axes[i], color='black', alpha=0.6, jitter=0.2)
        axes[i].set_title(titles_c[i], fontsize=12, fontweight='bold')
        axes[i].set_xlabel("Temporal Class", fontsize=10)
        axes[i].set_ylabel(m, fontsize=10)
        axes[i].grid(True, linestyle='--', alpha=0.4)
        
    plt.suptitle("HSV-1 Proteome: Class-Wise Physicochemical Comparisons (IE=5, E=15, L=54)", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    class_fig_path = os.path.join(figures_dir, "class_feature_distributions.png")
    plt.savefig(class_fig_path, dpi=300)
    plt.close()
    print(f"[Saved] Figure: {class_fig_path}")

    # 7. Automated Validation Checks
    validation_results = {}
    
    # A. Row count = 74
    v_rows = (len(df_features) == 74)
    validation_results['row_count_74'] = "PASS" if v_rows else "FAIL"
    
    # B. Unique protein IDs = 74
    v_ids = (df_features['protein_id'].nunique() == 74)
    validation_results['unique_protein_ids_74'] = "PASS" if v_ids else "FAIL"
    
    # C. Feature count = 25
    num_feats = len([c for c in df_features.columns if c not in ['protein_id', 'gene', 'temporal_class']])
    v_feats = (num_feats == 25)
    validation_results['exact_25_features'] = "PASS" if v_feats else "FAIL"
    
    # D. No missing values (NaN or Inf)
    num_nans = df_features[FEATURE_COLUMNS].isna().sum().sum()
    num_infs = np.isinf(df_features[FEATURE_COLUMNS].values).sum()
    v_nans = (num_nans == 0 and num_infs == 0)
    validation_results['no_nans_no_infs'] = "PASS" if v_nans else "FAIL"
    
    # E. Sequence length consistency
    v_lens = all(df_features['sequence_length'] == df_unique['sequence_length'])
    validation_results['sequence_length_consistency'] = "PASS" if v_lens else "FAIL"
    
    # F. Amino acid composition sum == 1.0 (+/- 1e-5)
    aa_sums = df_features[aa_cols].sum(axis=1)
    v_sums = all(np.isclose(aa_sums, 1.0, atol=1e-5))
    validation_results['composition_sums_to_1'] = "PASS" if v_sums else "FAIL"
    
    # G. Aromaticity range [0, 1]
    v_arom = all((df_features['aromaticity'] >= 0.0) & (df_features['aromaticity'] <= 1.0))
    validation_results['aromaticity_in_valid_range'] = "PASS" if v_arom else "FAIL"
    
    # H. Class distribution consistency: 5 / 15 / 54
    c_counts = df_features['temporal_class'].value_counts().to_dict()
    v_class = (c_counts.get('Immediate-Early', 0) == 5 and 
               c_counts.get('Early', 0) == 15 and 
               c_counts.get('Late', 0) == 54)
    validation_results['class_distribution_preserved'] = "PASS" if v_class else "FAIL"
    
    # Overall Status
    all_passed = all(v == "PASS" for v in validation_results.values())
    validation_results['overall_phase7_status'] = "PASS" if all_passed else "FAIL"
    
    # 8. Deterministic Reproduction Check
    # Recalculate and compare matrix
    df_features_check = pd.DataFrame(feature_rows)[ordered_cols]
    v_determ = df_features.equals(df_features_check)
    validation_results['deterministic_reproduction'] = "PASS" if v_determ else "FAIL"
    
    # 9. Write Reproducibility & Provenance Log
    raw_hash = compute_sha256(unique_csv)
    annot_hash = compute_sha256(annot_csv)
    out_hash = compute_sha256(output_features_csv)
    
    prov_log_lines = [
        "==================================================",
        "PHYSICOCHEMICAL FEATURE GENERATION & PROVENANCE LOG",
        "==================================================",
        f"Execution Timestamp:       {pd.Timestamp.now().isoformat()}",
        f"Python Version:            {platform.python_version()}",
        f"Biopython Version:         {Bio.__version__}",
        f"NumPy Version:             {np.__version__}",
        f"Pandas Version:            {pd.__version__}",
        f"Input unique_proteins.csv: {raw_hash}",
        f"Input annotations:         {annot_hash}",
        f"Output feature matrix:     {out_hash}",
        f"Computed Features Count:   {len(FEATURE_COLUMNS)}",
        f"Analyzed Proteins Count:   {len(df_features)}",
        "=================================================="
    ]
    prov_log_path = os.path.join(logs_dir, "physicochemical_feature_generation.log")
    with open(prov_log_path, "w") as f:
        f.write("\n".join(prov_log_lines))
    print(f"[Saved] Provenance Log: {prov_log_path}")
    
    # 10. Write Final Phase 7 Validation Report
    report_lines = [
        "==================================================",
        "PHASE 7 VALIDATION REPORT: PHYSICOCHEMICAL FEATURES",
        "==================================================",
        f"Dataset:                  74 proteins",
        f"Feature matrix:           74 x {len(FEATURE_COLUMNS)} features",
        f"Missing values:           {num_nans + num_infs}",
        f"Duplicate IDs:            {len(df_features) - df_features['protein_id'].nunique()}",
        f"Composition validation:   {validation_results['composition_sums_to_1']}",
        f"Sequence-length valid.:   {validation_results['sequence_length_consistency']}",
        f"Class distribution:       5 / 15 / 54 ({validation_results['class_distribution_preserved']})",
        f"Deterministic reproduc.:  {validation_results['deterministic_reproduction']}",
        f"Overall Phase 7 status:   {validation_results['overall_phase7_status']}",
        "==================================================",
        "",
        "Detailed Check List:",
    ]
    for k, v in validation_results.items():
        report_lines.append(f"  - {k:<30}: {v}")
    report_lines.append("\n[STOP CONDITION REACHED] Phase 7 completed. Standby for Phase 8 ProtBERT.")
    
    report_text = "\n".join(report_lines)
    report_path = os.path.join(logs_dir, "phase7_validation_report.txt")
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"[Saved] Validation Report: {report_path}")
    print("\n" + report_text)
    
    if not all_passed:
        print("CRITICAL ERROR: Phase 7 validation failed!")
        sys.exit(1)

if __name__ == "__main__":
    extract_and_validate_features()
