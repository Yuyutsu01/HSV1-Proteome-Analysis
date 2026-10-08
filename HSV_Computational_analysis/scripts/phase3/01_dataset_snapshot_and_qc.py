#!/usr/bin/env python3
"""
Phase 3 - Step 01: Dataset Snapshot, Sequence Quality & Length Analysis
=======================================================================

Author: Shiva & Antigravity IDE
Date: October 2026
Project: HSV_Computational_analysis

Scientific Purpose:
-------------------
1. Deterministically verify sequence normalization (whitespace/formatting cleanup)
   while ensuring zero residue modification (data/processed/phase3_sequence_manifest.csv).
2. Compute comprehensive dataset snapshot statistics directly from the frozen supervised cohort (N=16,657).
3. Evaluate sequence-level QC properties (amino acid frequencies, non-standard residues, entropy, hydrophobicity).
4. Perform rigorous length distribution statistics and visualizations across temporal classes (IE, Early, Late).
"""

import os
import sys
import yaml
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def run_step_01():
    print("=" * 80)
    print("PHASE 3 - STEP 01: DATASET SNAPSHOT, SEQUENCE QUALITY & LENGTH ANALYSIS")
    print("=" * 80)

    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    config_path = os.path.join(base_dir, "configs", "phase3_config.yaml")
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    # Output paths
    processed_dir = os.path.join(base_dir, cfg['paths']['processed_dir'])
    tables_dir = os.path.join(base_dir, cfg['paths']['tables_dir'])
    figures_dir = os.path.join(base_dir, cfg['paths']['figures_dir'])
    os.makedirs(processed_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    # 1. Load Frozen Supervised Dataset (N=16,657)
    sup_path = os.path.join(base_dir, cfg['dataset']['supervised_path'])
    df = pd.read_csv(sup_path)
    total_seqs = len(df)
    assert total_seqs == 16657, f"Expected 16,657 supervised sequences, found {total_seqs}"

    # 2. Sequence Normalization Manifest
    manifest_rows = []
    for _, row in df.iterrows():
        cid = row['canonical_id']
        orig_seq = row['sequence']
        norm_seq = orig_seq.strip().upper()
        orig_len = len(orig_seq)
        norm_len = len(norm_seq)
        changes = 'NONE' if orig_seq == norm_seq else 'WHITESPACE_TRIMMED'
        manifest_rows.append({
            'canonical_id': cid,
            'original_sequence': orig_seq,
            'normalized_sequence': norm_seq,
            'original_length': orig_len,
            'normalized_length': norm_len,
            'normalization_changes': changes
        })

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = os.path.join(processed_dir, "phase3_sequence_manifest.csv")
    manifest_df.to_csv(manifest_path, index=False)
    print(f"Saved Sequence Manifest: {manifest_path} (N={len(manifest_df)})")

    # 3. Dataset Snapshot Table (results/tables/phase3_dataset_snapshot.csv)
    c_ie = (df['temporal_class'] == 'IMMEDIATE_EARLY').sum()
    c_early = (df['temporal_class'] == 'EARLY').sum()
    c_late = (df['temporal_class'] == 'LATE').sum()
    c_hsv1 = (df['species'] == 'HSV-1').sum()
    c_hsv2 = (df['species'] == 'HSV-2').sum()
    c_genes = df['gene'].nunique()
    c_prots = df['protein'].nunique()

    lens = df['length']
    l_min = lens.min()
    l_max = lens.max()
    l_mean = lens.mean()
    l_med = lens.median()
    l_std = lens.std()
    l_q1 = lens.quantile(0.25)
    l_q3 = lens.quantile(0.75)

    snapshot_rows = [{
        'total_sequences': total_seqs,
        'IE_count': c_ie,
        'Early_count': c_early,
        'Late_count': c_late,
        'class_fraction': f"IE: {c_ie/total_seqs:.4f}, Early: {c_early/total_seqs:.4f}, Late: {c_late/total_seqs:.4f}",
        'HSV1_count': c_hsv1,
        'HSV2_count': c_hsv2,
        'unique_genes': c_genes,
        'unique_proteins': c_prots,
        'length_min': int(l_min),
        'length_max': int(l_max),
        'length_mean': round(float(l_mean), 2),
        'length_median': round(float(l_med), 2),
        'length_std': round(float(l_std), 2),
        'length_q1': round(float(l_q1), 2),
        'length_q3': round(float(l_q3), 2)
    }]
    snapshot_df = pd.DataFrame(snapshot_rows)
    snapshot_path = os.path.join(tables_dir, "phase3_dataset_snapshot.csv")
    snapshot_df.to_csv(snapshot_path, index=False)
    print(f"Saved Dataset Snapshot: {snapshot_path}")

    # 4. Length Statistics by Temporal Class (results/tables/phase3_length_statistics.csv)
    length_stat_rows = []
    for tc in ['IMMEDIATE_EARLY', 'EARLY', 'LATE', 'ALL']:
        sub_lens = df['length'] if tc == 'ALL' else df[df['temporal_class'] == tc]['length']
        length_stat_rows.append({
            'temporal_class': tc,
            'n': len(sub_lens),
            'mean': round(float(sub_lens.mean()), 2),
            'median': round(float(sub_lens.median()), 2),
            'std': round(float(sub_lens.std()), 2),
            'min': int(sub_lens.min()),
            'max': int(sub_lens.max()),
            'Q1': round(float(sub_lens.quantile(0.25)), 2),
            'Q3': round(float(sub_lens.quantile(0.75)), 2)
        })
    length_stat_df = pd.DataFrame(length_stat_rows)
    length_stat_path = os.path.join(tables_dir, "phase3_length_statistics.csv")
    length_stat_df.to_csv(length_stat_path, index=False)
    print(f"Saved Length Statistics Table: {length_stat_path}")

    # 5. Sequence Quality & Composition Analysis (results/tables/phase3_sequence_quality.csv)
    standard_aas = set(cfg['features']['standard_amino_acids'])
    hydrophobic_aas = {'A', 'I', 'L', 'M', 'F', 'W', 'V', 'P', 'G'}
    pos_charge_aas = {'R', 'K', 'H'}
    neg_charge_aas = {'D', 'E'}

    quality_rows = []
    for _, row in df.iterrows():
        seq = row['sequence'].upper()
        L = len(seq)
        if L == 0:
            continue
        
        # Composition and entropy
        counts = {}
        non_std_cnt = 0
        for aa in seq:
            counts[aa] = counts.get(aa, 0) + 1
            if aa not in standard_aas:
                non_std_cnt += 1
        
        # Shannon entropy
        entropy = 0.0
        for aa, cnt in counts.items():
            p = cnt / L
            if p > 0:
                entropy -= p * math.log2(p)

        # Hydrophobicity and charge fractions
        hydro_cnt = sum(counts.get(aa, 0) for aa in hydrophobic_aas)
        pos_cnt = sum(counts.get(aa, 0) for aa in pos_charge_aas)
        neg_cnt = sum(counts.get(aa, 0) for aa in neg_charge_aas)

        quality_rows.append({
            'canonical_id': row['canonical_id'],
            'length': L,
            'temporal_class': row['temporal_class'],
            'non_standard_residues': non_std_cnt,
            'non_standard_fraction': round(non_std_cnt / L, 6),
            'shannon_entropy': round(entropy, 4),
            'hydrophobic_fraction': round(hydro_cnt / L, 4),
            'positive_charge_fraction': round(pos_cnt / L, 4),
            'negative_charge_fraction': round(neg_cnt / L, 4),
            'net_charge': pos_cnt - neg_cnt
        })

    quality_df = pd.DataFrame(quality_rows)
    quality_path = os.path.join(tables_dir, "phase3_sequence_quality.csv")
    quality_df.to_csv(quality_path, index=False)
    print(f"Saved Sequence Quality Table: {quality_path} (N={len(quality_df)})")

    # 6. Generate Figures: Length Distributions
    sns.set_theme(style="whitegrid", palette="muted")
    
    # Figure 1: Overall Length Distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(df['length'], bins=60, kde=True, color='royalblue', edgecolor='black', alpha=0.7)
    plt.title("Sequence Length Distribution across Supervised Dataset (N=16,657)", fontsize=14, fontweight='bold')
    plt.xlabel("Sequence Length (Amino Acids)", fontsize=12)
    plt.ylabel("Sequence Count", fontsize=12)
    plt.yscale('log')
    plt.axvline(l_median := df['length'].median(), color='red', linestyle='--', label=f'Median Length ({int(l_median)} AA)')
    plt.legend(fontsize=11)
    plt.tight_layout()
    fig1_path = os.path.join(figures_dir, "phase3_sequence_length_distribution.png")
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"Saved Figure: {fig1_path}")

    # Figure 2: Length by Temporal Class Boxplot & Violin Plot
    plt.figure(figsize=(10, 6))
    order = ['IMMEDIATE_EARLY', 'EARLY', 'LATE']
    palette = {'IMMEDIATE_EARLY': '#E69F00', 'EARLY': '#56B4E9', 'LATE': '#009E73'}
    sns.boxplot(x='temporal_class', y='length', data=df, order=order, palette=palette, showfliers=False, width=0.5)
    sns.stripplot(x='temporal_class', y='length', data=df.sample(min(1500, len(df)), random_state=42),
                  order=order, color='black', alpha=0.15, jitter=0.2, size=3)
    plt.title("Protein Sequence Length Distribution by Temporal Kinetic Class", fontsize=14, fontweight='bold')
    plt.xlabel("Temporal Expression Class", fontsize=12)
    plt.ylabel("Sequence Length (Amino Acids)", fontsize=12)
    plt.tight_layout()
    fig2_path = os.path.join(figures_dir, "phase3_length_by_temporal_class.png")
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"Saved Figure: {fig2_path}")

    print("Step 01 completed successfully.\n")

if __name__ == "__main__":
    run_step_01()
