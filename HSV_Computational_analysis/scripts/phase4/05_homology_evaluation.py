#!/usr/bin/env python3
"""
Phase 4 — Homology-Aware Evaluation & Comparative Degradation
Script: scripts/phase4/05_homology_evaluation.py

Purpose:
Collects, verifies, and analyzes model performance under the Homology-Aware partition regime (CD-HIT 70% clusters).
Calculates performance degradation: Random Stratified -> Homology-Aware.
Ensures zero cluster leakage between partitions:
  Train ∩ Validation = 0
  Train ∩ Test = 0
  Validation ∩ Test = 0

Outputs:
  - results/tables/phase4_homology_evaluation.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(PROJECT_ROOT))

def main():
    print("=" * 70)
    print("Phase 4: Homology-Aware Partition Evaluation & Degradation Analysis")
    print("=" * 70)

    # 1. Verify split cluster isolation
    split_manifest_path = PROJECT_ROOT / "data" / "processed" / "phase3_split_manifest.csv"
    if not split_manifest_path.exists():
        raise FileNotFoundError(f"Missing split manifest: {split_manifest_path}")

    df_splits = pd.read_csv(split_manifest_path)
    
    # Check homology cluster disjointness
    hom_train_c = set(df_splits[df_splits["split_homology"] == "TRAIN"]["homology_cluster_id"])
    hom_val_c = set(df_splits[df_splits["split_homology"] == "VALIDATION"]["homology_cluster_id"])
    hom_test_c = set(df_splits[df_splits["split_homology"] == "TEST"]["homology_cluster_id"])

    overlap_train_val = len(hom_train_c.intersection(hom_val_c))
    overlap_train_test = len(hom_train_c.intersection(hom_test_c))
    overlap_val_test = len(hom_val_c.intersection(hom_test_c))

    print(f"Homology Cluster Partition Isolation:")
    print(f"  Train clusters: {len(hom_train_c)}")
    print(f"  Validation clusters: {len(hom_val_c)}")
    print(f"  Test clusters: {len(hom_test_c)}")
    print(f"  Train intersect Val Overlap: {overlap_train_val} (Expected: 0)")
    print(f"  Train intersect Test Overlap: {overlap_train_test} (Expected: 0)")
    print(f"  Val intersect Test Overlap: {overlap_val_test} (Expected: 0)")

    if overlap_train_val > 0 or overlap_train_test > 0 or overlap_val_test > 0:
        raise ValueError("CRITICAL LEAKAGE: Homology clusters overlap across partitions!")

    # 2. Collect evaluation tables
    tables_dir = PROJECT_ROOT / "results" / "tables" / "phase4"
    results = []

    # Baselines
    baseline_file = tables_dir / "baseline_models_results.csv"
    if baseline_file.exists():
        df_base = pd.read_csv(baseline_file)
        results.append(df_base)

    # Classical
    classical_file = tables_dir / "classical_models_results.csv"
    if classical_file.exists():
        df_class = pd.read_csv(classical_file)
        results.append(df_class)

    # ESM2
    esm2_file = tables_dir / "esm2_models_results.csv"
    if esm2_file.exists():
        df_esm = pd.read_csv(esm2_file)
        results.append(df_esm)

    if not results:
        print("No evaluation results files found yet. Run baseline, classical, and esm2 scripts first.")
        return

    df_all = pd.concat(results, ignore_index=True)
    
    # Filter homology-aware results
    df_hom = df_all[df_all["split_regime"] == "Homology-Aware"].copy()

    # Calculate degradation against random stratified where model and representation match
    df_rand = df_all[df_all["split_regime"] == "Random Stratified"].copy()
    
    comparison_rows = []
    for _, h_row in df_hom.iterrows():
        rep = h_row["representation"]
        mod = h_row["model"]
        
        # Find matching random row
        r_match = df_rand[(df_rand["representation"] == rep) & (df_rand["model"] == mod)]
        if not r_match.empty:
            r_row = r_match.iloc[0]
            rand_macro_f1 = r_row["macro_f1"]
            hom_macro_f1 = h_row["macro_f1"]
            delta_f1 = hom_macro_f1 - rand_macro_f1
            rel_drop = (delta_f1 / rand_macro_f1) * 100 if rand_macro_f1 > 0 else 0
            
            rand_mcc = r_row["mcc"]
            hom_mcc = h_row["mcc"]
            delta_mcc = hom_mcc - rand_mcc
            
            comparison_rows.append({
                "representation": rep,
                "model": mod,
                "random_macro_f1": rand_macro_f1,
                "homology_macro_f1": hom_macro_f1,
                "macro_f1_delta": delta_f1,
                "macro_f1_rel_change_pct": rel_drop,
                "random_mcc": rand_mcc,
                "homology_mcc": hom_mcc,
                "mcc_delta": delta_mcc,
                "homology_accuracy": h_row["accuracy"],
                "homology_balanced_accuracy": h_row["balanced_accuracy"],
                "homology_ie_f1": h_row["ie_f1"],
                "homology_early_f1": h_row["early_f1"],
                "homology_late_f1": h_row["late_f1"]
            })

    df_comp = pd.DataFrame(comparison_rows)
    out_path = PROJECT_ROOT / "results" / "tables" / "phase4_homology_evaluation.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_comp.to_csv(out_path, index=False)
    print(f"\nSaved homology-aware evaluation comparison to: {out_path}")
    print(df_comp[["representation", "model", "random_macro_f1", "homology_macro_f1", "macro_f1_rel_change_pct", "homology_mcc"]].to_string(index=False))

if __name__ == "__main__":
    main()
