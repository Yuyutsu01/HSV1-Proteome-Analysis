"""
Phase 4 - Step 06: Gene-Family Disjoint Evaluation & Performance Degradation Analysis.

Biological & Computational Concept:
1. Gene-Family Confound Experiment:
   In viral proteomics, proteins belonging to the same gene family (e.g. UL23 Thymidine Kinase,
   UL30 DNA Polymerase, ICP4 transcription factor) invariably belong to a single temporal class.
   If training and test sets contain sequences from the same gene families (even from different isolates),
   classifiers can achieve high accuracy by memorizing gene family boundaries rather than learning
   generalizable temporal regulatory sequence motifs.

2. Disjoint Evaluation Architecture:
   - Partitions whole gene/protein families into TRAIN, VALIDATION, and TEST sets.
   - Strict Group Invariant: Train Gene Families ∩ Test Gene Families = ∅.
   - Compares performance degradation across 3 increasingly strict evaluation regimes:
     Regime A: Random Stratified Split
     Regime B: Homology-Aware Split (70% sequence identity clusters)
     Regime C: Gene-Family-Disjoint Split (Whole gene families held out)

3. Scientific Interpretation:
   - If performance collapses under Gene-Family-Disjoint evaluation, it establishes that the predictive
     signal is strongly coupled to specific gene-family identity.
   - If performance retains discriminative power above random baselines, it provides evidence for
     family-independent temporal sequence signatures.

Creates: results/tables/phase4_gene_family_degradation.csv
"""

import os
import yaml
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from evaluation_utils import compute_comprehensive_metrics, CLASS_LABELS


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_gene_family_split(df_sup: pd.DataFrame, train_ratio: float = 0.70, test_ratio: float = 0.30, random_seed: int = 42) -> pd.Series:
    """
    Construct group-aware split where whole gene/protein families are assigned to either TRAIN or TEST.
    Guarantees zero gene-family overlap between train and test.
    """
    np.random.seed(random_seed)
    
    # Identify gene family identifier (use gene if available, else primary protein name prefix)
    genes = df_sup['gene'].fillna('UNKNOWN').values
    proteins = df_sup['protein'].fillna('UNKNOWN').values
    
    group_ids = []
    for g, p in zip(genes, proteins):
        if g != 'UNKNOWN':
            group_ids.append(f"GENE_{g}")
        else:
            # Group by protein name first token
            clean_p = p.split()[0].replace(',', '').upper()
            group_ids.append(f"PROT_{clean_p}")
            
    df_sup = df_sup.copy()
    df_sup['group_id'] = group_ids
    
    # Group class counts
    group_data = defaultdict(lambda: {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0, 'indices': []})
    for idx, row in df_sup.iterrows():
        grp = row['group_id']
        tc = row['temporal_class']
        group_data[grp]['indices'].append(idx)
        group_data[grp]['total'] += 1
        if tc == 'IMMEDIATE_EARLY':
            group_data[grp]['IE'] += 1
        elif tc == 'EARLY':
            group_data[grp]['Early'] += 1
        elif tc == 'LATE':
            group_data[grp]['Late'] += 1

    total_seqs = len(df_sup)
    target_train = int(total_seqs * train_ratio)
    target_test = total_seqs - target_train
    
    total_ie = sum(d['IE'] for d in group_data.values())
    total_early = sum(d['Early'] for d in group_data.values())
    total_late = sum(d['Late'] for d in group_data.values())
    
    target_dist = {
        'TRAIN': {'IE': total_ie * train_ratio, 'Early': total_early * train_ratio, 'Late': total_late * train_ratio, 'target_total': target_train},
        'TEST': {'IE': total_ie * test_ratio, 'Early': total_early * test_ratio, 'Late': total_late * test_ratio, 'target_total': target_test},
    }
    
    current_counts = {
        'TRAIN': {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0},
        'TEST': {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0},
    }
    
    # Sort groups by size descending
    sorted_groups = sorted(group_data.items(), key=lambda x: x[1]['total'], reverse=True)
    
    split_assignment = np.full(total_seqs, 'TRAIN', dtype=object)
    
    for grp_id, g_info in sorted_groups:
        best_split = None
        best_score = float('inf')
        
        for split in ['TRAIN', 'TEST']:
            current_tot = current_counts[split]['total']
            target_tot = target_dist[split]['target_total']
            capacity_ratio = (current_tot + g_info['total']) / (target_tot + 1e-5)
            
            new_ie = current_counts[split]['IE'] + g_info['IE']
            new_early = current_counts[split]['Early'] + g_info['Early']
            new_late = current_counts[split]['Late'] + g_info['Late']
            
            ie_def = max(0, target_dist[split]['IE'] - new_ie) / (target_dist[split]['IE'] + 1e-5)
            early_def = max(0, target_dist[split]['Early'] - new_early) / (target_dist[split]['Early'] + 1e-5)
            late_def = max(0, target_dist[split]['Late'] - new_late) / (target_dist[split]['Late'] + 1e-5)
            
            score = (capacity_ratio * 3.0) + (ie_def * 2.0) + (early_def * 1.5) + (late_def * 1.0)
            if score < best_score:
                best_score = score
                best_split = split
                
        current_counts[best_split]['IE'] += g_info['IE']
        current_counts[best_split]['Early'] += g_info['Early']
        current_counts[best_split]['Late'] += g_info['Late']
        current_counts[best_split]['total'] += g_info['total']
        
        for idx in g_info['indices']:
            split_assignment[idx] = best_split
            
    # Verification of zero group leakage
    train_groups = set(df_sup.loc[split_assignment == 'TRAIN', 'group_id'])
    test_groups = set(df_sup.loc[split_assignment == 'TEST', 'group_id'])
    overlap = train_groups.intersection(test_groups)
    assert len(overlap) == 0, f"Gene group leakage detected! Overlap: {overlap}"
    
    print(f"Gene-Family-Disjoint Split Created: Train={sum(split_assignment=='TRAIN')}, Test={sum(split_assignment=='TEST')} | Groups: Train={len(train_groups)}, Test={len(test_groups)}, Overlap={len(overlap)}")
    return pd.Series(split_assignment, index=df_sup.index)


def run_gene_family_evaluation():
    """Run gene family confound evaluation across representations."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 06: GENE-FAMILY DISJOINT & PERFORMANCE DEGRADATION ANALYSIS")
    print("=" * 80)

    # 1. Load Supervised Dataset and Split Manifest
    df_sup = pd.read_csv(config["dataset"]["supervised_path"])
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    
    # 2. Construct Gene-Family Disjoint Split
    df_sup['split_gene_disjoint'] = build_gene_family_split(df_sup, train_ratio=0.70, test_ratio=0.30, random_seed=config["random_seed"])
    
    # 3. Load Representations
    # A. AAC
    df_aac = pd.read_csv(config["features"]["aac_path"])
    aac_cols = [c for c in df_aac.columns if c not in ['canonical_id', 'temporal_class']]
    feats_aac = df_aac[aac_cols].values
    
    # B. Classical Descriptors
    df_class = pd.read_csv(config["features"]["classical_path"])
    class_cols = [c for c in df_class.columns if c not in ['canonical_id', 'temporal_class']]
    feats_class = df_class[class_cols].values
    
    # C. ESM-2
    esm_data = np.load(config["features"]["esm2_path"], allow_pickle=True)
    feats_esm2 = esm_data['embeddings']

    reps = {
        "AAC (20-dim)": feats_aac,
        "Classical (13-dim)": feats_class,
        "ESM-2 (320-dim)": feats_esm2
    }

    results = []

    # Evaluate each representation across Random -> Homology-Aware -> Gene-Family-Disjoint
    for rep_name, X_full in reps.items():
        print(f"\n--- Evaluating Degradation for {rep_name} ---")
        
        # Regime A: Random Stratified
        tr_rnd = (df_splits['split_random'] == 'TRAIN').values
        te_rnd = (df_splits['split_random'] == 'TEST').values
        scaler_rnd = StandardScaler()
        X_tr_rnd = scaler_rnd.fit_transform(X_full[tr_rnd])
        X_te_rnd = scaler_rnd.transform(X_full[te_rnd])
        y_tr_rnd = df_splits.loc[tr_rnd, 'temporal_class'].values
        y_te_rnd = df_splits.loc[te_rnd, 'temporal_class'].values
        
        lr_rnd = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
        lr_rnd.fit(X_tr_rnd, y_tr_rnd)
        y_pred_rnd = lr_rnd.predict(X_te_rnd)
        y_prob_rnd = lr_rnd.predict_proba(X_te_rnd)
        m_rnd = compute_comprehensive_metrics(y_te_rnd, y_pred_rnd, y_prob_rnd, rep_name, "Logistic Regression", "Random Stratified")
        results.append(m_rnd)

        # Regime B: Homology-Aware
        tr_hom = (df_splits['split_homology'] == 'TRAIN').values
        te_hom = (df_splits['split_homology'] == 'TEST').values
        scaler_hom = StandardScaler()
        X_tr_hom = scaler_hom.fit_transform(X_full[tr_hom])
        X_te_hom = scaler_hom.transform(X_full[te_hom])
        y_tr_hom = df_splits.loc[tr_hom, 'temporal_class'].values
        y_te_hom = df_splits.loc[te_hom, 'temporal_class'].values
        
        lr_hom = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
        lr_hom.fit(X_tr_hom, y_tr_hom)
        y_pred_hom = lr_hom.predict(X_te_hom)
        y_prob_hom = lr_hom.predict_proba(X_te_hom)
        m_hom = compute_comprehensive_metrics(y_te_hom, y_pred_hom, y_prob_hom, rep_name, "Logistic Regression", "Homology-Aware")
        results.append(m_hom)

        # Regime C: Gene-Family Disjoint
        tr_gd = (df_sup['split_gene_disjoint'] == 'TRAIN').values
        te_gd = (df_sup['split_gene_disjoint'] == 'TEST').values
        scaler_gd = StandardScaler()
        X_tr_gd = scaler_gd.fit_transform(X_full[tr_gd])
        X_te_gd = scaler_gd.transform(X_full[te_gd])
        y_tr_gd = df_sup.loc[tr_gd, 'temporal_class'].values
        y_te_gd = df_sup.loc[te_gd, 'temporal_class'].values
        
        lr_gd = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
        lr_gd.fit(X_tr_gd, y_tr_gd)
        y_pred_gd = lr_gd.predict(X_te_gd)
        y_prob_gd = lr_gd.predict_proba(X_te_gd)
        m_gd = compute_comprehensive_metrics(y_te_gd, y_pred_gd, y_prob_gd, rep_name, "Logistic Regression", "Gene-Family Disjoint")
        results.append(m_gd)

        print(f"  Random Stratified   -> MacroF1: {m_rnd['macro_f1']}, BalAcc: {m_rnd['balanced_accuracy']}, MCC: {m_rnd['mcc']}")
        print(f"  Homology-Aware      -> MacroF1: {m_hom['macro_f1']}, BalAcc: {m_hom['balanced_accuracy']}, MCC: {m_hom['mcc']}")
        print(f"  Gene-Family Disjoint -> MacroF1: {m_gd['macro_f1']}, BalAcc: {m_gd['balanced_accuracy']}, MCC: {m_gd['mcc']}")

    df_degradation = pd.DataFrame(results)
    out_path = config["paths"]["gene_family_degradation"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_degradation.to_csv(out_path, index=False)
    print(f"\nSaved Gene Family Degradation Table: {out_path}")
    return df_degradation


if __name__ == "__main__":
    run_gene_family_evaluation()
