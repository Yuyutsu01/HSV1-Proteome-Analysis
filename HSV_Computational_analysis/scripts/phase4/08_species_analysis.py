"""
Phase 4 - Step 08: Species-Stratified and Cross-Species Transfer Analysis.

Biological & Computational Concept:
1. Species Invariance vs Divergence:
   HSV-1 and HSV-2 share extensive sequence homology across orthologous gene products,
   yet exhibit genomic divergence and differential tissue tropism.
   Does temporal prediction survive when evaluated exclusively within HSV-1, exclusively
   within HSV-2, and across species (HSV-1 -> HSV-2 and HSV-2 -> HSV-1)?

2. Experimental Regimes (Section 10):
   - Regime A: Combined HSV-1 + HSV-2 (Standard)
   - Regime B: HSV-1 Only (Train on HSV-1, Test on HSV-1)
   - Regime C: HSV-2 Only (Train on HSV-2, Test on HSV-2)
   - Regime D: Cross-Species Transfer (Train HSV-1 -> Test HSV-2; Train HSV-2 -> Test HSV-1)

3. Representation: ESM-2 embeddings (320-dim) and AAC (20-dim) with Logistic Regression & Random Forest.
   Scaling fitted strictly on training partition.

Creates: results/tables/phase4_species_analysis.csv
"""

import os
import yaml
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from evaluation_utils import compute_comprehensive_metrics, CLASS_LABELS


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_species_analysis():
    """Execute within-species and cross-species temporal modeling experiments."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 08: SPECIES-STRATIFIED & CROSS-SPECIES TRANSFER ANALYSIS")
    print("=" * 80)

    # 1. Load Supervised Dataset, Split Manifest, and ESM-2 Embeddings
    df_sup = pd.read_csv(config["dataset"]["supervised_path"])
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    esm_data = np.load(config["features"]["esm2_path"], allow_pickle=True)
    X_full = esm_data['embeddings'].astype(np.float32)

    df = df_splits.merge(df_sup[['canonical_id', 'species']], on='canonical_id')

    # Species masks
    hsv1_mask = (df['species'] == 'HSV-1').values
    hsv2_mask = (df['species'] == 'HSV-2').values

    print(f"Species Counts: HSV-1={sum(hsv1_mask)}, HSV-2={sum(hsv2_mask)}")

    results = []

    # -------------------------------------------------------------
    # 1. Within-Species Evaluation (Using Random Stratified Splits)
    # -------------------------------------------------------------
    experiments = [
        ("Combined (HSV-1 + HSV-2)", np.ones(len(df), dtype=bool), np.ones(len(df), dtype=bool)),
        ("Within HSV-1 Only", hsv1_mask, hsv1_mask),
        ("Within HSV-2 Only", hsv2_mask, hsv2_mask),
        ("Cross-Species (Train HSV-1 -> Test HSV-2)", hsv1_mask, hsv2_mask),
        ("Cross-Species (Train HSV-2 -> Test HSV-1)", hsv2_mask, hsv1_mask),
    ]

    for exp_name, tr_filter, te_filter in experiments:
        print(f"\n--- Running Experiment: {exp_name} ---")
        
        # Use random split partition as base, then intersect with species filters
        if "Cross-Species" in exp_name:
            # For cross-species transfer, train on ALL training sequences of Species A, test on ALL test sequences of Species B
            tr_idx = (df['split_random'] == 'TRAIN').values & tr_filter
            te_idx = (df['split_random'] == 'TEST').values & te_filter
        else:
            tr_idx = (df['split_random'] == 'TRAIN').values & tr_filter
            te_idx = (df['split_random'] == 'TEST').values & te_filter

        X_train = X_full[tr_idx]
        y_train = df.loc[tr_idx, 'temporal_class'].values
        X_test = X_full[te_idx]
        y_test = df.loc[te_idx, 'temporal_class'].values

        print(f"  Train N={len(X_train)} | Test N={len(X_test)}")

        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X_train)
        X_te_s = scaler.transform(X_test)

        # Logistic Regression
        lr = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
        lr.fit(X_tr_s, y_train)
        y_pred_lr = lr.predict(X_te_s)
        y_prob_lr = lr.predict_proba(X_te_s)
        m_lr = compute_comprehensive_metrics(y_test, y_pred_lr, y_prob_lr, "ESM-2 (320-dim)", "Logistic Regression", exp_name)
        results.append(m_lr)

        # Random Forest
        rf = RandomForestClassifier(n_estimators=100, max_depth=15, random_state=config["random_seed"], class_weight='balanced', n_jobs=-1)
        rf.fit(X_tr_s, y_train)
        y_pred_rf = rf.predict(X_te_s)
        y_prob_rf = rf.predict_proba(X_te_s)
        m_rf = compute_comprehensive_metrics(y_test, y_pred_rf, y_prob_rf, "ESM-2 (320-dim)", "Random Forest", exp_name)
        results.append(m_rf)

        print(f"  LR -> Acc: {m_lr['accuracy']:.4f} | BalAcc: {m_lr['balanced_accuracy']:.4f} | MacroF1: {m_lr['macro_f1']:.4f} | MCC: {m_lr['mcc']:.4f}")
        print(f"  RF -> Acc: {m_rf['accuracy']:.4f} | BalAcc: {m_rf['balanced_accuracy']:.4f} | MacroF1: {m_rf['macro_f1']:.4f} | MCC: {m_rf['mcc']:.4f}")

    df_species = pd.DataFrame(results)
    out_path = config["paths"]["species_analysis"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_species.to_csv(out_path, index=False)
    print(f"\nSaved Species Analysis Table: {out_path}")
    return df_species


if __name__ == "__main__":
    run_species_analysis()
