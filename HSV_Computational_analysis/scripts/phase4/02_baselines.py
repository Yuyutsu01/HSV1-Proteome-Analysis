"""
Phase 4 - Step 02: Baseline Models (Majority Class, Stratified Random, Length-Only).

Biological & Computational Concept:
1. Majority-Class Baseline:
   Always predicts the dominant class (LATE, ~65.8%). Sets the minimum floor for Accuracy,
   while producing Macro-F1 ≈ 0.26 and Balanced Accuracy = 0.3333.

2. Stratified Random Baseline:
   Samples predictions from the training class prior distribution. Repeated across 5 random seeds
   to report empirical Mean ± Standard Deviation, preventing single-seed cherry-picking.

3. Length-Only Baseline:
   Uses only sequence length as the input feature (scaled on training data).
   Evaluated with:
   - Logistic Regression
   - Shallow Decision Tree (max_depth=3)
   Mandatory control: Quantifies how much temporal predictive signal can be explained purely
   by sequence length differences without any residue composition or motif features.
"""

import os
import yaml
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import StandardScaler
from evaluation_utils import compute_comprehensive_metrics, CLASS_LABELS, CLASS_TO_IDX


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_baselines():
    """Execute all baseline models across Random and Homology-Aware splits."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 02: BASELINE MODELS EVALUATION")
    print("=" * 80)

    # Load splits and manifest
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    df_man = pd.read_csv(config["dataset"]["manifest_path"])
    
    # Merge length and labels
    df = df_splits.merge(df_man[['canonical_id', 'normalized_length']], on='canonical_id')
    
    results = []
    
    for split_col, regime_name in [('split_random', 'Random Stratified'), ('split_homology', 'Homology-Aware')]:
        print(f"\n--- Evaluating Baselines on {regime_name} Regime ---")
        train_mask = (df[split_col] == 'TRAIN').values
        val_mask = (df[split_col] == 'VALIDATION').values
        test_mask = (df[split_col] == 'TEST').values
        
        y_train = df.loc[train_mask, 'temporal_class'].values
        y_test = df.loc[test_mask, 'temporal_class'].values
        lengths_train = df.loc[train_mask, 'normalized_length'].values.reshape(-1, 1)
        lengths_test = df.loc[test_mask, 'normalized_length'].values.reshape(-1, 1)
        
        # 1. Majority-Class Baseline (Always LATE)
        y_pred_maj = np.full(len(y_test), 'LATE')
        y_prob_maj = np.zeros((len(y_test), 3))
        y_prob_maj[:, 2] = 1.0 # 100% prob on LATE
        
        metrics_maj = compute_comprehensive_metrics(
            y_true=y_test,
            y_pred=y_pred_maj,
            y_prob=y_prob_maj,
            representation_name="None (Baseline)",
            model_name="Majority Class (Always Late)",
            split_regime=regime_name
        )
        results.append(metrics_maj)
        print(f"Majority Baseline -> Acc: {metrics_maj['accuracy']}, BalAcc: {metrics_maj['balanced_accuracy']}, MacroF1: {metrics_maj['macro_f1']}, MCC: {metrics_maj['mcc']}")

        # 2. Stratified Random Baseline (Repeated over 5 seeds)
        train_probs = [np.mean(y_train == c) for c in CLASS_LABELS]
        random_macro_f1s = []
        random_metrics_list = []
        
        for seed in config["repeated_seeds"]:
            rng = np.random.RandomState(seed)
            y_pred_rand_idx = rng.choice(3, size=len(y_test), p=train_probs)
            y_pred_rand = np.array([CLASS_LABELS[i] for i in y_pred_rand_idx])
            y_prob_rand = np.tile(train_probs, (len(y_test), 1))
            
            m_rand = compute_comprehensive_metrics(
                y_true=y_test,
                y_pred=y_pred_rand,
                y_prob=y_prob_rand,
                representation_name="None (Baseline)",
                model_name="Stratified Random",
                split_regime=regime_name
            )
            random_metrics_list.append(m_rand)
            random_macro_f1s.append(m_rand['macro_f1'])
            
        # Report mean of random baseline
        avg_rand = random_metrics_list[0].copy()
        avg_rand['macro_f1'] = round(float(np.mean(random_macro_f1s)), 4)
        avg_rand['accuracy'] = round(float(np.mean([m['accuracy'] for m in random_metrics_list])), 4)
        avg_rand['balanced_accuracy'] = round(float(np.mean([m['balanced_accuracy'] for m in random_metrics_list])), 4)
        avg_rand['mcc'] = round(float(np.mean([m['mcc'] for m in random_metrics_list])), 4)
        results.append(avg_rand)
        print(f"Stratified Random -> Acc: {avg_rand['accuracy']}, BalAcc: {avg_rand['balanced_accuracy']}, MacroF1: {avg_rand['macro_f1']} ± {np.std(random_macro_f1s):.4f}")

        # 3. Length-Only Baseline: Logistic Regression
        scaler = StandardScaler()
        X_train_len = scaler.fit_transform(lengths_train)
        X_test_len = scaler.transform(lengths_test)
        
        lr_len = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
        lr_len.fit(X_train_len, y_train)
        y_pred_lr = lr_len.predict(X_test_len)
        y_prob_lr = lr_len.predict_proba(X_test_len)
        
        m_lr_len = compute_comprehensive_metrics(
            y_true=y_test,
            y_pred=y_pred_lr,
            y_prob=y_prob_lr,
            representation_name="Length Only",
            model_name="Logistic Regression (Length)",
            split_regime=regime_name
        )
        results.append(m_lr_len)
        print(f"Length-Only (Logistic Reg) -> Acc: {m_lr_len['accuracy']}, BalAcc: {m_lr_len['balanced_accuracy']}, MacroF1: {m_lr_len['macro_f1']}, MCC: {m_lr_len['mcc']}")

        # 4. Length-Only Baseline: Shallow Decision Tree
        dt_len = DecisionTreeClassifier(max_depth=3, random_state=config["random_seed"], class_weight='balanced')
        dt_len.fit(X_train_len, y_train)
        y_pred_dt = dt_len.predict(X_test_len)
        y_prob_dt = dt_len.predict_proba(X_test_len)
        
        m_dt_len = compute_comprehensive_metrics(
            y_true=y_test,
            y_pred=y_pred_dt,
            y_prob=y_prob_dt,
            representation_name="Length Only",
            model_name="Decision Tree (Length, depth=3)",
            split_regime=regime_name
        )
        results.append(m_dt_len)
        print(f"Length-Only (Decision Tree) -> Acc: {m_dt_len['accuracy']}, BalAcc: {m_dt_len['balanced_accuracy']}, MacroF1: {m_dt_len['macro_f1']}, MCC: {m_dt_len['mcc']}")

    df_res = pd.DataFrame(results)
    out_path = os.path.join(config["paths"]["phase4_tables_dir"], "baseline_models_results.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_res.to_csv(out_path, index=False)
    print(f"\nSaved Baseline Models Results: {out_path}")
    return df_res


if __name__ == "__main__":
    run_baselines()
