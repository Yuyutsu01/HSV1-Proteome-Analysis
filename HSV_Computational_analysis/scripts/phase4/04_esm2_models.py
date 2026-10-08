"""
Phase 4 - Step 04: Pretrained Protein Language Model Models (ESM-2).

Biological & Computational Concept:
Evaluates whether high-dimensional contextual embeddings from ESM-2 (facebook/esm2_t6_8M_UR50D, D=320)
provide superior temporal-class discriminative power compared to classical frequency-based descriptors.

Classifiers evaluated:
- Logistic Regression (L2 regularized, class-weighted)
- Linear SVM (LinearSVC, class-weighted)
- Random Forest (100 trees, max_depth=15, class-weighted)
- HistGradientBoostingClassifier (Gradient-boosted decision trees)

Evaluated under both Random Stratified and Homology-Aware evaluation regimes.
"""

import os
import yaml
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from evaluation_utils import compute_comprehensive_metrics, CLASS_LABELS


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_esm2_models():
    """Run ESM-2 representation modeling experiments."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 04: ESM-2 PROTEIN LANGUAGE MODEL EVALUATION")
    print("=" * 80)

    # 1. Load splits and ESM-2 embeddings
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    esm_data = np.load(config["features"]["esm2_path"], allow_pickle=True)
    X_full = esm_data['embeddings']
    
    assert len(X_full) == len(df_splits), "Mismatch between embeddings and split manifest!"
    print(f"Loaded ESM-2 Embeddings: Shape={X_full.shape}")

    all_results = []

    for regime_col, regime_name in [('split_random', 'Random Stratified'), ('split_homology', 'Homology-Aware')]:
        print(f"\n==================== REGIME: {regime_name} ====================")
        train_mask = (df_splits[regime_col] == 'TRAIN').values
        test_mask = (df_splits[regime_col] == 'TEST').values
        
        y_train = df_splits.loc[train_mask, 'temporal_class'].values
        y_test = df_splits.loc[test_mask, 'temporal_class'].values

        X_train = X_full[train_mask]
        X_test = X_full[test_mask]

        scaler = StandardScaler()
        X_tr_scaled = scaler.fit_transform(X_train)
        X_te_scaled = scaler.transform(X_test)

        # 1. Logistic Regression
        print(f"  Training Logistic Regression (ESM-2 | {regime_name})...")
        lr = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced', solver='lbfgs')
        lr.fit(X_tr_scaled, y_train)
        y_pred_lr = lr.predict(X_te_scaled)
        y_prob_lr = lr.predict_proba(X_te_scaled)
        m_lr = compute_comprehensive_metrics(y_test, y_pred_lr, y_prob_lr, "ESM-2 (320-dim)", "Logistic Regression", regime_name)
        all_results.append(m_lr)
        print(f"    -> Acc: {m_lr['accuracy']}, BalAcc: {m_lr['balanced_accuracy']}, MacroF1: {m_lr['macro_f1']}, MCC: {m_lr['mcc']}")

        # 2. Linear SVM
        print(f"  Training Linear SVM (ESM-2 | {regime_name})...")
        svm = LinearSVC(C=1.0, max_iter=2000, random_state=config["random_seed"], class_weight='balanced', dual='auto', tol=1e-3)
        svm.fit(X_tr_scaled, y_train)
        y_pred_svm = svm.predict(X_te_scaled)
        dec_func = svm.decision_function(X_te_scaled)
        exp_dec = np.exp(dec_func - np.max(dec_func, axis=1, keepdims=True))
        y_prob_svm = exp_dec / np.sum(exp_dec, axis=1, keepdims=True)
        m_svm = compute_comprehensive_metrics(y_test, y_pred_svm, y_prob_svm, "ESM-2 (320-dim)", "Linear SVM", regime_name)
        all_results.append(m_svm)
        print(f"    -> Acc: {m_svm['accuracy']}, BalAcc: {m_svm['balanced_accuracy']}, MacroF1: {m_svm['macro_f1']}, MCC: {m_svm['mcc']}")

        # 3. Random Forest
        print(f"  Training Random Forest (ESM-2 | {regime_name})...")
        rf = RandomForestClassifier(n_estimators=100, max_depth=15, random_state=config["random_seed"], class_weight='balanced', n_jobs=-1)
        rf.fit(X_tr_scaled, y_train)
        y_pred_rf = rf.predict(X_te_scaled)
        y_prob_rf = rf.predict_proba(X_te_scaled)
        m_rf = compute_comprehensive_metrics(y_test, y_pred_rf, y_prob_rf, "ESM-2 (320-dim)", "Random Forest", regime_name)
        all_results.append(m_rf)
        print(f"    -> Acc: {m_rf['accuracy']}, BalAcc: {m_rf['balanced_accuracy']}, MacroF1: {m_rf['macro_f1']}, MCC: {m_rf['mcc']}")

        # 4. HistGradientBoostingClassifier
        print(f"  Training Gradient Boosting (ESM-2 | {regime_name})...")
        gb = HistGradientBoostingClassifier(max_iter=100, max_depth=6, random_state=config["random_seed"], class_weight='balanced')
        gb.fit(X_tr_scaled, y_train)
        y_pred_gb = gb.predict(X_te_scaled)
        y_prob_gb = gb.predict_proba(X_te_scaled)
        m_gb = compute_comprehensive_metrics(y_test, y_pred_gb, y_prob_gb, "ESM-2 (320-dim)", "Gradient Boosting", regime_name)
        all_results.append(m_gb)
        print(f"    -> Acc: {m_gb['accuracy']}, BalAcc: {m_gb['balanced_accuracy']}, MacroF1: {m_gb['macro_f1']}, MCC: {m_gb['mcc']}")

    df_res = pd.DataFrame(all_results)
    out_path = os.path.join(config["paths"]["phase4_tables_dir"], "esm2_models_results.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_res.to_csv(out_path, index=False)
    print(f"\nSaved ESM-2 Models Results: {out_path}")
    return df_res


if __name__ == "__main__":
    run_esm2_models()
