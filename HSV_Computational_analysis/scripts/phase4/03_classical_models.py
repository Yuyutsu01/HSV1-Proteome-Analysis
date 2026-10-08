"""
Phase 4 - Step 03: Classical Sequence Representation Models (AAC, Classical, 2-mer, 3-mer).

Biological & Computational Concept:
Evaluates whether interpretable, alignment-free sequence representations can predict viral temporal class:
1. Amino Acid Composition (AAC, 20-dim): Global residue frequency bias.
2. Classical Physicochemical Descriptors (13-dim): Biophysical properties (MW, hydrophobicity, aromaticity, charge, polarity, etc.).
3. 2-mer Frequencies (400-dim): Dipeptide motifs capturing local sequential dependencies.
4. 3-mer Frequencies (8,000-dim): Higher-order tripeptide motif composition.

Classifiers evaluated:
- Logistic Regression (L2 regularized, class-weighted)
- Linear SVM (LinearSVC, dual='auto', class-weighted)
- Random Forest (50 trees, max_depth=12, class-weighted)
- HistGradientBoostingClassifier (Gradient-boosted decision trees, max_iter=50)

Evaluated under both Random Stratified and Homology-Aware regimes. Preprocessing (StandardScaler)
is fitted strictly on the training partition.
"""

import os
import yaml
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from evaluation_utils import compute_comprehensive_metrics, CLASS_LABELS


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def train_and_eval_classifiers(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    rep_name: str,
    regime_name: str,
    random_seed: int = 42
) -> List[dict]:
    """Train standard classifiers on scaled training data and evaluate on test data."""
    scaler = StandardScaler()
    X_tr_scaled = scaler.fit_transform(X_train)
    X_te_scaled = scaler.transform(X_test)
    
    results = []
    
    # 1. Logistic Regression
    print(f"  Training Logistic Regression ({rep_name} | {regime_name})...", flush=True)
    lr = LogisticRegression(
        max_iter=1000, 
        random_state=random_seed, 
        class_weight='balanced',
        solver='lbfgs'
    )
    lr.fit(X_tr_scaled, y_train)
    y_pred_lr = lr.predict(X_te_scaled)
    y_prob_lr = lr.predict_proba(X_te_scaled)
    m_lr = compute_comprehensive_metrics(y_test, y_pred_lr, y_prob_lr, rep_name, "Logistic Regression", regime_name)
    results.append(m_lr)
    print(f"    -> Acc: {m_lr['accuracy']:.4f} | BalAcc: {m_lr['balanced_accuracy']:.4f} | MacroF1: {m_lr['macro_f1']:.4f} | MCC: {m_lr['mcc']:.4f}", flush=True)

    # 2. Linear SVM
    print(f"  Training Linear SVM ({rep_name} | {regime_name})...", flush=True)
    svm = LinearSVC(
        C=1.0, 
        max_iter=2000, 
        random_state=random_seed, 
        class_weight='balanced',
        dual='auto',
        tol=1e-3
    )
    svm.fit(X_tr_scaled, y_train)
    y_pred_svm = svm.predict(X_te_scaled)
    dec_func = svm.decision_function(X_te_scaled)
    exp_dec = np.exp(dec_func - np.max(dec_func, axis=1, keepdims=True))
    y_prob_svm = exp_dec / np.sum(exp_dec, axis=1, keepdims=True)
    m_svm = compute_comprehensive_metrics(y_test, y_pred_svm, y_prob_svm, rep_name, "Linear SVM", regime_name)
    results.append(m_svm)
    print(f"    -> Acc: {m_svm['accuracy']:.4f} | BalAcc: {m_svm['balanced_accuracy']:.4f} | MacroF1: {m_svm['macro_f1']:.4f} | MCC: {m_svm['mcc']:.4f}", flush=True)

    # 3. Random Forest
    print(f"  Training Random Forest ({rep_name} | {regime_name})...", flush=True)
    rf = RandomForestClassifier(
        n_estimators=50, 
        max_depth=12, 
        random_state=random_seed, 
        class_weight='balanced',
        n_jobs=-1
    )
    rf.fit(X_tr_scaled, y_train)
    y_pred_rf = rf.predict(X_te_scaled)
    y_prob_rf = rf.predict_proba(X_te_scaled)
    m_rf = compute_comprehensive_metrics(y_test, y_pred_rf, y_prob_rf, rep_name, "Random Forest", regime_name)
    results.append(m_rf)
    print(f"    -> Acc: {m_rf['accuracy']:.4f} | BalAcc: {m_rf['balanced_accuracy']:.4f} | MacroF1: {m_rf['macro_f1']:.4f} | MCC: {m_rf['mcc']:.4f}", flush=True)

    # 4. HistGradientBoostingClassifier (Fast on high-dim inputs)
    print(f"  Training Gradient Boosting ({rep_name} | {regime_name})...", flush=True)
    gb = HistGradientBoostingClassifier(
        max_iter=50, 
        max_depth=5, 
        random_state=random_seed, 
        class_weight='balanced'
    )
    gb.fit(X_tr_scaled, y_train)
    y_pred_gb = gb.predict(X_te_scaled)
    y_prob_gb = gb.predict_proba(X_te_scaled)
    m_gb = compute_comprehensive_metrics(y_test, y_pred_gb, y_prob_gb, rep_name, "Gradient Boosting", regime_name)
    results.append(m_gb)
    print(f"    -> Acc: {m_gb['accuracy']:.4f} | BalAcc: {m_gb['balanced_accuracy']:.4f} | MacroF1: {m_gb['macro_f1']:.4f} | MCC: {m_gb['mcc']:.4f}", flush=True)

    return results


def run_classical_models():
    """Run all classical representation modeling experiments."""
    config = load_config()
    print("=" * 80, flush=True)
    print("PHASE 4 - STEP 03: CLASSICAL SEQUENCE REPRESENTATIONS MODELING", flush=True)
    print("=" * 80, flush=True)

    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])

    print("Loading feature representations...", flush=True)
    # A. AAC
    df_aac = pd.read_csv(config["features"]["aac_path"])
    aac_cols = [c for c in df_aac.columns if c not in ['canonical_id', 'temporal_class']]
    feats_aac = df_aac[aac_cols].values.astype(np.float32)

    # B. Classical Physicochemical
    df_class = pd.read_csv(config["features"]["classical_path"])
    class_cols = [c for c in df_class.columns if c not in ['canonical_id', 'temporal_class']]
    feats_class = df_class[class_cols].values.astype(np.float32)

    # C. 2-mer
    k2_data = np.load(config["features"]["kmer_k2_path"], allow_pickle=True)
    feats_k2 = k2_data['features'].astype(np.float32)

    # D. 3-mer
    k3_data = np.load(config["features"]["kmer_k3_path"], allow_pickle=True)
    feats_k3 = k3_data['features'].astype(np.float32)

    feature_dict = {
        "AAC": feats_aac,
        "Classical Descriptors": feats_class,
        "2-mer (Dipeptides)": feats_k2,
        "3-mer (Tripeptides)": feats_k3
    }

    all_results = []

    for regime_col, regime_name in [('split_random', 'Random Stratified'), ('split_homology', 'Homology-Aware')]:
        print(f"\n==================== REGIME: {regime_name} ====================", flush=True)
        train_mask = (df_splits[regime_col] == 'TRAIN').values
        test_mask = (df_splits[regime_col] == 'TEST').values
        
        y_train = df_splits.loc[train_mask, 'temporal_class'].values
        y_test = df_splits.loc[test_mask, 'temporal_class'].values

        for rep_name, X_full in feature_dict.items():
            print(f"\nEvaluating {rep_name} (Shape: {X_full.shape}) under {regime_name}...", flush=True)
            X_train = X_full[train_mask]
            X_test = X_full[test_mask]
            
            res = train_and_eval_classifiers(
                X_train=X_train,
                y_train=y_train,
                X_test=X_test,
                y_test=y_test,
                rep_name=rep_name,
                regime_name=regime_name,
                random_seed=config["random_seed"]
            )
            all_results.extend(res)

    df_res = pd.DataFrame(all_results)
    out_path = os.path.join(config["paths"]["phase4_tables_dir"], "classical_models_results.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_res.to_csv(out_path, index=False)
    print(f"\nSaved Classical Models Results: {out_path}", flush=True)
    return df_res


if __name__ == "__main__":
    run_classical_models()
