"""
Phase 12: Supervised Machine Learning Classification Benchmarks (Strict In-Fold Scaling Pipeline)

Methodological Rigor & Data Leakage Prevention:
- To prevent any information leakage between training and testing folds, all feature standardization
  is executed strictly inside each fold using scikit-learn Pipeline([('scaler', StandardScaler()), ('classifier', clf)]).
  Means and standard deviations are computed purely from training folds and projected onto test folds.
- Target: Viral Temporal Expression Class (Immediate-Early vs Early vs Late).
- Benchmark Models:
  1. Logistic Regression (L2 regularized)
  2. Support Vector Machine (SVM with RBF kernel)
  3. Random Forest Classifier
  4. XGBoost Classifier
- Input Representations:
  A. Physicochemical (25-dim raw)
  B. ProtBERT (1024-dim raw)
  C. Combined Multimodal (1049-dim raw)
- Validation: Repeated Stratified 5-Fold Cross-Validation (5 splits x 5 repeats = 25 evaluations).
- Metrics reported: Mean & Std for Accuracy, Macro-Precision, Macro-Recall, Macro-F1.

Author: Computational Biology Pipeline
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix

def run_classification_benchmarks(physico_path="processed_data/physicochemical_features.csv",
                                  protbert_path="embeddings/protbert_embeddings.npy",
                                  output_csv="results/classification_results.csv",
                                  output_cm_csv="results/classification_confusion_matrices.csv"):
    """
    Run multi-model, multi-representation stratified cross-validation with strict in-fold pipeline scaling.
    """
    df_phys = pd.read_csv(physico_path)
    feature_cols = [c for c in df_phys.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']]
    
    # Raw unscaled matrices
    X_phys_raw = df_phys[feature_cols].values.astype(np.float32)
    X_bert_raw = np.load(protbert_path).astype(np.float32)
    X_comb_raw = np.hstack([X_phys_raw, X_bert_raw])
    
    # Target: Temporal Class
    y_raw = df_phys['temporal_class'].values
    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    class_names = le.classes_
    
    representations = {
        'Physicochemical': X_phys_raw,
        'ProtBERT': X_bert_raw,
        'Combined': X_comb_raw
    }
    
    models = {
        'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
        'SVM (RBF)': SVC(kernel='rbf', random_state=42),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'XGBoost': XGBClassifier(n_estimators=50, max_depth=3, eval_metric='mlogloss', random_state=42)
    }
    
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
    
    results = []
    cm_records = []
    
    for rep_name, X_mat in representations.items():
        for model_name, base_clf in models.items():
            acc_list, prec_list, rec_list, f1_list = [], [], [], []
            y_true_all, y_pred_all = [], []
            
            for train_idx, test_idx in rskf.split(X_mat, y):
                X_train, X_test = X_mat[train_idx], X_mat[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]
                
                # Enforce strict in-fold scaling via Pipeline
                pipe = Pipeline([
                    ('scaler', StandardScaler()),
                    ('clf', base_clf)
                ])
                
                pipe.fit(X_train, y_train)
                y_pred = pipe.predict(X_test)
                
                acc_list.append(accuracy_score(y_test, y_pred))
                prec_list.append(precision_score(y_test, y_pred, average='macro', zero_division=0))
                rec_list.append(recall_score(y_test, y_pred, average='macro', zero_division=0))
                f1_list.append(f1_score(y_test, y_pred, average='macro', zero_division=0))
                
                y_true_all.extend(y_test)
                y_pred_all.extend(y_pred)
                
            results.append({
                'Representation': rep_name,
                'Classifier': model_name,
                'Accuracy_Mean': np.mean(acc_list),
                'Accuracy_Std': np.std(acc_list),
                'Macro_Precision_Mean': np.mean(prec_list),
                'Macro_Precision_Std': np.std(prec_list),
                'Macro_Recall_Mean': np.mean(rec_list),
                'Macro_Recall_Std': np.std(rec_list),
                'Macro_F1_Mean': np.mean(f1_list),
                'Macro_F1_Std': np.std(f1_list)
            })
            
            cm = confusion_matrix(y_true_all, y_pred_all)
            cm_records.append({
                'Representation': rep_name,
                'Classifier': model_name,
                'Confusion_Matrix': cm.tolist()
            })
            
    res_df = pd.DataFrame(results)
    res_df.to_csv(output_csv, index=False)
    
    cm_df = pd.DataFrame(cm_records)
    cm_df.to_csv(output_cm_csv, index=False)
    
    print("\n[LEAK-FREE] Supervised Classification Performance Summary (Repeated 5-Fold CV):")
    print(res_df.sort_values(by='Macro_F1_Mean', ascending=False).to_string(index=False))
    
    return res_df, cm_df

if __name__ == "__main__":
    run_classification_benchmarks()
