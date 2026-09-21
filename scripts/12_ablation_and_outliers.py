"""
Phase 13 & 14: Systematic Feature Ablation Study and Candidate Outlier Profiling (Strict Pipeline Scaling)

Methodology:
1. Feature Ablation Experiments:
   - Exp 1: Global Properties Only (5 features: Length, MW, pI, Instability Index, Aromaticity)
   - Exp 2: Amino Acid Composition Only (20 features: aa_A ... aa_Y)
   - Exp 3: All Physicochemical Descriptors (25 features)
   - Exp 4: ProtBERT Embeddings Only (1024 features)
   - Exp 5: Combined Multimodal Representation (1049 features)
   - Strict in-fold Pipeline scaling across Repeated Stratified 5-Fold Cross-Validation.

2. Candidate Outlier & Disagreement Profiling:
   - Evaluates per-protein predictions across models and representations.
   - Identifies proteins exhibiting prediction conflict or low confidence.

Author: Computational Biology Pipeline
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, silhouette_score, adjusted_rand_score

def run_ablation_and_outliers(physico_path="processed_data/physicochemical_features.csv",
                              protbert_path="embeddings/protbert_embeddings.npy",
                              output_ablation_csv="results/ablation_results.csv",
                              output_outliers_csv="results/candidate_outliers.csv"):
    """
    Run systematic feature ablation and candidate outlier detection using strict Pipeline scaling.
    """
    df = pd.read_csv(physico_path)
    
    global_cols = ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']
    comp_cols = [c for c in df.columns if c.startswith('aa_')]
    all_phys_cols = global_cols + comp_cols
    
    X_global_raw = df[global_cols].values.astype(np.float32)
    X_comp_raw = df[comp_cols].values.astype(np.float32)
    X_phys_raw = df[all_phys_cols].values.astype(np.float32)
    X_bert_raw = np.load(protbert_path).astype(np.float32)
    X_comb_raw = np.hstack([X_phys_raw, X_bert_raw])
    
    y_raw = df['temporal_class'].values
    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    
    ablation_subsets = {
        'Global Descriptors Only (5-dim)': X_global_raw,
        'AA Composition Only (20-dim)': X_comp_raw,
        'All Physicochemical (25-dim)': X_phys_raw,
        'ProtBERT Embeddings (1024-dim)': X_bert_raw,
        'Combined Multimodal (1049-dim)': X_comb_raw
    }
    
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
    ablation_results = []
    
    # Baseline for delta: All Physicochemical
    baseline_f1 = None
    baseline_acc = None
    
    for name, X_mat in ablation_subsets.items():
        acc_list, f1_list = [], []
        for train_idx, test_idx in rskf.split(X_mat, y):
            pipe = Pipeline([
                ('scaler', StandardScaler()),
                ('clf', RandomForestClassifier(n_estimators=100, random_state=42))
            ])
            pipe.fit(X_mat[train_idx], y[train_idx])
            preds = pipe.predict(X_mat[test_idx])
            acc_list.append(accuracy_score(y[test_idx], preds))
            f1_list.append(f1_score(y[test_idx], preds, average='macro', zero_division=0))
            
        # Unsupervised clustering metrics on full standardized matrix for reference
        scaler_full = StandardScaler()
        X_mat_scaled = scaler_full.fit_transform(X_mat)
        km = KMeans(n_clusters=4, random_state=42, n_init=20)
        km_labels = km.fit_predict(X_mat_scaled)
        sil = silhouette_score(X_mat_scaled, km_labels)
        ari = adjusted_rand_score(y, km_labels)
        
        mean_acc = np.mean(acc_list)
        mean_f1 = np.mean(f1_list)
        
        if name == 'All Physicochemical (25-dim)':
            baseline_f1 = mean_f1
            baseline_acc = mean_acc
            
        ablation_results.append({
            'Representation_Subset': name,
            'Feature_Count': X_mat.shape[1],
            'Accuracy_Mean': mean_acc,
            'Accuracy_Std': np.std(acc_list),
            'Macro_F1_Mean': mean_f1,
            'Macro_F1_Std': np.std(f1_list),
            'Silhouette_Score_k4': sil,
            'Temporal_ARI_k4': ari
        })
        
    ablation_df = pd.DataFrame(ablation_results)
    ablation_df['Delta_Macro_F1_vs_Phys'] = ablation_df['Macro_F1_Mean'] - baseline_f1
    ablation_df['Delta_Accuracy_vs_Phys'] = ablation_df['Accuracy_Mean'] - baseline_acc
    ablation_df.to_csv(output_ablation_csv, index=False)
    
    print("\n[LEAK-FREE] Ablation Study Results:")
    print(ablation_df[['Representation_Subset', 'Feature_Count', 'Accuracy_Mean', 'Macro_F1_Mean', 'Delta_Accuracy_vs_Phys', 'Delta_Macro_F1_vs_Phys', 'Temporal_ARI_k4']].to_string(index=False))
    
    # 2. Outlier & Disagreement Profiling
    pipe_phys = Pipeline([('scaler', StandardScaler()), ('clf', RandomForestClassifier(n_estimators=100, random_state=42))])
    pipe_phys.fit(X_phys_raw, y)
    phys_probs = pipe_phys.predict_proba(X_phys_raw)
    phys_preds = le.inverse_transform(np.argmax(phys_probs, axis=1))
    phys_conf = np.max(phys_probs, axis=1)
    
    pipe_bert = Pipeline([('scaler', StandardScaler()), ('clf', RandomForestClassifier(n_estimators=100, random_state=42))])
    pipe_bert.fit(X_bert_raw, y)
    bert_probs = pipe_bert.predict_proba(X_bert_raw)
    bert_preds = le.inverse_transform(np.argmax(bert_probs, axis=1))
    bert_conf = np.max(bert_probs, axis=1)
    
    km_phys = KMeans(n_clusters=4, random_state=42, n_init=20).fit_predict(StandardScaler().fit_transform(X_phys_raw))
    km_bert = KMeans(n_clusters=4, random_state=42, n_init=20).fit_predict(StandardScaler().fit_transform(X_bert_raw))
    
    outliers = []
    for i, row in df.iterrows():
        prot_id = row['protein_id']
        gene = row['gene']
        true_temp = row['temporal_class']
        func_cat = row['functional_category']
        
        is_disagreement = (phys_preds[i] != bert_preds[i])
        is_phys_error = (phys_preds[i] != true_temp)
        is_bert_error = (bert_preds[i] != true_temp)
        is_low_conf = (phys_conf[i] < 0.65 or bert_conf[i] < 0.65)
        
        if is_disagreement or is_phys_error or is_bert_error or is_low_conf:
            outliers.append({
                'protein_id': prot_id,
                'gene': gene,
                'product': row['product'],
                'true_temporal_class': true_temp,
                'functional_category': func_cat,
                'phys_pred_class': phys_preds[i],
                'phys_confidence': phys_conf[i],
                'bert_pred_class': bert_preds[i],
                'bert_confidence': bert_conf[i],
                'phys_kmeans_cluster': km_phys[i],
                'bert_kmeans_cluster': km_bert[i],
                'length': row['length'],
                'mw_kda': row['mw_kda'],
                'isoelectric_point': row['isoelectric_point'],
                'instability_index': row['instability_index'],
                'aromaticity': row['aromaticity'],
                'disagreement_flag': 'Phys_vs_BERT_Disagreement' if is_disagreement else ('Phys_Error' if is_phys_error else 'Low_Confidence')
            })
            
    outliers_df = pd.DataFrame(outliers)
    outliers_df.to_csv(output_outliers_csv, index=False)
    
    print(f"\nIdentified {len(outliers_df)} candidate outlier/disagreement proteins. Saved to {output_outliers_csv}")
    
    return ablation_df, outliers_df

if __name__ == "__main__":
    run_ablation_and_outliers()
