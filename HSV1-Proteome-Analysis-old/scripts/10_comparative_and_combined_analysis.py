"""
Phase 10 & 11: Physicochemical vs ProtBERT Comparative Analysis and Multimodal Combination

Methodology:
- Loads standardized 25-feature Physicochemical matrix and 1024-dim ProtBERT embeddings.
- Creates Combined Multimodal Representation:
  StandardScaler(Physicochemical) concatenated with StandardScaler(ProtBERT).
- Performs unified benchmarking across all three representations:
  1. Physicochemical Descriptors Only (25-dim)
  2. ProtBERT Embeddings Only (1024-dim)
  3. Combined Representation (1049-dim)
- Evaluates:
  * Optimal K selection across K=2..10
  * Best K Silhouette Score, Calinski-Harabasz, Davies-Bouldin
  * Post-hoc External Biological Validation (Temporal & Functional ARI, NMI, Purity)
- Saves representation_comparison.csv and combined_features.npy

Author: Computational Biology Pipeline
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.cluster import contingency_matrix

def compute_purity(y_true, y_pred):
    matrix = contingency_matrix(y_true, y_pred)
    return np.sum(np.amax(matrix, axis=0)) / np.sum(matrix)

def run_comparative_analysis(physico_path="processed_data/physicochemical_features.csv",
                             protbert_path="embeddings/protbert_embeddings.npy",
                             output_combined_npy="embeddings/combined_features.npy",
                             output_comparison_csv="results/representation_comparison.csv"):
    """
    Compare representations head-to-head and evaluate combined multimodal fusion.
    """
    df_phys = pd.read_csv(physico_path)
    feature_cols = [c for c in df_phys.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']]
    X_phys = df_phys[feature_cols].values
    
    scaler_phys = StandardScaler()
    X_phys_scaled = scaler_phys.fit_transform(X_phys)
    
    X_bert = np.load(protbert_path)
    scaler_bert = StandardScaler()
    X_bert_scaled = scaler_bert.fit_transform(X_bert)
    
    # Combined representation: [X_phys_scaled, X_bert_scaled]
    X_combined = np.hstack([X_phys_scaled, X_bert_scaled])
    np.save(output_combined_npy, X_combined)
    print(f"Saved combined representation matrix: {X_combined.shape} to {output_combined_npy}")
    
    y_temp = df_phys['temporal_class'].values
    y_func = df_phys['functional_category'].values
    
    representations = {
        'Physicochemical (25-dim)': X_phys_scaled,
        'ProtBERT (1024-dim)': X_bert_scaled,
        'Combined (1049-dim)': X_combined
    }
    
    comparison_rows = []
    
    for rep_name, X_mat in representations.items():
        # Evaluate for K=2..10 to find best K by Silhouette
        best_sil = -1
        best_k = 2
        for k in range(2, 11):
            km = KMeans(n_clusters=k, random_state=42, n_init=20)
            lbls = km.fit_predict(X_mat)
            sil = silhouette_score(X_mat, lbls)
            if sil > best_sil:
                best_sil = sil
                best_k = k
                
        # Evaluate at K=4 (paper baseline) and Best K
        for eval_k, k_tag in [(4, "K=4 (Baseline)"), (best_k, f"Best K (K={best_k})")]:
            km = KMeans(n_clusters=eval_k, random_state=42, n_init=20)
            labels = km.fit_predict(X_mat)
            
            sil = silhouette_score(X_mat, labels)
            ch = calinski_harabasz_score(X_mat, labels)
            db = davies_bouldin_score(X_mat, labels)
            
            ari_temp = adjusted_rand_score(y_temp, labels)
            nmi_temp = normalized_mutual_info_score(y_temp, labels)
            pur_temp = compute_purity(y_temp, labels)
            
            ari_func = adjusted_rand_score(y_func, labels)
            nmi_func = normalized_mutual_info_score(y_func, labels)
            pur_func = compute_purity(y_func, labels)
            
            comparison_rows.append({
                'Representation': rep_name,
                'Evaluation_Setting': k_tag,
                'K': eval_k,
                'Silhouette_Score': sil,
                'Calinski_Harabasz': ch,
                'Davies_Bouldin': db,
                'Temporal_ARI': ari_temp,
                'Temporal_NMI': nmi_temp,
                'Temporal_Purity': pur_temp,
                'Functional_ARI': ari_func,
                'Functional_NMI': nmi_func,
                'Functional_Purity': pur_func
            })
            
    comp_df = pd.DataFrame(comparison_rows)
    comp_df.to_csv(output_comparison_csv, index=False)
    
    print("\nHead-to-Head Representation Comparison:")
    print(comp_df[['Representation', 'Evaluation_Setting', 'Silhouette_Score', 'Temporal_ARI', 'Temporal_Purity', 'Functional_ARI', 'Functional_Purity']].to_string(index=False))
    
    return comp_df

if __name__ == "__main__":
    run_comparative_analysis()
