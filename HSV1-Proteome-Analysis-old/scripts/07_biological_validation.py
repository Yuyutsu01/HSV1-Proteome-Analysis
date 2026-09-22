"""
Phase 7: Biological Validation of Unsupervised Clusters

Methodology & Ethical Separation:
- Unsupervised discovery was performed purely on feature distances without access to biological labels.
- Biological validation is performed POST-HOC by comparing unsupervised clusters against independently curated
  functional annotations (Temporal classes: Immediate-Early, Early, Late; Functional classes: Capsid/Structural, Tegument, Glycoprotein, Replication/Enzymes, Regulatory).
- Metrics computed:
  * Adjusted Rand Index (ARI): Corrected-for-chance measure of cluster overlap [-1, 1]
  * Normalized Mutual Information (NMI): Information-theoretic overlap [0, 1]
  * Cluster Purity: Fraction of dominant class samples in each cluster [0, 1]
  * Contingency Matrix: Full cross-tabulation table

Author: Computational Biology Pipeline
"""

import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.cluster import contingency_matrix

def compute_cluster_purity(y_true, y_pred):
    """
    Calculate cluster purity: sum_k(max_j |w_k cap c_j|) / N
    """
    matrix = contingency_matrix(y_true, y_pred)
    return np.sum(np.amax(matrix, axis=0)) / np.sum(matrix)

def run_biological_validation(clusters_path="results/clustering_results.csv",
                              output_val_path="results/biological_validation_metrics.csv",
                              output_contingency_path="results/contingency_matrix_kmeans.csv"):
    """
    Evaluate alignment between unsupervised clusters and ground-truth biological annotations.
    """
    df = pd.read_csv(clusters_path)
    
    # Filter for proteins with assigned labels
    valid_df = df[df['temporal_class'] != 'Unassigned'].copy()
    
    # Validation against Temporal Class
    y_temporal = valid_df['temporal_class'].values
    y_func = valid_df['functional_category'].values
    
    cluster_cols = ['KMeans_k4_cluster', 'Hierarchical_Ward_k4', 'Hierarchical_Average_k4']
    
    val_results = []
    
    for ccol in cluster_cols:
        pred_labels = valid_df[ccol].values
        
        # Temporal metrics
        ari_temp = adjusted_rand_score(y_temporal, pred_labels)
        nmi_temp = normalized_mutual_info_score(y_temporal, pred_labels)
        purity_temp = compute_cluster_purity(y_temporal, pred_labels)
        
        # Functional metrics
        ari_func = adjusted_rand_score(y_func, pred_labels)
        nmi_func = normalized_mutual_info_score(y_func, pred_labels)
        purity_func = compute_cluster_purity(y_func, pred_labels)
        
        val_results.append({
            'Method': ccol,
            'Temporal_ARI': ari_temp,
            'Temporal_NMI': nmi_temp,
            'Temporal_Purity': purity_temp,
            'Functional_ARI': ari_func,
            'Functional_NMI': nmi_func,
            'Functional_Purity': purity_func
        })
        
    val_df = pd.DataFrame(val_results)
    val_df.to_csv(output_val_path, index=False)
    
    print("\nBiological Validation Metrics (Unsupervised vs Curated Ground Truth):")
    print(val_df.to_string(index=False))
    
    # Generate contingency cross-tabulation table for K-Means k=4 vs Functional & Temporal
    cm_temp = pd.crosstab(valid_df['temporal_class'], valid_df['KMeans_k4_cluster'], rownames=['Temporal Class'], colnames=['KMeans Cluster'])
    cm_func = pd.crosstab(valid_df['functional_category'], valid_df['KMeans_k4_cluster'], rownames=['Functional Category'], colnames=['KMeans Cluster'])
    
    print("\nContingency Matrix (Temporal Class vs K-Means k=4):")
    print(cm_temp)
    
    print("\nContingency Matrix (Functional Category vs K-Means k=4):")
    print(cm_func)
    
    cm_func.to_csv(output_contingency_path)
    return val_df, cm_temp, cm_func

if __name__ == "__main__":
    run_biological_validation()
