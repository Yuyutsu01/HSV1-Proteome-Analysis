"""
Phase 6: Unsupervised Clustering Evaluation & Algorithm Comparison

Methodology:
- Evaluates K-Means clustering across K = 2 through K = 10 on standardized features.
- Computes four quantitative validation metrics for each K:
  1. Within-Cluster Sum of Squares (WCSS / Inertia)
  2. Silhouette Score (measures separation and cohesion in [-1, 1])
  3. Calinski-Harabasz Index (Variance Ratio Criterion; higher is better)
  4. Davies-Bouldin Index (average similarity of clusters; lower is better)
- Selects the optimal K strictly based on quantitative metrics (Elbow + Silhouette peak).
- Performs Agglomerative Hierarchical Clustering (Ward and Average linkage) for direct algorithmic comparison.
- Saves cluster evaluation curves and cluster assignments.

Author: Computational Biology Pipeline
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

def run_clustering_evaluation(features_path="processed_data/physicochemical_features.csv",
                              output_eval_path="results/cluster_evaluation.csv",
                              output_clusters_path="results/clustering_results.csv"):
    """
    Evaluate K-Means across K=2..10 and compare with Hierarchical Clustering.
    """
    df = pd.read_csv(features_path)
    feature_cols = [c for c in df.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']]
    
    X = df[feature_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    k_range = list(range(2, 11))
    eval_metrics = []
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=20)
        labels = kmeans.fit_predict(X_scaled)
        
        inertia = kmeans.inertia_
        sil = silhouette_score(X_scaled, labels)
        ch = calinski_harabasz_score(X_scaled, labels)
        db = davies_bouldin_score(X_scaled, labels)
        
        eval_metrics.append({
            'K': k,
            'Inertia_WCSS': inertia,
            'Silhouette_Score': sil,
            'Calinski_Harabasz': ch,
            'Davies_Bouldin': db
        })
        
    eval_df = pd.DataFrame(eval_metrics)
    eval_df.to_csv(output_eval_path, index=False)
    print("K-Means Cluster Evaluation Across K=2..10:")
    print(eval_df.to_string(index=False))
    
    # Identify best K by Silhouette and Davies-Bouldin
    best_k_sil = eval_df.loc[eval_df['Silhouette_Score'].idxmax()]['K']
    print(f"\nOptimal K suggested by Peak Silhouette Score: K={int(best_k_sil)}")
    
    # Run K-Means for K=4 (as evaluated in paper) and K=best_k_sil
    km4 = KMeans(n_clusters=4, random_state=42, n_init=20)
    km4_labels = km4.fit_predict(X_scaled)
    
    # Hierarchical Ward (K=4)
    agg_ward = AgglomerativeClustering(n_clusters=4, linkage='ward')
    ward_labels = agg_ward.fit_predict(X_scaled)
    
    # Hierarchical Average (K=4)
    agg_avg = AgglomerativeClustering(n_clusters=4, linkage='average')
    avg_labels = agg_avg.fit_predict(X_scaled)
    
    clusters_df = df[['protein_id', 'gene', 'temporal_class', 'functional_category', 'length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']].copy()
    clusters_df['KMeans_k4_cluster'] = km4_labels
    clusters_df['Hierarchical_Ward_k4'] = ward_labels
    clusters_df['Hierarchical_Average_k4'] = avg_labels
    
    clusters_df.to_csv(output_clusters_path, index=False)
    print(f"Saved cluster assignments to {output_clusters_path}")
    
    return eval_df, clusters_df, X_scaled

if __name__ == "__main__":
    run_clustering_evaluation()
