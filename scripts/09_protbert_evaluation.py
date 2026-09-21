"""
Phase 9: ProtBERT Representation Analysis, Manifold Reductions, and Clustering

Methodology:
- Loads the 1024-dimensional ProtBERT embeddings matrix.
- Performs PCA on ProtBERT embeddings: explained variance, scree analysis.
- Computes t-SNE and UMAP 2D projections of ProtBERT embeddings.
- Evaluates K-Means clustering across K = 2 through K = 10 on ProtBERT embeddings:
  * WCSS (Inertia)
  * Silhouette Score
  * Calinski-Harabasz Index
  * Davies-Bouldin Index
- Evaluates post-hoc biological validation metrics:
  * Temporal ARI, NMI, Purity
  * Functional ARI, NMI, Purity
- Compares K-Means with Agglomerative Hierarchical clustering on ProtBERT space.
- Saves results to results/protbert_clustering_results.csv and results/protbert_evaluation_metrics.csv.

Author: Computational Biology Pipeline
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.cluster import contingency_matrix

def compute_cluster_purity(y_true, y_pred):
    matrix = contingency_matrix(y_true, y_pred)
    return np.sum(np.amax(matrix, axis=0)) / np.sum(matrix)

def evaluate_protbert(embeddings_path="embeddings/protbert_embeddings.npy",
                      metadata_path="embeddings/protbert_metadata.csv",
                      curated_path="processed_data/curated_annotations.csv",
                      output_eval_path="results/protbert_clustering_evaluation.csv",
                      output_clusters_path="results/protbert_clustering_results.csv",
                      output_projections_path="results/protbert_manifold_embeddings.csv"):
    """
    Evaluate manifold projections and clustering on 1024-dim ProtBERT embeddings.
    """
    X_embed = np.load(embeddings_path)
    meta_df = pd.read_csv(metadata_path)
    curated_df = pd.read_csv(curated_path)
    
    # Merge metadata with curated labels
    merged_df = meta_df.merge(curated_df[['protein_id', 'temporal_class', 'functional_category']], on='protein_id', how='left')
    
    # Standardize or L2 normalize embeddings for distance-based ML
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_embed)
    
    # 1. PCA on ProtBERT
    pca = PCA(n_components=min(50, X_embed.shape[0]), random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    print(f"ProtBERT PCA PC1 explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%")
    print(f"ProtBERT PCA PC2 explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%")
    print(f"ProtBERT PCA Top 10 PCs cumulative variance: {np.sum(pca.explained_variance_ratio_[:10])*100:.2f}%")
    
    # 2. t-SNE & UMAP
    print("Computing t-SNE on ProtBERT embeddings...")
    tsne = TSNE(n_components=2, perplexity=15, random_state=42, n_iter_without_progress=300)
    X_tsne = tsne.fit_transform(X_scaled)
    
    print("Computing UMAP on ProtBERT embeddings...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='cosine', random_state=42)
    X_umap = reducer.fit_transform(X_scaled)
    
    proj_df = merged_df[['protein_id', 'gene', 'temporal_class', 'functional_category']].copy()
    proj_df['ProtBERT_PC1'] = X_pca[:, 0]
    proj_df['ProtBERT_PC2'] = X_pca[:, 1]
    proj_df['ProtBERT_tSNE1'] = X_tsne[:, 0]
    proj_df['ProtBERT_tSNE2'] = X_tsne[:, 1]
    proj_df['ProtBERT_UMAP1'] = X_umap[:, 0]
    proj_df['ProtBERT_UMAP2'] = X_umap[:, 1]
    proj_df.to_csv(output_projections_path, index=False)
    
    # 3. K-Means Evaluation across K=2..10
    k_range = list(range(2, 11))
    eval_metrics = []
    
    y_temp = merged_df['temporal_class'].values
    y_func = merged_df['functional_category'].values
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=20)
        labels = kmeans.fit_predict(X_scaled)
        
        sil = silhouette_score(X_scaled, labels)
        ch = calinski_harabasz_score(X_scaled, labels)
        db = davies_bouldin_score(X_scaled, labels)
        
        ari_temp = adjusted_rand_score(y_temp, labels)
        nmi_temp = normalized_mutual_info_score(y_temp, labels)
        purity_temp = compute_cluster_purity(y_temp, labels)
        
        ari_func = adjusted_rand_score(y_func, labels)
        nmi_func = normalized_mutual_info_score(y_func, labels)
        purity_func = compute_cluster_purity(y_func, labels)
        
        eval_metrics.append({
            'K': k,
            'Inertia_WCSS': kmeans.inertia_,
            'Silhouette_Score': sil,
            'Calinski_Harabasz': ch,
            'Davies_Bouldin': db,
            'Temporal_ARI': ari_temp,
            'Temporal_NMI': nmi_temp,
            'Temporal_Purity': purity_temp,
            'Functional_ARI': ari_func,
            'Functional_NMI': nmi_func,
            'Functional_Purity': purity_func
        })
        
    eval_df = pd.DataFrame(eval_metrics)
    eval_df.to_csv(output_eval_path, index=False)
    print("\nProtBERT K-Means Evaluation (K=2..10):")
    print(eval_df[['K', 'Silhouette_Score', 'Calinski_Harabasz', 'Davies_Bouldin', 'Temporal_ARI', 'Functional_ARI']].to_string(index=False))
    
    # Best K-Means and Hierarchical at K=4 for direct comparison
    km4 = KMeans(n_clusters=4, random_state=42, n_init=20)
    km4_labels = km4.fit_predict(X_scaled)
    
    agg_ward = AgglomerativeClustering(n_clusters=4, linkage='ward')
    ward_labels = agg_ward.fit_predict(X_scaled)
    
    clusters_df = merged_df[['protein_id', 'gene', 'temporal_class', 'functional_category']].copy()
    clusters_df['ProtBERT_KMeans_k4'] = km4_labels
    clusters_df['ProtBERT_Hierarchical_k4'] = ward_labels
    clusters_df.to_csv(output_clusters_path, index=False)
    
    print(f"\nSaved ProtBERT evaluation metrics to {output_eval_path}")
    print(f"Saved ProtBERT cluster assignments to {output_clusters_path}")
    
    return eval_df, clusters_df, proj_df

if __name__ == "__main__":
    evaluate_protbert()
