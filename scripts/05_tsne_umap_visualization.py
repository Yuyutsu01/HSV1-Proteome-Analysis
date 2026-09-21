"""
Phase 5: t-SNE and UMAP Non-Linear Manifold Projections

Methodology:
- Applies t-SNE and UMAP strictly as non-linear visual exploratory methods.
- Fixed random_state=42 for exact reproducibility.
- Hyperparameters documented:
  * t-SNE: n_components=2, perplexity=15, early_exaggeration=12.0, learning_rate='auto', max_iter=1000, random_state=42, metric='euclidean'
  * UMAP: n_components=2, n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42
- Saves 2D coordinates to results/manifold_embeddings.csv

Author: Computational Biology Pipeline
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
import umap

def run_tsne_umap(features_path="processed_data/physicochemical_features.csv",
                  output_path="results/physicochemical_manifold_embeddings.csv"):
    """
    Generate t-SNE and UMAP 2D coordinates from standardized physicochemical features.
    """
    df = pd.read_csv(features_path)
    feature_cols = [c for c in df.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']]
    
    X = df[feature_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # t-SNE
    print("Computing t-SNE (perplexity=15, random_state=42)...")
    tsne = TSNE(n_components=2, perplexity=15, random_state=42, n_iter_without_progress=300)
    X_tsne = tsne.fit_transform(X_scaled)
    
    # UMAP
    print("Computing UMAP (n_neighbors=15, min_dist=0.1, random_state=42)...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42)
    X_umap = reducer.fit_transform(X_scaled)
    
    results_df = df[['protein_id', 'gene', 'temporal_class', 'functional_category']].copy()
    results_df['tSNE_1'] = X_tsne[:, 0]
    results_df['tSNE_2'] = X_tsne[:, 1]
    results_df['UMAP_1'] = X_umap[:, 0]
    results_df['UMAP_2'] = X_umap[:, 1]
    
    results_df.to_csv(output_path, index=False)
    print(f"Saved manifold projections to {output_path}")
    return results_df

if __name__ == "__main__":
    run_tsne_umap()
