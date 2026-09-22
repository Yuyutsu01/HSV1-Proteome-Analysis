"""
Phase 4: Principal Component Analysis (PCA) on Physicochemical Features

Methodology:
- Standardizes the 25-feature physicochemical matrix (Z-score normalization with StandardScaler)
- Performs full PCA decomposition
- Computes explained variance ratio and cumulative variance across all components
- Computes PC1 and PC2 loadings and identifies top contributing features
- Generates coordinates for PC1, PC2, PC3 and saves results to results/pca_results.csv

Author: Computational Biology Pipeline
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def run_pca(features_path="processed_data/physicochemical_features.csv",
            output_pca_path="results/pca_results.csv",
            output_loadings_path="results/pca_loadings.csv",
            output_variance_path="results/pca_variance_explained.csv"):
    """
    Standardize feature matrix and compute PCA with detailed loadings and variance statistics.
    """
    df = pd.read_csv(features_path)
    feature_cols = [c for c in df.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']]
    
    X = df[feature_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=X.shape[1], random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    
    # Explained variance
    explained_var = pca.explained_variance_ratio_
    cum_var = np.cumsum(explained_var)
    
    var_df = pd.DataFrame({
        'Component': [f"PC{i+1}" for i in range(len(explained_var))],
        'Explained_Variance_Ratio': explained_var,
        'Cumulative_Explained_Variance': cum_var
    })
    var_df.to_csv(output_variance_path, index=False)
    
    print("\nPCA Variance Explained (Top 10 Components):")
    print(var_df.head(10).to_string(index=False))
    
    # Loadings (eigenvectors * sqrt(eigenvalues))
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    loadings_df = pd.DataFrame(loadings, index=feature_cols, columns=[f"PC{i+1}" for i in range(len(explained_var))])
    loadings_df.to_csv(output_loadings_path)
    
    print("\nTop 5 Absolute Loadings for PC1:")
    print(loadings_df['PC1'].abs().sort_values(ascending=False).head(5))
    
    print("\nTop 5 Absolute Loadings for PC2:")
    print(loadings_df['PC2'].abs().sort_values(ascending=False).head(5))
    
    # Save transformed coordinates per protein
    pca_results_df = df[['protein_id', 'gene', 'temporal_class', 'functional_category']].copy()
    for i in range(min(5, X_pca.shape[1])):
        pca_results_df[f'PC{i+1}'] = X_pca[:, i]
        
    pca_results_df.to_csv(output_pca_path, index=False)
    print(f"\nSaved PCA results to {output_pca_path}")
    
    return pca, X_scaled, X_pca, var_df, loadings_df

if __name__ == "__main__":
    run_pca()
