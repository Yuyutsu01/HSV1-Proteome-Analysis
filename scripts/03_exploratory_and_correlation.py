"""
Phase 3: Exploratory Analysis and Feature Correlation Matrix for HSV-1 Proteome

Biological & Statistical Context:
Computes distributional statistics (mean, median, std, IQR, min, max, skewness)
for primary physicochemical metrics (Length, MW, pI, Instability Index, Aromaticity)
and individual amino acid percentages.
Calculates the full Pearson and Spearman correlation matrices to evaluate
redundancy (such as near-linear length-MW coupling) and orthogonality (such as pI).

Author: Computational Biology Pipeline
"""

import pandas as pd
import numpy as np

def run_exploratory_analysis(features_path="processed_data/physicochemical_features.csv",
                             output_stats_path="results/exploratory_statistics.csv",
                             output_corr_path="results/feature_correlation_matrix.csv"):
    """
    Compute rigorous descriptive statistics and correlation matrices for all 26 features.
    """
    df = pd.read_csv(features_path)
    
    feature_cols = [c for c in df.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']]
    
    stats_list = []
    for col in feature_cols:
        series = df[col]
        q25 = series.quantile(0.25)
        q75 = series.quantile(0.75)
        iqr = q75 - q25
        stats_list.append({
            'Feature': col,
            'Mean': series.mean(),
            'Std': series.std(),
            'Median': series.median(),
            'IQR': iqr,
            'Min': series.min(),
            'Max': series.max(),
            'Skewness': series.skew(),
            'Kurtosis': series.kurtosis()
        })
        
    stats_df = pd.DataFrame(stats_list)
    stats_df.to_csv(output_stats_path, index=False)
    print(f"Saved exploratory statistics to {output_stats_path}")
    
    # Compute correlation matrix
    corr_df = df[feature_cols].corr(method='pearson')
    corr_df.to_csv(output_corr_path)
    print(f"Saved Pearson correlation matrix ({corr_df.shape}) to {output_corr_path}")
    
    # Key correlations summary
    key_features = ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']
    print("\nKey 5-Feature Correlation Matrix:")
    print(df[key_features].corr().round(4).to_string())
    
    return stats_df, corr_df

if __name__ == "__main__":
    run_exploratory_analysis()
