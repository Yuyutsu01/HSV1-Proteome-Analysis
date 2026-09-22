"""
Complete Publication-Quality Figure Generation Pipeline for HSV-1 Proteome Analysis

Generates all 15 publication figures (300 DPI, modern visual palettes, consistent typography).

Author: Computational Biology Pipeline
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import ast

# Global styling
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

os.makedirs('figures', exist_ok=True)

PALETTE_TEMPORAL = {'Immediate-Early': '#e74c3c', 'Early': '#e67e22', 'Late': '#2980b9', 'Latency': '#8e44ad', 'Unassigned': '#7f8c8d'}

def generate_all_figures():
    print("Executing complete 15-figure generation suite...")
    
    # Load base data
    df_phys = pd.read_csv("processed_data/physicochemical_features.csv")
    df_curated = pd.read_csv("processed_data/curated_annotations.csv")
    
    # -------------------------------------------------------------
    # Figure 1: Dataset Composition & Functional Breakdown
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    temp_counts = df_curated['temporal_class'].value_counts()
    colors_temp = [PALETTE_TEMPORAL.get(c, '#95a5a6') for c in temp_counts.index]
    axes[0].pie(temp_counts.values, labels=[f"{k}\n(n={v})" for k, v in temp_counts.items()],
                autopct='%1.1f%%', colors=colors_temp, startangle=140,
                wedgeprops={'edgecolor': 'white', 'linewidth': 1.5})
    axes[0].set_title("A. Temporal Expression Class Distribution (N=74)", weight='bold')
    
    func_counts = df_curated['functional_category'].value_counts().head(8)
    sns.barplot(x=func_counts.values, y=func_counts.index, ax=axes[1], palette='Blues_r')
    axes[1].set_xlabel("Number of Proteins")
    axes[1].set_title("B. Top Functional Categories in HSV-1 Proteome", weight='bold')
    for i, v in enumerate(func_counts.values):
        axes[1].text(v + 0.2, i, str(v), va='center', fontweight='bold', fontsize=9)
    plt.tight_layout()
    plt.savefig("figures/fig01_dataset_composition.png")
    plt.close()
    print("Saved fig01_dataset_composition.png")
    
    # -------------------------------------------------------------
    # Figure 2: Sequence Length and Molecular Weight Distribution
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.histplot(df_phys['length'], bins=20, kde=True, ax=axes[0], color='#2980b9', edgecolor='black')
    axes[0].axvline(df_phys['length'].mean(), color='#e74c3c', linestyle='--', linewidth=2, label=f"Mean: {df_phys['length'].mean():.1f} aa")
    axes[0].axvline(df_phys['length'].median(), color='#27ae60', linestyle=':', linewidth=2, label=f"Median: {df_phys['length'].median():.1f} aa")
    axes[0].set_xlabel("Sequence Length (amino acids)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("A. Sequence Length Distribution", weight='bold')
    axes[0].legend()
    
    sns.histplot(df_phys['mw_kda'], bins=20, kde=True, ax=axes[1], color='#8e44ad', edgecolor='black')
    axes[1].axvline(df_phys['mw_kda'].mean(), color='#e74c3c', linestyle='--', linewidth=2, label=f"Mean: {df_phys['mw_kda'].mean():.1f} kDa")
    axes[1].axvline(df_phys['mw_kda'].median(), color='#27ae60', linestyle=':', linewidth=2, label=f"Median: {df_phys['mw_kda'].median():.1f} kDa")
    axes[1].set_xlabel("Molecular Weight (kDa)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("B. Molecular Weight Distribution", weight='bold')
    axes[1].legend()
    plt.tight_layout()
    plt.savefig("figures/fig02_length_mw_distribution.png")
    plt.close()
    print("Saved fig02_length_mw_distribution.png")
    
    # -------------------------------------------------------------
    # Figure 3: Isoelectric Point, Instability Index & Aromaticity
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sns.histplot(df_phys['isoelectric_point'], bins=18, kde=True, ax=axes[0], color='#16a085', edgecolor='black')
    axes[0].axvline(7.0, color='black', linestyle='--', linewidth=1.5, label="Neutral pH 7.0")
    axes[0].set_xlabel("Theoretical Isoelectric Point (pI)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("A. Bimodal Isoelectric Point (pI)", weight='bold')
    axes[0].legend()
    
    sns.histplot(df_phys['instability_index'], bins=18, kde=True, ax=axes[1], color='#d35400', edgecolor='black')
    axes[1].axvline(40.0, color='#c0392b', linestyle='--', linewidth=2, label="Threshold (40.0)")
    axes[1].set_xlabel("Guruprasad Instability Index")
    axes[1].set_ylabel("Count")
    axes[1].set_title("B. In Vitro Instability Index", weight='bold')
    axes[1].legend()
    
    sns.histplot(df_phys['aromaticity'], bins=18, kde=True, ax=axes[2], color='#2c3e50', edgecolor='black')
    axes[2].axvline(df_phys['aromaticity'].mean(), color='#e74c3c', linestyle='--', linewidth=2, label=f"Mean: {df_phys['aromaticity'].mean():.3f}")
    axes[2].set_xlabel("Aromaticity (Phe + Trp + Tyr fraction)")
    axes[2].set_ylabel("Count")
    axes[2].set_title("C. Aromaticity Distribution", weight='bold')
    axes[2].legend()
    plt.tight_layout()
    plt.savefig("figures/fig03_pi_instability_aromaticity.png")
    plt.close()
    print("Saved fig03_pi_instability_aromaticity.png")
    
    # -------------------------------------------------------------
    # Figure 4: Amino Acid Composition Profile
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 5))
    aa_cols = [c for c in df_phys.columns if c.startswith('aa_')]
    aa_names = [c.replace('aa_', '') for c in aa_cols]
    mean_aa = df_phys[aa_cols].mean() * 100
    std_aa = df_phys[aa_cols].std() * 100
    bars = ax.bar(aa_names, mean_aa, yerr=std_aa, capsize=4, color='#34495e', edgecolor='black', alpha=0.85)
    ax.set_xlabel("Amino Acid Residue")
    ax.set_ylabel("Mean Molar Percentage (%)")
    ax.set_title("Molar Amino Acid Composition of HSV-1 Non-Redundant Proteome (N=74)", weight='bold')
    max_idx = np.argmax(mean_aa.values)
    bars[max_idx].set_color('#e74c3c')
    plt.tight_layout()
    plt.savefig("figures/fig04_aa_composition_profile.png")
    plt.close()
    print("Saved fig04_aa_composition_profile.png")
    
    # -------------------------------------------------------------
    # Figure 5: Feature Correlation Heatmap
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6.5))
    key_features = ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']
    key_labels = ['Length (aa)', 'MW (kDa)', 'pI', 'Instability Index', 'Aromaticity']
    corr_mat = df_phys[key_features].corr()
    sns.heatmap(corr_mat, annot=True, fmt='.3f', cmap='coolwarm', vmin=-1, vmax=1,
                xticklabels=key_labels, yticklabels=key_labels, ax=ax, cbar_kws={'label': 'Pearson r'})
    ax.set_title("Pearson Correlation Heatmap of Primary Physicochemical Descriptors", weight='bold')
    plt.tight_layout()
    plt.savefig("figures/fig05_correlation_heatmap.png")
    plt.close()
    print("Saved fig05_correlation_heatmap.png")
    
    # -------------------------------------------------------------
    # Figure 6: PCA Scree Plot and Biplot
    # -------------------------------------------------------------
    if os.path.exists("results/pca_variance_explained.csv") and os.path.exists("results/pca_results.csv"):
        var_df = pd.read_csv("results/pca_variance_explained.csv")
        pca_df = pd.read_csv("results/pca_results.csv")
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        axes[0].bar(range(1, 11), var_df['Explained_Variance_Ratio'][:10] * 100, color='#3498db', edgecolor='black', alpha=0.8, label='Individual')
        axes[0].step(range(1, 11), var_df['Cumulative_Explained_Variance'][:10] * 100, where='mid', color='#e74c3c', linewidth=2, label='Cumulative')
        axes[0].set_xlabel("Principal Component")
        axes[0].set_ylabel("Variance Explained (%)")
        axes[0].set_title("A. PCA Scree Plot (Top 10 Components)", weight='bold')
        axes[0].set_xticks(range(1, 11))
        axes[0].legend()
        
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='temporal_class', palette=PALETTE_TEMPORAL, s=80, edgecolor='black', ax=axes[1])
        axes[1].set_xlabel(f"PC1 ({var_df['Explained_Variance_Ratio'][0]*100:.1f}% variance)")
        axes[1].set_ylabel(f"PC2 ({var_df['Explained_Variance_Ratio'][1]*100:.1f}% variance)")
        axes[1].set_title("B. PCA Projection (PC1 vs PC2) by Temporal Class", weight='bold')
        axes[1].legend(title="Temporal Class")
        plt.tight_layout()
        plt.savefig("figures/fig06_pca_scree_and_biplot.png")
        plt.close()
        print("Saved fig06_pca_scree_and_biplot.png")
        
    # -------------------------------------------------------------
    # Figure 7: Physicochemical t-SNE and UMAP Projections
    # -------------------------------------------------------------
    if os.path.exists("results/physicochemical_manifold_embeddings.csv"):
        mani_df = pd.read_csv("results/physicochemical_manifold_embeddings.csv")
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        sns.scatterplot(data=mani_df, x='tSNE_1', y='tSNE_2', hue='temporal_class', palette=PALETTE_TEMPORAL, s=80, edgecolor='black', ax=axes[0])
        axes[0].set_title("A. t-SNE of Physicochemical Descriptors", weight='bold')
        axes[0].set_xlabel("t-SNE Dimension 1")
        axes[0].set_ylabel("t-SNE Dimension 2")
        axes[0].legend(title="Temporal Class")
        
        sns.scatterplot(data=mani_df, x='UMAP_1', y='UMAP_2', hue='temporal_class', palette=PALETTE_TEMPORAL, s=80, edgecolor='black', ax=axes[1])
        axes[1].set_title("B. UMAP of Physicochemical Descriptors", weight='bold')
        axes[1].set_xlabel("UMAP Dimension 1")
        axes[1].set_ylabel("UMAP Dimension 2")
        axes[1].legend(title="Temporal Class")
        plt.tight_layout()
        plt.savefig("figures/fig07_tsne_umap_physicochemical.png")
        plt.close()
        print("Saved fig07_tsne_umap_physicochemical.png")
        
    # -------------------------------------------------------------
    # Figure 8: K-Means Evaluation Curves Across K=2..10
    # -------------------------------------------------------------
    if os.path.exists("results/cluster_evaluation.csv"):
        eval_df = pd.read_csv("results/cluster_evaluation.csv")
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        axes[0, 0].plot(eval_df['K'], eval_df['Inertia_WCSS'], marker='o', color='#2980b9', linewidth=2)
        axes[0, 0].set_title("A. Within-Cluster Sum of Squares (Elbow)", weight='bold')
        axes[0, 0].set_xlabel("Number of Clusters (K)")
        axes[0, 0].set_ylabel("Inertia / WCSS")
        axes[0, 0].set_xticks(eval_df['K'])
        
        axes[0, 1].plot(eval_df['K'], eval_df['Silhouette_Score'], marker='s', color='#27ae60', linewidth=2)
        axes[0, 1].set_title("B. Silhouette Score vs K", weight='bold')
        axes[0, 1].set_xlabel("Number of Clusters (K)")
        axes[0, 1].set_ylabel("Silhouette Score")
        axes[0, 1].set_xticks(eval_df['K'])
        
        axes[1, 0].plot(eval_df['K'], eval_df['Calinski_Harabasz'], marker='^', color='#8e44ad', linewidth=2)
        axes[1, 0].set_title("C. Calinski-Harabasz Index (Higher is better)", weight='bold')
        axes[1, 0].set_xlabel("Number of Clusters (K)")
        axes[1, 0].set_ylabel("Calinski-Harabasz Score")
        axes[1, 0].set_xticks(eval_df['K'])
        
        axes[1, 1].plot(eval_df['K'], eval_df['Davies_Bouldin'], marker='d', color='#d35400', linewidth=2)
        axes[1, 1].set_title("D. Davies-Bouldin Index (Lower is better)", weight='bold')
        axes[1, 1].set_xlabel("Number of Clusters (K)")
        axes[1, 1].set_ylabel("Davies-Bouldin Score")
        axes[1, 1].set_xticks(eval_df['K'])
        plt.tight_layout()
        plt.savefig("figures/fig08_kmeans_evaluation_curves.png")
        plt.close()
        print("Saved fig08_kmeans_evaluation_curves.png")

    # -------------------------------------------------------------
    # Figure 9: Clustering Comparison (K-Means vs Hierarchical)
    # -------------------------------------------------------------
    if os.path.exists("results/clustering_results.csv") and os.path.exists("results/pca_results.csv"):
        clust_df = pd.read_csv("results/clustering_results.csv")
        pca_df = pd.read_csv("results/pca_results.csv")
        merged_c = clust_df.merge(pca_df[['protein_id', 'PC1', 'PC2']], on='protein_id')
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        sns.scatterplot(data=merged_c, x='PC1', y='PC2', hue='KMeans_k4_cluster', palette='Set1', s=85, edgecolor='black', ax=axes[0])
        axes[0].set_title("A. K-Means (K=4) Cluster Partition in PCA Space", weight='bold')
        axes[0].set_xlabel("PC1")
        axes[0].set_ylabel("PC2")
        axes[0].legend(title="KMeans Cluster")
        
        sns.scatterplot(data=merged_c, x='PC1', y='PC2', hue='Hierarchical_Ward_k4', palette='Set2', s=85, edgecolor='black', ax=axes[1])
        axes[1].set_title("B. Hierarchical Ward (K=4) Cluster Partition", weight='bold')
        axes[1].set_xlabel("PC1")
        axes[1].set_ylabel("PC2")
        axes[1].legend(title="Ward Cluster")
        plt.tight_layout()
        plt.savefig("figures/fig09_clustering_comparison_hierarchical.png")
        plt.close()
        print("Saved fig09_clustering_comparison_hierarchical.png")

    # -------------------------------------------------------------
    # Figure 10: ProtBERT Manifold Projections
    # -------------------------------------------------------------
    if os.path.exists("results/protbert_manifold_embeddings.csv"):
        pbert_df = pd.read_csv("results/protbert_manifold_embeddings.csv")
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
        
        sns.scatterplot(data=pbert_df, x='ProtBERT_PC1', y='ProtBERT_PC2', hue='temporal_class', palette=PALETTE_TEMPORAL, s=85, edgecolor='black', ax=axes[0])
        axes[0].set_title("A. ProtBERT PCA Projection", weight='bold')
        axes[0].set_xlabel("ProtBERT PC1")
        axes[0].set_ylabel("ProtBERT PC2")
        axes[0].legend(title="Temporal Class")
        
        sns.scatterplot(data=pbert_df, x='ProtBERT_tSNE1', y='ProtBERT_tSNE2', hue='temporal_class', palette=PALETTE_TEMPORAL, s=85, edgecolor='black', ax=axes[1])
        axes[1].set_title("B. ProtBERT t-SNE Projection", weight='bold')
        axes[1].set_xlabel("ProtBERT t-SNE 1")
        axes[1].set_ylabel("ProtBERT t-SNE 2")
        axes[1].legend(title="Temporal Class")
        
        sns.scatterplot(data=pbert_df, x='ProtBERT_UMAP1', y='ProtBERT_UMAP2', hue='temporal_class', palette=PALETTE_TEMPORAL, s=85, edgecolor='black', ax=axes[2])
        axes[2].set_title("C. ProtBERT UMAP Projection", weight='bold')
        axes[2].set_xlabel("ProtBERT UMAP 1")
        axes[2].set_ylabel("ProtBERT UMAP 2")
        axes[2].legend(title="Temporal Class")
        
        plt.tight_layout()
        plt.savefig("figures/fig10_protbert_manifold_projections.png")
        plt.close()
        print("Saved fig10_protbert_manifold_projections.png")

    # -------------------------------------------------------------
    # Figure 11: ProtBERT vs Physicochemical Clustering Benchmarks
    # -------------------------------------------------------------
    if os.path.exists("results/representation_comparison.csv"):
        comp_df = pd.read_csv("results/representation_comparison.csv")
        k4_df = comp_df[comp_df['Evaluation_Setting'].str.contains('K=4')].copy()
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        
        metrics_plot = ['Silhouette_Score', 'Temporal_ARI', 'Temporal_Purity', 'Functional_ARI']
        labels_plot = ['Silhouette', 'Temporal ARI', 'Temporal Purity', 'Functional ARI']
        
        df_melt = pd.melt(k4_df, id_vars=['Representation'], value_vars=metrics_plot, var_name='Metric', value_name='Score')
        df_melt['Metric'] = df_melt['Metric'].map(dict(zip(metrics_plot, labels_plot)))
        
        sns.barplot(data=df_melt, x='Metric', y='Score', hue='Representation', palette='Set2', ax=axes[0], edgecolor='black')
        axes[0].set_title("A. Clustering Validation Metrics at K=4", weight='bold')
        axes[0].set_ylabel("Metric Value")
        axes[0].set_xlabel("")
        axes[0].legend(title="Representation")
        
        # Best K comparison
        best_df = comp_df[comp_df['Evaluation_Setting'].str.contains('Best K')].copy()
        df_melt_best = pd.melt(best_df, id_vars=['Representation'], value_vars=metrics_plot, var_name='Metric', value_name='Score')
        df_melt_best['Metric'] = df_melt_best['Metric'].map(dict(zip(metrics_plot, labels_plot)))
        
        sns.barplot(data=df_melt_best, x='Metric', y='Score', hue='Representation', palette='Set2', ax=axes[1], edgecolor='black')
        axes[1].set_title("B. Clustering Validation Metrics at Optimal K", weight='bold')
        axes[1].set_ylabel("Metric Value")
        axes[1].set_xlabel("")
        axes[1].legend(title="Representation")
        
        plt.tight_layout()
        plt.savefig("figures/fig11_protbert_vs_physicochem_clustering.png")
        plt.close()
        print("Saved fig11_protbert_vs_physicochem_clustering.png")

    # -------------------------------------------------------------
    # Figure 12: Supervised Classification Benchmarks
    # -------------------------------------------------------------
    if os.path.exists("results/classification_results.csv"):
        clf_df = pd.read_csv("results/classification_results.csv")
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        
        sns.barplot(data=clf_df, x='Classifier', y='Accuracy_Mean', hue='Representation', palette='Set1', ax=axes[0], edgecolor='black')
        axes[0].set_title("A. Supervised Classification Accuracy (5-Fold CV)", weight='bold')
        axes[0].set_ylabel("Mean Accuracy")
        axes[0].set_ylim(0.4, 1.0)
        axes[0].legend(title="Representation")
        
        sns.barplot(data=clf_df, x='Classifier', y='Macro_F1_Mean', hue='Representation', palette='Set1', ax=axes[1], edgecolor='black')
        axes[1].set_title("B. Supervised Classification Macro-F1 (5-Fold CV)", weight='bold')
        axes[1].set_ylabel("Mean Macro-F1")
        axes[1].set_ylim(0.3, 1.0)
        axes[1].legend(title="Representation")
        
        plt.tight_layout()
        plt.savefig("figures/fig12_supervised_classification_benchmarks.png")
        plt.close()
        print("Saved fig12_supervised_classification_benchmarks.png")

    # -------------------------------------------------------------
    # Figure 13: Feature Ablation Delta Analysis
    # -------------------------------------------------------------
    if os.path.exists("results/ablation_results.csv"):
        abl_df = pd.read_csv("results/ablation_results.csv")
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        
        sns.barplot(data=abl_df, x='Representation_Subset', y='Macro_F1_Mean', palette='viridis', ax=axes[0], edgecolor='black')
        axes[0].set_title("A. Feature Ablation Macro-F1 Scores", weight='bold')
        axes[0].set_ylabel("Macro-F1")
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
        
        # Delta plot
        colors_delta = ['#27ae60' if d >= 0 else '#c0392b' for d in abl_df['Delta_Macro_F1_vs_Phys']]
        axes[1].bar(range(len(abl_df)), abl_df['Delta_Macro_F1_vs_Phys'], color=colors_delta, edgecolor='black')
        axes[1].axhline(0, color='black', linestyle='--', linewidth=1)
        axes[1].set_xticks(range(len(abl_df)))
        axes[1].set_xticklabels(abl_df['Representation_Subset'], rotation=30, ha='right')
        axes[1].set_title("B. Performance Delta vs Baseline Physicochemical (25-dim)", weight='bold')
        axes[1].set_ylabel("Delta Macro-F1")
        
        plt.tight_layout()
        plt.savefig("figures/fig13_ablation_delta_analysis.png")
        plt.close()
        print("Saved fig13_ablation_delta_analysis.png")

    # -------------------------------------------------------------
    # Figure 14: Confusion Matrices
    # -------------------------------------------------------------
    if os.path.exists("results/classification_confusion_matrices.csv"):
        cm_df = pd.read_csv("results/classification_confusion_matrices.csv")
        classes = ['Early', 'Immediate-Early', 'Late']
        
        # Plot confusion matrix for Random Forest across the 3 representations
        rf_records = cm_df[cm_df['Classifier'] == 'Random Forest']
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        
        for idx, (_, row) in enumerate(rf_records.iterrows()):
            cm_data = np.array(ast.literal_eval(row['Confusion_Matrix']))
            sns.heatmap(cm_data, annot=True, fmt='d', cmap='Blues', xticklabels=classes, yticklabels=classes, ax=axes[idx], cbar=False)
            axes[idx].set_title(f"{row['Representation']}\n(Random Forest)", weight='bold')
            axes[idx].set_xlabel("Predicted Label")
            axes[idx].set_ylabel("True Label")
            
        plt.tight_layout()
        plt.savefig("figures/fig14_confusion_matrices.png")
        plt.close()
        print("Saved fig14_confusion_matrices.png")

    # -------------------------------------------------------------
    # Figure 15: Outlier & Disagreement Profiles
    # -------------------------------------------------------------
    if os.path.exists("results/candidate_outliers.csv"):
        out_df = pd.read_csv("results/candidate_outliers.csv")
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
        
        disagree_counts = out_df['disagreement_flag'].value_counts()
        sns.barplot(x=disagree_counts.values, y=disagree_counts.index, ax=axes[0], palette='Reds_r', edgecolor='black')
        axes[0].set_title("A. Distribution of Disagreement / Uncertainty Flags", weight='bold')
        axes[0].set_xlabel("Protein Count")
        
        sns.scatterplot(data=out_df, x='instability_index', y='mw_kda', hue='true_temporal_class',
                        style='disagreement_flag', s=100, palette=PALETTE_TEMPORAL, ax=axes[1], edgecolor='black')
        axes[1].axvline(40, color='red', linestyle='--', label='Instability Threshold')
        axes[1].set_title("B. Candidate Outliers: Instability Index vs Molecular Weight", weight='bold')
        axes[1].set_xlabel("Instability Index")
        axes[1].set_ylabel("Molecular Weight (kDa)")
        axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig("figures/fig15_outlier_disagreement_profiles.png")
        plt.close()
        print("Saved fig15_outlier_disagreement_profiles.png")

    print("All 15 figures generation completed.")

if __name__ == "__main__":
    generate_all_figures()
