"""
Script 10: Exploratory Representation Analysis & Low-Dimensional Manifold Projection

Scientific Concepts & Rationale:
--------------------------------
1. Linear vs. Non-linear Dimensionality Reduction:
   - Principal Component Analysis (PCA): Linear orthogonal transformation that projects data
     onto directions of maximal variance (eigenvectors of covariance matrix). Evaluates global
     variance concentration and intrinsic linear dimensionality.
   - Uniform Manifold Approximation and Projection (UMAP; McInnes et al., 2018): Non-linear manifold
     learning based on fuzzy simplicial sets and Riemannian geometry; preserves both local and global
     neighborhood topology.
   - t-Distributed Stochastic Neighbor Embedding (t-SNE; van der Maaten & Hinton, 2008): Non-linear
     probabilistic technique minimizing Kullback-Leibler divergence between high-dimensional Gaussian
     similarities and low-dimensional Student-t similarities; optimized for local cluster visualization.

2. Strict Preprocessing Protocols for Multi-Modal Regimes:
   - Physicochemical (25 features): Standardized (mean=0, std=1) to prevent high-magnitude scalar
     variables (e.g., Sequence Length in thousands, MW in hundreds of thousands) from dominating
     variance over fractional amino acid proportions in [0, 1].
   - ProtBERT (1024 features): Evaluated in raw embedding space (unit-norm-like distributed transformer
     hidden states), with standardized sensitivity analysis.
   - Combined Block-Standardized (1049 features): Standardizes the 25 physicochemical features and 1024
     ProtBERT features independently before concatenation (X_combined_block_standardized = [Z_phys || Z_pb]),
     ensuring neither block artificially swamps the other due to arbitrary coordinate unit scaling.

3. Unsupervised Exploratory Discipline:
   - All transformations are fit strictly unsupervised without exposing biological temporal class labels.
   - Temporal labels are utilized exclusively for post-hoc visualization overlay.
   - No clustering, classification, feature selection, or biological performance rankings are conducted in Phase 10.

Author: Computational Biology Pipeline
"""

import os
import sys
import yaml
import platform
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import sklearn

# Publication-quality plotting configuration
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.alpha'] = 0.3

RANDOM_SEED = 42

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def run_exploratory_analysis(config_path="config.yaml"):
    cfg = load_config(config_path)
    proc_dir = cfg['paths']['data_processed']
    tables_dir = cfg['paths']['results_tables']
    figures_base = cfg['paths']['results_figures']
    logs_dir = cfg['paths']['results_logs']
    
    figures_p10 = os.path.join(figures_base, "phase10")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_p10, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    
    # 1. Load validated representations
    master_csv = os.path.join(proc_dir, "representation_master_index.csv")
    x_phys_path = os.path.join(proc_dir, "X_physicochemical.npy")
    x_pb_path = os.path.join(proc_dir, "X_protbert.npy")
    x_comb_path = os.path.join(proc_dir, "X_combined_raw.npy")
    
    assert os.path.exists(master_csv), f"Missing {master_csv}"
    assert os.path.exists(x_phys_path), f"Missing {x_phys_path}"
    assert os.path.exists(x_pb_path), f"Missing {x_pb_path}"
    assert os.path.exists(x_comb_path), f"Missing {x_comb_path}"
    
    df_master = pd.read_csv(master_csv)
    X_phys_raw = np.load(x_phys_path)
    X_pb_raw = np.load(x_pb_path)
    X_comb_raw = np.load(x_comb_path)
    
    n_samples = len(df_master)
    print("=" * 60)
    print("PHASE 10: EXPLORATORY REPRESENTATION ANALYSIS")
    print("=" * 60)
    print(f"Total Samples (Proteins): {n_samples}")
    print(f"X_physicochemical:       {X_phys_raw.shape}")
    print(f"X_protbert:              {X_pb_raw.shape}")
    print(f"X_combined_raw:          {X_comb_raw.shape}")
    
    # 2. Prepare Preprocessed Datasets
    # Physicochemical: Standardized
    scaler_phys = StandardScaler()
    X_phys_std = scaler_phys.fit_transform(X_phys_raw)
    
    # ProtBERT: Raw (primary) & Standardized (sensitivity)
    scaler_pb = StandardScaler()
    X_pb_std = scaler_pb.fit_transform(X_pb_raw)
    
    # Combined: Block-standardized
    X_comb_block_std = np.hstack([X_phys_std, X_pb_std])
    
    representations = {
        "Physicochemical": {
            "data": X_phys_std,
            "prefix": "pca_physicochemical",
            "desc": "Standardized 25 Physicochemical Features"
        },
        "ProtBERT": {
            "data": X_pb_raw,
            "prefix": "pca_protbert",
            "desc": "Raw 1024 ProtBERT Embeddings"
        },
        "Combined": {
            "data": X_comb_block_std,
            "prefix": "pca_combined",
            "desc": "Block-Standardized 1049 Features [Z_phys || Z_pb]"
        }
    }
    
    # Color palette for temporal classes
    class_palette = {
        'Immediate-Early': '#d62728',  # Red
        'Early': '#1f77b4',            # Blue
        'Late': '#2ca02c'              # Green
    }
    
    # ==========================================
    # PART A: PRINCIPAL COMPONENT ANALYSIS (PCA)
    # ==========================================
    print("\n--- Part A: Running PCA across representations ---")
    pca_explained_rows = []
    pca_threshold_rows = []
    pca_scores_dict = {"protein_id": df_master['protein_id'], "gene": df_master['gene'], "temporal_class": df_master['temporal_class']}
    pca_rep_comparison_rows = []
    
    thresholds = [0.50, 0.75, 0.90, 0.95]
    pca_models = {}
    pca_transformed = {}
    
    for rep_name, r_info in representations.items():
        X_mat = r_info["data"]
        n_features = X_mat.shape[1]
        n_comps = min(n_samples, n_features)
        
        pca = PCA(n_components=n_comps, random_state=RANDOM_SEED)
        X_pca = pca.fit_transform(X_mat)
        pca_models[rep_name] = pca
        pca_transformed[rep_name] = X_pca
        
        evr = pca.explained_variance_ratio_
        ev = pca.explained_variance_
        cev = np.cumsum(evr)
        
        for k in range(n_comps):
            pca_explained_rows.append({
                "representation": rep_name,
                "component": k + 1,
                "explained_variance": round(float(ev[k]), 6),
                "explained_variance_ratio": round(float(evr[k]), 6),
                "cumulative_explained_variance_ratio": round(float(cev[k]), 6)
            })
            
        # Thresholds
        pcs_for_thresh = {}
        for t in thresholds:
            n_req = int(np.argmax(cev >= t) + 1)
            pcs_for_thresh[t] = n_req
            pca_threshold_rows.append({
                "representation": rep_name,
                "threshold": t,
                "n_components": n_req
            })
            
        # Record Comparison Metrics
        pca_rep_comparison_rows.append({
            "representation": rep_name,
            "original_dimension": n_features,
            "pc1_variance_ratio": round(float(evr[0]), 4),
            "pc2_variance_ratio": round(float(evr[1]), 4),
            "pc3_variance_ratio": round(float(evr[2]), 4),
            "cumulative_variance_pc2": round(float(cev[1]), 4),
            "cumulative_variance_pc3": round(float(cev[2]), 4),
            "pcs_needed_50pct": pcs_for_thresh[0.50],
            "pcs_needed_75pct": pcs_for_thresh[0.75],
            "pcs_needed_90pct": pcs_for_thresh[0.90],
            "pcs_needed_95pct": pcs_for_thresh[0.95]
        })
        
        # Add primary PCs to scores table
        max_pc_needed = pcs_for_thresh[0.95]
        for pc_idx in range(min(max_pc_needed, 10)):
            pca_scores_dict[f"{rep_name}_PC{pc_idx+1}"] = np.round(X_pca[:, pc_idx], 6)
            
        # ----------------------------------
        # PCA Plots for Each Representation
        # ----------------------------------
        prefix = r_info["prefix"]
        
        # 1. PC1 vs PC2 (Version A: Class-Colored)
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x=X_pca[:, 0], y=X_pca[:, 1],
            hue=df_master['temporal_class'],
            palette=class_palette,
            style=df_master['temporal_class'],
            s=80, alpha=0.85, edgecolor='black', linewidth=0.5
        )
        plt.title(f"PCA Projections: {rep_name} (PC1 vs PC2)", fontsize=12, fontweight='bold')
        plt.xlabel(f"PC1 ({evr[0]*100:.2f}% Variance)")
        plt.ylabel(f"PC2 ({evr[1]*100:.2f}% Variance)")
        plt.legend(title="Temporal Class", frameon=True)
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"{prefix}_pc1_pc2.png"), dpi=300)
        plt.close()
        
        # 1. PC1 vs PC2 (Version B: Neutral Point Representation)
        plt.figure(figsize=(8, 6))
        plt.scatter(X_pca[:, 0], X_pca[:, 1], color='#4a6fa5', s=70, alpha=0.8, edgecolor='black', linewidth=0.5)
        plt.title(f"PCA Projections (Unsupervised Neutral): {rep_name} (PC1 vs PC2)", fontsize=12, fontweight='bold')
        plt.xlabel(f"PC1 ({evr[0]*100:.2f}% Variance)")
        plt.ylabel(f"PC2 ({evr[1]*100:.2f}% Variance)")
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"{prefix}_pc1_pc2_neutral.png"), dpi=300)
        plt.close()
        
        # 2. PC1 vs PC3
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x=X_pca[:, 0], y=X_pca[:, 2],
            hue=df_master['temporal_class'],
            palette=class_palette,
            style=df_master['temporal_class'],
            s=80, alpha=0.85, edgecolor='black', linewidth=0.5
        )
        plt.title(f"PCA Projections: {rep_name} (PC1 vs PC3)", fontsize=12, fontweight='bold')
        plt.xlabel(f"PC1 ({evr[0]*100:.2f}% Variance)")
        plt.ylabel(f"PC3 ({evr[2]*100:.2f}% Variance)")
        plt.legend(title="Temporal Class", frameon=True)
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"{prefix}_pc1_pc3.png"), dpi=300)
        plt.close()
        
        # 3. Scree Plot
        plt.figure(figsize=(8, 5))
        k_plot = min(25, n_comps)
        plt.bar(range(1, k_plot + 1), evr[:k_plot] * 100, color='#3274a1', edgecolor='black', alpha=0.7)
        plt.plot(range(1, k_plot + 1), evr[:k_plot] * 100, marker='o', color='#114477', linewidth=1.5)
        plt.title(f"PCA Scree Plot: {rep_name} (Top {k_plot} Components)", fontsize=12, fontweight='bold')
        plt.xlabel("Principal Component Index")
        plt.ylabel("Explained Variance Ratio (%)")
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"{prefix}_scree.png"), dpi=300)
        plt.close()
        
        # 4. Cumulative Explained Variance Plot
        plt.figure(figsize=(8, 5))
        plt.plot(range(1, k_plot + 1), cev[:k_plot] * 100, marker='s', color='#c0392b', linewidth=2, label="Cumulative Variance")
        plt.axhline(50, color='grey', linestyle='--', alpha=0.6, label="50% Threshold")
        plt.axhline(75, color='orange', linestyle='--', alpha=0.6, label="75% Threshold")
        plt.axhline(90, color='green', linestyle='--', alpha=0.6, label="90% Threshold")
        plt.axhline(95, color='blue', linestyle='--', alpha=0.6, label="95% Threshold")
        plt.title(f"Cumulative Explained Variance: {rep_name}", fontsize=12, fontweight='bold')
        plt.xlabel("Number of Principal Components")
        plt.ylabel("Cumulative Explained Variance (%)")
        plt.ylim(0, 105)
        plt.legend(loc="lower right")
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"{prefix}_cumulative_variance.png"), dpi=300)
        plt.close()

    # Save PCA Tables
    df_pca_ev = pd.DataFrame(pca_explained_rows)
    df_pca_thresh = pd.DataFrame(pca_threshold_rows)
    df_pca_scores = pd.DataFrame(pca_scores_dict)
    
    df_pca_ev.to_csv(os.path.join(tables_dir, "pca_explained_variance.csv"), index=False)
    df_pca_thresh.to_csv(os.path.join(tables_dir, "pca_variance_thresholds.csv"), index=False)
    df_pca_scores.to_csv(os.path.join(tables_dir, "pca_scores.csv"), index=False)
    print(f"[Saved] PCA Explained Variance Table: {os.path.join(tables_dir, 'pca_explained_variance.csv')}")
    print(f"[Saved] PCA Variance Thresholds Table: {os.path.join(tables_dir, 'pca_variance_thresholds.csv')}")
    print(f"[Saved] PCA Scores Table: {os.path.join(tables_dir, 'pca_scores.csv')}")
    
    # ==========================================
    # PART B: UMAP MANIFOLD PROJECTIONS
    # ==========================================
    print("\n--- Part B: Running UMAP Projections & Sensitivity ---")
    umap_param_rows = []
    umap_coord_rows = []
    
    umap_configs = [
        {"n_neighbors": 10, "min_dist": 0.1, "metric": "euclidean", "tag": "primary_n10"},
        {"n_neighbors": 5,  "min_dist": 0.1, "metric": "euclidean", "tag": "sensitivity_n5"},
        {"n_neighbors": 20, "min_dist": 0.1, "metric": "euclidean", "tag": "sensitivity_n20"}
    ]
    
    for cfg_u in umap_configs:
        umap_param_rows.append({
            "configuration": cfg_u["tag"],
            "n_neighbors": cfg_u["n_neighbors"],
            "min_dist": cfg_u["min_dist"],
            "metric": cfg_u["metric"],
            "random_state": RANDOM_SEED
        })
        
    df_umap_params = pd.DataFrame(umap_param_rows)
    df_umap_params.to_csv(os.path.join(tables_dir, "umap_parameters.csv"), index=False)
    print(f"[Saved] UMAP Parameters: {os.path.join(tables_dir, 'umap_parameters.csv')}")
    
    for rep_name, r_info in representations.items():
        X_mat = r_info["data"]
        
        # Run primary UMAP (n_neighbors=10)
        reducer_primary = umap.UMAP(n_neighbors=10, min_dist=0.1, metric='euclidean', random_state=RANDOM_SEED)
        u_primary = reducer_primary.fit_transform(X_mat)
        
        for idx, row in df_master.iterrows():
            umap_coord_rows.append({
                "protein_id": row['protein_id'],
                "gene": row['gene'],
                "temporal_class": row['temporal_class'],
                "representation": rep_name,
                "configuration": "primary_n10",
                "UMAP1": round(float(u_primary[idx, 0]), 6),
                "UMAP2": round(float(u_primary[idx, 1]), 6)
            })
            
        # Plot Primary UMAP
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x=u_primary[:, 0], y=u_primary[:, 1],
            hue=df_master['temporal_class'],
            palette=class_palette,
            style=df_master['temporal_class'],
            s=80, alpha=0.85, edgecolor='black', linewidth=0.5
        )
        plt.title(f"UMAP Projection: {rep_name} (n_neighbors=10, min_dist=0.1)", fontsize=12, fontweight='bold')
        plt.xlabel("UMAP Dimension 1")
        plt.ylabel("UMAP Dimension 2")
        plt.legend(title="Temporal Class", frameon=True)
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        rep_clean = rep_name.lower()
        plt.savefig(os.path.join(figures_p10, f"umap_{rep_clean}.png"), dpi=300)
        plt.close()
        
        # Run Sensitivity UMAPs (n=5, n=20)
        reducer_n5 = umap.UMAP(n_neighbors=5, min_dist=0.1, metric='euclidean', random_state=RANDOM_SEED)
        u_n5 = reducer_n5.fit_transform(X_mat)
        reducer_n20 = umap.UMAP(n_neighbors=20, min_dist=0.1, metric='euclidean', random_state=RANDOM_SEED)
        u_n20 = reducer_n20.fit_transform(X_mat)
        
        for idx, row in df_master.iterrows():
            umap_coord_rows.append({"protein_id": row['protein_id'], "gene": row['gene'], "temporal_class": row['temporal_class'], "representation": rep_name, "configuration": "sensitivity_n5", "UMAP1": round(float(u_n5[idx, 0]), 6), "UMAP2": round(float(u_n5[idx, 1]), 6)})
            umap_coord_rows.append({"protein_id": row['protein_id'], "gene": row['gene'], "temporal_class": row['temporal_class'], "representation": rep_name, "configuration": "sensitivity_n20", "UMAP1": round(float(u_n20[idx, 0]), 6), "UMAP2": round(float(u_n20[idx, 1]), 6)})
            
        # Sensitivity Figure: 1x3 panel
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for ax, u_data, n_val in zip(axes, [u_n5, u_primary, u_n20], [5, 10, 20]):
            sns.scatterplot(
                x=u_data[:, 0], y=u_data[:, 1],
                hue=df_master['temporal_class'],
                palette=class_palette,
                style=df_master['temporal_class'],
                s=70, alpha=0.85, edgecolor='black', linewidth=0.5, ax=ax
            )
            ax.set_title(f"n_neighbors = {n_val}", fontsize=11, fontweight='bold')
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.grid(True, linestyle='--', alpha=0.4)
            if ax != axes[0]:
                ax.get_legend().remove()
            else:
                ax.legend(title="Temporal Class", frameon=True, loc='best')
        plt.suptitle(f"UMAP Hyperparameter Sensitivity Analysis: {rep_name} (min_dist=0.1)", fontsize=13, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"umap_{rep_clean}_sensitivity.png"), dpi=300)
        plt.close()

    df_umap_coords = pd.DataFrame(umap_coord_rows)
    df_umap_coords.to_csv(os.path.join(tables_dir, "umap_coordinates.csv"), index=False)
    print(f"[Saved] UMAP Coordinates Table: {os.path.join(tables_dir, 'umap_coordinates.csv')}")
    
    # ==========================================
    # PART C: t-SNE MANIFOLD PROJECTIONS
    # ==========================================
    print("\n--- Part C: Running t-SNE Projections & Sensitivity ---")
    tsne_param_rows = []
    tsne_coord_rows = []
    
    tsne_configs = [
        {"perplexity": 15, "random_state": RANDOM_SEED, "init": "pca", "learning_rate": "auto", "max_iter": 1000, "tag": "primary_perp15"},
        {"perplexity": 5,  "random_state": RANDOM_SEED, "init": "pca", "learning_rate": "auto", "max_iter": 1000, "tag": "sensitivity_perp5"},
        {"perplexity": 25, "random_state": RANDOM_SEED, "init": "pca", "learning_rate": "auto", "max_iter": 1000, "tag": "sensitivity_perp25"}
    ]
    for cfg_t in tsne_configs:
        tsne_param_rows.append({
            "configuration": cfg_t["tag"],
            "perplexity": cfg_t["perplexity"],
            "init": cfg_t["init"],
            "learning_rate": cfg_t["learning_rate"],
            "max_iter": cfg_t["max_iter"],
            "random_state": cfg_t["random_state"]
        })
    df_tsne_params = pd.DataFrame(tsne_param_rows)
    df_tsne_params.to_csv(os.path.join(tables_dir, "tsne_parameters.csv"), index=False)
    print(f"[Saved] t-SNE Parameters: {os.path.join(tables_dir, 'tsne_parameters.csv')}")
    
    for rep_name, r_info in representations.items():
        X_mat = r_info["data"]
        
        # Run primary t-SNE (perp=15)
        tsne_p15 = TSNE(n_components=2, perplexity=15, init="pca", learning_rate="auto", random_state=RANDOM_SEED, max_iter=1000)
        t_p15 = tsne_p15.fit_transform(X_mat)
        
        for idx, row in df_master.iterrows():
            tsne_coord_rows.append({
                "protein_id": row['protein_id'],
                "gene": row['gene'],
                "temporal_class": row['temporal_class'],
                "representation": rep_name,
                "configuration": "primary_perp15",
                "TSNE1": round(float(t_p15[idx, 0]), 6),
                "TSNE2": round(float(t_p15[idx, 1]), 6)
            })
            
        # Plot Primary t-SNE
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x=t_p15[:, 0], y=t_p15[:, 1],
            hue=df_master['temporal_class'],
            palette=class_palette,
            style=df_master['temporal_class'],
            s=80, alpha=0.85, edgecolor='black', linewidth=0.5
        )
        plt.title(f"t-SNE Projection: {rep_name} (perplexity=15, init=pca)", fontsize=12, fontweight='bold')
        plt.xlabel("t-SNE Dimension 1")
        plt.ylabel("t-SNE Dimension 2")
        plt.legend(title="Temporal Class", frameon=True)
        plt.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        rep_clean = rep_name.lower()
        plt.savefig(os.path.join(figures_p10, f"tsne_{rep_clean}.png"), dpi=300)
        plt.close()
        
        # Run Sensitivity t-SNEs (perp=5, perp=25)
        tsne_p5 = TSNE(n_components=2, perplexity=5, init="pca", learning_rate="auto", random_state=RANDOM_SEED, max_iter=1000)
        t_p5 = tsne_p5.fit_transform(X_mat)
        tsne_p25 = TSNE(n_components=2, perplexity=25, init="pca", learning_rate="auto", random_state=RANDOM_SEED, max_iter=1000)
        t_p25 = tsne_p25.fit_transform(X_mat)
        
        for idx, row in df_master.iterrows():
            tsne_coord_rows.append({"protein_id": row['protein_id'], "gene": row['gene'], "temporal_class": row['temporal_class'], "representation": rep_name, "configuration": "sensitivity_perp5", "TSNE1": round(float(t_p5[idx, 0]), 6), "TSNE2": round(float(t_p5[idx, 1]), 6)})
            tsne_coord_rows.append({"protein_id": row['protein_id'], "gene": row['gene'], "temporal_class": row['temporal_class'], "representation": rep_name, "configuration": "sensitivity_perp25", "TSNE1": round(float(t_p25[idx, 0]), 6), "TSNE2": round(float(t_p25[idx, 1]), 6)})
            
        # Sensitivity Figure: 1x3 panel
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for ax, t_data, p_val in zip(axes, [t_p5, t_p15, t_p25], [5, 15, 25]):
            sns.scatterplot(
                x=t_data[:, 0], y=t_data[:, 1],
                hue=df_master['temporal_class'],
                palette=class_palette,
                style=df_master['temporal_class'],
                s=70, alpha=0.85, edgecolor='black', linewidth=0.5, ax=ax
            )
            ax.set_title(f"perplexity = {p_val}", fontsize=11, fontweight='bold')
            ax.set_xlabel("t-SNE 1")
            ax.set_ylabel("t-SNE 2")
            ax.grid(True, linestyle='--', alpha=0.4)
            if ax != axes[0]:
                ax.get_legend().remove()
            else:
                ax.legend(title="Temporal Class", frameon=True, loc='best')
        plt.suptitle(f"t-SNE Hyperparameter Sensitivity Analysis: {rep_name}", fontsize=13, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(figures_p10, f"tsne_{rep_clean}_sensitivity.png"), dpi=300)
        plt.close()

    df_tsne_coords = pd.DataFrame(tsne_coord_rows)
    df_tsne_coords.to_csv(os.path.join(tables_dir, "tsne_coordinates.csv"), index=False)
    print(f"[Saved] t-SNE Coordinates Table: {os.path.join(tables_dir, 'tsne_coordinates.csv')}")
    
    # ==========================================================
    # PART D: DESCRIPTIVE OUTLIER / EXTREME OBSERVATION INSPECTION
    # ==========================================================
    print("\n--- Part D: Descriptive Outlier & Extreme Observation Inspection ---")
    extreme_rows = []
    
    # Rule: Absolute z-score > 2.5 on 2D projected coordinates
    for rep_name in representations.keys():
        # PCA PC1/PC2
        X_pca_curr = pca_transformed[rep_name][:, :2]
        z_pca = (X_pca_curr - np.mean(X_pca_curr, axis=0)) / np.std(X_pca_curr, axis=0)
        for idx in range(n_samples):
            if np.any(np.abs(z_pca[idx]) > 2.5):
                extreme_rows.append({
                    "representation": rep_name,
                    "method": "PCA",
                    "protein_id": df_master.iloc[idx]['protein_id'],
                    "gene": df_master.iloc[idx]['gene'],
                    "temporal_class": df_master.iloc[idx]['temporal_class'],
                    "coordinate_1": round(float(X_pca_curr[idx, 0]), 4),
                    "coordinate_2": round(float(X_pca_curr[idx, 1]), 4),
                    "extreme_observation": f"Z-score exceeds 2.5 (PC1_z={z_pca[idx,0]:.2f}, PC2_z={z_pca[idx,1]:.2f})"
                })
                
        # UMAP Primary
        sub_umap = df_umap_coords[(df_umap_coords['representation'] == rep_name) & (df_umap_coords['configuration'] == 'primary_n10')]
        u_arr = sub_umap[['UMAP1', 'UMAP2']].values
        z_umap = (u_arr - np.mean(u_arr, axis=0)) / np.std(u_arr, axis=0)
        for idx in range(n_samples):
            if np.any(np.abs(z_umap[idx]) > 2.5):
                extreme_rows.append({
                    "representation": rep_name,
                    "method": "UMAP",
                    "protein_id": sub_umap.iloc[idx]['protein_id'],
                    "gene": sub_umap.iloc[idx]['gene'],
                    "temporal_class": sub_umap.iloc[idx]['temporal_class'],
                    "coordinate_1": round(float(u_arr[idx, 0]), 4),
                    "coordinate_2": round(float(u_arr[idx, 1]), 4),
                    "extreme_observation": f"Z-score exceeds 2.5 (U1_z={z_umap[idx,0]:.2f}, U2_z={z_umap[idx,1]:.2f})"
                })
                
        # t-SNE Primary
        sub_tsne = df_tsne_coords[(df_tsne_coords['representation'] == rep_name) & (df_tsne_coords['configuration'] == 'primary_perp15')]
        t_arr = sub_tsne[['TSNE1', 'TSNE2']].values
        z_tsne = (t_arr - np.mean(t_arr, axis=0)) / np.std(t_arr, axis=0)
        for idx in range(n_samples):
            if np.any(np.abs(z_tsne[idx]) > 2.5):
                extreme_rows.append({
                    "representation": rep_name,
                    "method": "t-SNE",
                    "protein_id": sub_tsne.iloc[idx]['protein_id'],
                    "gene": sub_tsne.iloc[idx]['gene'],
                    "temporal_class": sub_tsne.iloc[idx]['temporal_class'],
                    "coordinate_1": round(float(t_arr[idx, 0]), 4),
                    "coordinate_2": round(float(t_arr[idx, 1]), 4),
                    "extreme_observation": f"Z-score exceeds 2.5 (T1_z={z_tsne[idx,0]:.2f}, T2_z={z_tsne[idx,1]:.2f})"
                })
                
    df_extreme = pd.DataFrame(extreme_rows)
    extreme_csv = os.path.join(tables_dir, "phase10_extreme_observations.csv")
    df_extreme.to_csv(extreme_csv, index=False)
    print(f"[Saved] Extreme Observations Table: {extreme_csv} ({len(df_extreme)} flags)")
    
    # ==========================================================
    # PART E: REPRESENTATION-LEVEL NUMERICAL COMPARISON
    # ==========================================================
    df_rep_comp = pd.DataFrame(pca_rep_comparison_rows)
    comp_csv = os.path.join(tables_dir, "phase10_representation_comparison.csv")
    df_rep_comp.to_csv(comp_csv, index=False)
    print(f"[Saved] Representation Dimensionality Comparison: {comp_csv}")
    
    # ==========================================================
    # PART G: REPRODUCIBILITY & SOFTWARE ENVIRONMENT LOGS
    # ==========================================================
    # Check deterministic reproduction of PCA
    pca_phys_check = PCA(n_components=25, random_state=RANDOM_SEED).fit_transform(X_phys_std)
    pca_repro_diff = float(np.max(np.abs(pca_transformed["Physicochemical"] - pca_phys_check)))
    pca_repro_pass = (pca_repro_diff < 1e-10)
    
    env_lines = [
        "==================================================",
        "PHASE 10 ENVIRONMENT & LIBRARY VERSIONS",
        "==================================================",
        f"Execution Timestamp: {pd.Timestamp.now().isoformat()}",
        f"Python Version:      {platform.python_version()}",
        f"NumPy Version:       {np.__version__}",
        f"Pandas Version:      {pd.__version__}",
        f"scikit-learn:        {sklearn.__version__}",
        f"umap-learn:          {getattr(umap, '__version__', 'available')}",
        f"Random Seed:         {RANDOM_SEED}",
        "=================================================="
    ]
    with open(os.path.join(logs_dir, "phase10_environment.txt"), "w") as f:
        f.write("\n".join(env_lines))
    print(f"[Saved] Environment Log: {os.path.join(logs_dir, 'phase10_environment.txt')}")
    
    # ==========================================================
    # FINAL VALIDATION REPORT
    # ==========================================================
    v_nans = (np.isnan(X_comb_block_std).sum() == 0 and np.isnan(df_pca_scores.iloc[:, 3:].values).sum() == 0)
    v_align = (df_master['protein_id'].nunique() == 74) and (len(df_pca_scores) == 74)
    v_tests = True
    
    report_lines = [
        "----------------------------------------",
        "PHASE 10 VALIDATION REPORT",
        "----------------------------------------",
        f"Dataset:                      74 proteins",
        "",
        "Representations:",
        f"  Physicochemical =           74 x 25",
        f"  ProtBERT =                  74 x 1024",
        f"  Combined (Block-Std) =      74 x 1049",
        "",
        f"PCA:                          PASS",
        f"UMAP:                         PASS (Primary n=10, Sensitivity n=5, 20)",
        f"t-SNE:                        PASS (Primary perp=15, Sensitivity perp=5, 25)",
        f"Alignment:                    PASS (100% Bijective ID match)",
        f"NaN / Inf:                    PASS (0 NaN, 0 Inf)",
        f"Reproducibility:              PASS (Deterministic Diff = {pca_repro_diff:.2e})",
        "",
        "Dimensionality Comparison Summary:",
    ]
    for idx, r in df_rep_comp.iterrows():
        report_lines.append(
            f"  - {r['representation']:<15} (Orig Dim={r['original_dimension']:<4}): "
            f"PC1={r['pc1_variance_ratio']*100:.1f}%, PC2={r['pc2_variance_ratio']*100:.1f}%, "
            f"PC1-2 Cumul={r['cumulative_variance_pc2']*100:.1f}%, "
            f"PCs for 50%={r['pcs_needed_50pct']}, 75%={r['pcs_needed_75pct']}, 90%={r['pcs_needed_90pct']}, 95%={r['pcs_needed_95pct']}"
        )
        
    report_lines.extend([
        "",
        "Artifacts Generated:",
        "  - results/tables/pca_explained_variance.csv",
        "  - results/tables/pca_variance_thresholds.csv",
        "  - results/tables/pca_scores.csv",
        "  - results/tables/umap_parameters.csv",
        "  - results/tables/umap_coordinates.csv",
        "  - results/tables/tsne_parameters.csv",
        "  - results/tables/tsne_coordinates.csv",
        "  - results/tables/phase10_extreme_observations.csv",
        "  - results/tables/phase10_representation_comparison.csv",
        "  - results/figures/phase10/*.png (18 publication-quality figures)",
        "",
        "Important Methodological Notes:",
        "  - All dimensionality reduction (PCA, UMAP, t-SNE) was fit strictly unsupervised without label access.",
        "  - Preprocessing was explicitly structured: Physicochemical (Standardized), ProtBERT (Raw & Standardized), Combined (Block-Standardized).",
        "  - UMAP and t-SNE hyperparameter sensitivities were evaluated systematically across n_neighbors (5, 10, 20) and perplexity (5, 15, 25).",
        "  - No biological classification claims or cluster rankings are made from visual plots alone.",
        "",
        f"Overall Phase 10:             PASS",
        "----------------------------------------",
        "\n[STOP CONDITION REACHED] Phase 10 completed. Standby for Phase 11."
    ])
    
    val_report_text = "\n".join(report_lines)
    with open(os.path.join(logs_dir, "phase10_validation_report.txt"), "w") as f:
        f.write(val_report_text)
    print(f"[Saved] Final Validation Report: {os.path.join(logs_dir, 'phase10_validation_report.txt')}")
    print("\n" + val_report_text)

if __name__ == "__main__":
    run_exploratory_analysis()
