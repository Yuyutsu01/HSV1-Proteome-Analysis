"""
Phase 3 - Step 05: Representation Space Analysis, Dimensionality Reduction,
and Proxy / Confound Analysis.

Biological & Computational Concept:
1. Representation Space Exploration:
   To understand what structure is captured by each representation:
   - Amino Acid Composition (AAC, 20-dim)
   - Classical Physicochemical Descriptors (13-dim)
   - K-mer Frequencies (k=2, 400-dim)
   - ESM-2 Pretrained Protein Language Model Embeddings (320-dim)
   We apply Principal Component Analysis (PCA) to evaluate the primary axes of variance.

2. Confound / Proxy Risk in Viral Proteomics:
   A fundamental scientific risk in sequence-based temporal classification:
   Does a representation exhibit temporal structure because it genuinely captures temporal
   regulatory signals, or is temporal class merely confounded by:
   - Specific Gene / Protein Families (e.g. all UL23 TK sequences are Early)?
   - HSV Species divergence (HSV-1 vs HSV-2)?
   - Sequence Length (IE proteins are systematically longer on average)?
   - Homology clusters / isolate duplication?
   
   To answer this rigorously, we compute quantitative statistical associations:
   - One-way ANOVA F-statistic & Explained Variance (Eta-squared) across classes, genes, and species.
   - Spearman rank correlation with sequence length.
   - Silhouette scores measuring cluster tightness by temporal class, gene family, and species.
"""

import os
import sys
import yaml
import platform
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score


# Color palette for temporal classes
CLASS_PALETTE = {
    'IMMEDIATE_EARLY': '#E64B35', # Vermilion / Red
    'EARLY': '#4DBBD5',           # Blue
    'LATE': '#00A087'             # Teal / Green
}


def load_config(config_path: str = "configs/phase3_config.yaml") -> dict:
    """Load Phase 3 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def compute_eta_squared(values: np.ndarray, categories: pd.Series) -> float:
    """
    Compute Eta-squared (proportion of total variance explained by categorical factor).
    One-way ANOVA effect size: SS_between / SS_total.
    """
    df_temp = pd.DataFrame({'val': values, 'cat': categories}).dropna()
    if df_temp['cat'].nunique() <= 1:
        return 0.0
    
    grand_mean = df_temp['val'].mean()
    ss_total = ((df_temp['val'] - grand_mean) ** 2).sum()
    if ss_total == 0:
        return 0.0
    
    group_means = df_temp.groupby('cat')['val'].mean()
    group_sizes = df_temp.groupby('cat')['val'].count()
    ss_between = (group_sizes * ((group_means - grand_mean) ** 2)).sum()
    
    return round(float(ss_between / ss_total), 4)


def plot_pca(
    pca_df: pd.DataFrame, 
    rep_name: str, 
    var_exp: Tuple[float, float], 
    out_path: str
):
    """Generate and save publication-quality PCA scatter plot."""
    plt.figure(figsize=(9, 7), dpi=300)
    sns.set_theme(style="whitegrid", font_scale=1.1)
    
    order = ['IMMEDIATE_EARLY', 'EARLY', 'LATE']
    
    ax = sns.scatterplot(
        data=pca_df,
        x='PC1',
        y='PC2',
        hue='temporal_class',
        hue_order=order,
        palette=CLASS_PALETTE,
        alpha=0.6,
        s=30,
        edgecolor='none'
    )
    
    plt.title(f"Phase 3: PCA Representation Space — {rep_name}", fontsize=14, weight='bold', pad=12)
    plt.xlabel(f"PC1 ({var_exp[0]*100:.1f}% Variance Explained)", fontsize=12)
    plt.ylabel(f"PC2 ({var_exp[1]*100:.1f}% Variance Explained)", fontsize=12)
    plt.legend(title="Temporal Class", frameon=True, loc='best')
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved PCA Figure: {out_path}")


def analyze_representation_space(
    name: str,
    feats: np.ndarray,
    df_meta: pd.DataFrame,
    figures_dir: str,
    random_seed: int = 42
) -> Tuple[dict, dict, pd.DataFrame]:
    """
    Perform PCA, evaluate variance, calculate silhouette scores, and compute confound statistics.
    """
    print(f"\nAnalyzing representation space: {name} (Shape: {feats.shape})...")
    
    # 1. Standardize features
    scaler = StandardScaler()
    feats_scaled = scaler.fit_transform(feats)
    
    # 2. PCA
    pca = PCA(n_components=min(10, feats.shape[1]), random_state=random_seed)
    pca_coords = pca.fit_transform(feats_scaled)
    var_exp = pca.explained_variance_ratio_
    
    pc1_var = float(var_exp[0])
    pc2_var = float(var_exp[1]) if len(var_exp) > 1 else 0.0
    cum_var_2 = pc1_var + pc2_var
    
    # Build plotting DataFrame
    pca_df = pd.DataFrame({
        'canonical_id': df_meta['canonical_id'].values,
        'temporal_class': df_meta['temporal_class'].values,
        'species': df_meta['species'].values,
        'gene': df_meta['gene'].values,
        'length': df_meta['length'].values,
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1]
    })
    
    # Plot PCA
    plot_file = f"phase3_pca_{name.lower().replace(' ', '_').replace('-', '_')}.png"
    plot_path = os.path.join(figures_dir, plot_file)
    plot_pca(pca_df, name, (pc1_var, pc2_var), plot_path)
    
    # 3. Silhouette Score (Sampled for computational efficiency)
    sample_size = min(3000, len(feats_scaled))
    rng = np.random.RandomState(random_seed)
    sample_idx = rng.choice(len(feats_scaled), size=sample_size, replace=False)
    
    feats_sample = feats_scaled[sample_idx]
    tc_sample = df_meta['temporal_class'].values[sample_idx]
    gene_sample = df_meta['gene'].values[sample_idx]
    species_sample = df_meta['species'].values[sample_idx]
    
    sil_tc = round(float(silhouette_score(feats_sample, tc_sample)), 4)
    sil_gene = round(float(silhouette_score(feats_sample, gene_sample)), 4)
    sil_species = round(float(silhouette_score(feats_sample, species_sample)), 4)
    
    rep_summary = {
        'representation': name,
        'feature_dimension': feats.shape[1],
        'pc1_explained_variance': round(pc1_var, 4),
        'pc2_explained_variance': round(pc2_var, 4),
        'pc1_pc2_cumulative_variance': round(cum_var_2, 4),
        'silhouette_temporal_class': sil_tc,
        'silhouette_gene_identity': sil_gene,
        'silhouette_species': sil_species
    }
    
    # 4. Proxy / Confound Analysis: Variance explained in PC1 and PC2
    eta_tc_pc1 = compute_eta_squared(pca_coords[:, 0], df_meta['temporal_class'])
    eta_gene_pc1 = compute_eta_squared(pca_coords[:, 0], df_meta['gene'])
    eta_species_pc1 = compute_eta_squared(pca_coords[:, 0], df_meta['species'])
    
    rho_len_pc1, _ = stats.spearmanr(pca_coords[:, 0], df_meta['length'])
    rho_len_pc2, _ = stats.spearmanr(pca_coords[:, 1], df_meta['length'])
    r2_len_pc1 = round(float(rho_len_pc1 ** 2), 4)
    
    confound_row = {
        'representation': name,
        'pc1_eta2_temporal_class': eta_tc_pc1,
        'pc1_eta2_gene_identity': eta_gene_pc1,
        'pc1_eta2_species': eta_species_pc1,
        'pc1_r2_sequence_length': r2_len_pc1,
        'pc1_spearman_length': round(float(rho_len_pc1), 4),
        'pc2_spearman_length': round(float(rho_len_pc2), 4)
    }
    
    return rep_summary, confound_row, pca_df


def log_environment_metadata(out_path: str, config: dict):
    """Record execution environment metadata for full scientific reproducibility."""
    import torch
    import transformers
    import sklearn
    import scipy
    
    env_text = f"""================================================================================
PHASE 3 ENVIRONMENT & REPRODUCIBILITY MANIFEST
================================================================================
Timestamp: {pd.Timestamp.now().isoformat()}
Operating System: {platform.system()} {platform.release()} ({platform.version()})
Platform: {platform.platform()}
Processor: {platform.processor()}
Python Version: {sys.version}

Key Dependencies:
- PyTorch: {torch.__version__} (CUDA Available: {torch.cuda.is_available()})
- Transformers: {transformers.__version__}
- Scikit-learn: {sklearn.__version__}
- SciPy: {scipy.__version__}
- Pandas: {pd.__version__}
- NumPy: {np.__version__}
- Matplotlib: {plt.matplotlib.__version__}

Pretrained Models & Checkpoints:
- ESM-2 Checkpoint: {config['embeddings']['esm2']['model_name']}
- Embedding Dimension: {config['embeddings']['esm2']['embedding_dim']}
- Window Size: {config['embeddings']['esm2']['max_chunk_length']}
- Window Overlap: {config['embeddings']['esm2']['chunk_overlap']}

Homology Clustering & Partitioning:
- Thresholds Evaluated: {config['homology_clustering']['thresholds']}
- Primary Homology Threshold: {config['homology_clustering']['primary_threshold']}
- Word Size: {config['homology_clustering']['word_size']}
- Random Seed: {config['random_seed']}
- Split Proportions: Train={config['splitting']['train_ratio']}, Val={config['splitting']['val_ratio']}, Test={config['splitting']['test_ratio']}
================================================================================
"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(env_text)
    print(f"\nSaved Environment Log: {out_path}")


def run_representation_and_confound_pipeline():
    """Main execution function for Step 05."""
    config = load_config()
    
    print("=" * 80)
    print("PHASE 3 - STEP 05: REPRESENTATION SPACE & PROXY/CONFOUND ANALYSIS")
    print("=" * 80)
    
    # 1. Load Metadata & Sequences
    manifest_path = config["paths"]["sequence_manifest"]
    supervised_path = config["dataset"]["supervised_path"]
    
    df_manifest = pd.read_csv(manifest_path)
    df_sup = pd.read_csv(supervised_path)
    
    df_meta = df_manifest.merge(
        df_sup[['canonical_id', 'temporal_class', 'species', 'gene', 'protein']],
        on='canonical_id',
        how='inner'
    )
    df_meta['length'] = df_meta['normalized_sequence'].apply(len)
    
    # 2. Load Representations
    # A. Classical Physicochemical Features
    df_classical = pd.read_csv(config["paths"]["classical_features"])
    classical_cols = [c for c in df_classical.columns if c not in ['canonical_id', 'temporal_class']]
    feats_classical = df_classical[classical_cols].values
    
    # B. Amino Acid Composition (AAC)
    df_aac = pd.read_csv(config["paths"]["aac_features"])
    aac_cols = [c for c in df_aac.columns if c not in ['canonical_id', 'temporal_class']]
    feats_aac = df_aac[aac_cols].values
    
    # C. K-mer Features (k=2)
    kmer_path = os.path.join(config["paths"]["kmer_dir"], "kmer_k2_features.npz")
    kmer_data = np.load(kmer_path, allow_pickle=True)
    feats_kmer = kmer_data['features']
    
    # D. ESM-2 Embeddings
    esm_path = os.path.join(config["paths"]["embeddings_dir"], "esm2", "esm2_embeddings.npz")
    if not os.path.exists(esm_path):
        raise FileNotFoundError(f"ESM-2 embeddings not found at: {esm_path}. Run Step 04 first.")
    esm_data = np.load(esm_path, allow_pickle=True)
    feats_esm2 = esm_data['embeddings']
    
    representations = [
        ("Classical Features", feats_classical),
        ("AAC", feats_aac),
        ("K-mer (k=2)", feats_kmer),
        ("ESM-2", feats_esm2)
    ]
    
    figures_dir = config["paths"]["figures_dir"]
    rep_summaries = []
    confound_rows = []
    
    for name, feats in representations:
        summary, conf, _ = analyze_representation_space(
            name=name,
            feats=feats,
            df_meta=df_meta,
            figures_dir=figures_dir,
            random_seed=config["random_seed"]
        )
        rep_summaries.append(summary)
        confound_rows.append(conf)
        
    # Save Representation Summary Table
    df_rep_summary = pd.DataFrame(rep_summaries)
    rep_summary_path = config["paths"]["representation_summary"]
    os.makedirs(os.path.dirname(rep_summary_path), exist_ok=True)
    df_rep_summary.to_csv(rep_summary_path, index=False)
    print(f"\nSaved Representation Summary Table: {rep_summary_path}")
    print(df_rep_summary.to_string(index=False))
    
    # Save Proxy / Confound Analysis Table
    df_confounds = pd.DataFrame(confound_rows)
    confounds_path = config["paths"]["proxy_confound_analysis"]
    df_confounds.to_csv(confounds_path, index=False)
    print(f"\nSaved Proxy / Confound Analysis Table: {confounds_path}")
    print(df_confounds.to_string(index=False))
    
    # Save Environment Manifest
    env_path = config["paths"]["environment_log"]
    log_environment_metadata(env_path, config)
    
    print("\nStep 05 completed successfully.")


if __name__ == "__main__":
    run_representation_and_confound_pipeline()
