"""
Phase 3 - Step 03: Sequence Redundancy, Homology Clustering, Cluster Purity Analysis,
and Leakage-Free Train/Validation/Test Partitioning.

Scientific & Computational Concept:
1. Sequence Homology & Data Leakage:
   In viral proteomic datasets, multiple clinical isolates and strain variants of the same
   viral proteins exist (e.g. DNA polymerase UL30, thymidine kinase UL23, envelope glycoproteins).
   If sequences sharing high sequence identity are randomly partitioned into training and test
   sets, predictive models can achieve deceptively inflated accuracy through trivial sequence
   memorization (homology leakage) rather than learning generalizable temporal signatures.

2. Homology Clustering Architecture:
   - Deterministic CD-HIT style greedy incremental clustering on length-sorted sequences.
   - Evaluated at multiple sequence-identity thresholds: 90%, 70%, and 50% sequence identity.
   - Vectorized dipeptide/tripeptide compositional representation with length constraints.
   - Sequences sorted by length descending; longest sequence initiates a cluster representative.
   - For query sequences, candidates are filtered by length-ratio bounds and similarity.
   - Primary threshold (70% identity) defines the homology partitions for modeling.

3. Homology Cluster Purity & Biological Mixed Clusters:
   - For every cluster, calculate class proportions, Shannon entropy, and purity:
     purity = max(IE_count, Early_count, Late_count) / cluster_size
     class_entropy = - sum(p_c * log2(p_c))
   - Clusters containing multiple temporal classes are biologically informative: they reveal
     whether homologous proteins can exhibit distinct temporal regulation or functional divergence.

4. Dual Splitting Regimes:
   - Regime A (Random Stratified Split): Baseline random split preserving global temporal class ratio (70/15/15).
   - Regime B (Homology-Aware Split): Entire homology clusters are assigned to Train, Val, or Test.
     Zero cluster overlap between splits is mathematically guaranteed.
"""

import os
import math
import yaml
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from collections import Counter, defaultdict
from sklearn.preprocessing import normalize


def load_config(config_path: str = "configs/phase3_config.yaml") -> dict:
    """Load Phase 3 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def calculate_entropy(counts: List[int]) -> float:
    """Calculate Shannon entropy (base 2) for a distribution of class counts."""
    total = sum(counts)
    if total <= 0:
        return 0.0
    entropy = 0.0
    for count in counts:
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def cluster_sequences_vectorized(
    canonical_ids: np.ndarray,
    lengths: np.ndarray,
    feats: np.ndarray,
    threshold: float
) -> Tuple[Dict[str, str], List[dict]]:
    """
    Vectorized greedy CD-HIT clustering.
    Sorts sequences by length descending.
    Matches sequences against cluster representatives within length and similarity bounds.
    """
    # Normalize features for fast cosine similarity via dot product
    feats_norm = normalize(feats, norm='l2', axis=1)
    
    # Sort sequences by length descending
    sorted_order = np.argsort(-lengths)
    lengths_sorted = lengths[sorted_order]
    feats_sorted = feats_norm[sorted_order]
    cids_sorted = canonical_ids[sorted_order]
    
    rep_indices = [] # Indices relative to sorted_order
    cluster_member_cids = defaultdict(list)
    mapping = {}
    
    for i in range(len(sorted_order)):
        cid = cids_sorted[i]
        q_len = lengths_sorted[i]
        
        if len(rep_indices) == 0:
            rep_indices.append(i)
            c_id = f"CLUST_{threshold:.2f}_{len(rep_indices):05d}"
            mapping[cid] = c_id
            cluster_member_cids[c_id].append(cid)
            continue
            
        rep_idx_arr = np.array(rep_indices)
        rep_lens = lengths_sorted[rep_idx_arr]
        
        # Length ratio constraint: shorter sequence must be at least threshold * 0.70 of rep length
        len_mask = (q_len / rep_lens) >= (threshold * 0.70)
        
        assigned = False
        if np.any(len_mask):
            valid_rep_indices = rep_idx_arr[len_mask]
            sims = np.dot(feats_sorted[valid_rep_indices], feats_sorted[i])
            max_sim_idx = np.argmax(sims)
            max_sim = sims[max_sim_idx]
            
            if max_sim >= threshold:
                best_rep_sorted_idx = valid_rep_indices[max_sim_idx]
                rep_rank = rep_indices.index(best_rep_sorted_idx) + 1
                c_id = f"CLUST_{threshold:.2f}_{rep_rank:05d}"
                mapping[cid] = c_id
                cluster_member_cids[c_id].append(cid)
                assigned = True
                
        if not assigned:
            rep_indices.append(i)
            c_id = f"CLUST_{threshold:.2f}_{len(rep_indices):05d}"
            mapping[cid] = c_id
            cluster_member_cids[c_id].append(cid)
            
    # Assemble cluster summary metadata
    cluster_summaries = []
    for rank, rep_i in enumerate(rep_indices):
        c_id = f"CLUST_{threshold:.2f}_{rank+1:05d}"
        rep_cid = cids_sorted[rep_i]
        members = cluster_member_cids[c_id]
        cluster_summaries.append({
            'cluster_id': c_id,
            'representative_id': rep_cid,
            'members': members,
            'size': len(members)
        })
        
    return mapping, cluster_summaries


def partition_homology_aware(
    df: pd.DataFrame, 
    cluster_mapping: Dict[str, str],
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    random_seed: int = 42
) -> Dict[str, str]:
    """
    Partition dataset into TRAIN, VALIDATION, and TEST such that entire homology clusters
    are assigned to a single partition (guaranteeing zero cluster leakage).
    Uses stratified cluster distribution to balance class proportions across TRAIN, VALIDATION, and TEST.
    """
    np.random.seed(random_seed)
    
    cluster_data = defaultdict(lambda: {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0, 'members': []})
    for _, row in df.iterrows():
        cid = row['canonical_id']
        tclass = row['temporal_class']
        clust_id = cluster_mapping[cid]
        cluster_data[clust_id]['members'].append(cid)
        cluster_data[clust_id]['total'] += 1
        if tclass == 'IMMEDIATE_EARLY':
            cluster_data[clust_id]['IE'] += 1
        elif tclass == 'EARLY':
            cluster_data[clust_id]['Early'] += 1
        elif tclass == 'LATE':
            cluster_data[clust_id]['Late'] += 1
            
    total_seqs = len(df)
    target_train = int(total_seqs * train_ratio)
    target_val = int(total_seqs * val_ratio)
    target_test = total_seqs - target_train - target_val
    
    total_ie = sum(d['IE'] for d in cluster_data.values())
    total_early = sum(d['Early'] for d in cluster_data.values())
    total_late = sum(d['Late'] for d in cluster_data.values())
    
    target_dist = {
        'TRAIN': {'IE': total_ie * train_ratio, 'Early': total_early * train_ratio, 'Late': total_late * train_ratio, 'target_total': target_train},
        'VALIDATION': {'IE': total_ie * val_ratio, 'Early': total_early * val_ratio, 'Late': total_late * val_ratio, 'target_total': target_val},
        'TEST': {'IE': total_ie * test_ratio, 'Early': total_early * test_ratio, 'Late': total_late * test_ratio, 'target_total': target_test},
    }
    
    current_counts = {
        'TRAIN': {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0},
        'VALIDATION': {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0},
        'TEST': {'IE': 0, 'Early': 0, 'Late': 0, 'total': 0},
    }
    
    # Sort clusters by size descending (largest clusters allocated first)
    sorted_clusters = sorted(cluster_data.items(), key=lambda x: x[1]['total'], reverse=True)
    
    split_assignment = {}
    
    for clust_id, c_info in sorted_clusters:
        # Determine candidate partitions that still have capacity or minimal relative overshoot
        best_split = None
        best_score = float('inf')
        
        for split in ['TRAIN', 'VALIDATION', 'TEST']:
            current_tot = current_counts[split]['total']
            target_tot = target_dist[split]['target_total']
            
            # If a split is severely over capacity, penalize heavily
            capacity_ratio = (current_tot + c_info['total']) / (target_tot + 1e-5)
            
            new_ie = current_counts[split]['IE'] + c_info['IE']
            new_early = current_counts[split]['Early'] + c_info['Early']
            new_late = current_counts[split]['Late'] + c_info['Late']
            
            ie_deficit = max(0, target_dist[split]['IE'] - new_ie) / (target_dist[split]['IE'] + 1e-5)
            early_deficit = max(0, target_dist[split]['Early'] - new_early) / (target_dist[split]['Early'] + 1e-5)
            late_deficit = max(0, target_dist[split]['Late'] - new_late) / (target_dist[split]['Late'] + 1e-5)
            
            # Score balances filling deficit and staying close to target total
            score = (capacity_ratio * 3.0) + (ie_deficit * 2.0) + (early_deficit * 1.5) + (late_deficit * 1.0)
            
            if score < best_score:
                best_score = score
                best_split = split
                
        current_counts[best_split]['IE'] += c_info['IE']
        current_counts[best_split]['Early'] += c_info['Early']
        current_counts[best_split]['Late'] += c_info['Late']
        current_counts[best_split]['total'] += c_info['total']
        
        for cid in c_info['members']:
            split_assignment[cid] = best_split
            
    return split_assignment


def partition_random_stratified(
    df: pd.DataFrame, 
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    random_seed: int = 42
) -> Dict[str, str]:
    """
    Standard stratified random partitioning across temporal classes.
    """
    np.random.seed(random_seed)
    split_assignment = {}
    
    for tclass, group in df.groupby('temporal_class'):
        cids = group['canonical_id'].tolist()
        np.random.shuffle(cids)
        
        n = len(cids)
        n_train = int(n * train_ratio)
        n_val = int(n * val_ratio)
        
        train_ids = cids[:n_train]
        val_ids = cids[n_train:n_train + n_val]
        test_ids = cids[n_train + n_val:]
        
        for cid in train_ids:
            split_assignment[cid] = 'TRAIN'
        for cid in val_ids:
            split_assignment[cid] = 'VALIDATION'
        for cid in test_ids:
            split_assignment[cid] = 'TEST'
            
    return split_assignment


def run_homology_and_split_pipeline():
    """Main execution function for Step 03."""
    config = load_config()
    
    print("=" * 80)
    print("PHASE 3 - STEP 03: HOMOLOGY REDUNDANCY & DUAL SPLIT PIPELINE")
    print("=" * 80)
    
    # 1. Load Manifest & Supervised Dataset
    manifest_path = config["paths"]["sequence_manifest"]
    supervised_path = config["dataset"]["supervised_path"]
    
    df_manifest = pd.read_csv(manifest_path)
    df_sup = pd.read_csv(supervised_path)
    
    df = df_manifest.merge(
        df_sup[['canonical_id', 'temporal_class', 'species', 'gene', 'protein']],
        on='canonical_id',
        how='inner'
    )
    df['length'] = df['normalized_sequence'].apply(len)
    
    print(f"Dataset loaded: N={len(df)}")
    
    # 2. Load Precomputed K-mer Features (k=2) for fast similarity calculation
    kmer_k2_path = os.path.join(config["paths"]["kmer_dir"], "kmer_k2_features.npz")
    kmer_data = np.load(kmer_k2_path, allow_pickle=True)
    canonical_ids = kmer_data['canonical_ids']
    feats = kmer_data['features']
    lengths = df['length'].values
    
    # Ensure matching order
    assert np.array_equal(canonical_ids, df['canonical_id'].values), "Canonical ID order mismatch between manifest and kmer features!"
    
    # -------------------------------------------------------------
    # 3. Multi-Threshold Redundancy Analysis (90%, 70%, 50%)
    # -------------------------------------------------------------
    thresholds = config["homology_clustering"]["thresholds"]
    summary_rows = []
    cluster_mappings = {}
    
    for thresh in thresholds:
        print(f"\nClustering sequences at identity threshold = {thresh:.2f}...")
        mapping, clusters = cluster_sequences_vectorized(canonical_ids, lengths, feats, threshold=thresh)
        cluster_mappings[thresh] = (mapping, clusters)
        
        cluster_sizes = [c['size'] for c in clusters]
        singletons = sum(1 for s in cluster_sizes if s == 1)
        large_clusters = sum(1 for s in cluster_sizes if s >= 10)
        
        summary_rows.append({
            'identity_threshold': thresh,
            'total_clusters': len(clusters),
            'singleton_clusters': singletons,
            'singleton_fraction': round(singletons / len(clusters), 4),
            'large_clusters_ge_10': large_clusters,
            'max_cluster_size': max(cluster_sizes),
            'mean_cluster_size': round(float(np.mean(cluster_sizes)), 2),
            'median_cluster_size': float(np.median(cluster_sizes)),
            'total_sequences_clustered': len(df)
        })
        
    df_summary = pd.DataFrame(summary_rows)
    summary_path = config["paths"]["homology_cluster_summary"]
    os.makedirs(os.path.dirname(summary_path), exist_ok=True)
    df_summary.to_csv(summary_path, index=False)
    print(f"\nSaved Homology Cluster Summary: {summary_path}")
    print(df_summary.to_string(index=False))
    
    # -------------------------------------------------------------
    # 4. Cluster Purity Analysis at Primary Threshold (70%)
    # -------------------------------------------------------------
    primary_thresh = config["homology_clustering"]["primary_threshold"]
    primary_mapping, primary_clusters = cluster_mappings[primary_thresh]
    
    print(f"\nAnalyzing Cluster Purity and Class Composition at {primary_thresh*100:.0f}% identity...")
    composition_rows = []
    
    for c in primary_clusters:
        c_id = c['cluster_id']
        members = c['members']
        m_df = df[df['canonical_id'].isin(members)]
        
        counts = m_df['temporal_class'].value_counts()
        ie_cnt = counts.get('IMMEDIATE_EARLY', 0)
        early_cnt = counts.get('EARLY', 0)
        late_cnt = counts.get('LATE', 0)
        c_size = len(members)
        
        majority_class = counts.idxmax() if not counts.empty else 'UNKNOWN'
        majority_count = counts.max() if not counts.empty else 0
        purity = round(majority_count / c_size, 4)
        entropy = calculate_entropy([ie_cnt, early_cnt, late_cnt])
        is_mixed = (sum(1 for cnt in [ie_cnt, early_cnt, late_cnt] if cnt > 0) > 1)
        
        composition_rows.append({
            'cluster_id': c_id,
            'representative_id': c['representative_id'],
            'cluster_size': c_size,
            'IE_count': ie_cnt,
            'Early_count': early_cnt,
            'Late_count': late_cnt,
            'majority_class': majority_class,
            'purity': purity,
            'class_entropy': entropy,
            'is_mixed': is_mixed
        })
        
    df_comp = pd.DataFrame(composition_rows)
    comp_path = config["paths"]["homology_cluster_class_composition"]
    df_comp.to_csv(comp_path, index=False)
    print(f"Saved Cluster Composition & Purity Table: {comp_path}")
    
    mixed_count = df_comp['is_mixed'].sum()
    print(f"Primary Clusters: {len(df_comp)} | Pure Clusters: {len(df_comp) - mixed_count} | Mixed Clusters: {mixed_count} ({mixed_count/len(df_comp)*100:.2f}%)")
    
    # -------------------------------------------------------------
    # 5. Dual Splitting Manifest Generation
    # -------------------------------------------------------------
    print("\nGenerating Dual Splitting Manifest (Random Stratified vs Homology-Aware)...")
    split_random = partition_random_stratified(
        df, 
        train_ratio=config["splitting"]["train_ratio"],
        val_ratio=config["splitting"]["val_ratio"],
        test_ratio=config["splitting"]["test_ratio"],
        random_seed=config["random_seed"]
    )
    
    split_homology = partition_homology_aware(
        df, 
        primary_mapping,
        train_ratio=config["splitting"]["train_ratio"],
        val_ratio=config["splitting"]["val_ratio"],
        test_ratio=config["splitting"]["test_ratio"],
        random_seed=config["random_seed"]
    )
    
    manifest_rows = []
    for cid in df['canonical_id']:
        manifest_rows.append({
            'canonical_id': cid,
            'split_random': split_random[cid],
            'split_homology': split_homology[cid],
            'homology_cluster_id': primary_mapping[cid],
            'temporal_class': df.loc[df['canonical_id'] == cid, 'temporal_class'].values[0]
        })
        
    df_splits = pd.DataFrame(manifest_rows)
    split_manifest_path = config["paths"]["split_manifest"]
    df_splits.to_csv(split_manifest_path, index=False)
    print(f"Saved Dual Split Manifest: {split_manifest_path} (N={len(df_splits)})")
    
    # Invariant Verification
    print("\n--- SPLIT INVARIANT VERIFICATION ---")
    print(f"Total records in split manifest: {len(df_splits)} (Expected: {len(df)})")
    assert len(df_splits) == len(df), "Record count mismatch in split manifest!"
    assert df_splits['canonical_id'].nunique() == len(df), "Duplicate canonical IDs found!"
    
    # Verify zero cluster leakage in homology split
    train_clusters = set(df_splits[df_splits['split_homology'] == 'TRAIN']['homology_cluster_id'])
    val_clusters = set(df_splits[df_splits['split_homology'] == 'VALIDATION']['homology_cluster_id'])
    test_clusters = set(df_splits[df_splits['split_homology'] == 'TEST']['homology_cluster_id'])
    
    train_val_overlap = train_clusters.intersection(val_clusters)
    train_test_overlap = train_clusters.intersection(test_clusters)
    val_test_overlap = val_clusters.intersection(test_clusters)
    
    print(f"Homology Cluster Overlap - Train/Val: {len(train_val_overlap)}, Train/Test: {len(train_test_overlap)}, Val/Test: {len(val_test_overlap)}")
    assert len(train_val_overlap) == 0, f"Train/Val cluster leakage detected! {train_val_overlap}"
    assert len(train_test_overlap) == 0, f"Train/Test cluster leakage detected! {train_test_overlap}"
    assert len(val_test_overlap) == 0, f"Val/Test cluster leakage detected! {val_test_overlap}"
    print("Homology partition leakage check: PASSED (Zero cluster overlap).")
    
    print("\nRandom Stratified Split Counts:")
    print(pd.crosstab(df_splits['split_random'], df_splits['temporal_class'], margins=True))
    
    print("\nHomology-Aware Split Counts:")
    print(pd.crosstab(df_splits['split_homology'], df_splits['temporal_class'], margins=True))
    
    print("\nStep 03 completed successfully.")


if __name__ == "__main__":
    run_homology_and_split_pipeline()
