"""
Phase 3 - Step 02: Amino Acid Composition (AAC), Classical Physicochemical Features,
and K-mer Frequency Representations.

Biological & Computational Concept:
1. Amino Acid Composition (AAC):
   Measures the relative frequency of each of the 20 standard amino acids in a protein sequence.
   Formula: freq(aa) = count(aa) / length
   Provides baseline composition information reflecting codon usage, structural class, and proteome-wide biophysical bias.

2. Classical Physicochemical Features:
   Biologically interpretable descriptors characterizing gross biochemical properties:
   - Molecular Weight: sum of standard residue masses minus water loss for peptide bonds.
   - Aromaticity (F, W, Y): tendency for hydrophobic core packing and UV absorbance.
   - Hydrophobicity (A, V, I, L, M, F, W, P): Kyte-Doolittle / Eisenberg classification.
   - Charge & Isoelectric tendencies:
     - Positive / Basic (K, R, H)
     - Negative / Acidic (D, E)
     - Charged (K, R, H, D, E)
     - Basic-to-Acidic ratio: (K + R + H) / (D + E + 1e-5) (smoothed to avoid division by zero)
   - Polarity (S, T, N, Q, Y, C): uncharged polar residues often found on protein surfaces.
   - Aliphatic (A, V, I, L): hydrophobic branched-chain stability.
   - Structural markers: Glycine (conformational flexibility), Proline (helix breaking / rigidity), Cysteine (disulfide bonding).

3. K-mer Frequency Representations (k=2, k=3):
   Captures local sequence context, short motifs, and dipeptide/tripeptide compositional bias without alignment.
   - k=2: 20^2 = 400 dipeptide features.
   - k=3: 20^3 = 8,000 tripeptide features.
   Stored efficiently in compressed sparse / dense NumPy (.npz) and Parquet formats with explicit vocabulary mappings.
"""

import os
import json
import yaml
import itertools
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from collections import Counter


# Standard 20 Amino Acids in canonical alphabetical order
STANDARD_AMINO_ACIDS = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
]

# Monoisotopic residue weights (average Da)
AA_WEIGHTS = {
    'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12, 'F': 147.18,
    'G': 57.05, 'H': 137.14, 'I': 113.16, 'K': 128.17, 'L': 113.16,
    'M': 131.20, 'N': 114.10, 'P': 97.12,  'Q': 128.13, 'R': 156.19,
    'S': 87.08,  'T': 101.11, 'V': 99.13,  'W': 186.21, 'Y': 163.18
}

# Physicochemical grouping sets
AROMATIC_AA = {'F', 'W', 'Y'}
HYDROPHOBIC_AA = {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'}
POSITIVE_AA = {'K', 'R', 'H'}
NEGATIVE_AA = {'D', 'E'}
CHARGED_AA = {'K', 'R', 'H', 'D', 'E'}
POLAR_AA = {'S', 'T', 'N', 'Q', 'Y', 'C'}
ALIPHATIC_AA = {'A', 'V', 'I', 'L'}


def load_config(config_path: str = "configs/phase3_config.yaml") -> dict:
    """Load Phase 3 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def compute_aac(sequence: str) -> Dict[str, float]:
    """Compute 20 standard amino-acid frequencies for a sequence."""
    seq_len = len(sequence)
    if seq_len == 0:
        return {aa: 0.0 for aa in STANDARD_AMINO_ACIDS}
    
    counts = Counter(sequence)
    # Frequencies normalized across sequence length
    return {aa: counts.get(aa, 0) / seq_len for aa in STANDARD_AMINO_ACIDS}


def compute_classical_features(sequence: str) -> Dict[str, float]:
    """
    Compute biologically interpretable sequence descriptors.
    Formulas are explicitly documented and deterministic.
    """
    seq_len = len(sequence)
    if seq_len == 0:
        return {
            'length': 0,
            'molecular_weight_estimate': 0.0,
            'aromatic_fraction': 0.0,
            'hydrophobic_fraction': 0.0,
            'charged_fraction': 0.0,
            'positive_fraction': 0.0,
            'negative_fraction': 0.0,
            'polar_fraction': 0.0,
            'aliphatic_fraction': 0.0,
            'glycine_fraction': 0.0,
            'proline_fraction': 0.0,
            'cysteine_fraction': 0.0,
            'basic_to_acidic_ratio': 0.0
        }
    
    counts = Counter(sequence)
    
    # Molecular weight: sum of residue masses + terminal H2O (18.015 Da)
    mw = sum(counts.get(aa, 0) * AA_WEIGHTS.get(aa, 110.0) for aa in counts) + 18.015
    
    # Fractions
    aromatic_cnt = sum(counts.get(aa, 0) for aa in AROMATIC_AA)
    hydrophobic_cnt = sum(counts.get(aa, 0) for aa in HYDROPHOBIC_AA)
    positive_cnt = sum(counts.get(aa, 0) for aa in POSITIVE_AA)
    negative_cnt = sum(counts.get(aa, 0) for aa in NEGATIVE_AA)
    charged_cnt = sum(counts.get(aa, 0) for aa in CHARGED_AA)
    polar_cnt = sum(counts.get(aa, 0) for aa in POLAR_AA)
    aliphatic_cnt = sum(counts.get(aa, 0) for aa in ALIPHATIC_AA)
    glycine_cnt = counts.get('G', 0)
    proline_cnt = counts.get('P', 0)
    cysteine_cnt = counts.get('C', 0)
    
    # Basic to acidic ratio with small epsilon smoothing
    basic_acidic_ratio = positive_cnt / (negative_cnt + 1e-5)
    
    return {
        'length': seq_len,
        'molecular_weight_estimate': round(mw, 2),
        'aromatic_fraction': round(aromatic_cnt / seq_len, 6),
        'hydrophobic_fraction': round(hydrophobic_cnt / seq_len, 6),
        'charged_fraction': round(charged_cnt / seq_len, 6),
        'positive_fraction': round(positive_cnt / seq_len, 6),
        'negative_fraction': round(negative_cnt / seq_len, 6),
        'polar_fraction': round(polar_cnt / seq_len, 6),
        'aliphatic_fraction': round(aliphatic_cnt / seq_len, 6),
        'glycine_fraction': round(glycine_cnt / seq_len, 6),
        'proline_fraction': round(proline_cnt / seq_len, 6),
        'cysteine_fraction': round(cysteine_cnt / seq_len, 6),
        'basic_to_acidic_ratio': round(basic_acidic_ratio, 6)
    }


def generate_kmer_vocab(k: int) -> List[str]:
    """Generate all 20^k standard amino-acid k-mers in lexicographical order."""
    return [''.join(p) for p in itertools.product(STANDARD_AMINO_ACIDS, repeat=k)]


def compute_kmer_frequencies(sequences: List[str], k: int, vocab: List[str]) -> np.ndarray:
    """
    Compute normalized k-mer frequencies for a list of sequences.
    Matrix shape: (N_sequences, 20^k).
    Each row sums to ~1.0 (for sequences of length >= k containing standard residues).
    """
    vocab_to_idx = {kmer: idx for idx, kmer in enumerate(vocab)}
    n_seqs = len(sequences)
    n_vocab = len(vocab)
    matrix = np.zeros((n_seqs, n_vocab), dtype=np.float32)
    
    for i, seq in enumerate(sequences):
        seq_len = len(seq)
        if seq_len < k:
            continue
        
        num_kmers = seq_len - k + 1
        kmer_counts = Counter(seq[j:j+k] for j in range(num_kmers))
        
        for kmer, count in kmer_counts.items():
            if kmer in vocab_to_idx:
                matrix[i, vocab_to_idx[kmer]] = count / num_kmers
                
    return matrix


def run_classical_and_kmer_feature_engineering():
    """Main execution function for Step 02."""
    config = load_config()
    
    print("=" * 80)
    print("PHASE 3 - STEP 02: AAC, CLASSICAL PHYSICOCHEMICAL & K-MER FEATURES")
    print("=" * 80)
    
    # 1. Load Normalized Sequence Manifest (from Step 01)
    manifest_path = config["paths"]["sequence_manifest"]
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Sequence manifest not found at: {manifest_path}. Run Step 01 first.")
    
    df_manifest = pd.read_csv(manifest_path)
    print(f"Loaded Sequence Manifest: {manifest_path} (N={len(df_manifest)})")
    
    # 2. Load Supervised Dataset to obtain ground-truth temporal class metadata
    supervised_path = config["dataset"]["supervised_path"]
    df_sup = pd.read_csv(supervised_path)
    
    # Merge temporal class and metadata
    df_merged = df_manifest.merge(
        df_sup[['canonical_id', 'temporal_class', 'species', 'gene', 'protein']], 
        on='canonical_id', 
        how='inner'
    )
    assert len(df_merged) == len(df_manifest), f"Mismatch in manifest merge: {len(df_merged)} vs {len(df_manifest)}"
    
    # -------------------------------------------------------------
    # 3. Compute Amino Acid Composition (AAC)
    # -------------------------------------------------------------
    print("\nComputing 20-standard Amino Acid Composition (AAC)...")
    aac_rows = []
    for idx, row in df_merged.iterrows():
        seq = row['normalized_sequence']
        aac_dict = compute_aac(seq)
        aac_dict['canonical_id'] = row['canonical_id']
        aac_dict['temporal_class'] = row['temporal_class']
        aac_rows.append(aac_dict)
        
    df_aac = pd.DataFrame(aac_rows)
    # Reorder columns: canonical_id, temporal_class, A, C, D, ...
    aac_cols = ['canonical_id', 'temporal_class'] + STANDARD_AMINO_ACIDS
    df_aac = df_aac[aac_cols]
    
    aac_out_path = config["paths"]["aac_features"]
    os.makedirs(os.path.dirname(aac_out_path), exist_ok=True)
    df_aac.to_csv(aac_out_path, index=False)
    print(f"Saved AAC Features: {aac_out_path} (Shape: {df_aac.shape})")
    
    # Validation check: Row sum approx 1.0
    aac_sums = df_aac[STANDARD_AMINO_ACIDS].sum(axis=1)
    print(f"AAC Row Sums - Min: {aac_sums.min():.4f}, Max: {aac_sums.max():.4f}, Mean: {aac_sums.mean():.4f}")
    
    # -------------------------------------------------------------
    # 4. Compute Additional Classical Physicochemical Descriptors
    # -------------------------------------------------------------
    print("\nComputing Biologically Interpretable Physicochemical Descriptors...")
    classical_rows = []
    for idx, row in df_merged.iterrows():
        seq = row['normalized_sequence']
        feat_dict = compute_classical_features(seq)
        feat_dict['canonical_id'] = row['canonical_id']
        feat_dict['temporal_class'] = row['temporal_class']
        classical_rows.append(feat_dict)
        
    df_classical = pd.DataFrame(classical_rows)
    classical_cols = [
        'canonical_id', 'temporal_class', 'length', 'molecular_weight_estimate',
        'aromatic_fraction', 'hydrophobic_fraction', 'charged_fraction',
        'positive_fraction', 'negative_fraction', 'polar_fraction',
        'aliphatic_fraction', 'glycine_fraction', 'proline_fraction',
        'cysteine_fraction', 'basic_to_acidic_ratio'
    ]
    df_classical = df_classical[classical_cols]
    
    classical_out_path = config["paths"]["classical_features"]
    df_classical.to_csv(classical_out_path, index=False)
    print(f"Saved Classical Physicochemical Features: {classical_out_path} (Shape: {df_classical.shape})")
    
    # -------------------------------------------------------------
    # 5. Compute K-mer Representations (k=2, k=3)
    # -------------------------------------------------------------
    kmer_dir = config["paths"]["kmer_dir"]
    os.makedirs(kmer_dir, exist_ok=True)
    
    k_values = config["features"].get("kmer_sizes", [2, 3])
    sequences = df_merged['normalized_sequence'].tolist()
    canonical_ids = df_merged['canonical_id'].values
    
    for k in k_values:
        print(f"\nComputing k-mer representations for k={k} ({20**k} features)...")
        vocab = generate_kmer_vocab(k)
        matrix = compute_kmer_frequencies(sequences, k, vocab)
        
        # Save vocabulary metadata
        vocab_path = os.path.join(kmer_dir, f"kmer_k{k}_vocab.json")
        with open(vocab_path, "w", encoding="utf-8") as f:
            json.dump({
                "k": k,
                "dimension": len(vocab),
                "vocabulary": vocab
            }, f, indent=2)
            
        # Save compressed binary .npz with canonical IDs and matrix
        matrix_path = os.path.join(kmer_dir, f"kmer_k{k}_features.npz")
        np.savez_compressed(
            matrix_path,
            canonical_ids=canonical_ids,
            features=matrix,
            k=k,
            feature_names=vocab
        )
        print(f"Saved k={k} k-mer matrix: {matrix_path} (Shape: {matrix.shape}, Size: {os.path.getsize(matrix_path)/(1024*1024):.2f} MB)")
        
    print("\nStep 02 feature engineering completed successfully.")


if __name__ == "__main__":
    run_classical_and_kmer_feature_engineering()
