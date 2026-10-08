"""
Phase 4 - Step 07: Length-Controlled Ablation Study (A1 - A9).

Biological & Computational Concept:
1. Length-Controlled Hypothesis:
   Phase 3 established that Immediate-Early and Early proteins have longer average sequence lengths
   than Late structural proteins. Does the addition of sequence composition/embeddings provide
   significant predictive signal beyond length alone, or does combining sequence features with length
   yield synergistic gains?

2. Required Ablation Models (Section 14):
   A1. Length only
   A2. AAC
   A3. AAC + Length
   A4. Classical Descriptors
   A5. Classical Descriptors + Length
   A6. 2-mer (Dipeptides)
   A7. 3-mer (Tripeptides)
   A8. ESM-2
   A9. ESM-2 + Length

3. Incremental Analysis:
   - Delta(Representation - Length)
   - Delta(Representation + Length - Representation)
   Evaluated with Logistic Regression (class-weighted) on Random and Homology-Aware partitions.
   Scaling (StandardScaler) fitted strictly on the training partition.

Creates: results/tables/phase4_length_control_analysis.csv
"""

import os
import yaml
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from evaluation_utils import compute_comprehensive_metrics


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_length_control_analysis():
    """Execute complete length-controlled ablation study across Random and Homology splits."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 07: LENGTH-CONTROLLED ABLATION STUDY (A1 - A9)")
    print("=" * 80)

    # 1. Load splits and sequence manifest for length
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    df_man = pd.read_csv(config["dataset"]["manifest_path"])
    
    df = df_splits.merge(df_man[['canonical_id', 'normalized_length']], on='canonical_id')
    lengths = df['normalized_length'].values.reshape(-1, 1).astype(np.float32)

    # 2. Load Feature Representations
    # AAC
    df_aac = pd.read_csv(config["features"]["aac_path"])
    aac_cols = [c for c in df_aac.columns if c not in ['canonical_id', 'temporal_class']]
    feats_aac = df_aac[aac_cols].values.astype(np.float32)

    # Classical
    df_class = pd.read_csv(config["features"]["classical_path"])
    class_cols = [c for c in df_class.columns if c not in ['canonical_id', 'temporal_class', 'length']]
    feats_class = df_class[class_cols].values.astype(np.float32)

    # 2-mer
    k2_data = np.load(config["features"]["kmer_k2_path"], allow_pickle=True)
    feats_k2 = k2_data['features'].astype(np.float32)

    # 3-mer
    k3_data = np.load(config["features"]["kmer_k3_path"], allow_pickle=True)
    feats_k3 = k3_data['features'].astype(np.float32)

    # ESM-2
    esm_data = np.load(config["features"]["esm2_path"], allow_pickle=True)
    feats_esm2 = esm_data['embeddings'].astype(np.float32)

    # Define the 9 Ablation Configurations
    ablations = [
        ("A1: Length Only", lengths),
        ("A2: AAC Only", feats_aac),
        ("A3: AAC + Length", np.hstack([feats_aac, lengths])),
        ("A4: Classical Descriptors Only", feats_class),
        ("A5: Classical Descriptors + Length", np.hstack([feats_class, lengths])),
        ("A6: 2-mer (Dipeptides)", feats_k2),
        ("A7: 3-mer (Tripeptides)", feats_k3),
        ("A8: ESM-2 Only", feats_esm2),
        ("A9: ESM-2 + Length", np.hstack([feats_esm2, lengths]))
    ]

    all_results = []

    for split_col, regime_name in [('split_random', 'Random Stratified'), ('split_homology', 'Homology-Aware')]:
        print(f"\n==================== REGIME: {regime_name} ====================")
        train_mask = (df[split_col] == 'TRAIN').values
        test_mask = (df[split_col] == 'TEST').values
        
        y_train = df.loc[train_mask, 'temporal_class'].values
        y_test = df.loc[test_mask, 'temporal_class'].values

        for abl_name, X_full in ablations:
            X_train = X_full[train_mask]
            X_test = X_full[test_mask]

            # Fit scaler strictly on training partition
            scaler = StandardScaler()
            X_tr_s = scaler.fit_transform(X_train)
            X_te_s = scaler.transform(X_test)

            lr = LogisticRegression(
                max_iter=1000, 
                random_state=config["random_seed"], 
                class_weight='balanced',
                solver='lbfgs'
            )
            lr.fit(X_tr_s, y_train)
            y_pred = lr.predict(X_te_s)
            y_prob = lr.predict_proba(X_te_s)

            metrics = compute_comprehensive_metrics(
                y_true=y_test,
                y_pred=y_pred,
                y_prob=y_prob,
                representation_name=abl_name,
                model_name="Logistic Regression",
                split_regime=regime_name
            )
            metrics['feature_dimension'] = X_full.shape[1]
            all_results.append(metrics)
            print(f"  {abl_name:<35} -> Dim={X_full.shape[1]:<4} | Acc={metrics['accuracy']:.4f} | BalAcc={metrics['balanced_accuracy']:.4f} | MacroF1={metrics['macro_f1']:.4f} | MCC={metrics['mcc']:.4f}")

    df_res = pd.DataFrame(all_results)
    out_path = config["paths"]["length_control"]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_res.to_csv(out_path, index=False)
    print(f"\nSaved Length-Controlled Analysis Table: {out_path}")
    return df_res


if __name__ == "__main__":
    run_length_control_analysis()
