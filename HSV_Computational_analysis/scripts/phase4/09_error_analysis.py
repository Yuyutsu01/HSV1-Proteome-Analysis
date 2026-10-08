"""
Phase 4 - Step 09: Error Stratification & Biological Misclassification Analysis.

Biological & Computational Concept:
1. Error Analysis Objective (Section 16):
   To understand what biological features cause model confusion without manually altering
   or cherry-picking the dataset:
   - False Immediate-Early (IE predicted, True Early/Late)
   - False Early (Early predicted, True IE/Late)
   - False Late (Late predicted, True IE/Early)

2. Stratification Dimensions:
   - Gene Symbol / Protein annotation
   - Virus Species (HSV-1 vs HSV-2)
   - Sequence Length quartile bins
   - Homology Cluster ID
   - Prediction Confidence (max softmax probability)

Creates: results/tables/phase4_error_analysis.csv
"""

import os
import yaml
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from evaluation_utils import CLASS_LABELS


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_error_analysis():
    """Execute detailed error stratification on test predictions."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 09: ERROR STRATIFICATION & MISCLASSIFICATION ANALYSIS")
    print("=" * 80)

    # 1. Load Data
    df_sup = pd.read_csv(config["dataset"]["supervised_path"])
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    esm_data = np.load(config["features"]["esm2_path"], allow_pickle=True)
    X_full = esm_data['embeddings'].astype(np.float32)

    df = df_splits.merge(
        df_sup[['canonical_id', 'species', 'gene', 'protein', 'length']], 
        on='canonical_id'
    )

    # Evaluate on Homology-Aware test partition
    train_mask = (df['split_homology'] == 'TRAIN').values
    test_mask = (df['split_homology'] == 'TEST').values

    X_train = X_full[train_mask]
    y_train = df.loc[train_mask, 'temporal_class'].values
    X_test = X_full[test_mask]
    y_test = df.loc[test_mask, 'temporal_class'].values
    df_test = df.loc[test_mask].copy().reset_index(drop=True)

    # Fit Logistic Regression
    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_train)
    X_te_s = scaler.transform(X_test)

    lr = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
    lr.fit(X_tr_s, y_train)
    
    y_pred = lr.predict(X_te_s)
    y_prob = lr.predict_proba(X_te_s)

    df_test['predicted_class'] = y_pred
    df_test['prediction_confidence'] = np.max(y_prob, axis=1)
    df_test['is_correct'] = (df_test['temporal_class'] == df_test['predicted_class'])
    df_test['error_type'] = np.where(
        df_test['is_correct'], 
        'CORRECT', 
        df_test['temporal_class'] + '_PREDICTED_AS_' + df_test['predicted_class']
    )

    # Length bins
    df_test['length_bin'] = pd.qcut(df_test['length'], q=4, labels=['Q1_Short', 'Q2_Medium_Low', 'Q3_Medium_High', 'Q4_Long'])

    # Error summary by Gene
    gene_error = df_test.groupby('gene').agg(
        total=('canonical_id', 'count'),
        correct=('is_correct', 'sum'),
        accuracy=('is_correct', 'mean'),
        avg_confidence=('prediction_confidence', 'mean')
    ).reset_index()
    gene_error['error_rate'] = 1.0 - gene_error['accuracy']
    gene_error = gene_error.sort_values(by='total', ascending=False)

    # Error summary by Species
    species_error = df_test.groupby('species').agg(
        total=('canonical_id', 'count'),
        correct=('is_correct', 'sum'),
        accuracy=('is_correct', 'mean'),
        avg_confidence=('prediction_confidence', 'mean')
    ).reset_index()

    # Error summary by Length Bin
    len_error = df_test.groupby('length_bin', observed=False).agg(
        total=('canonical_id', 'count'),
        correct=('is_correct', 'sum'),
        accuracy=('is_correct', 'mean'),
        avg_confidence=('prediction_confidence', 'mean')
    ).reset_index()

    # Confusion breakdown table
    confusion_breakdown = df_test['error_type'].value_counts().reset_index()
    confusion_breakdown.columns = ['error_category', 'count']
    confusion_breakdown['percentage'] = (confusion_breakdown['count'] / len(df_test)) * 100

    print("\nError Summary by Species:")
    print(species_error.to_string(index=False))

    print("\nError Summary by Length Bin:")
    print(len_error.to_string(index=False))

    print("\nTop 10 Misclassification Categories:")
    print(confusion_breakdown.head(10).to_string(index=False))

    # Save complete record-level and summary tables
    out_dir = config["paths"]["phase4_tables_dir"]
    os.makedirs(out_dir, exist_ok=True)
    
    record_path = config["paths"]["error_analysis"]
    df_test.to_csv(record_path, index=False)
    print(f"\nSaved Detailed Error Analysis Table: {record_path} (N={len(df_test)})")

    gene_err_path = os.path.join(out_dir, "error_summary_by_gene.csv")
    gene_error.to_csv(gene_err_path, index=False)
    
    return df_test


if __name__ == "__main__":
    run_error_analysis()
