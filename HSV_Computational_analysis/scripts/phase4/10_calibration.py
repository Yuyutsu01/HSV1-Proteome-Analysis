"""
Phase 4 - Step 10: Model Calibration & Reliability Analysis.

Biological & Computational Concept:
1. Probability Calibration:
   In biological classification, predicted probability scores are often interpreted as confidence.
   However, raw uncalibrated outputs (e.g. from tree models or unregularized classifiers) can be
   overconfident.
   We compute:
   - Multiclass Brier Score (lower is better, measuring mean squared error of probability vectors)
   - Expected Calibration Error (ECE) per class
   - Reliability diagrams plotting predicted probability bins vs observed empirical frequency.

Creates:
- results/tables/phase4_calibration_summary.csv
- results/figures/phase4/phase4_calibration_curve.png
"""

import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import calibration_curve
from sklearn.preprocessing import StandardScaler
from evaluation_utils import CLASS_LABELS


def load_config(config_path: str = "configs/phase4_config.yaml") -> dict:
    """Load Phase 4 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def compute_multiclass_brier(y_true_indices: np.ndarray, y_prob: np.ndarray) -> float:
    """Compute overall multiclass Brier score."""
    n_samples = len(y_true_indices)
    one_hot = np.zeros_like(y_prob)
    for i, idx in enumerate(y_true_indices):
        one_hot[i, idx] = 1.0
    return float(np.mean(np.sum((y_prob - one_hot) ** 2, axis=1)))


def compute_ece(y_true_binary: np.ndarray, y_prob_binary: np.ndarray, n_bins: int = 10) -> float:
    """Compute Expected Calibration Error (ECE) for binary 1-vs-rest."""
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    total_samples = len(y_true_binary)
    
    for i in range(n_bins):
        bin_mask = (y_prob_binary >= bin_edges[i]) & (y_prob_binary < bin_edges[i+1])
        bin_count = np.sum(bin_mask)
        if bin_count > 0:
            bin_acc = np.mean(y_true_binary[bin_mask])
            bin_conf = np.mean(y_prob_binary[bin_mask])
            ece += (bin_count / total_samples) * abs(bin_acc - bin_conf)
            
    return float(ece)


def run_calibration_analysis():
    """Execute model calibration and reliability analysis."""
    config = load_config()
    print("=" * 80)
    print("PHASE 4 - STEP 10: MODEL CALIBRATION & RELIABILITY ANALYSIS")
    print("=" * 80)

    # 1. Load splits and ESM-2 embeddings
    df_splits = pd.read_csv(config["dataset"]["split_manifest_path"])
    esm_data = np.load(config["features"]["esm2_path"], allow_pickle=True)
    X_full = esm_data['embeddings'].astype(np.float32)

    # Use Homology-Aware partition
    train_mask = (df_splits['split_homology'] == 'TRAIN').values
    test_mask = (df_splits['split_homology'] == 'TEST').values

    X_train = X_full[train_mask]
    y_train = df_splits.loc[train_mask, 'temporal_class'].values
    X_test = X_full[test_mask]
    y_test = df_splits.loc[test_mask, 'temporal_class'].values

    class_to_idx = {c: i for i, c in enumerate(CLASS_LABELS)}
    y_test_idx = np.array([class_to_idx[c] for c in y_test])

    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_train)
    X_te_s = scaler.transform(X_test)

    # Train Logistic Regression and Random Forest
    lr = LogisticRegression(max_iter=1000, random_state=config["random_seed"], class_weight='balanced')
    lr.fit(X_tr_s, y_train)
    y_prob_lr = lr.predict_proba(X_te_s)

    rf = RandomForestClassifier(n_estimators=100, max_depth=15, random_state=config["random_seed"], class_weight='balanced', n_jobs=-1)
    rf.fit(X_tr_s, y_train)
    y_prob_rf = rf.predict_proba(X_te_s)

    # Calibration Summary Records
    brier_lr = compute_multiclass_brier(y_test_idx, y_prob_lr)
    brier_rf = compute_multiclass_brier(y_test_idx, y_prob_rf)

    calib_rows = []
    for model_name, y_probs, brier in [("Logistic Regression (ESM-2)", y_prob_lr, brier_lr), ("Random Forest (ESM-2)", y_prob_rf, brier_rf)]:
        for c_idx, c_name in enumerate(CLASS_LABELS):
            y_bin = (y_test_idx == c_idx).astype(int)
            prob_bin = y_probs[:, c_idx]
            ece = compute_ece(y_bin, prob_bin)
            calib_rows.append({
                'model': model_name,
                'temporal_class': c_name,
                'overall_brier_score': round(brier, 4),
                'expected_calibration_error_ece': round(ece, 4)
            })

    df_calib = pd.DataFrame(calib_rows)
    out_table = config["paths"]["calibration"]
    os.makedirs(os.path.dirname(out_table), exist_ok=True)
    df_calib.to_csv(out_table, index=False)
    print(f"\nSaved Calibration Summary Table: {out_table}")
    print(df_calib.to_string(index=False))

    # Plot Calibration Curves (Reliability Diagram)
    plt.figure(figsize=(10, 7), dpi=300)
    sns.set_theme(style="whitegrid", font_scale=1.1)
    
    colors = {'IMMEDIATE_EARLY': '#E64B35', 'EARLY': '#4DBBD5', 'LATE': '#00A087'}
    
    for c_idx, c_name in enumerate(CLASS_LABELS):
        y_bin = (y_test_idx == c_idx).astype(int)
        prob_bin = y_prob_lr[:, c_idx]
        fraction_of_positives, mean_predicted_value = calibration_curve(y_bin, prob_bin, n_bins=10)
        plt.plot(mean_predicted_value, fraction_of_positives, "s-", label=f"LR - {c_name}", color=colors[c_name], linewidth=2)

    plt.plot([0, 1], [0, 1], "k--", label="Perfect Calibration", alpha=0.7)
    plt.xlabel("Mean Predicted Probability", fontsize=12)
    plt.ylabel("Observed Class Fraction", fontsize=12)
    plt.title("Phase 4: Multi-Class Probability Calibration Curves (Homology-Aware Split)", fontsize=14, weight='bold', pad=12)
    plt.legend(frameon=True, loc='best')
    plt.tight_layout()

    fig_dir = config["paths"]["figures_dir"]
    os.makedirs(fig_dir, exist_ok=True)
    fig_path = os.path.join(fig_dir, "phase4_calibration_curve.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"Saved Calibration Figure: {fig_path}")

    return df_calib


if __name__ == "__main__":
    run_calibration_analysis()
