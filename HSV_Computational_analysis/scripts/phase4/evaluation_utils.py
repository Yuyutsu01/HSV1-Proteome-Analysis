"""
Phase 4 Shared Evaluation Utilities.

Biological & Computational Concept:
Provides standard, reproducible metric calculation functions across all Phase 4 models:
- Multiclass Accuracy, Balanced Accuracy, Macro-F1, Weighted-F1, Matthews Correlation Coefficient (MCC).
- Per-class Precision, Recall, and F1 for IMMEDIATE_EARLY, EARLY, and LATE.
- Confusion matrix computation and Brier score calibration metric.
- Leakage-safe scaling and preprocessing transformations.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Any
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    precision_score, recall_score, matthews_corrcoef,
    confusion_matrix, brier_score_loss
)


CLASS_LABELS = ["IMMEDIATE_EARLY", "EARLY", "LATE"]
CLASS_TO_IDX = {c: i for i, c in enumerate(CLASS_LABELS)}
IDX_TO_CLASS = {i: c for i, c in enumerate(CLASS_LABELS)}


def compute_comprehensive_metrics(
    y_true: np.ndarray, 
    y_pred: np.ndarray, 
    y_prob: np.ndarray = None,
    representation_name: str = "",
    model_name: str = "",
    split_regime: str = ""
) -> Dict[str, Any]:
    """
    Compute full suite of classification and calibration metrics for 3-class temporal modeling.
    """
    # Ensure numeric labels for calculation
    if isinstance(y_true[0], str):
        y_true_num = np.array([CLASS_TO_IDX[c] for c in y_true])
        y_pred_num = np.array([CLASS_TO_IDX[c] for c in y_pred])
    else:
        y_true_num = y_true
        y_pred_num = y_pred

    acc = accuracy_score(y_true_num, y_pred_num)
    bal_acc = balanced_accuracy_score(y_true_num, y_pred_num)
    macro_f1 = f1_score(y_true_num, y_pred_num, average="macro", zero_division=0)
    weighted_f1 = f1_score(y_true_num, y_pred_num, average="weighted", zero_division=0)
    mcc = matthews_corrcoef(y_true_num, y_pred_num)

    # Per-class metrics (labels 0: IE, 1: Early, 2: Late)
    prec_per_class = precision_score(y_true_num, y_pred_num, labels=[0, 1, 2], average=None, zero_division=0)
    rec_per_class = recall_score(y_true_num, y_pred_num, labels=[0, 1, 2], average=None, zero_division=0)
    f1_per_class = f1_score(y_true_num, y_pred_num, labels=[0, 1, 2], average=None, zero_division=0)

    # Multi-class Brier Score (if probabilities provided)
    brier = 0.0
    if y_prob is not None and y_prob.shape[1] == 3:
        # Sum of squared differences between probability and one-hot true vector
        one_hot_true = np.zeros((len(y_true_num), 3))
        for idx, val in enumerate(y_true_num):
            one_hot_true[idx, val] = 1.0
        brier = float(np.mean(np.sum((y_prob - one_hot_true) ** 2, axis=1)))

    cm = confusion_matrix(y_true_num, y_pred_num, labels=[0, 1, 2])

    metrics = {
        'representation': representation_name,
        'model': model_name,
        'split_regime': split_regime,
        'accuracy': round(float(acc), 4),
        'balanced_accuracy': round(float(bal_acc), 4),
        'macro_f1': round(float(macro_f1), 4),
        'weighted_f1': round(float(weighted_f1), 4),
        'mcc': round(float(mcc), 4),
        'ie_precision': round(float(prec_per_class[0]), 4),
        'ie_recall': round(float(rec_per_class[0]), 4),
        'ie_f1': round(float(f1_per_class[0]), 4),
        'early_precision': round(float(prec_per_class[1]), 4),
        'early_recall': round(float(rec_per_class[1]), 4),
        'early_f1': round(float(f1_per_class[1]), 4),
        'late_precision': round(float(prec_per_class[2]), 4),
        'late_recall': round(float(rec_per_class[2]), 4),
        'late_f1': round(float(f1_per_class[2]), 4),
        'brier_score': round(float(brier), 4),
        'confusion_matrix': cm.tolist()
    }
    return metrics
