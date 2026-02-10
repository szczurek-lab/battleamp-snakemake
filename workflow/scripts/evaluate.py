"""
Evaluate model predictions against ground truth labels.

Supports classification and regression tasks. Metrics are written to a JSON
file for downstream aggregation.

Usage: called by Snakemake via script: directive
"""

import json
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
    average_precision_score,
    confusion_matrix,
    mean_squared_error,
    r2_score,
)
from scipy.stats import spearmanr


# ---------------------------------------------------------------------------
# Classification metrics
# ---------------------------------------------------------------------------


def compute_classification_metrics(y_true, y_pred, y_prob, requested_metrics):
    """Compute classification metrics.

    Args:
        y_true: binary ground truth (1 = positive, 0 = negative)
        y_pred: binary predictions (1 = positive, 0 = negative)
        y_prob: predicted probability of positive class
        requested_metrics: list of metric names to compute

    Returns:
        dict of metric_name -> value
    """
    results = {}
    n = len(y_true)

    if n == 0:
        return {m: None for m in requested_metrics}

    # Confusion matrix components
    tn, fp, fn, tp = 0, 0, 0, 0
    if len(set(y_true)) > 1 or len(set(y_pred)) > 1:
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()

    metric_funcs = {
        "accuracy": lambda: accuracy_score(y_true, y_pred),
        "mcc": lambda: matthews_corrcoef(y_true, y_pred),
        "f1": lambda: f1_score(y_true, y_pred, zero_division=0),
        "fpr": lambda: fp / (fp + tn) if (fp + tn) > 0 else None,
        "tpr": lambda: tp / (tp + fn) if (tp + fn) > 0 else None,
        "tnr": lambda: tn / (tn + fp) if (tn + fp) > 0 else None,
        "lr_plus": lambda: _lr_plus(tp, fn, fp, tn),
        "auroc": lambda: _safe_auroc(y_true, y_prob),
        "auprc": lambda: _safe_auprc(y_true, y_prob),
    }

    for metric_name in requested_metrics:
        if metric_name in metric_funcs:
            try:
                results[metric_name] = metric_funcs[metric_name]()
            except Exception as e:
                results[metric_name] = None
                print(
                    f"Warning: could not compute {metric_name}: {e}",
                    file=sys.stderr,
                )
        else:
            results[metric_name] = None
            print(
                f"Warning: unknown metric '{metric_name}'", file=sys.stderr
            )

    # Always include sample counts
    results["n_samples"] = n
    results["n_positive"] = int(y_true.sum())
    results["n_negative"] = int(n - y_true.sum())
    results["n_predicted_positive"] = int(y_pred.sum())

    return results


def _lr_plus(tp, fn, fp, tn):
    """Positive likelihood ratio = TPR / FPR."""
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    if fpr == 0:
        return float("inf") if tpr > 0 else None
    return tpr / fpr


def _safe_auroc(y_true, y_prob):
    """AUROC with safety check for single-class ground truth."""
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return roc_auc_score(y_true, y_prob)


def _safe_auprc(y_true, y_prob):
    """AUPRC with safety check."""
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return average_precision_score(y_true, y_prob)


# ---------------------------------------------------------------------------
# Regression metrics
# ---------------------------------------------------------------------------


def compute_regression_metrics(y_true, y_pred, requested_metrics):
    """Compute regression metrics.

    Args:
        y_true: ground truth MIC values (numeric, in target units)
        y_pred: predicted MIC values (numeric, in target units)
        requested_metrics: list of metric names

    Returns:
        dict of metric_name -> value
    """
    results = {}
    n = len(y_true)

    if n == 0:
        return {m: None for m in requested_metrics}

    metric_funcs = {
        "msle_ln": lambda: _msle_ln(y_true, y_pred),
        "spearman": lambda: spearmanr(y_true, y_pred).statistic,
        "r2": lambda: r2_score(y_true, y_pred),
        "rmse": lambda: np.sqrt(mean_squared_error(y_true, y_pred)),
        "mae": lambda: np.mean(np.abs(y_true - y_pred)),
    }

    for metric_name in requested_metrics:
        if metric_name in metric_funcs:
            try:
                results[metric_name] = metric_funcs[metric_name]()
            except Exception as e:
                results[metric_name] = None
                print(
                    f"Warning: could not compute {metric_name}: {e}",
                    file=sys.stderr,
                )
        else:
            results[metric_name] = None
            print(
                f"Warning: unknown metric '{metric_name}'", file=sys.stderr
            )

    results["n_samples"] = n
    return results


def _msle_ln(y_true, y_pred):
    """Mean squared log error using natural log.
    MSLE = mean((ln(y_true + 1) - ln(y_pred + 1))^2)
    """
    # Clip to avoid log of negative
    y_true_safe = np.maximum(y_true, 0)
    y_pred_safe = np.maximum(y_pred, 0)
    return np.mean(
        (np.log1p(y_true_safe) - np.log1p(y_pred_safe)) ** 2
    )


# ---------------------------------------------------------------------------
# MIC unit conversion
# ---------------------------------------------------------------------------

# Default average molecular weight for AMPs (rough estimate)
# Used when per-peptide MW is not available
DEFAULT_AMP_MW = 2500.0  # Da


def convert_mic_um_to_ugml(mic_um, molecular_weight=None):
    """Convert MIC from uM to ug/ml.
    ug/ml = uM * MW / 1000
    """
    mw = molecular_weight if molecular_weight is not None else DEFAULT_AMP_MW
    return mic_um * mw / 1000.0


def convert_mic_ugml_to_um(mic_ugml, molecular_weight=None):
    """Convert MIC from ug/ml to uM.
    uM = ug/ml * 1000 / MW
    """
    mw = molecular_weight if molecular_weight is not None else DEFAULT_AMP_MW
    return mic_ugml * 1000.0 / mw


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(snakemake):
    predictions_path = snakemake.input.predictions
    labels_path = snakemake.input.labels
    output_path = snakemake.output.metrics

    task_config = snakemake.params.task_config
    seq_col = snakemake.params.sequence_column
    task_type = task_config["type"]

    # Load data
    pred_df = pd.read_csv(predictions_path, sep="\t")
    label_df = pd.read_csv(labels_path, sep="\t")

    # Merge on sequence
    merged = pd.merge(
        pred_df, label_df, on=seq_col, how="inner", suffixes=("_pred", "_label")
    )

    report = {
        "task": task_config.get("description", ""),
        "task_type": task_type,
        "n_predictions": len(pred_df),
        "n_labels": len(label_df),
        "n_matched": len(merged),
        "n_unmatched_predictions": len(pred_df) - len(merged),
        "n_unmatched_labels": len(label_df) - len(merged),
    }

    if len(merged) == 0:
        report["metrics"] = {}
        report["error"] = "No matching sequences between predictions and labels"
        print(
            f"ERROR: No matching sequences between predictions and labels",
            file=sys.stderr,
        )
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)
        return

    requested_metrics = task_config.get("metrics", [])

    if task_type == "classification":
        label_col = task_config["label_column"]
        positive_label = task_config["positive_label"]

        # Binary encode
        y_true = (merged[label_col] == positive_label).astype(int).values

        # Handle prediction column
        if "Prediction" in merged.columns:
            y_pred = (merged["Prediction"] == "AMP").astype(int).values
        else:
            y_pred = np.zeros(len(merged), dtype=int)

        # Probability scores
        y_prob = None
        if "Probability_score" in merged.columns:
            y_prob = pd.to_numeric(
                merged["Probability_score"], errors="coerce"
            ).values

        report["metrics"] = compute_classification_metrics(
            y_true, y_pred, y_prob, requested_metrics
        )

    elif task_type == "regression":
        mic_col = task_config.get("mic_column", "MIC")
        target_unit = task_config.get("target_unit", "ug/ml")

        # Get predicted MIC
        y_pred = pd.to_numeric(merged["MIC"], errors="coerce").values

        # Handle unit conversion for predictions
        if "MIC_unit" in merged.columns:
            pred_units = merged["MIC_unit"].values
            for i in range(len(y_pred)):
                if pred_units[i] == "uM" and target_unit == "ug/ml":
                    y_pred[i] = convert_mic_um_to_ugml(y_pred[i])
                elif pred_units[i] == "ug/ml" and target_unit == "uM":
                    y_pred[i] = convert_mic_ugml_to_um(y_pred[i])

        # Get ground truth MIC
        y_true = pd.to_numeric(merged[mic_col], errors="coerce").values

        # Drop NaN pairs
        valid = ~(np.isnan(y_true) | np.isnan(y_pred))
        y_true = y_true[valid]
        y_pred = y_pred[valid]

        report["metrics"] = compute_regression_metrics(
            y_true, y_pred, requested_metrics
        )

    else:
        report["error"] = f"Unknown task type: {task_type}"
        report["metrics"] = {}

    # Write report
    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)

    print(
        f"Evaluation complete: {len(merged)} matched samples, "
        f"metrics: {list(report['metrics'].keys())}",
        file=sys.stderr,
    )


main(snakemake)
