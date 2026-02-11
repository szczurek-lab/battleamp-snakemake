"""
Evaluate model predictions against ground truth labels.

Supports classification and regression tasks. Metrics are written to a JSON
file for downstream aggregation.

Classification labels TSV format:
    sequence    label
    ACDEF...    AMP
    GHIJK...    non-AMP

Regression labels TSV format:
    sequence    MIC         MIC_unit
    ACDEF...    64.0        ug/ml

Usage: called by Snakemake via script: directive
"""

import json
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
    average_precision_score,
    confusion_matrix,
    mean_squared_error,
    mean_squared_log_error,
    r2_score,
)
from scipy.stats import spearmanr


# ---------------------------------------------------------------------------
# Default metrics (matching original BattleAMP benchmark)
# ---------------------------------------------------------------------------

DEFAULT_CLASSIFICATION_METRICS = [
    "accuracy", "fpr", "tpr", "tnr", "informedness", "mcc",
    "precision", "recall", "f1", "pos_preds",
    "auroc", "auprc",
]

DEFAULT_REGRESSION_METRICS = [
    "r2", "mse", "rmse", "msle", "rmsle", "msl2e", "rmsl2e",
    "spearman",
]


# ---------------------------------------------------------------------------
# Classification metrics
# ---------------------------------------------------------------------------


def compute_classification_metrics(y_true, y_pred, y_prob, requested_metrics):
    results = {}
    n = len(y_true)

    if n == 0:
        return {m: None for m in requested_metrics}

    # Confusion matrix
    tn, fp, fn, tp = 0, 0, 0, 0
    if len(set(y_true)) > 1 or len(set(y_pred)) > 1:
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()

    def safe_div(num, denom):
        return float(num) / denom if denom > 0 else 0.0

    metric_funcs = {
        "accuracy": lambda: accuracy_score(y_true, y_pred),
        "mcc": lambda: matthews_corrcoef(y_true, y_pred),
        "f1": lambda: f1_score(y_true, y_pred, zero_division=0),
        "precision": lambda: safe_div(tp, tp + fp),
        "recall": lambda: safe_div(tp, tp + fn),
        "fpr": lambda: safe_div(fp, fp + tn),
        "tpr": lambda: safe_div(tp, tp + fn),
        "tnr": lambda: safe_div(tn, tn + fp),
        "informedness": lambda: balanced_accuracy_score(y_true, y_pred, adjusted=True),
        "pos_preds": lambda: safe_div(tp + fp, tp + tn + fp + fn),
        "auroc": lambda: _safe_auroc(y_true, y_prob),
        "auprc": lambda: _safe_auprc(y_true, y_prob),
    }

    for name in requested_metrics:
        if name in metric_funcs:
            try:
                results[name] = metric_funcs[name]()
            except Exception as e:
                results[name] = None
                print(f"Warning: could not compute {name}: {e}", file=sys.stderr)
        else:
            results[name] = None

    results["n_samples"] = n
    results["n_positive"] = int(y_true.sum())
    results["n_negative"] = int(n - y_true.sum())
    results["n_predicted_positive"] = int(y_pred.sum())
    results["errors"] = 0

    return results


def _safe_auroc(y_true, y_prob):
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return float(roc_auc_score(y_true, y_prob))


def _safe_auprc(y_true, y_prob):
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return float(average_precision_score(y_true, y_prob))


# ---------------------------------------------------------------------------
# Regression metrics
# ---------------------------------------------------------------------------


def compute_regression_metrics(y_true, y_pred, requested_metrics):
    results = {}
    n = len(y_true)

    if n == 0:
        return {m: None for m in requested_metrics}

    metric_funcs = {
        "r2": lambda: r2_score(y_true, y_pred),
        "mse": lambda: float(mean_squared_error(y_true, y_pred)),
        "rmse": lambda: float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "msle": lambda: _safe_msle(y_true, y_pred),
        "rmsle": lambda: _safe_rmsle(y_true, y_pred),
        "msl2e": lambda: float(np.mean(np.square(
            np.log2(1 + np.maximum(y_true, 0)) - np.log2(1 + np.maximum(y_pred, 0))
        ))),
        "rmsl2e": lambda: float(np.sqrt(np.mean(np.square(
            np.log2(1 + np.maximum(y_true, 0)) - np.log2(1 + np.maximum(y_pred, 0))
        )))),
        "spearman": lambda: float(spearmanr(y_true, y_pred).statistic),
        "mae": lambda: float(np.mean(np.abs(y_true - y_pred))),
    }

    for name in requested_metrics:
        if name in metric_funcs:
            try:
                results[name] = metric_funcs[name]()
            except Exception as e:
                results[name] = None
                print(f"Warning: could not compute {name}: {e}", file=sys.stderr)
        else:
            results[name] = None

    results["n_samples"] = n
    results["errors"] = 0
    return results


def _safe_msle(y_true, y_pred):
    yt = np.maximum(y_true, 0)
    yp = np.maximum(y_pred, 0)
    return float(mean_squared_log_error(yt, yp))


def _safe_rmsle(y_true, y_pred):
    return float(np.sqrt(_safe_msle(y_true, y_pred)))


# ---------------------------------------------------------------------------
# MIC unit conversion (per-peptide molecular weight)
# ---------------------------------------------------------------------------

# Average molecular weights of the 20 standard amino acids (Da).
AA_MW = {
    "A":  71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G":  57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P":  97.05276, "S":  87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V":  99.06841,
}

WATER_MW = 18.01056  # Da, lost per peptide bond formation + terminal groups


def peptide_mw(sequence):
    """Compute molecular weight (Da) from amino acid sequence.

    Sum of residue weights plus one water molecule (N-terminal H and
    C-terminal OH).  Non-standard residues are skipped with a warning.
    """
    mw = WATER_MW
    for aa in sequence.upper():
        mw += AA_MW.get(aa, 0.0)
    return mw


def convert_mic_um_to_ugml(mic_um, mw):
    """uM -> ug/mL: C[ug/mL] = C[uM] * MW[Da] / 1000"""
    return mic_um * mw / 1000.0


def convert_mic_ugml_to_um(mic_ugml, mw):
    """ug/mL -> uM: C[uM] = C[ug/mL] * 1000 / MW[Da]"""
    return mic_ugml * 1000.0 / mw


def harmonize_mic_units(values, units, sequences, target_unit):
    """Convert MIC values to the target unit using per-peptide MW.

    Parameters
    ----------
    values : np.ndarray       MIC values (float)
    units : np.ndarray        Unit strings per row ("uM" or "ug/ml")
    sequences : np.ndarray    Amino acid sequences (for MW calculation)
    target_unit : str         "uM" or "ug/ml"

    Returns
    -------
    np.ndarray  MIC values in target_unit
    int         Number of conversions performed
    """
    out = values.copy()
    n_converted = 0
    for i in range(len(out)):
        src = units[i]
        if src == target_unit:
            continue
        mw = peptide_mw(sequences[i])
        if src == "uM" and target_unit == "ug/ml":
            out[i] = convert_mic_um_to_ugml(out[i], mw)
        elif src == "ug/ml" and target_unit == "uM":
            out[i] = convert_mic_ugml_to_um(out[i], mw)
        n_converted += 1
    return out, n_converted


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
            f"ERROR: No matching sequences found (predictions={len(pred_df)}, "
            f"labels={len(label_df)})",
            file=sys.stderr,
        )
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)
        return

    # Get requested metrics or use defaults
    requested_metrics = task_config.get("metrics", [])

    if task_type == "classification":
        if not requested_metrics:
            requested_metrics = DEFAULT_CLASSIFICATION_METRICS

        # Convention: label column is "label", positive is "AMP"
        label_col = task_config.get("label_column", "label")
        positive_label = task_config.get("positive_label", "AMP")

        y_true = (merged[label_col] == positive_label).astype(int).values

        # Prediction column
        if "Prediction" in merged.columns:
            y_pred = (merged["Prediction"] == "AMP").astype(int).values
        else:
            y_pred = np.zeros(len(merged), dtype=int)
            print("Warning: no 'Prediction' column found, using all-negative",
                  file=sys.stderr)

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
        if not requested_metrics:
            requested_metrics = DEFAULT_REGRESSION_METRICS

        benchmark_unit = snakemake.params.benchmark_unit
        sequences = merged[seq_col].values

        # --- Predicted MIC ---
        mic_pred_col = "MIC_pred" if "MIC_pred" in merged.columns else "MIC"
        y_pred = pd.to_numeric(merged[mic_pred_col], errors="coerce").values

        pred_unit_col = (
            "MIC_unit_pred" if "MIC_unit_pred" in merged.columns
            else "MIC_unit" if "MIC_unit" in merged.columns
            else None
        )
        if pred_unit_col:
            y_pred, n_pred_conv = harmonize_mic_units(
                y_pred, merged[pred_unit_col].values, sequences, benchmark_unit
            )
        else:
            n_pred_conv = 0

        # --- Ground truth MIC ---
        mic_col = task_config.get("mic_column", "MIC")
        if mic_col in merged.columns:
            y_true = pd.to_numeric(merged[mic_col], errors="coerce").values
            true_unit_col = "MIC_unit" if "MIC_unit" in merged.columns else None
        elif "MIC_label" in merged.columns:
            y_true = pd.to_numeric(merged["MIC_label"], errors="coerce").values
            true_unit_col = (
                "MIC_unit_label" if "MIC_unit_label" in merged.columns else None
            )
        else:
            report["error"] = f"No MIC column found in merged data"
            report["metrics"] = {}
            with open(output_path, "w") as f:
                json.dump(report, f, indent=2)
            return

        if true_unit_col:
            y_true, n_true_conv = harmonize_mic_units(
                y_true, merged[true_unit_col].values, sequences, benchmark_unit
            )
        else:
            n_true_conv = 0

        # Drop NaN pairs
        valid = ~(np.isnan(y_true) | np.isnan(y_pred))
        y_true = y_true[valid]
        y_pred = y_pred[valid]

        report["benchmark_unit"] = benchmark_unit
        report["n_pred_converted"] = n_pred_conv
        report["n_label_converted"] = n_true_conv

        report["metrics"] = compute_regression_metrics(
            y_true, y_pred, requested_metrics
        )

    else:
        report["error"] = f"Unknown task type: {task_type}"
        report["metrics"] = {}

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)

    n_matched = report["n_matched"]
    metric_names = list(report["metrics"].keys())
    print(
        f"Evaluation complete: {n_matched} matched samples, "
        f"metrics: {metric_names}",
        file=sys.stderr,
    )


main(snakemake)
