"""
Evaluate model predictions against ground truth labels.

Supports classification and regression tasks, including cross-evaluation
where a regressor is evaluated on a classification task (MIC predictions
are thresholded into binary labels) or a classifier is evaluated on a
regression task (probability scores are reported but regression metrics
are skipped).

Regressor-on-classification evaluation:
    MIC predictions are converted to the threshold unit and binarized:
      - active (1):   MIC <= active_threshold
      - inactive (0): MIC >= inactive_threshold
      - grey zone:    active_threshold < MIC < inactive_threshold (excluded)
    The probability proxy for AUROC/AUPRC is -MIC (lower MIC = more active).

Regression evaluation:
    Predictions and ground truth are optionally clamped to a configurable
    assay range (mic_clamp) after unit conversion.  This prevents extreme
    outlier predictions from dominating error metrics.  Both linear-scale
    and log2-scale metrics are computed.

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
# Default metrics
# ---------------------------------------------------------------------------

DEFAULT_CLASSIFICATION_METRICS = [
    "mcc", "f1", "precision", "fpr", "tpr", "tnr",
    "auroc", "auprc",
    "precision_at_k",
    "pauroc_01", "pauroc_001",
]

DEFAULT_REGRESSION_METRICS = [
    "r2_log2",
    "msl2e", "rmsl2e",
    "spearman",
]

# Default k for precision@k.  Can be overridden per task via
# task_config["precision_at_k"]["k"].
DEFAULT_PRECISION_AT_K = 100


# ---------------------------------------------------------------------------
# Classification metrics
# ---------------------------------------------------------------------------


def compute_classification_metrics(
    y_true, y_pred, y_prob, requested_metrics, task_config=None,
):
    """Compute requested classification metrics.

    Parameters
    ----------
    y_true : np.ndarray (int)       Ground-truth binary labels (0/1).
    y_pred : np.ndarray (int)       Predicted binary labels (0/1).
    y_prob : np.ndarray or None     Continuous scores for ranking.
        For classifiers this is the predicted probability of the positive
        class.  For regressors evaluated on a classification task this is
        -MIC (lower MIC = higher score).
    requested_metrics : list[str]   Which metrics to compute.
    task_config : dict or None      Full task configuration.  Used to read
        per-task overrides such as ``precision_at_k.k``.

    Returns
    -------
    dict  Metric name -> value (or None if the metric could not be computed).
    """
    if task_config is None:
        task_config = {}

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

    # Read the k value for precision@k from the task config, falling back
    # to the module-level default.
    pak_config = task_config.get("precision_at_k", {})
    k = pak_config.get("k", DEFAULT_PRECISION_AT_K)

    metric_funcs = {
        "accuracy": lambda: accuracy_score(y_true, y_pred),
        "mcc": lambda: matthews_corrcoef(y_true, y_pred),
        "f1": lambda: f1_score(y_true, y_pred, zero_division=0),
        "precision": lambda: safe_div(tp, tp + fp),
        "fpr": lambda: safe_div(fp, fp + tn),
        "tpr": lambda: safe_div(tp, tp + fn),
        "tnr": lambda: safe_div(tn, tn + fp),
        "informedness": lambda: balanced_accuracy_score(y_true, y_pred, adjusted=True),
        "pos_preds": lambda: safe_div(tp + fp, tp + tn + fp + fn),
        "auroc": lambda: _safe_auroc(y_true, y_prob),
        "auprc": lambda: _safe_auprc(y_true, y_prob),
        "precision_at_k": lambda: _precision_at_k(y_true, y_prob, k),
        "pauroc_01": lambda: _safe_pauroc(y_true, y_prob, max_fpr=0.1),
        "pauroc_001": lambda: _safe_pauroc(y_true, y_prob, max_fpr=0.01),
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

    # If precision_at_k was computed, record the k that was used so it is
    # unambiguous in the output JSON.
    if "precision_at_k" in requested_metrics:
        results["precision_at_k_k"] = k

    results["n_samples"] = n
    results["n_positive"] = int(y_true.sum())
    results["n_negative"] = int(n - y_true.sum())
    results["n_predicted_positive"] = int(y_pred.sum())
    results["errors"] = 0

    return results


# ---------------------------------------------------------------------------
# Helpers for score-based classification metrics
# ---------------------------------------------------------------------------


def _safe_auroc(y_true, y_prob):
    """Full AUROC.  Returns None when undefined (single class or no scores)."""
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return float(roc_auc_score(y_true, y_prob))


def _safe_pauroc(y_true, y_prob, max_fpr):
    """Partial AUROC restricted to FPR <= max_fpr (McClish standardization).

    sklearn's ``roc_auc_score(max_fpr=...)`` computes the area under the ROC
    curve only up to the given false-positive rate and then applies the
    McClish standardization so the result is in [0, 1] (0.5 = random).

    Returns None when the metric is undefined.
    """
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return float(roc_auc_score(y_true, y_prob, max_fpr=max_fpr))


def _safe_auprc(y_true, y_prob):
    """Area under the precision-recall curve.  Returns None when undefined."""
    if len(set(y_true)) < 2:
        return None
    if y_prob is None or len(y_prob) == 0:
        return None
    return float(average_precision_score(y_true, y_prob))


def _precision_at_k(y_true, y_prob, k):
    """Precision among the top-k highest-scoring predictions.

    Sequences are ranked by ``y_prob`` in descending order (highest score
    first).  For classifiers ``y_prob`` is the predicted probability; for
    regressors evaluated on a classification task it is ``-MIC``, so lower
    MIC values rank first.

    If fewer than k scored samples are available, precision is computed
    over all of them and a warning is printed.

    Returns None if scores are unavailable.
    """
    if y_prob is None or len(y_prob) == 0:
        return None
    if k <= 0:
        return None

    # Only consider samples that have a finite score.
    valid = np.isfinite(y_prob)
    y_true_v = y_true[valid]
    y_prob_v = y_prob[valid]

    if len(y_prob_v) == 0:
        return None

    effective_k = min(k, len(y_prob_v))
    if effective_k < k:
        print(
            f"Warning: precision_at_k requested k={k} but only "
            f"{len(y_prob_v)} scored samples available; using k={effective_k}",
            file=sys.stderr,
        )

    top_idx = np.argsort(y_prob_v)[::-1][:effective_k]
    return float(y_true_v[top_idx].sum()) / effective_k


# ---------------------------------------------------------------------------
# Regression metrics
# ---------------------------------------------------------------------------


def compute_regression_metrics(y_true, y_pred, requested_metrics):
    results = {}
    n = len(y_true)

    if n == 0:
        return {m: None for m in requested_metrics}

    # Pre-compute log2 transforms for log-scale metrics.
    # log2(1 + x) avoids log(0) and is monotonic for x >= 0.
    yt_log2 = np.log2(1 + np.maximum(y_true, 0))
    yp_log2 = np.log2(1 + np.maximum(y_pred, 0))

    metric_funcs = {
        "r2": lambda: r2_score(y_true, y_pred),
        "r2_log2": lambda: r2_score(yt_log2, yp_log2),
        "mse": lambda: float(mean_squared_error(y_true, y_pred)),
        "rmse": lambda: float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "msle": lambda: _safe_msle(y_true, y_pred),
        "rmsle": lambda: _safe_rmsle(y_true, y_pred),
        "msl2e": lambda: float(np.mean(np.square(yt_log2 - yp_log2))),
        "rmsl2e": lambda: float(np.sqrt(np.mean(np.square(yt_log2 - yp_log2)))),
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

# Monoisotopic residue weights of the 20 standard amino acids (Da).
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


def compute_mw_array(sequences):
    """Compute molecular weights for an array of sequences."""
    return np.array([peptide_mw(s) for s in sequences])


def convert_mic_um_to_ugml(mic_um, mw):
    """uM -> ug/mL: C[ug/mL] = C[uM] * MW[Da] / 1000"""
    return mic_um * mw / 1000.0


def convert_mic_ugml_to_um(mic_ugml, mw):
    """ug/mL -> uM: C[uM] = C[ug/mL] * 1000 / MW[Da]"""
    return mic_ugml * 1000.0 / mw


def convert_mic_values(values, mw_array, src_unit, dst_unit):
    """Convert MIC values between units using per-peptide MW (vectorized).

    Parameters
    ----------
    values : np.ndarray       MIC values (float)
    mw_array : np.ndarray     Per-peptide molecular weights (Da)
    src_unit : str            Source unit ("uM" or "ug/ml")
    dst_unit : str            Destination unit ("uM" or "ug/ml")

    Returns
    -------
    np.ndarray  Converted MIC values
    """
    if src_unit == dst_unit:
        return values.copy()
    if src_unit == "uM" and dst_unit == "ug/ml":
        return values * mw_array / 1000.0
    if src_unit == "ug/ml" and dst_unit == "uM":
        return values * 1000.0 / mw_array
    raise ValueError(f"Unknown unit conversion: {src_unit} -> {dst_unit}")


def harmonize_mic_units(values, units, sequences, target_unit):
    """Convert MIC values to the target unit using per-peptide MW.

    Handles mixed-unit columns where each row may have a different unit.

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
# MIC clamping
# ---------------------------------------------------------------------------


def clamp_mic_values(values, sequences, current_unit, clamp_config):
    """Clamp MIC values to an assay-realistic range.

    The clamp range is defined in clamp_config["unit"]. If this differs
    from current_unit, thresholds are converted per-peptide using MW.

    Parameters
    ----------
    values : np.ndarray       MIC values in current_unit
    sequences : np.ndarray    Amino acid sequences (for MW if conversion needed)
    current_unit : str        Unit of values ("uM" or "ug/ml")
    clamp_config : dict       {"min": float, "max": float, "unit": str}

    Returns
    -------
    np.ndarray  Clamped MIC values in current_unit
    int         Number of values clamped
    """
    if clamp_config is None:
        return values.copy(), 0

    clamp_min = clamp_config["min"]
    clamp_max = clamp_config["max"]
    clamp_unit = clamp_config["unit"]

    is_valid = ~np.isnan(values)

    if current_unit == clamp_unit:
        # Direct comparison, same unit for all peptides.
        lo = np.full(len(values), clamp_min)
        hi = np.full(len(values), clamp_max)
    else:
        # Convert clamp bounds to current_unit per peptide.
        # E.g. clamp is in ug/ml but values are in uM: the uM equivalent
        # of 512 ug/ml depends on each peptide's MW.
        mw_array = compute_mw_array(sequences)
        lo = convert_mic_values(
            np.full(len(values), clamp_min), mw_array,
            clamp_unit, current_unit,
        )
        hi = convert_mic_values(
            np.full(len(values), clamp_max), mw_array,
            clamp_unit, current_unit,
        )

    out = values.copy()
    clamped_low = is_valid & (out < lo)
    clamped_high = is_valid & (out > hi)
    out[clamped_low] = lo[clamped_low]
    out[clamped_high] = hi[clamped_high]

    n_clamped = int(clamped_low.sum() + clamped_high.sum())
    return out, n_clamped


# ---------------------------------------------------------------------------
# Regressor -> classification conversion
# ---------------------------------------------------------------------------


def extract_mic_predictions(merged, seq_col, benchmark_unit):
    """Read MIC predictions from merged dataframe and convert to benchmark unit.

    Returns
    -------
    mic_values : np.ndarray   MIC in benchmark_unit (NaN where missing)
    n_converted : int         Number of unit conversions performed
    """
    sequences = merged[seq_col].values

    mic_pred_col = "MIC_pred" if "MIC_pred" in merged.columns else "MIC"
    if mic_pred_col not in merged.columns:
        return None, 0

    mic_values = pd.to_numeric(merged[mic_pred_col], errors="coerce").values

    pred_unit_col = (
        "MIC_unit_pred" if "MIC_unit_pred" in merged.columns
        else "MIC_unit" if "MIC_unit" in merged.columns
        else None
    )
    if pred_unit_col:
        mic_values, n_converted = harmonize_mic_units(
            mic_values, merged[pred_unit_col].values, sequences, benchmark_unit
        )
    else:
        n_converted = 0

    return mic_values, n_converted


def regressor_to_classification(
    mic_values, sequences, benchmark_unit,
    active_threshold, inactive_threshold, threshold_unit,
):
    """Convert MIC predictions to binary classification with grey zone exclusion.

    Thresholds are defined in threshold_unit. If this differs from
    benchmark_unit, predictions are converted to threshold_unit for
    comparison (per-peptide MW is used for the conversion). This matters
    because a fixed ug/ml threshold corresponds to a different uM value
    for each peptide depending on molecular weight.

    Parameters
    ----------
    mic_values : np.ndarray
        Predicted MIC in benchmark_unit.
    sequences : np.ndarray
        Amino acid sequences (needed for MW if unit conversion required).
    benchmark_unit : str
        Unit of mic_values ("uM" or "ug/ml").
    active_threshold : float
        MIC at or below this value is active. In threshold_unit.
    inactive_threshold : float
        MIC at or above this value is inactive. In threshold_unit.
    threshold_unit : str
        Unit the thresholds are expressed in ("uM" or "ug/ml").

    Returns
    -------
    y_pred : np.ndarray (int)
        Binary predictions: 1=active, 0=inactive. Grey zone gets 0.
    y_prob : np.ndarray (float)
        Probability proxy for ranking. -MIC so lower MIC = higher score.
        NaN entries get the worst possible score.
    keep_mask : np.ndarray (bool)
        True for samples outside the grey zone (to keep for evaluation).
        False for grey zone samples and NaN predictions.
    n_grey : int
        Number of samples excluded due to the grey zone.
    """
    is_valid = ~np.isnan(mic_values)

    # Convert predictions to threshold unit for comparison.
    if benchmark_unit == threshold_unit:
        mic_for_threshold = mic_values.copy()
    else:
        mw_array = compute_mw_array(sequences)
        mic_for_threshold = np.full_like(mic_values, np.nan)
        mic_for_threshold[is_valid] = convert_mic_values(
            mic_values[is_valid], mw_array[is_valid],
            benchmark_unit, threshold_unit,
        )

    # Classify: active / grey zone / inactive
    is_active = is_valid & (mic_for_threshold <= active_threshold)
    is_inactive = is_valid & (mic_for_threshold >= inactive_threshold)
    is_grey = is_valid & ~is_active & ~is_inactive

    # Samples to keep: active or inactive (not grey, not NaN)
    keep_mask = is_active | is_inactive

    y_pred = np.zeros(len(mic_values), dtype=int)
    y_pred[is_active] = 1

    # Probability proxy: -MIC (in benchmark_unit) so lower MIC = higher rank.
    # NaN and grey zone entries get the worst score.
    worst_score = -1e6
    y_prob = np.full(len(mic_values), worst_score)
    y_prob[is_valid] = -mic_values[is_valid]

    n_grey = int(is_grey.sum())

    return y_pred, y_prob, keep_mask, n_grey


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(snakemake):
    predictions_path = snakemake.input.predictions
    labels_path = snakemake.input.labels
    output_path = snakemake.output.metrics

    task_config = snakemake.params.task_config
    seq_col = snakemake.params.sequence_column
    benchmark_unit = snakemake.params.benchmark_unit

    # Activity thresholds for regressor-on-classification evaluation
    activity_thresholds = snakemake.params.activity_thresholds
    active_threshold = activity_thresholds["active"]
    inactive_threshold = activity_thresholds["inactive"]
    threshold_unit = activity_thresholds["unit"]

    # MIC clamping config (may be None to disable)
    mic_clamp = snakemake.params.mic_clamp

    task_type = task_config["type"]           # what the task expects
    variant_type = snakemake.params.variant_type  # what the model produces

    # Load data
    pred_df = pd.read_csv(predictions_path, sep="\t")
    label_df = pd.read_csv(labels_path, sep="\t")

    # Merge on sequence
    merged = pd.merge(
        pred_df, label_df, on=seq_col, how="inner", suffixes=("_pred", "_label")
    )

    report = {
        "task_type": task_type,
        "variant_type": variant_type,
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

    # -----------------------------------------------------------------
    # Classification task
    # -----------------------------------------------------------------
    if task_type == "classification":
        if not requested_metrics:
            requested_metrics = DEFAULT_CLASSIFICATION_METRICS

        # Ground truth: binary labels from the classification task
        label_col = task_config.get("label_column", "label")
        positive_label = task_config.get("positive_label", "AMP")
        y_true = (merged[label_col] == positive_label).astype(int).values

        if variant_type == "regressor":
            # Regressor on a classification task: threshold MIC predictions.
            mic_values, n_conv = extract_mic_predictions(
                merged, seq_col, benchmark_unit
            )
            if mic_values is None:
                report["error"] = (
                    "Regressor variant has no MIC column in predictions"
                )
                report["metrics"] = {m: None for m in requested_metrics}
                with open(output_path, "w") as f:
                    json.dump(report, f, indent=2)
                return

            sequences = merged[seq_col].values
            y_pred, y_prob, keep_mask, n_grey = regressor_to_classification(
                mic_values, sequences, benchmark_unit,
                active_threshold, inactive_threshold, threshold_unit,
            )

            # Exclude grey zone from both predictions and ground truth.
            y_true = y_true[keep_mask]
            y_pred = y_pred[keep_mask]
            y_prob = y_prob[keep_mask]

            n_valid_mic = int(np.sum(~np.isnan(mic_values)))
            n_kept = int(keep_mask.sum())

            report["regressor_on_classification"] = True
            report["active_threshold"] = active_threshold
            report["inactive_threshold"] = inactive_threshold
            report["threshold_unit"] = threshold_unit
            report["benchmark_unit"] = benchmark_unit
            report["n_mic_converted"] = n_conv
            report["n_valid_mic_predictions"] = n_valid_mic
            report["n_grey_zone_excluded"] = n_grey
            report["n_evaluated"] = n_kept

            print(
                f"Regressor on classification task: "
                f"active<={active_threshold} {threshold_unit}, "
                f"inactive>={inactive_threshold} {threshold_unit}, "
                f"{n_valid_mic}/{len(mic_values)} valid MIC, "
                f"{n_grey} grey zone excluded, "
                f"{n_kept} evaluated, "
                f"{int(y_pred.sum())} predicted active",
                file=sys.stderr,
            )

        else:
            # Standard classifier on classification task.
            if "Prediction" in merged.columns:
                y_pred = (merged["Prediction"] == "AMP").astype(int).values
            else:
                y_pred = np.zeros(len(merged), dtype=int)
                print(
                    "Warning: no 'Prediction' column found, using all-negative",
                    file=sys.stderr,
                )

            y_prob = None
            if "Probability_score" in merged.columns:
                y_prob = pd.to_numeric(
                    merged["Probability_score"], errors="coerce"
                ).values

        report["metrics"] = compute_classification_metrics(
            y_true, y_pred, y_prob, requested_metrics,
            task_config=task_config,
        )

    # -----------------------------------------------------------------
    # Regression task
    # -----------------------------------------------------------------
    elif task_type == "regression":
        if not requested_metrics:
            requested_metrics = DEFAULT_REGRESSION_METRICS

        sequences = merged[seq_col].values

        if variant_type == "classifier":
            # Classifier on a regression task: no meaningful regression
            # metrics can be computed from probability scores.
            report["error"] = (
                "Classifier variant cannot produce regression metrics. "
                "Skipped."
            )
            report["metrics"] = {m: None for m in requested_metrics}
            report["metrics"]["n_samples"] = len(merged)
            report["metrics"]["errors"] = 0
            print(
                "Warning: classifier variant on regression task, "
                "no regression metrics computed",
                file=sys.stderr,
            )
            with open(output_path, "w") as f:
                json.dump(report, f, indent=2)
            return

        # Regressor on regression task (standard path).

        # --- Predicted MIC ---
        mic_values, n_pred_conv = extract_mic_predictions(
            merged, seq_col, benchmark_unit
        )
        if mic_values is None:
            report["error"] = "No MIC column found in predictions"
            report["metrics"] = {}
            with open(output_path, "w") as f:
                json.dump(report, f, indent=2)
            return
        y_pred = mic_values

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
            report["error"] = "No MIC column found in merged data"
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
        sequences_valid = sequences[valid]

        # Clamp both predictions and ground truth to assay range.
        # Ground truth is clamped too because experimental MIC values
        # beyond the assay range are themselves unreliable.
        n_clamped_pred = 0
        n_clamped_true = 0
        if mic_clamp is not None:
            y_pred, n_clamped_pred = clamp_mic_values(
                y_pred, sequences_valid, benchmark_unit, mic_clamp
            )
            y_true, n_clamped_true = clamp_mic_values(
                y_true, sequences_valid, benchmark_unit, mic_clamp
            )

        report["benchmark_unit"] = benchmark_unit
        report["n_pred_converted"] = n_pred_conv
        report["n_label_converted"] = n_true_conv
        if mic_clamp is not None:
            report["mic_clamp"] = mic_clamp
            report["n_clamped_pred"] = n_clamped_pred
            report["n_clamped_true"] = n_clamped_true

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