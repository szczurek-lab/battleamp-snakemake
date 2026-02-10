"""
Validate model predictions against reference outputs.

Checks:
  - All reference sequences are present in predictions
  - Predictions (AMP/non-AMP) match exactly
  - Probability scores match within tolerance
  - For regressors: MIC values match within tolerance

Produces a JSON report with per-sequence comparison and overall pass/fail.

Usage: called by Snakemake via script: directive
"""

import json
import sys

import numpy as np
import pandas as pd


def validate_classifier(pred_df, ref_df, seq_col, tolerance):
    """Compare classifier outputs against reference."""
    results = {
        "prediction_mismatches": [],
        "probability_mismatches": [],
        "missing_sequences": [],
        "extra_sequences": [],
    }

    ref_seqs = set(ref_df[seq_col])
    pred_seqs = set(pred_df[seq_col])

    # Check coverage
    missing = ref_seqs - pred_seqs
    extra = pred_seqs - ref_seqs
    results["missing_sequences"] = sorted(list(missing))[:20]  # cap at 20
    results["extra_sequences"] = sorted(list(extra))[:20]
    results["n_missing"] = len(missing)
    results["n_extra"] = len(extra)

    # Merge on sequence for comparison
    merged = pd.merge(
        ref_df, pred_df, on=seq_col, suffixes=("_ref", "_pred"), how="inner"
    )
    results["n_compared"] = len(merged)

    if len(merged) == 0:
        results["passed"] = False
        results["error"] = "No matching sequences between predictions and reference"
        return results

    # Compare predictions (exact match)
    pred_col_ref = "Prediction_ref"
    pred_col_pred = "Prediction_pred"
    if pred_col_ref in merged.columns and pred_col_pred in merged.columns:
        mismatches = merged[merged[pred_col_ref] != merged[pred_col_pred]]
        for _, row in mismatches.iterrows():
            results["prediction_mismatches"].append({
                "sequence": row[seq_col][:30] + "..." if len(row[seq_col]) > 30 else row[seq_col],
                "reference": row[pred_col_ref],
                "predicted": row[pred_col_pred],
            })

    # Compare probability scores (within tolerance)
    prob_col_ref = "Probability_score_ref"
    prob_col_pred = "Probability_score_pred"
    if prob_col_ref in merged.columns and prob_col_pred in merged.columns:
        ref_prob = pd.to_numeric(merged[prob_col_ref], errors="coerce").values
        pred_prob = pd.to_numeric(merged[prob_col_pred], errors="coerce").values
        diffs = np.abs(ref_prob - pred_prob)

        over_tolerance = np.where(diffs > tolerance)[0]
        for idx in over_tolerance[:20]:  # cap at 20
            results["probability_mismatches"].append({
                "sequence": merged.iloc[idx][seq_col][:30] + "...",
                "reference": float(ref_prob[idx]),
                "predicted": float(pred_prob[idx]),
                "difference": float(diffs[idx]),
            })

        results["max_probability_diff"] = float(np.nanmax(diffs)) if len(diffs) > 0 else 0.0
        results["mean_probability_diff"] = float(np.nanmean(diffs)) if len(diffs) > 0 else 0.0
        results["n_probability_over_tolerance"] = int(len(over_tolerance))

    results["n_prediction_mismatches"] = len(results["prediction_mismatches"])
    results["passed"] = (
        results["n_prediction_mismatches"] == 0
        and results.get("n_probability_over_tolerance", 0) == 0
        and results["n_missing"] == 0
    )

    return results


def validate_regressor(pred_df, ref_df, seq_col, tolerance):
    """Compare regressor outputs against reference."""
    results = {
        "mic_mismatches": [],
        "missing_sequences": [],
        "extra_sequences": [],
    }

    ref_seqs = set(ref_df[seq_col])
    pred_seqs = set(pred_df[seq_col])

    missing = ref_seqs - pred_seqs
    extra = pred_seqs - ref_seqs
    results["missing_sequences"] = sorted(list(missing))[:20]
    results["extra_sequences"] = sorted(list(extra))[:20]
    results["n_missing"] = len(missing)
    results["n_extra"] = len(extra)

    merged = pd.merge(
        ref_df, pred_df, on=seq_col, suffixes=("_ref", "_pred"), how="inner"
    )
    results["n_compared"] = len(merged)

    if len(merged) == 0:
        results["passed"] = False
        results["error"] = "No matching sequences between predictions and reference"
        return results

    # Compare MIC values (within tolerance, using relative difference)
    mic_col_ref = "MIC_ref"
    mic_col_pred = "MIC_pred"
    if mic_col_ref in merged.columns and mic_col_pred in merged.columns:
        ref_mic = pd.to_numeric(merged[mic_col_ref], errors="coerce").values
        pred_mic = pd.to_numeric(merged[mic_col_pred], errors="coerce").values

        # Relative difference: |ref - pred| / max(|ref|, epsilon)
        epsilon = 1e-10
        rel_diffs = np.abs(ref_mic - pred_mic) / np.maximum(np.abs(ref_mic), epsilon)

        over_tolerance = np.where(rel_diffs > tolerance)[0]
        for idx in over_tolerance[:20]:
            results["mic_mismatches"].append({
                "sequence": merged.iloc[idx][seq_col][:30] + "...",
                "reference": float(ref_mic[idx]),
                "predicted": float(pred_mic[idx]),
                "relative_difference": float(rel_diffs[idx]),
            })

        results["max_mic_relative_diff"] = float(np.nanmax(rel_diffs)) if len(rel_diffs) > 0 else 0.0
        results["mean_mic_relative_diff"] = float(np.nanmean(rel_diffs)) if len(rel_diffs) > 0 else 0.0
        results["n_mic_over_tolerance"] = int(len(over_tolerance))

    results["n_mic_mismatches"] = len(results["mic_mismatches"])
    results["passed"] = (
        results["n_mic_mismatches"] == 0
        and results["n_missing"] == 0
    )

    return results


def main(snakemake):
    pred_path = snakemake.input.predictions
    ref_path = snakemake.input.reference
    report_path = snakemake.output.report

    model_type = snakemake.params.model_type
    tolerance = float(snakemake.params.tolerance)
    seq_col = snakemake.params.sequence_column

    pred_df = pd.read_csv(pred_path, sep="\t")
    ref_df = pd.read_csv(ref_path, sep="\t")

    report = {
        "model_type": model_type,
        "tolerance": tolerance,
        "n_reference": len(ref_df),
        "n_predictions": len(pred_df),
    }

    if model_type == "classifier":
        result = validate_classifier(pred_df, ref_df, seq_col, tolerance)
    elif model_type == "regressor":
        result = validate_regressor(pred_df, ref_df, seq_col, tolerance)
    else:
        result = {"passed": False, "error": f"Unknown model type: {model_type}"}

    report.update(result)

    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    status = "PASSED" if report["passed"] else "FAILED"
    print(f"Validation {status} for {model_type}", file=sys.stderr)

    if not report["passed"]:
        if report.get("n_prediction_mismatches", 0) > 0:
            print(
                f"  {report['n_prediction_mismatches']} prediction mismatches",
                file=sys.stderr,
            )
        if report.get("n_probability_over_tolerance", 0) > 0:
            print(
                f"  {report['n_probability_over_tolerance']} probability scores "
                f"exceed tolerance {tolerance} "
                f"(max diff: {report.get('max_probability_diff', '?')})",
                file=sys.stderr,
            )
        if report.get("n_mic_over_tolerance", 0) > 0:
            print(
                f"  {report['n_mic_over_tolerance']} MIC values "
                f"exceed tolerance {tolerance}",
                file=sys.stderr,
            )
        if report.get("n_missing", 0) > 0:
            print(
                f"  {report['n_missing']} reference sequences missing from predictions",
                file=sys.stderr,
            )


main(snakemake)
