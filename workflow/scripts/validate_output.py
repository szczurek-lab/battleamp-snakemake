"""
Validate model inference output against the expected schema.

Checks:
  - File exists and is non-empty
  - Required columns present (case-sensitive)
  - No null/empty values in required columns
  - Probability_score in [0, 1] for classifiers
  - Prediction values are valid for classifiers
  - MIC values are numeric and positive for regressors

Produces a JSON validation report.

Usage: called by Snakemake via script: directive
"""

import json
import sys
import pandas as pd


def validate_classifier(df, seq_col):
    """Validate classifier output."""
    errors = []
    warnings = []

    required_cols = [seq_col, "Prediction", "Probability_score"]
    for col in required_cols:
        if col not in df.columns:
            errors.append(f"Missing required column: '{col}'")

    if errors:
        return {"valid": False, "errors": errors, "warnings": warnings}

    # Check for nulls
    for col in required_cols:
        null_count = df[col].isna().sum()
        if null_count > 0:
            errors.append(f"Column '{col}' has {null_count} null values")

    # Check probability range
    prob = pd.to_numeric(df["Probability_score"], errors="coerce")
    out_of_range = ((prob < 0) | (prob > 1)).sum()
    if out_of_range > 0:
        errors.append(
            f"Probability_score has {out_of_range} values outside [0, 1]"
        )

    coercion_failures = prob.isna().sum() - df["Probability_score"].isna().sum()
    if coercion_failures > 0:
        errors.append(
            f"Probability_score has {coercion_failures} non-numeric values"
        )

    # Check prediction values
    valid_predictions = {"AMP", "non-AMP"}
    actual_predictions = set(df["Prediction"].dropna().unique())
    unexpected = actual_predictions - valid_predictions
    if unexpected:
        warnings.append(
            f"Unexpected Prediction values: {unexpected}. "
            f"Expected: {valid_predictions}"
        )

    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
    }


def validate_regressor(df, seq_col):
    """Validate regressor output."""
    errors = []
    warnings = []

    required_cols = [seq_col, "MIC", "MIC_unit"]
    for col in required_cols:
        if col not in df.columns:
            errors.append(f"Missing required column: '{col}'")

    if errors:
        return {"valid": False, "errors": errors, "warnings": warnings}

    # Check for nulls
    for col in required_cols:
        null_count = df[col].isna().sum()
        if null_count > 0:
            errors.append(f"Column '{col}' has {null_count} null values")

    # Check MIC is numeric and positive
    mic = pd.to_numeric(df["MIC"], errors="coerce")
    negative_count = (mic < 0).sum()
    if negative_count > 0:
        warnings.append(f"MIC has {negative_count} negative values")

    coercion_failures = mic.isna().sum() - df["MIC"].isna().sum()
    if coercion_failures > 0:
        errors.append(f"MIC has {coercion_failures} non-numeric values")

    # Check MIC_unit values
    valid_units = {"ug/ml", "uM"}
    actual_units = set(df["MIC_unit"].dropna().unique())
    unexpected = actual_units - valid_units
    if unexpected:
        warnings.append(
            f"Unexpected MIC_unit values: {unexpected}. "
            f"Expected: {valid_units}"
        )

    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
    }


def main(snakemake):
    tsv_path = snakemake.input.tsv
    report_path = snakemake.output.report
    model_type = snakemake.params.model_type
    seq_col = snakemake.params.sequence_column

    report = {
        "file": str(tsv_path),
        "model_type": model_type,
    }

    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        report["valid"] = False
        report["errors"] = [f"Failed to read TSV: {e}"]
        report["warnings"] = []
        report["row_count"] = 0
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
        print(f"VALIDATION FAILED: {e}", file=sys.stderr)
        return

    report["row_count"] = len(df)
    report["columns"] = list(df.columns)

    if model_type == "classifier":
        result = validate_classifier(df, seq_col)
    elif model_type == "regressor":
        result = validate_regressor(df, seq_col)
    else:
        result = {
            "valid": False,
            "errors": [f"Unknown model type: {model_type}"],
            "warnings": [],
        }

    report.update(result)

    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    status = "PASSED" if report["valid"] else "FAILED"
    print(
        f"Validation {status}: {len(df)} rows, "
        f"{len(report['errors'])} errors, {len(report['warnings'])} warnings",
        file=sys.stderr,
    )

    if not report["valid"]:
        for err in report["errors"]:
            print(f"  ERROR: {err}", file=sys.stderr)
        # Don't fail the pipeline -- let evaluation handle missing/bad data
        # But the report is available for inspection.


main(snakemake)
