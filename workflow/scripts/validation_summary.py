"""
Aggregate all per-model validation reports into a summary table.

Usage: called by Snakemake via script: directive
"""

import json
import sys
from pathlib import Path

import pandas as pd


def main(snakemake):
    report_files = snakemake.input.reports
    summary_path = snakemake.output.summary

    rows = []
    for report_file in report_files:
        model = Path(report_file).parts[-2]

        with open(report_file) as f:
            data = json.load(f)

        row = {
            "model": model,
            "passed": data.get("passed", False),
            "model_type": data.get("model_type", "unknown"),
            "n_reference": data.get("n_reference", 0),
            "n_compared": data.get("n_compared", 0),
            "n_missing": data.get("n_missing", 0),
        }

        if data.get("model_type") == "classifier":
            row["n_prediction_mismatches"] = data.get("n_prediction_mismatches", 0)
            row["n_probability_over_tolerance"] = data.get(
                "n_probability_over_tolerance", 0
            )
            row["max_probability_diff"] = data.get("max_probability_diff", None)
        elif data.get("model_type") == "regressor":
            row["n_mic_mismatches"] = data.get("n_mic_mismatches", 0)
            row["max_mic_relative_diff"] = data.get("max_mic_relative_diff", None)

        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.sort_values("model").reset_index(drop=True)
    df.to_csv(summary_path, sep="\t", index=False)

    n_passed = df["passed"].sum()
    n_total = len(df)
    print(
        f"Validation summary: {n_passed}/{n_total} models passed",
        file=sys.stderr,
    )

    failed = df[~df["passed"]]
    if len(failed) > 0:
        print("FAILED models:", file=sys.stderr)
        for _, row in failed.iterrows():
            print(f"  - {row['model']}", file=sys.stderr)


main(snakemake)
