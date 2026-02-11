"""
Aggregate all per-variant per-task metric files into summary tables.

Produces:
  - summary.tsv: all metrics in one table (variant x task x metric)
  - classification_results.tsv: classification tasks only
  - regression_results.tsv: regression tasks only

Usage: called by Snakemake via script: directive
"""

import json
import sys
from pathlib import Path

import pandas as pd


def main(snakemake):
    metric_files = snakemake.input.metrics
    summary_path = snakemake.output.summary
    classification_path = snakemake.output.classification
    regression_path = snakemake.output.regression

    rows = []

    for metric_file in metric_files:
        # Extract variant and task from path
        # Path pattern: .../evaluation/{variant}/{task}/metrics.json
        parts = Path(metric_file).parts
        task = parts[-2]
        variant = parts[-3]

        with open(metric_file) as f:
            data = json.load(f)

        task_type = data.get("task_type", "unknown")
        metrics = data.get("metrics", {})
        n_matched = data.get("n_matched", 0)

        row = {
            "variant": variant,
            "task": task,
            "task_type": task_type,
            "n_matched": n_matched,
        }
        row.update(metrics)
        rows.append(row)

    if not rows:
        # Write empty files
        pd.DataFrame().to_csv(summary_path, sep="\t", index=False)
        pd.DataFrame().to_csv(classification_path, sep="\t", index=False)
        pd.DataFrame().to_csv(regression_path, sep="\t", index=False)
        print("Warning: no metric files found", file=sys.stderr)
        return

    df = pd.DataFrame(rows)

    # Sort by task, then variant
    df = df.sort_values(["task", "variant"]).reset_index(drop=True)

    # Add coverage: fraction of task sequences the model could evaluate.
    # task_size = max n_matched across all models for that task (at least
    # one model with no length constraints covers every sequence).
    task_size = df.groupby("task")["n_matched"].transform("max")
    col_pos = df.columns.get_loc("n_matched") + 1
    df.insert(col_pos, "task_size", task_size.astype(int))
    df.insert(col_pos + 1, "coverage", df["n_matched"] / df["task_size"])

    # Write full summary
    df.to_csv(summary_path, sep="\t", index=False)

    # Split by task type
    clf_df = df[df["task_type"] == "classification"]
    reg_df = df[df["task_type"] == "regression"]

    clf_df.to_csv(classification_path, sep="\t", index=False)
    reg_df.to_csv(regression_path, sep="\t", index=False)

    print(
        f"Aggregated {len(df)} results: "
        f"{len(clf_df)} classification, {len(reg_df)} regression",
        file=sys.stderr,
    )


main(snakemake)
