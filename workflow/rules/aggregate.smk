"""
Rule: aggregation
Combines all per-variant per-task metric files into summary tables.
"""


rule aggregate_results:
    """Combine all metrics into summary tables."""
    input:
        metrics = ALL_METRICS,
    output:
        summary = f"{OUTPUT_DIR}/aggregated/summary.tsv",
        classification = f"{OUTPUT_DIR}/aggregated/classification_results.tsv",
        regression = f"{OUTPUT_DIR}/aggregated/regression_results.tsv",
    log:
        f"{OUTPUT_DIR}/logs/aggregate.log",
    script:
        "../scripts/aggregate_results.py"
