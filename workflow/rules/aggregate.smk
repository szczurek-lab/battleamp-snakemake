"""
Rule: aggregation
Combines all per-variant per-task metric files into summary tables,
and aggregates Snakemake benchmark logs into resource usage tables.
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


rule aggregate_resources:
    """Aggregate Snakemake benchmark logs into resource usage tables."""
    input:
        metrics = ALL_METRICS,  # ensures inference is complete
    output:
        usage = f"{OUTPUT_DIR}/aggregated/resource_usage.tsv",
        summary = f"{OUTPUT_DIR}/aggregated/resource_summary.tsv",
    params:
        benchmark_dir = f"{OUTPUT_DIR}/benchmarks",
    log:
        f"{OUTPUT_DIR}/logs/aggregate_resources.log",
    script:
        "../scripts/aggregate_resources.py"


rule extract_md_predictions:
    """Extract predictions for MD validation datasets from battleamp-all inference.

    Pulls control (safe) and activity cliff predictions from existing
    battleamp-all inference output and formats them for Figure 5.

    Core data files live under core-data/data/md/ in the repo root.
    """
    input:
        metrics = ALL_METRICS,  # ensures inference is complete
        control_labels = str(REPO_ROOT / "md" / "control-md.tsv"),
        cliff_labels = str(REPO_ROOT  / "md" / "activity-cliffs.tsv"),
    output:
        control = f"{OUTPUT_DIR}/aggregated/control-md_predictions.tsv",
        cliffs = f"{OUTPUT_DIR}/aggregated/activity-cliffs_predictions.tsv",
    params:
        inference_dir = f"{OUTPUT_DIR}/inference",
        dataset = "battleamp-all",
    log:
        f"{OUTPUT_DIR}/logs/extract_md_predictions.log",
    script:
        "../scripts/extract_md_predictions.py"