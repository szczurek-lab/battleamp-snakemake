"""
Rule: evaluation
Computes metrics for each variant x task combination by comparing
model predictions against ground truth labels.

The label file is specified per-task in config.yaml (not per-dataset),
so multiple tasks can evaluate the same inference output differently.

When a regressor variant is evaluated on a classification task, its MIC
predictions are thresholded using activity_thresholds from config:
  - active:   MIC <= active threshold
  - inactive: MIC >= inactive threshold
  - grey zone (between thresholds) is excluded from evaluation
Thresholds are specified with an explicit unit; per-peptide molecular
weight is used if conversion is needed.

For regression tasks, predictions and ground truth are clamped to an
assay-realistic range (mic_clamp) before computing metrics.  This
prevents extreme outlier predictions from dominating error metrics.
"""


rule evaluate_task:
    """Compute metrics for one variant on one task."""
    input:
        predictions = lambda wc: (
            f"{OUTPUT_DIR}/inference/{wc.variant}/"
            f"{TASKS[wc.task]['dataset']}/predictions.tsv"
        ),
        validation = lambda wc: (
            f"{OUTPUT_DIR}/inference/{wc.variant}/"
            f"{TASKS[wc.task]['dataset']}/validation_report.json"
        ),
        labels = lambda wc: get_task_labels(wc.task),
    output:
        metrics = f"{OUTPUT_DIR}/evaluation/{{variant}}/{{task}}/metrics.json",
    params:
        task_config = lambda wc: TASKS[wc.task],
        variant_type = lambda wc: VARIANTS[wc.variant]["type"],
        sequence_column = SEQ_COL,
        benchmark_unit = config.get("benchmark_unit", "ug/ml"),
        activity_thresholds = config.get("activity_thresholds", {
            "active": 32,
            "inactive": 128,
            "unit": "ug/ml",
        }),
        mic_clamp = config.get("mic_clamp", None),
    log:
        f"{OUTPUT_DIR}/logs/evaluation/{{variant}}/{{task}}.log",
    script:
        "../scripts/evaluate.py"