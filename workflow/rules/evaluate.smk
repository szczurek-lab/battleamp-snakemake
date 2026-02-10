"""
Rule: evaluation
Computes metrics for each variant x task combination by comparing
model predictions against ground truth labels.
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
        labels = lambda wc: str(
            DATASETS_DIR / TASKS[wc.task]["dataset"] / "labels.tsv"
        ),
    output:
        metrics = f"{OUTPUT_DIR}/evaluation/{{variant}}/{{task}}/metrics.json",
    params:
        task_config = lambda wc: TASKS[wc.task],
        sequence_column = SEQ_COL,
    log:
        f"{OUTPUT_DIR}/logs/evaluation/{{variant}}/{{task}}.log",
    script:
        "../scripts/evaluate.py"
