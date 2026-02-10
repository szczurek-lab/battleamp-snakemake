"""
Rule: model inference
Runs inference.sh inside the model's conda environment on a pre-filtered dataset.
The conda env is managed by Snakemake (--use-conda) from environment.yaml.
inference.sh should NOT activate any conda env -- it is already active.

Validation runs in the pipeline's own environment (no conda: directive).
"""


rule run_inference:
    """Run a model variant on a dataset."""
    input:
        setup_done = lambda wc: (
            f"{OUTPUT_DIR}/setup/{VARIANTS[wc.variant]['model']}/.setup_done"
        ),
        fasta = f"{OUTPUT_DIR}/prefiltered/{{variant}}/{{dataset}}/sequences.fasta",
    output:
        tsv = f"{OUTPUT_DIR}/inference/{{variant}}/{{dataset}}/predictions.tsv",
    params:
        model_dir = lambda wc: str(MODELS_DIR / VARIANTS[wc.variant]["model"]),
        extra_args = lambda wc: " ".join(
            str(a) for a in VARIANTS[wc.variant]["args"]
        ),
    log:
        f"{OUTPUT_DIR}/logs/inference/{{variant}}/{{dataset}}.log",
    resources:
        gpu = lambda wc: 1 if VARIANTS[wc.variant]["gpu"] else 0,
    conda:
        lambda wc: VARIANTS[wc.variant]["env"]
    shell:
        """
        INPUT_ABS=$(realpath {input.fasta})
        OUTPUT_ABS=$(realpath -m {output.tsv})
        LOG_ABS=$(realpath -m {log})
        mkdir -p $(dirname "$OUTPUT_ABS")
        mkdir -p $(dirname "$LOG_ABS")
        cd {params.model_dir} && \
        bash inference.sh \
            "$INPUT_ABS" \
            "$OUTPUT_ABS" \
            {params.extra_args} \
            >> "$LOG_ABS" 2>&1
        """


rule validate_output:
    """Validate that inference output conforms to the expected schema."""
    input:
        tsv = f"{OUTPUT_DIR}/inference/{{variant}}/{{dataset}}/predictions.tsv",
    output:
        report = f"{OUTPUT_DIR}/inference/{{variant}}/{{dataset}}/validation_report.json",
    params:
        model_type = lambda wc: VARIANTS[wc.variant]["type"],
        sequence_column = SEQ_COL,
    log:
        f"{OUTPUT_DIR}/logs/validation/{{variant}}/{{dataset}}.log",
    script:
        "../scripts/validate_output.py"
