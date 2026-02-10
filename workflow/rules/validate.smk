"""
Rule: model validation
Runs each model on its reference inputs and compares predictions against
stored reference outputs. This verifies that code modifications (hotfixes,
adapter scripts, environment changes) did not alter model behavior.

Each model that has a validation/ directory with reference_input.fasta and
reference_output.tsv will be validated. Models without validation data are
skipped with a warning.
"""


def get_models_with_validation():
    """Return list of model names that have validation data."""
    validated = []
    for model_name in MODELS:
        ref_input = MODELS_DIR / model_name / "validation" / "reference_input.fasta"
        ref_output = MODELS_DIR / model_name / "validation" / "reference_output.tsv"
        if ref_input.exists() and ref_output.exists():
            validated.append(model_name)
    return validated


VALIDATED_MODELS = get_models_with_validation()


rule run_validation_inference:
    """Run model on its reference input sequences."""
    input:
        setup_done = f"{OUTPUT_DIR}/setup/{{model}}/.setup_done",
        fasta = lambda wc: str(
            MODELS_DIR / wc.model / "validation" / "reference_input.fasta"
        ),
    output:
        tsv = f"{OUTPUT_DIR}/validation/{{model}}/predictions.tsv",
    params:
        model_dir = lambda wc: str(MODELS_DIR / wc.model),
    log:
        f"{OUTPUT_DIR}/logs/validation_inference/{{model}}.log",
    resources:
        gpu = lambda wc: 1 if MODELS[wc.model].get("gpu_required", False) else 0,
    conda:
        lambda wc: MODEL_ENV_PATHS[wc.model]
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
            >> "$LOG_ABS" 2>&1
        """


rule compare_validation_output:
    """Compare model predictions against reference outputs."""
    input:
        predictions = f"{OUTPUT_DIR}/validation/{{model}}/predictions.tsv",
        reference = lambda wc: str(
            MODELS_DIR / wc.model / "validation" / "reference_output.tsv"
        ),
    output:
        report = f"{OUTPUT_DIR}/validation/{{model}}/validation_report.json",
    params:
        model_type = lambda wc: MODELS[wc.model]["type"],
        tolerance = lambda wc: MODELS[wc.model].get("validation_tolerance", 1e-4),
        sequence_column = SEQ_COL,
    log:
        f"{OUTPUT_DIR}/logs/validation_compare/{{model}}.log",
    conda:
        PIPELINE_ENV
    script:
        "../scripts/validate_model.py"


rule validation_summary:
    """Aggregate all validation reports into a single summary."""
    input:
        reports = expand(
            f"{OUTPUT_DIR}/validation/{{model}}/validation_report.json",
            model=VALIDATED_MODELS,
        ),
    output:
        summary = f"{OUTPUT_DIR}/validation/validation_summary.tsv",
    conda:
        PIPELINE_ENV
    script:
        "../scripts/validation_summary.py"
