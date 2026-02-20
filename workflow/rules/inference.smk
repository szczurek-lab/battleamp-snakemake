"""
Rule: model inference
Runs inference.sh inside the model's conda environment on a pre-filtered dataset.
The conda env is managed by Snakemake (--use-conda) from environment.yaml.
inference.sh should NOT activate any conda env -- it is already active.

Multi-output models
-------------------
Some models (e.g. APEX) produce predictions for multiple variants in a single
inference pass. These declare ``multioutput: true`` in model.yaml. For them:

1. ``run_multioutput_inference`` runs once per model and writes all variant
   prediction files into a shared directory.
2. ``resolve_multioutput_predictions`` copies each variant's specific file
   to the standard per-variant path so downstream rules work unchanged.

Wildcard constraints keep the two inference paths mutually exclusive.
"""


# ---------------------------------------------------------------------------
# Standard inference (non-multioutput variants)
# ---------------------------------------------------------------------------

rule run_inference:
    """Run a model variant on a dataset."""
    wildcard_constraints:
        variant=STANDARD_VARIANT_RE,
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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/inference/{{variant}}/{{dataset}}.tsv"
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
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        export CUDA_VISIBLE_DEVICES=${{CUDA_VISIBLE_DEVICES:-0}}
        cd {params.model_dir} && \
        bash inference.sh \
            "$INPUT_ABS" \
            "$OUTPUT_ABS" \
            {params.extra_args} \
            >> "$LOG_ABS" 2>&1
        """


# ---------------------------------------------------------------------------
# Multi-output inference (runs once per model, produces all variant files)
# ---------------------------------------------------------------------------

rule run_multioutput_inference:
    """Run a multi-output model once; produces prediction files for all its variants."""
    wildcard_constraints:
        model=MULTIOUTPUT_MODEL_RE,
    input:
        setup_done = f"{OUTPUT_DIR}/setup/{{model}}/.setup_done",
        # Reuse the primary variant's prefiltered FASTA (all variants of a
        # multioutput model share the same length constraints).
        fasta = lambda wc: (
            f"{OUTPUT_DIR}/prefiltered/"
            f"{MULTIOUTPUT_MODELS[wc.model]['primary']}/{wc.dataset}/sequences.fasta"
        ),
    output:
        done = touch(
            f"{OUTPUT_DIR}/inference/{{model}}/{{dataset}}/.multioutput_done"
        ),
        tsv = f"{OUTPUT_DIR}/inference/{{model}}/{{dataset}}/predictions.tsv",
    params:
        model_dir = lambda wc: str(MODELS_DIR / wc.model),
    log:
        f"{OUTPUT_DIR}/logs/inference/{{model}}/{{dataset}}.log",
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/inference/{{model}}/{{dataset}}.tsv"
    resources:
        gpu = lambda wc: 1 if MODELS[wc.model].get("gpu_required", False) else 0,
    conda:
        lambda wc: MODEL_ENV_PATHS[wc.model]
    shell:
        """
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        export CUDA_VISIBLE_DEVICES=${{CUDA_VISIBLE_DEVICES:-0}}
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


rule resolve_multioutput_predictions:
    """Copy a multioutput variant's predictions to the standard per-variant path.

    After the model-level inference finishes, each variant has a specific file
    (e.g. predictions-min.tsv) in the model's output directory.  This rule
    copies it to results/inference/{variant}/{dataset}/predictions.tsv so that
    downstream rules (validate, evaluate) work unchanged.
    """
    wildcard_constraints:
        variant=MULTIOUTPUT_VARIANT_RE,
    input:
        done = lambda wc: (
            f"{OUTPUT_DIR}/inference/"
            f"{VARIANTS[wc.variant]['model']}/{wc.dataset}/.multioutput_done"
        ),
    output:
        tsv = f"{OUTPUT_DIR}/inference/{{variant}}/{{dataset}}/predictions.tsv",
    params:
        src = lambda wc: (
            f"{OUTPUT_DIR}/inference/"
            f"{VARIANTS[wc.variant]['model']}/{wc.dataset}/"
            f"{VARIANTS[wc.variant]['predictions_file']}"
        ),
    run:
        import shutil
        shutil.copy2(params.src, output.tsv)


# ---------------------------------------------------------------------------
# Output validation (shared by both paths)
# ---------------------------------------------------------------------------

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