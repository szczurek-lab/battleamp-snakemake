"""
Rule: model setup
Runs setup.sh inside the model's conda environment to download weights,
compile extensions, etc. The conda env itself is created automatically
by Snakemake from the model's environment.yaml (via --use-conda).
setup.sh should NOT create or activate any conda env.
"""


rule model_setup:
    """Run setup.sh for a model (once per model, not per variant)."""
    output:
        done = touch(f"{OUTPUT_DIR}/setup/{{model}}/.setup_done"),
    params:
        model_dir = lambda wc: str(MODELS_DIR / wc.model),
        setup_script = lambda wc: str(MODELS_DIR / wc.model / "setup.sh"),
    log:
        f"{OUTPUT_DIR}/logs/setup/{{model}}.log",
    conda:
        lambda wc: MODEL_ENV_PATHS[wc.model]
    shell:
        """
        LOG_ABS=$(realpath -m {log})
        mkdir -p $(dirname "$LOG_ABS")
        if [ -f {params.setup_script} ]; then
            echo "Running setup.sh for {wildcards.model}..." >> "$LOG_ABS"
            cd {params.model_dir} && bash setup.sh >> "$LOG_ABS" 2>&1
        else
            echo "No setup.sh found for {wildcards.model}, skipping." >> "$LOG_ABS"
        fi
        """
