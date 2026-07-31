"""
Rule: model setup
Runs setup.sh inside the model's conda environment to download weights,
compile extensions, etc. The conda env itself is created automatically
by Snakemake from the model's environment.yaml (via --use-conda).
setup.sh should NOT create or activate any conda env.
"""


def _setup_deps(model):
    """Files whose change must invalidate the model's setup marker.

    setup.sh is optional (example-model has none), so only existing paths are
    returned -- a missing input that no rule produces is a DAG build error.
    """
    deps = [MODEL_ENV_PATHS[model]]
    setup_script = MODELS_DIR / model / "setup.sh"
    if setup_script.exists():
        deps.append(str(setup_script))
    return deps


rule model_setup:
    """Run setup.sh for a model (once per model, not per variant).

    The marker is deliberately mtime-coupled to setup.sh and environment.yaml.
    Without that, bumping a model submodule to fix its setup script leaves the
    old marker in place and the broken environment is reused forever.
    """
    input:
        lambda wc: _setup_deps(wc.model),
    output:
        done = touch(f"{OUTPUT_DIR}/setup/{{model}}/.setup_done"),
        env = f"{OUTPUT_DIR}/setup/{{model}}/env.txt",
    params:
        model_dir = lambda wc: str(MODELS_DIR / wc.model),
        setup_script = lambda wc: str(MODELS_DIR / wc.model / "setup.sh"),
    log:
        f"{OUTPUT_DIR}/logs/setup/{{model}}.log",
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/setup/{{model}}.tsv"
    conda:
        lambda wc: MODEL_ENV_PATHS[wc.model]
    shell:
        """
        LOG_ABS=$(realpath -m {log})
        ENV_ABS=$(realpath -m {output.env})
        mkdir -p $(dirname "$LOG_ABS") $(dirname "$ENV_ABS")
        # Truncate: setup now re-runs when setup.sh changes, and a log that
        # mixes the old and new attempt is worse than no log.
        : > "$LOG_ABS"
        if [ -f {params.setup_script} ]; then
            echo "Running setup.sh for {wildcards.model}..." >> "$LOG_ABS"
            cd {params.model_dir} && bash setup.sh >> "$LOG_ABS" 2>&1
        else
            echo "No setup.sh found for {wildcards.model}, skipping." >> "$LOG_ABS"
        fi
        # Recorded so the Snakefile can spot a rebuilt env and re-run setup.
        echo "${{CONDA_PREFIX:-}}" > "$ENV_ABS"
        """