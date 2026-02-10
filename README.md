# battleamp-snakemake

A general-purpose Snakemake pipeline for benchmarking antimicrobial peptide (AMP)
prediction models. Model-agnostic, dataset-agnostic. Bring your own models and data,
get standardized evaluation.

## Prerequisites

- Python >= 3.10
- conda (Miniconda or Anaconda) for model environment management
- Git

## Installation

### Option A: pip (recommended, works on clusters without admin access)

```bash
# Create a virtual environment in your home directory
python3 -m venv ~/.venvs/battleamp-snakemake
source ~/.venvs/battleamp-snakemake/bin/activate

# Install Snakemake and pipeline dependencies
pip install snakemake pulp pandas scikit-learn scipy numpy pyyaml

# Verify
snakemake --version
```

On HPC clusters, you may need to load a Python module first:
```bash
module avail python          # see what's available
module load python/3.10      # load it (version may vary)
```

### Option B: conda

```bash
conda create -n battleamp-snakemake -c conda-forge -c bioconda \
    snakemake pandas scikit-learn scipy numpy pyyaml
conda activate battleamp-snakemake
```

### Clone the repository

```bash
git clone https://github.com/szczurek-lab/battleamp-snakemake.git
cd battleamp-snakemake
```

## Quick Start

```bash
# Activate the environment (pick whichever you installed above)
source ~/.venvs/battleamp-snakemake/bin/activate   # if pip
# conda activate battleamp-snakemake               # if conda

# Run the pipeline on the included example (no --use-conda needed for example)
snakemake --snakefile workflow/Snakefile --cores 4

# Check results
cat results/aggregated/summary.tsv
```

When running with real models, add `--use-conda` so Snakemake creates each
model's isolated conda environment from its `environment.yaml`:

```bash
snakemake --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1
```

On first run, Snakemake will automatically create a separate conda environment
for each model. These are cached in `.snakemake/conda/` and reused on
subsequent runs. You do not need to create or manage any model environments
yourself.

## How It Works

The pipeline runs five stages:

1. **Setup** -- downloads model weights and compiles extensions (inside each
   model's auto-created conda env)
2. **Pre-filter** -- generates model-specific FASTA files filtered by each model's
   declared sequence length constraints
3. **Inference** -- runs each model variant on each dataset
4. **Evaluate** -- compares predictions to ground truth and computes metrics
5. **Aggregate** -- combines all metrics into summary tables

All stages are orchestrated by Snakemake, which handles dependency resolution,
parallelism, and caching (won't re-run steps whose outputs are already up to date).

## Running on an HPC Cluster

### SLURM

```bash
# Basic: run on a single node
snakemake --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1

# With a SLURM profile (submits each rule as a separate job):
snakemake --snakefile workflow/Snakefile --use-conda --profile slurm/
```

### GPU Scheduling

Models declare `gpu_required: true/false` in their `model.yaml`. Control GPU
scheduling with the `--resources` flag:

```bash
# Single GPU workstation: GPU models queue, CPU models parallelize
snakemake --use-conda --cores 8 --resources gpu=1

# Two GPUs available:
snakemake --use-conda --cores 16 --resources gpu=2

# No GPU (CPU models only; GPU models will fail):
snakemake --use-conda --cores 8 --resources gpu=0
```

## Adding Your Own Model

See [docs/adding_a_model.md](docs/adding_a_model.md) for full instructions.

In short, create a directory under `models/` with:

```
models/your-model/
    model.yaml         # Metadata: name, type, length limits, GPU requirement
    environment.yaml   # Conda env spec (Snakemake creates the env from this)
    inference.sh       # Prediction script: $1=input.fasta, $2=output.tsv
    setup.sh           # (Optional) Download weights, compile extensions
```

Then add `your-model` to the `models:` list in `config/config.yaml`.

## Adding Your Own Dataset

See [docs/adding_a_dataset.md](docs/adding_a_dataset.md) for full instructions.

Create a directory under `datasets/` with:

```
datasets/your-dataset/
    sequences.fasta    # Input sequences
    labels.tsv         # Ground truth: sequence + label/MIC columns
    dataset.yaml       # (Optional) Metadata
```

Then add a task definition in `config/config.yaml` that references it.

## Model Validation

The pipeline includes a validation stage that verifies code modifications did
not alter model predictions. Each model can store reference inputs and expected
outputs; the pipeline re-runs the model and compares.

```bash
# Run validation for all models that have reference data
snakemake --snakefile workflow/Snakefile validate --use-conda --cores 4

# Check results
cat results/validation/validation_summary.tsv
```

See [docs/validation.md](docs/validation.md) for details on generating reference
data and interpreting results.

## Output Structure

```
results/
    setup/{model}/.setup_done
    prefiltered/{variant}/{dataset}/sequences.fasta
    inference/{variant}/{dataset}/predictions.tsv
    inference/{variant}/{dataset}/validation_report.json
    evaluation/{variant}/{task}/metrics.json
    aggregated/
        summary.tsv
        classification_results.tsv
        regression_results.tsv
    validation/
        {model}/predictions.tsv
        {model}/validation_report.json
        validation_summary.tsv
    logs/
        setup/{model}.log
        prefilter/{variant}/{dataset}.log
        inference/{variant}/{dataset}.log
        evaluation/{variant}/{task}.log
        validation_inference/{model}.log
        validation_compare/{model}.log
```

## Model Interface Contract

### Output TSV format

**Classifiers** must produce:
```
sequence    Prediction    Probability_score
KFLQ...     AMP           0.87
GIKL...     non-AMP       0.23
```

**Regressors** must produce:
```
sequence    MIC    MIC_unit
KFLQ...     4.2    ug/ml
GIKL...     16.0   uM
```

Column names are case-sensitive. Sequences the model cannot handle must be skipped
(not crash). Exit code 0 on success.

## Supported Metrics

### Classification
accuracy, mcc, f1, fpr, tpr, tnr, lr_plus, auroc, auprc

### Regression
msle_ln, spearman, r2, rmse, mae

## Troubleshooting

### conda create fails with HTTP 000 CONNECTION FAILED
Your conda installation may have SSL or channel configuration issues. Try:
```bash
conda config --set ssl_verify true
conda config --show channels
```
If conda channels are unreachable but pip works, use the pip installation method
above for Snakemake. The `--use-conda` flag will still use conda for model
environments, so conda itself must be functional. Check with your cluster
administrator if conda.anaconda.org is accessible.

### python3 not found
On HPC clusters, Python is often provided as a module:
```bash
module avail python
module load python/3.10
```

### Permission denied errors
The pipeline does not require admin/sudo access. All files are created in the
repository directory and `~/.venvs/` (if using pip). If you encounter permission
issues, make sure you are working within your home directory.

## License

MIT
