# battleamp-snakemake

A Snakemake pipeline for benchmarking antimicrobial peptide (AMP) prediction
models against a curated set of classification and regression tasks.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [I want to score my own peptides](#i-want-to-score-my-own-peptides)
- [I want to reproduce the benchmark](#i-want-to-reproduce-the-benchmark)
- [I want to add my own model](#i-want-to-add-my-own-model)
- [Reference](#reference)
  - [Available Models](#available-models)
  - [Evaluation Metrics](#evaluation-metrics)
  - [Configuration](#configuration)
  - [Output Structure](#output-structure)
  - [Running on an HPC Cluster](#running-on-an-hpc-cluster)
  - [Model Validation](#model-validation)


## Installation

**Prerequisites:** Python >= 3.10, conda >= 23.x, Git.

Install Snakemake and dependencies into an isolated environment:

```bash
# Option A: pip (recommended on clusters without admin access)
python3 -m venv ~/.venvs/battleamp-snakemake
source ~/.venvs/battleamp-snakemake/bin/activate
pip install snakemake pulp pandas scikit-learn scipy numpy pyyaml biopython

# Option B: conda
conda create -n battleamp-snakemake -c conda-forge -c bioconda \
    snakemake pandas scikit-learn scipy numpy pyyaml biopython
conda activate battleamp-snakemake
```

Clone the repository:

```bash
git clone https://github.com/szczurek-lab/battleamp-snakemake.git
cd battleamp-snakemake
git submodule update --init --recursive
```

> **Conda solver.** conda >= 23.x uses libmamba by default, which is required
> to reliably solve the conflicting dependencies across models. If you are on an
> older version, run `conda update -n base conda` first.


## Quick Start

The repository ships with `example-dataset` (20 peptides) and `example-model`
(a random baseline) for verifying that the installation works.

```bash
snakemake --profile profile/ score \
    --config fasta="datasets/example-dataset/sequences.fasta" run_models="example-model"
```

On the first run Snakemake creates conda environments for each model and caches
them in `.snakemake/conda/`. Subsequent runs reuse them.

Predictions are written to:

```
results/inference/example-model/example-dataset/predictions.tsv
```

> **Always use `--profile profile/`.** It sets `--use-conda`, `--keep-going`,
> and other required flags. Running snakemake directly without it will fail on
> any rule that has a `conda:` directive.


## I want to score my own peptides

Run all integrated models on your sequences and get AMP/non-AMP predictions
(classifiers) and estimated MIC values (regressors).

### Quickest path: pass a FASTA directly

No config changes needed. The dataset name is taken from the filename stem.

```bash
# All models
snakemake --profile profile/ score --config fasta="/path/to/my_peptides.fasta"

# Selected models only
snakemake --profile profile/ score \
    --config fasta="/path/to/my_peptides.fasta" run_models="ampeppy,amplify,ampredmfa"
```

Predictions land in `results/inference/{model_variant}/my_peptides/predictions.tsv`.

Two constraints apply: the filename stem must not clash with a dataset name
already defined in `config/config.yaml`, and evaluation against ground-truth
labels is not available via this path. Use the full registration workflow below
if you need metrics.

### Full registration: dataset + evaluation

**1. Add a directory under `datasets/`:**

```
datasets/my-peptides/
    sequences.fasta
    labels.tsv        # optional, needed only for evaluation
```

FASTA format:

```
>peptide_1
KFLQSARKILGK
>peptide_2
GIKLSARKVFPA
```

Classification labels TSV (`sequence`, `label`):

```
sequence        label
KFLQSARKILGK   AMP
GIKLSARKVFPA   non-AMP
```

Regression labels TSV (`sequence`, `MIC`, `MIC_unit`):

```
sequence        MIC    MIC_unit
KFLQSARKILGK   4.2    ug/ml
GIKLSARKVFPA   16.0   uM
```

Mixed units within the same file are supported; each row is converted
individually to the `benchmark_unit` set in `config/config.yaml`.

**2. Register in `config/config.yaml`:**

```yaml
datasets:
  my-peptides:
    sequences: datasets/my-peptides/sequences.fasta

tasks:
  my-classification:
    type: classification
    dataset: my-peptides
    labels: datasets/my-peptides/labels.tsv
    model_type: [classifier, regressor]   # regressors auto-converted via activity_thresholds

  my-regression:
    type: regression
    dataset: my-peptides
    labels: datasets/my-peptides/labels.tsv
    model_type: regressor
```

**3. Run:**

```bash
snakemake --profile profile/ score      # inference only
snakemake --profile profile/            # inference + evaluation
cat results/aggregated/summary.tsv
```

To restrict to specific models without editing the config, append
`--config run_models="ampeppy,amplify"`. Model names must match entries in the
`models:` list in `config/config.yaml`; an unrecognised name errors immediately
with a list of valid options.

### Output format

**Classifiers** produce a TSV with columns `sequence`, `Prediction`,
`Probability_score`. **Regressors** produce `sequence`, `MIC`, `MIC_unit`.
Sequences outside a model's accepted length range are silently skipped.
Evaluation metrics appear in `results/evaluation/{variant}/{task}/metrics.json`
and are aggregated into `results/aggregated/summary.tsv`.


## I want to reproduce the benchmark

The full benchmark runs all models in `config/config.yaml` across all defined
datasets and tasks.

```bash
snakemake --profile profile/
cat results/aggregated/summary.tsv
```

Preview what will run without executing anything:

```bash
snakemake --profile profile/ --dryrun
```

Run only a subset of models without editing the config:

```bash
snakemake --profile profile/ --config run_models="ampeppy,amplify,ampredmfa"
```

The pipeline runs five stages automatically: **setup** (download weights,
compile extensions), **pre-filter** (remove sequences outside each model's
length limits), **inference**, **evaluate** (compute metrics against labels),
and **aggregate** (combine all metrics into summary tables). Completed steps are
cached and not re-run.


## I want to add my own model

Each model in the benchmark lives in its own Git repository and is tracked as a
**Git submodule** under `models/`. This keeps model code, weights, and version
history independent from the pipeline while allowing reproducible pinning to
specific commits.

### Short version

```bash
# 1. Fork or create a repo for your model, then add it as a submodule

git submodule add https://github.com/your-org/your-model.git models/your-model

# 2. Ensure these interface files exist inside models/your-model/
#    model.yaml         - metadata (name, type, length limits, GPU)
#    environment.yaml   - conda env spec
#    inference.sh       - reads input FASTA ($1), writes output TSV ($2)
#    setup.sh           - optional: download weights, compile extensions

# 3. Register and test
#    Add your-model to models: in config/config.yaml
#    Add an entry to models/registry.yaml, then:
python scripts/generate_model_table.py
snakemake --profile profile/ score \
    --config score_datasets="[example-dataset]" run_models="your-model"
```

`inference.sh` must write a TSV in the classifier or regressor format described
in [Output format](#output-format) above. The pipeline pre-filters the input
FASTA to the length range declared in `model.yaml`, so `inference.sh` does not
need to handle out-of-range sequences.

See [docs/adding_a_model.md](docs/adding_a_model.md) for a full walkthrough
covering multi-variant models, validation reference data, and common pitfalls.


## Reference

### Available Models

<!-- MODEL_TABLE_START -->
| Model | Variants | Type | Framework | Accepted lengths | GPU |
|-------|----------|------|-----------|------------------|-----|
| AMPscannerv2 | (single) | classifier | CNN, LSTM | 10 to 200 | yes |
| AmPEPpy | (single) | classifier | RF | Unlimited | no |
| AMPlify | (single) | classifier | LSTM, ATT | 1 to 199 | yes |
| AMPred-MFA | (single) | classifier | LSTM, CNN, ATT | >= 3 | yes |
| HydrAMP | HydrAMP-AMP, HydrAMP-MIC | classifier | LSTM | 1 to 25 | yes |
| MolE | MolE-max | classifier | XGBoost | Unlimited | no |
| sAMP-VGG16 | (single) | classifier | CNN | Unlimited | yes |
| SenseXAMP | SenseXAMP-classifier | classifier | ESM, ATT | 6 to 25 | yes |
| AMPredictor | (single) | classifier | GCN, ESM | 1 to 65 | yes |
| APEX | APEX-Ecoli, APEX-Saureus, APEX-Kpneumoniae, APEX-min | regressor | ATT, RNN | 1 to 52 | yes |
| Deep-AMP | Deep-AMP-CNN-Gram+, Deep-AMP-CNN-Gram-, Deep-AMP-LSTM-Gram+, Deep-AMP-LSTM-Gram- | regressor | CNN, LSTM | 1 to 49 | yes |
| MBC-Attention | (single) | regressor | CNN, ATT | 5 to 60 | yes |
| SenseXAMP | SenseXAMP-Saureus, SenseXAMP-Ecoli | regressor | ESM, ATT | 6 to 25 | yes |

Total: 9 classifiers, 4 regressors, 21 variants
<!-- MODEL_TABLE_END -->


### Evaluation Metrics

#### Classification

| Metric | Description |
|--------|-------------|
| accuracy | Fraction of correctly classified sequences |
| mcc | Matthews correlation coefficient (-1 to +1) |
| f1 | Harmonic mean of precision and recall |
| precision | TP / (TP + FP) |
| recall / tpr | TP / (TP + FN) |
| fpr | FP / (FP + TN) |
| tnr | TN / (TN + FP) |
| informedness | TPR + TNR - 1; 0 = random, 1 = perfect |
| pos_preds | Fraction of all predictions that are positive |
| auroc | Area under the ROC curve |
| auprc | Area under the precision-recall curve |
| precision_at_k | Precision among top-k predictions (default k=100) |
| pauroc_01 | Partial AUROC at FPR <= 0.1 (McClish-normalized) |
| pauroc_001 | Partial AUROC at FPR <= 0.01 |

When a regressor is evaluated on a classification task, MIC predictions are
binarized using `activity_thresholds`. Sequences in the grey zone between
`active` and `inactive` thresholds are excluded. The negative MIC is used as
the ranking score for AUROC and precision@k.

#### Regression

| Metric | Description |
|--------|-------------|
| r2 | Coefficient of determination (linear scale) |
| r2_log2 | Coefficient of determination (log2 scale) |
| rmse | Root mean squared error (linear scale) |
| rmsle | Root mean squared log error (natural log) |
| rmsl2e | Root mean squared log2 error |
| spearman | Spearman rank correlation |
| mae | Mean absolute error (not in default set; add explicitly) |

Both predicted and ground-truth MIC values are clamped to the range in
`mic_clamp` before computing regression metrics, reflecting assay detection
limits. Remove the key from the config to disable clamping.


### Configuration

All settings live in `config/config.yaml`.

**`benchmark_unit`** (`"ug/ml"` or `"uM"`). Unit all MIC values are converted
to before evaluation. Per-peptide molecular weight is used for conversion.

**`activity_thresholds`**. Controls binarization of regressor output for
classification tasks. `active`: MIC at or below this is positive. `inactive`:
MIC at or above this is negative. Sequences between the two thresholds are
excluded from evaluation (grey zone).

**`mic_clamp`**. Clamps MIC values to an assay-realistic range before
regression metrics are computed.

**`models`**. List of model directory names to include. Override at runtime
without editing the file: `--config run_models="model1,model2"`.

**`datasets`**. Each entry points to a FASTA file. Models run inference once
per dataset. A dataset can also be supplied at runtime via
`--config fasta="/path/to/file.fasta"`.

**`tasks`**. Each task pairs a dataset with ground-truth labels and a metric
type (`classification` or `regression`). `model_type` controls which model
types are evaluated on that task.


### Output Structure

```
results/
    inference/{variant}/{dataset}/predictions.tsv
    evaluation/{variant}/{task}/metrics.json
    aggregated/
        summary.tsv
        classification_results.tsv
        regression_results.tsv
    validation/
        {model}/validation_report.json
        validation_summary.tsv
    logs/
        setup/{model}.log
        inference/{variant}/{dataset}.log
        evaluation/{variant}/{task}.log
```


### Running on an HPC Cluster

#### Single-node (interactive or batch allocation)

The bundled `profile/` works on a single node with direct GPU access. The
simplest way to use GPUs on a SLURM cluster without the executor plugin is to
request an interactive allocation and run snakemake inside it:

```bash
srun --partition=<your-partition> --gres=gpu:1 \
     --mem=64G --cpus-per-task=8 --time=4:00:00 --pty bash

# Once inside the allocation:
source ~/.venvs/battleamp-snakemake/bin/activate
snakemake --profile profile/ score \
    --config score_datasets="[<your-dataset>]"
```

This is the most robust option when the SLURM executor plugin is unavailable
or the cluster has strict job submission policies.

#### Multi-job SLURM execution

To submit each pipeline rule as a separate SLURM job, install the executor
plugin into the pipeline venv and use the bundled `slurm/` profile:

```bash
pip install snakemake-executor-plugin-slurm
snakemake --profile slurm/ score \
    --config score_datasets="[<your-dataset>]"
```

**Creating `slurm/config.yaml`.** The repository ships a template; copy and
adapt it to your cluster before first use:

```yaml
executor: slurm
jobs: 10                      # keep below your QOS job limit

use-conda: true
keep-going: true
conda-prefix: .snakemake/conda

default-resources:
  slurm_partition: <your-partition>
  slurm_account: <your-account>   # required; omitting this stalls submission silently
  mem_mb: 16000
  runtime: 240                    # minutes
  cpus_per_task: 4

set-resources:
  run_inference:
    mem_mb: 64000
    runtime: 720
    cpus_per_task: 8
    slurm_extra: "'--gres=gpu:1'"
  run_multioutput_inference:
    mem_mb: 64000
    runtime: 720
    cpus_per_task: 8
    slurm_extra: "'--gres=gpu:1'"
  model_setup:
    mem_mb: 16000
    runtime: 120
    cpus_per_task: 4
```

`slurm_account` must be set explicitly. Without it, snakemake attempts to
guess the account from `sacctmgr`, which can stall or pick the wrong account
silently.

`jobs` must stay at or below your cluster QOS limit. If you are unsure:

```bash
sacctmgr show qos format=name,maxjobspu
```

**Verifying submission.** Before running the full benchmark, test with the
example dataset to confirm jobs actually appear in the queue:

```bash
snakemake --profile slurm/ score \
    --config fasta="datasets/example-dataset/sequences.fasta" \
    run_models="example-model"

# In another terminal:
squeue -u $USER
```

#### Slow internet nodes

On nodes without direct internet access, conda environment creation can time
out when pulling packages. Set a longer timeout before running:

```bash
export PIP_DEFAULT_TIMEOUT=300
```

If conda repodata fails to download (error: `An error occurred when loading
cached repodata`), clear the index cache and retry:

```bash
conda clean --index-cache
```


### Model Validation

The validation stage checks that code changes have not altered model predictions
by re-running inference on stored reference inputs and comparing outputs.

```bash
snakemake --profile profile/ validate
cat results/validation/validation_summary.tsv
```

See [docs/validation.md](docs/validation.md) for details.


## License

MIT
