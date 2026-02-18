# battleamp-snakemake

A Snakemake pipeline for benchmarking antimicrobial peptide (AMP) prediction
models. It ships with a curated set of classifiers and regressors that you can
run on your own peptide sequences to get standardized predictions and evaluation
metrics. The pipeline is model-agnostic and dataset-agnostic: bring your own
data, bring your own models, or use the ones already integrated.

## Table of Contents

- [Installation](#installation)
- [The Benchmark](#the-benchmark)
- [Configuration](#configuration)
- [Available Models](#available-models)
- [Evaluation Metrics](#evaluation-metrics)
- [Scoring Your Own Peptides](#scoring-your-own-peptides)
- [Adding Your Own Model](#adding-your-own-model)
- [Running on an HPC Cluster](#running-on-an-hpc-cluster)
- [Model Validation](#model-validation)
- [Output Structure](#output-structure)
- [License](#license)


## Installation

### Prerequisites

- Python >= 3.10
- conda (Miniconda or Anaconda) for model environment management
- Git

### Option A: pip (recommended, works on clusters without admin access)

```bash
python3 -m venv ~/.venvs/battleamp-snakemake
source ~/.venvs/battleamp-snakemake/bin/activate
pip install snakemake pulp pandas scikit-learn scipy numpy pyyaml biopython
```

On HPC clusters you may need to load a Python module first
(`module load python/3.10` or similar).

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
git submodule update --init --recursive
```


## The Benchmark

The pipeline runs five stages, all orchestrated by Snakemake. Snakemake handles
dependency resolution, parallelism, and caching (steps whose outputs are already
up to date are not re-run).

1. **Setup** downloads model weights and compiles extensions inside each model's
   auto-created conda environment.
2. **Pre-filter** generates model-specific FASTA files, removing sequences that
   fall outside each model's declared length constraints.
3. **Inference** runs every selected model variant on every dataset.
4. **Evaluate** compares predictions to ground-truth labels and computes metrics
   (skipped if no labels are provided).
5. **Aggregate** combines all metrics into summary tables.

### Quick start with the included example

The repository ships with a small `example-dataset` (20 peptides) and an
`example-model` for testing. To score just this dataset without running the
full benchmark:

```bash
source ~/.venvs/battleamp-snakemake/bin/activate   # or: conda activate battleamp-snakemake
snakemake score --snakefile workflow/Snakefile --use-conda --cores 4 \
    --config score_datasets="[example-dataset]"
```

Predictions are written to
`results/inference/example-model/example-dataset/predictions.tsv`.

To run the full benchmark (all models, all datasets, with evaluation against
ground-truth labels):

```bash
snakemake --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1
cat results/aggregated/summary.tsv
```

On the first run Snakemake will build a separate conda environment for each model.
These are cached in `.snakemake/conda/` and reused afterwards. You do not need to
create or manage any model environments yourself.

### Model output format

**Classifiers** produce a TSV with columns (case-sensitive):

```
sequence    Prediction    Probability_score
KFLQ...     AMP           0.87
GIKL...     non-AMP       0.23
```

**Regressors** produce:

```
sequence    MIC    MIC_unit
KFLQ...     4.2    ug/ml
GIKL...     16.0   uM
```

Sequences a model cannot handle are silently skipped (not reported as errors).
Exit code 0 on success.


## Configuration

All pipeline settings live in `config/config.yaml`. The file has four sections:
global evaluation parameters, models, datasets, and tasks.

### Full annotated example

```yaml
# ── Global evaluation parameters ──────────────────────────────────────────

# All MIC values (predictions and ground truth) are converted to this unit
# before any comparison.  Per-peptide molecular weight is computed from the
# amino acid sequence for accurate uM <-> ug/ml conversion.
benchmark_unit: "ug/ml"     # options: "ug/ml" or "uM"

# When a regressor is evaluated on a classification task, its MIC predictions
# are binarized using these thresholds (expressed in the unit given below):
#   active:   MIC <= 32  ug/ml  ->  predicted positive  (AMP)
#   inactive: MIC >= 128 ug/ml  ->  predicted negative  (non-AMP)
#   grey zone: 32 < MIC < 128   ->  excluded from evaluation
# The grey zone avoids penalizing models for borderline cases where the
# ground truth itself is ambiguous.
activity_thresholds:
  active: 32
  inactive: 128
  unit: "ug/ml"

# MIC predictions and ground-truth values are clamped to this range before
# computing regression metrics.  Values outside the range are set to the
# nearest boundary.  This prevents extreme outlier predictions from
# dominating error metrics and reflects the realistic detection limits of
# MIC assays.
mic_clamp:
  min: 0.25
  max: 512
  unit: "ug/ml"

# ── Models ────────────────────────────────────────────────────────────────

# Each entry must have a corresponding directory under models/ containing
# at least model.yaml and inference.sh.
models:
  - ampeppy
  - amplify
  - ampscanner
  - ampredmfa
  - hydramp-amp-classifier
  - hydramp-mic-classifier
  - apex
  - mole-amp
  - mbc-attention
  - sensexamp
  - ampredictor
  - deep-amp

# ── Datasets ──────────────────────────────────────────────────────────────

# A dataset is a FASTA file that models run inference on.  Models run once
# per dataset.  Multiple tasks can share the same dataset.
datasets:
  battleamp-all:
    sequences: datasets/battleamp-all/sequences.fasta
  deep-amp:
    sequences: datasets/deep-amp/all_sequences.fasta

# ── Tasks ─────────────────────────────────────────────────────────────────

# A task pairs a dataset with ground-truth labels and a metric type.  Both
# classifiers and regressors can be evaluated on classification tasks (the
# pipeline converts MIC predictions to binary labels automatically).

tasks:
  amp:
    type: classification
    dataset: battleamp-all
    labels: tasks/amp/labels.tsv
    model_type:             # which model types this task applies to
      - classifier
      - regressor

  broad_activity:
    type: classification
    dataset: battleamp-all
    labels: tasks/broad_activity/labels.tsv
    model_type:
      - classifier
      - regressor

  regression_ecoli25922:
    type: regression
    dataset: battleamp-all
    labels: tasks/regression_ecoli25922/labels.tsv
    model_type: regressor   # regression tasks only apply to regressors
```

### Parameter reference

**benchmark_unit** (`"ug/ml"` or `"uM"`, default `"ug/ml"`). The unit all MIC
values are converted to before evaluation. Conversion uses each peptide's
molecular weight, computed from its amino acid sequence.

**activity_thresholds**. Controls how regressor predictions are turned into
binary labels for classification tasks. Three fields:

- `active` (float): MIC at or below this value is considered active (positive).
- `inactive` (float): MIC at or above this value is considered inactive
  (negative).
- `unit` (string): the unit the thresholds are expressed in. If this differs
  from `benchmark_unit`, the pipeline converts per-peptide using molecular
  weight.

Predictions that fall between the two thresholds (the "grey zone") are excluded
from classification evaluation entirely, because at these concentrations the
biological activity is ambiguous.

**mic_clamp** (optional, omit to disable). Clamps both predicted and ground-truth
MIC values to an assay-realistic range before computing regression metrics.
Three fields: `min`, `max`, and `unit`. Ground truth is clamped too because
experimental MIC values beyond the assay detection limits are themselves
unreliable. Set to `null` or remove the key entirely to disable clamping.

**datasets**. Each entry has a name (used as key) and a `sequences` field
pointing to a FASTA file. Models run inference once per dataset; the results
are then reused across all tasks that reference that dataset.

**tasks**. Each entry defines an evaluation and has the following fields:

- `type`: `classification` or `regression`.
- `dataset`: name of a dataset defined above.
- `labels`: path to a TSV file with ground-truth labels.
- `model_type`: which model types this task applies to. Can be a single string
  (`"regressor"`) or a list (`["classifier", "regressor"]`).
- `metrics` (optional): list of metric names to compute. If omitted, all
  default metrics for the task type are computed.
- `precision_at_k` (optional): override the default k for the precision@k
  metric with `precision_at_k: { k: 50 }`.

### Cross-evaluation: regressors on classification tasks

When a regressor is evaluated on a classification task, the pipeline
automatically converts its MIC predictions into binary labels using the
`activity_thresholds`. The probability proxy used for ranking-based metrics
(AUROC, AUPRC, precision@k, partial AUROC) is the negative MIC value, so that
peptides with lower predicted MIC (stronger predicted activity) rank higher.
Samples in the grey zone are excluded from evaluation. The evaluation report
records the number of grey-zone exclusions and all threshold settings for full
transparency.

Classifiers on regression tasks cannot produce meaningful regression metrics and
are skipped with a note in the output.


## Available Models

The table below is rendered from `models/registry.yaml`, which is the single
source of truth for all integrated models. When you add or change a model,
update that file and run:

```bash
python scripts/generate_model_table.py
```

This replaces the table in this README in place (between the marker comments).
You can also use `--check` in CI to verify the table is up to date, or
`--stdout` to print it without modifying any files.

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

The pre-filter step automatically removes sequences that fall outside a model's
accepted range, so you do not need to trim your input manually. Keep in mind that
very short or very long peptides will not receive predictions from every model.


## Evaluation Metrics

When ground-truth labels are provided, the pipeline computes metrics
appropriate to the task type. You can override the default set per task by
adding a `metrics:` list to the task definition in `config/config.yaml`.

### Classification metrics

All of these are computed by default for classification tasks. Threshold-based
metrics (accuracy through pos_preds) use the binary predicted label.
Score-based metrics (auroc through pauroc_001) use the continuous score: the
predicted probability for classifiers, or the negative MIC for regressors
evaluated on a classification task.

| Metric | Description |
|--------|-------------|
| accuracy | Fraction of correctly classified sequences |
| mcc | Matthews correlation coefficient; balanced measure that accounts for all four confusion matrix categories, ranging from -1 (total disagreement) to +1 (perfect prediction) |
| f1 | Harmonic mean of precision and recall |
| precision | Fraction of positive predictions that are true positives (TP / (TP + FP)) |
| recall | Fraction of actual positives correctly identified (TP / (TP + FN)); equivalent to tpr |
| fpr | False positive rate; fraction of negatives incorrectly predicted as positive (FP / (FP + TN)) |
| tpr | True positive rate (sensitivity); fraction of positives correctly identified (TP / (TP + FN)) |
| tnr | True negative rate (specificity); fraction of negatives correctly identified (TN / (TN + FP)) |
| informedness | Balanced accuracy adjusted for chance (TPR + TNR - 1); 0 means random, 1 means perfect |
| pos_preds | Fraction of all predictions that are positive ((TP + FP) / N); useful for checking whether a model predicts nearly everything as positive or negative |
| auroc | Area under the receiver operating characteristic curve; measures discrimination across all thresholds |
| auprc | Area under the precision-recall curve; especially informative when classes are imbalanced |
| precision_at_k | Precision among the top-k highest-scoring predictions (default k = 100); measures how many of the top-ranked candidates are true positives, directly relevant for prioritizing peptides for experimental validation. Override k per task with `precision_at_k: { k: 50 }` |
| pauroc_01 | Partial AUROC restricted to FPR <= 0.1 (McClish-standardized to [0, 1]); focuses evaluation on the low-false-positive regime where screening applications typically operate |
| pauroc_001 | Partial AUROC restricted to FPR <= 0.01; stricter variant for scenarios where even a 10% false positive rate is too costly |

### Regression metrics

All of these except `mae` are computed by default. Add `mae` to the task's
`metrics:` list if you want it.

| Metric | Description |
|--------|-------------|
| r2 | Coefficient of determination on the linear (raw MIC) scale |
| r2_log2 | Coefficient of determination on the log2(1 + MIC) scale; less sensitive to outliers and more appropriate when MIC values span several orders of magnitude |
| mse | Mean squared error (linear scale) |
| rmse | Root mean squared error (linear scale); in the same unit as MIC |
| msle | Mean squared logarithmic error (natural log); computed as mean((ln(1+y_true) - ln(1+y_pred))^2), so it penalizes under-predictions of low-MIC (potent) peptides more than over-predictions of high-MIC ones |
| rmsle | Root of msle; in log-space units |
| msl2e | Mean squared log2 error; same idea as msle but on the log2 scale, which aligns with the conventional doubling dilution series used in MIC assays |
| rmsl2e | Root of msl2e |
| spearman | Spearman rank correlation; measures how well predicted and observed MIC values agree in rank order, without assuming a linear relationship. Values near +1 mean the model preserves the activity ranking even if the absolute MIC values are off |
| mae | Mean absolute error (linear scale); average absolute deviation between predicted and observed MIC. Not in the default set; add explicitly if needed |

### Regression metrics and MIC clamping

When `mic_clamp` is set in the config, both predictions and ground-truth values
are clamped to the specified range before computing regression metrics. This is
important because MIC assays have physical detection limits (typically 0.25 to
512 ug/ml in a standard two-fold dilution setup). A model predicting 0.001 ug/ml
for a peptide whose true MIC is 0.25 ug/ml should not be heavily penalized,
since 0.25 is itself the assay floor. Clamping ensures that metrics reflect
real predictive accuracy within the assay's operational range. The number of
clamped values is recorded in each evaluation report.


## Scoring Your Own Peptides

This section walks you through adding your own dataset so that the integrated
models score your peptide sequences. By default the pipeline runs both
classifiers and regressors, so you get AMP/non-AMP predictions together with
estimated MIC values from a single run.

### 1. Prepare your input files

Create a directory under `datasets/`:

```
datasets/my-peptides/
    sequences.fasta    # your peptide sequences
    labels.tsv         # ground truth (optional, needed only for evaluation)
```

**sequences.fasta** contains your peptides in standard FASTA format:

```
>peptide_1
KFLQSARKILGK
>peptide_2
GIKLSARKVFPA
>peptide_3
RWKIFKKIEKMGRNIRDGIVKAGPAIEVLGSAKAIGK
```

**labels.tsv** is required only if you want the pipeline to evaluate predictions
against known ground truth. For a classification task, provide a tab-separated
file with two columns:

```
sequence	label
KFLQSARKILGK	AMP
GIKLSARKVFPA	non-AMP
```

For a regression task (MIC prediction):

```
sequence	MIC	MIC_unit
KFLQSARKILGK	4.2	ug/ml
GIKLSARKVFPA	16.0	uM
```

The pipeline handles mixed units within a single file: each row's `MIC_unit` is
used for per-peptide conversion to the `benchmark_unit`.

If you only need raw predictions with no evaluation, omit `labels.tsv` and do
not define any tasks for this dataset.

### 2. Register your dataset and tasks in the config

Edit `config/config.yaml`. Add your dataset under `datasets:`, then optionally
define tasks under `tasks:`:

```yaml
datasets:
  my-peptides:
    sequences: datasets/my-peptides/sequences.fasta

tasks:
  # Classification: is this peptide an AMP?
  my-classification:
    type: classification
    dataset: my-peptides
    labels: datasets/my-peptides/labels.tsv
    model_type:
      - classifier
      - regressor          # regressors are auto-converted via activity_thresholds

  # Regression: what is the predicted MIC?
  my-regression:
    type: regression
    dataset: my-peptides
    labels: datasets/my-peptides/labels_mic.tsv
    model_type: regressor  # only regressors can produce MIC predictions
```

The `model_type` field controls which models are evaluated on each task. Setting
it to `["classifier", "regressor"]` on a classification task means both model
types are evaluated. Regressors are automatically cross-evaluated using the
global `activity_thresholds`.

The `models:` list in the config already contains all integrated models by
default. Comment out any you want to skip.

### 3. Run the pipeline

To score your peptides without running the full benchmark evaluation:

    snakemake score --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1

This runs all models on all datasets listed in the config. To restrict
scoring to specific datasets:

    snakemake score --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1 \
        --config score_datasets="[my-peptides]"

To run the full benchmark (with evaluation against ground-truth labels):

    snakemake --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1

### 4. Collect your results

Predictions for each model variant and dataset are written to:

```
results/inference/{model_variant}/{dataset}/predictions.tsv
```

Classifiers will have produced `Prediction` and `Probability_score` columns;
regressors will have produced `MIC` and `MIC_unit` columns. If you defined
evaluation tasks, metrics appear in:

```
results/evaluation/{model_variant}/{task}/metrics.json
results/aggregated/summary.tsv
```

Each `metrics.json` file also records metadata about the evaluation: number of
matched samples, number of grey-zone exclusions (for cross-evaluated
regressors), number of clamped values, and the threshold/clamping settings
that were in effect.


## Adding Your Own Model

Create a directory under `models/` with the following structure:

```
models/your-model/
    model.yaml         # metadata: name, type, length limits, GPU requirement
    environment.yaml   # conda environment spec (Snakemake creates the env from this)
    inference.sh       # prediction script: $1 = input FASTA, $2 = output TSV
    setup.sh           # (optional) download weights, compile extensions
```

**model.yaml** declares the model's properties so the pipeline can pre-filter
sequences and schedule resources correctly. The fields should match the schema
used in `models/registry.yaml`.

**environment.yaml** is a standard conda environment file. Snakemake will create
and cache the environment automatically; you never need to activate it yourself.

**inference.sh** receives two positional arguments: the path to a (pre-filtered)
FASTA file and the path where the output TSV should be written. It must exit
with code 0 on success. The output TSV must follow the model output format
described above (classifier or regressor columns).

After creating the directory, do two things:

1. Add `your-model` to the `models:` list in `config/config.yaml`.
2. Add an entry to `models/registry.yaml` and regenerate the table:
   ```bash
   python scripts/generate_model_table.py
   ```

See [docs/adding_a_model.md](docs/adding_a_model.md) for a detailed walkthrough
with examples.


## Running on an HPC Cluster

### SLURM

```bash
# Single node
snakemake --snakefile workflow/Snakefile --use-conda --cores 8 --resources gpu=1

# Submit each rule as a separate SLURM job
snakemake --snakefile workflow/Snakefile --use-conda --profile slurm/
```

### GPU scheduling

Models declare `gpu_required: true/false` in their `model.yaml`. Control GPU
scheduling with the `--resources` flag:

```bash
snakemake --use-conda --cores 8  --resources gpu=1   # single GPU
snakemake --use-conda --cores 16 --resources gpu=2   # two GPUs
snakemake --use-conda --cores 8  --resources gpu=0   # CPU only (GPU models will fail)
```

If model environment creation fails with pip timeout errors on nodes with
slow internet, set a longer timeout before running:

```export PIP_DEFAULT_TIMEOUT=300```

## Model Validation

The pipeline includes a validation stage that checks whether code modifications
have altered model predictions. Each model can store reference inputs and expected
outputs; the pipeline re-runs inference and compares.

```bash
snakemake --snakefile workflow/Snakefile validate --use-conda --cores 4
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
```


## License

MIT