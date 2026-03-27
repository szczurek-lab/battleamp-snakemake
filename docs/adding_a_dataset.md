# Adding a Dataset

This guide explains how to register a new benchmark dataset in the pipeline so
that all models run inference on it and, optionally, predictions are evaluated
against ground-truth labels.

## Step 1: Create the dataset directory

```bash
mkdir datasets/my-dataset
```

## Step 2: Add sequences.fasta

Standard FASTA format with unique headers:

```
>peptide_001
GIGKFLHSAKKFGKAFVGEIMNS
>peptide_002
KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK
```

## Step 3: Add labels.tsv (optional)

Required only if you want evaluation metrics. The `sequence` column must
contain the actual amino acid sequences (not the FASTA headers).

Classification:

```
sequence	label
GIGKFLHSAKKFGKAFVGEIMNS	AMP
KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK	AMP
GIVECCCTSICSLYQLENYC	non-AMP
```

Regression (MIC):

```
sequence	MIC	MIC_unit
GIGKFLHSAKKFGKAFVGEIMNS	4.2	ug/ml
KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK	16.0	uM
```

Mixed units within the same file are supported. Additional columns (species,
source, etc.) are ignored by the pipeline.

## Step 4: Register in config/config.yaml

Add your dataset under `datasets:` and define any tasks you want to evaluate:

```yaml
datasets:
  my-dataset:
    sequences: datasets/my-dataset/sequences.fasta

tasks:
  my-classification:
    type: classification
    dataset: my-dataset
    labels: datasets/my-dataset/labels.tsv
    model_type: [classifier, regressor]

  my-regression:
    type: regression
    dataset: my-dataset
    labels: datasets/my-dataset/labels.tsv
    model_type: regressor
```

`model_type` controls which model types are evaluated on the task. Setting it
to `[classifier, regressor]` on a classification task means regressors are
automatically cross-evaluated by binarizing their MIC predictions using the
`activity_thresholds` in the config. Omit `labels` and skip the task
definition entirely if you only need raw predictions with no evaluation.

## Step 5: Run

```bash
# Inference only
snakemake --profile profile/ score

# Inference + evaluation
snakemake --profile profile/
cat results/aggregated/summary.tsv
```

To score only your new dataset without re-running existing ones:

```bash
snakemake --profile profile/ score --config score_datasets="[my-dataset]"
```
