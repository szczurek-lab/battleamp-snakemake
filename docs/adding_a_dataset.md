# Adding a Dataset

This guide explains how to add a new benchmark dataset to the pipeline.

## Step 1: Create the dataset directory

```bash
mkdir datasets/my-dataset
```

## Step 2: Create sequences.fasta

Standard FASTA format. Each sequence gets a unique header:

```
>peptide_001
GIGKFLHSAKKFGKAFVGEIMNS
>peptide_002
KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK
```

## Step 3: Create labels.tsv

Tab-separated file with ground truth. Must include a `sequence` column that
matches the actual amino acid sequences in the FASTA (not the headers).

### For classification tasks

```
sequence	label
GIGKFLHSAKKFGKAFVGEIMNS	AMP
KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK	AMP
GIVEQCCTSICSLYQLENYCN	non-AMP
```

### For regression tasks

```
sequence	MIC	MIC_unit	species
GIGKFLHSAKKFGKAFVGEIMNS	4.2	ug/ml	E. coli
KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK	16.0	ug/ml	E. coli
```

You can include additional columns (species, source, etc.) -- the pipeline
will only read the columns referenced in the task definition.

## Step 4: Create dataset.yaml (optional)

Metadata for documentation:

```yaml
name: my-dataset
description: "AMP/non-AMP test set from DBAASP v3"
n_sequences: 5000
source: "https://dbaasp.org"
label_column: label
positive_label: AMP
```

## Step 5: Define a task in config.yaml

Add a task that references your dataset:

```yaml
tasks:
  my_classification_task:
    dataset: my-dataset
    type: classification
    label_column: label
    positive_label: AMP
    metrics:
      - accuracy
      - mcc
      - f1
      - auroc
    models: all
    description: "My custom AMP classification benchmark"

  my_regression_task:
    dataset: my-regression-dataset
    type: regression
    mic_column: MIC
    mic_unit_column: MIC_unit
    target_unit: ug/ml
    metrics:
      - msle_ln
      - spearman
      - r2
    models:
      - apex-ecoli
      - mbc-attention
    description: "E. coli MIC regression"
```

The `models` field controls which model variants participate. Use `all` to
include every variant of the matching type (classifiers for classification
tasks, regressors for regression tasks), or list specific variant names.
