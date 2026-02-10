# Adding a Model

This guide explains how to add a new AMP prediction model to the pipeline.

## Step 1: Create the model directory

Create a directory under `models/` with your model name:

```bash
mkdir models/my-model
```

## Step 2: Create model.yaml

This file declares your model's metadata. The pipeline reads it to determine
how to handle the model (type, length constraints, GPU needs, variants).

```yaml
# models/my-model/model.yaml

name: my-model
version: 1.0.0
type: classifier           # "classifier" or "regressor"
framework: transformer     # free text, for documentation
length_min: 5              # minimum sequence length in amino acids, null = no limit
length_max: 200            # maximum sequence length, null = no limit
gpu_required: true         # whether the model needs a GPU
conda_env: my-model        # name of the conda env (must match environment.yaml)
```

### Multi-variant models

If your model produces multiple benchmark entries (e.g. species-specific predictions),
declare variants:

```yaml
name: my-model
version: 1.0.0
type: regressor
framework: transformer
length_min: null
length_max: null
gpu_required: true
conda_env: my-model

variants:
  - name: my-model-ecoli
    args: ["ecoli"]           # passed to inference.sh as $3, $4, ...
  - name: my-model-saureus
    args: ["saureus"]
  - name: my-model-min
    args: ["min"]
```

Each variant becomes a separate row in the benchmark results.

## Step 3: Create environment.yaml

Standard conda environment specification. Snakemake creates and manages this
environment automatically when the user runs with `--use-conda`. Nobody needs
to manually create or activate any conda env -- Snakemake does it all.

```yaml
# models/my-model/environment.yaml

name: my-model
channels:
  - conda-forge
  - pytorch
  - defaults
dependencies:
  - python=3.10
  - pytorch=2.0
  - numpy
  - pandas
  - pip:
    - fair-esm==2.0.0
```

## Step 4: Create inference.sh

This is the core script. The pipeline calls it with:

```bash
bash inference.sh <input.fasta> <output.tsv> [extra_args...]
```

### Classifier output format

```
sequence    Prediction    Probability_score
KFLQ...     AMP           0.87
GIKL...     non-AMP       0.23
```

### Regressor output format

```
sequence    MIC    MIC_unit
KFLQ...     4.2    ug/ml
GIKL...     16.0   uM
```

### Rules

- Column names are **case-sensitive** and must match exactly
- Sequences the model cannot process must be **skipped** with a warning to stderr
- Output order does not need to match input order
- Exit code 0 on success, non-zero on failure
- The conda env from environment.yaml is **already active** when this script runs.
  Do NOT run `conda activate` inside inference.sh.

### Example inference.sh

```bash
#!/bin/bash
INPUT_FASTA="$1"
OUTPUT_TSV="$2"
VARIANT_ARG="${3:-default}"

python predict.py \
    --input "$INPUT_FASTA" \
    --output "$OUTPUT_TSV" \
    --variant "$VARIANT_ARG"
```

## Step 5: Create setup.sh (optional)

If your model needs to download weights, compile extensions, or do other one-time
setup beyond what conda handles, put it in `setup.sh`. This runs inside the
model's conda env (already active). Do NOT create or activate any conda env here.

```bash
#!/bin/bash
# Download model weights
if [ ! -f weights/model.pt ]; then
    wget https://example.com/model.pt -O weights/model.pt
fi
```

This runs once before inference.

## Step 6: Generate validation reference data (recommended)

Before making any code changes, capture the original model's predictions
as a reference baseline. See [docs/validation.md](validation.md) for full details.

```bash
# Option A: From existing cached predictions
python workflow/scripts/generate_reference.py \
    --predictions /path/to/original_predictions.tsv \
    --fasta /path/to/input.fasta \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200

# Option B: Run the original code fresh
python workflow/scripts/generate_reference.py \
    --run-inference \
    --fasta /path/to/input.fasta \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200
```

This creates `models/my-model/validation/reference_input.fasta` and
`reference_output.tsv`. After adapting the model, run `snakemake validate`
to confirm your changes didn't alter predictions.

## Step 7: Register in config.yaml

Add your model to the `models:` list in `config/config.yaml`:

```yaml
models:
  - example-model
  - my-model        # <-- add this
```

If your model has variants, make sure the relevant tasks either use `models: all`
or list the variant names explicitly.

## Testing your model

Run the pipeline with just your model on the example dataset:

```bash
# Edit config to include only your model
snakemake --snakefile workflow/Snakefile --use-conda --cores 1

# Check the output
cat results/inference/my-model/example-dataset/predictions.tsv
cat results/inference/my-model/example-dataset/validation_report.json
```

The validation report will tell you if your output format is correct.
