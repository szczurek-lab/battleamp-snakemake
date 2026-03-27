# Adding a Model

This guide explains how to integrate a new AMP prediction model into the
pipeline.

## Step 1: Create the model directory

```bash
mkdir models/my-model
```

## Step 2: Create model.yaml

Declares the model's metadata. The pipeline reads this to pre-filter sequences
and schedule GPU resources.

```yaml
name: my-model
version: 1.0.0
type: classifier      # "classifier" or "regressor"
framework: transformer
length_min: 5         # null = no lower limit
length_max: 200       # null = no upper limit
gpu_required: true
```

### Multi-variant models

If your model produces multiple benchmark entries (e.g. species-specific
variants), declare them:

```yaml
name: my-model
version: 1.0.0
type: regressor
framework: transformer
length_min: null
length_max: null
gpu_required: true

variants:
  - name: my-model-ecoli
    args: ["ecoli"]       # passed to inference.sh as $3, $4, ...
  - name: my-model-saureus
    args: ["saureus"]
```

Each variant becomes a separate row in the benchmark results.

## Step 3: Create environment.yaml

Standard conda environment spec. Snakemake creates and caches this
automatically when running with `--use-conda` (enforced via the profile).
Do not create or activate it manually.

```yaml
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

The pipeline calls this as:

```bash
bash inference.sh <input.fasta> <output.tsv> [extra_args...]
```

The conda environment from `environment.yaml` is already active when this
script runs. Do NOT call `conda activate` inside it.

Output format for classifiers (`sequence`, `Prediction`, `Probability_score`):

```
sequence    Prediction    Probability_score
KFLQ...     AMP           0.87
GIKL...     non-AMP       0.23
```

Output format for regressors (`sequence`, `MIC`, `MIC_unit`):

```
sequence    MIC    MIC_unit
KFLQ...     4.2    ug/ml
GIKL...     16.0   uM
```

Rules:
- Column names are case-sensitive and must match exactly.
- Sequences the model cannot process must be skipped silently (warning to
  stderr is fine).
- Output order does not need to match input order.
- Exit code 0 on success, non-zero on failure.

Example:

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

For one-time setup beyond what conda handles (downloading weights, compiling
extensions, etc.). Runs inside the model's conda env before inference. Do NOT
call `conda activate` here.

```bash
#!/bin/bash
if [ ! -f weights/model.pt ]; then
    wget https://example.com/model.pt -O weights/model.pt
fi
```

## Step 6: Generate validation reference data (recommended)

Before adapting any model code, capture the original predictions as a
reference baseline. See [docs/validation.md](validation.md) for full details.

```bash
# From existing cached predictions
python workflow/scripts/generate_reference.py \
    --predictions /path/to/original_predictions.tsv \
    --fasta /path/to/input.fasta \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200

# Or by running the original code fresh
python workflow/scripts/generate_reference.py \
    --run-inference \
    --fasta /path/to/input.fasta \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200
```

After adapting the model, run `snakemake --profile profile/ validate` to
confirm predictions are unchanged.

## Step 7: Register in config/config.yaml

```yaml
models:
  - my-model
```

## Step 8: Test

Run your model on the example dataset to verify the output format:

```bash
snakemake --profile profile/ score \
    --config score_datasets="[example-dataset]" run_models="my-model"

cat results/inference/my-model/example-dataset/predictions.tsv
cat results/inference/my-model/example-dataset/validation_report.json
```

The validation report will flag any output format issues.
