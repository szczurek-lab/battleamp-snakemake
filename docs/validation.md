# Model Validation

This pipeline includes a validation stage that verifies code modifications
(hotfixes, adapter scripts, environment changes) did not alter model behavior.

## Why Validate?

Many AMP prediction models required hotfixes to run in a standardized
benchmark setting: fixing broken imports, adapting input/output formats,
updating deprecated API calls, etc. A reviewer will ask: "How do you know
these changes didn't affect model predictions?" The validation stage answers
this question with empirical evidence.

## How It Works

1. Each model stores **reference data** in `models/{model}/validation/`:
   - `reference_input.fasta`: a set of sequences (ideally 100-200)
   - `reference_output.tsv`: the model's predictions on those sequences
     from the ORIGINAL (pre-hotfix) code

2. The pipeline runs the adapted model on `reference_input.fasta` and
   compares the new predictions against `reference_output.tsv`

3. For classifiers, it checks:
   - Predictions (AMP/non-AMP) match exactly
   - Probability scores match within tolerance (default 1e-4)
   - All reference sequences are present in the output

4. For regressors, it checks:
   - MIC values match within relative tolerance (default 1e-4)
   - All reference sequences are present

5. A summary report lists which models passed and which failed

## Running Validation

```bash
# Validate all models that have reference data
snakemake --snakefile workflow/Snakefile validate --use-conda --cores 4

# Check results
cat results/validation/validation_summary.tsv

# Detailed per-model report
cat results/validation/ampeppy/validation_report.json
```

## Generating Reference Data

Reference data must be generated BEFORE adapting the model code. There are
two approaches:

### Approach A: From existing cached predictions

If the model repo already contains cached predictions (e.g., from previous
benchmark runs), you can extract reference data from those:

```bash
python workflow/scripts/generate_reference.py \
    --predictions /path/to/cached_predictions.tsv \
    --fasta /path/to/input_sequences.fasta \
    --model-dir models/ampeppy \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

If the cached predictions have sequence IDs but not actual sequences (common),
provide the original FASTA and the ID column name:

```bash
python workflow/scripts/generate_reference.py \
    --predictions /path/to/cached_predictions.tsv \
    --fasta /path/to/input_sequences.fasta \
    --seq-id-col seq_id \
    --model-dir models/ampeppy \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

### Approach B: Run the original model fresh

If you have the original (unmodified) model working:

```bash
# Activate the model's original environment
conda activate ampeppy-original

# Generate reference by running inference
python workflow/scripts/generate_reference.py \
    --run-inference \
    --fasta /path/to/benchmark_sequences.fasta \
    --model-dir models/ampeppy \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

This samples 200 sequences from the FASTA, runs the model's `inference.sh`,
and saves both the input and output as reference.

### Choosing sequences for validation

For the strongest evidence, include sequences from both:
- **Training data**: verifies the model was loaded/trained correctly
- **Test data**: verifies inference works correctly on unseen data
- **Both AMP and non-AMP**: ensures both classes are validated

The `generate_reference.py` script handles stratified sampling automatically.
With `--n-samples 200`, you get approximately 100 AMPs and 100 non-AMPs.

## Adjusting Tolerance

Some models produce slightly different float values across environments
(different numpy/sklearn versions, CPU vs GPU, etc.). The tolerance is
configurable per model in `model.yaml`:

```yaml
# model.yaml
validation_tolerance: 1e-4    # absolute difference for probabilities
                               # relative difference for MIC values
```

If validation fails only due to tiny float differences, increase the
tolerance. If predictions or classifications change, something is wrong.

## File Formats

### reference_input.fasta

Standard FASTA. Headers are used for logging only; matching is done by
sequence content.

```
>DBAASP_3951
FLGKVFKLASKVFPAVFGKV
>DBAASP_14675
GRFKRFRKKFKKLFKKLSPVIPLLHLG
```

### reference_output.tsv (classifier)

```
sequence	Prediction	Probability_score
FLGKVFKLASKVFPAVFGKV	AMP	0.9667
GRFKRFRKKFKKLFKKLSPVIPLLHLG	AMP	0.9031
```

### reference_output.tsv (regressor)

```
sequence	MIC	MIC_unit
FLGKVFKLASKVFPAVFGKV	4.2	ug/ml
GRFKRFRKKFKKLFKKLSPVIPLLHLG	16.0	uM
```
