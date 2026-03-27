# Model Validation

The validation stage checks that code modifications (hotfixes, adapter
scripts, environment updates) have not altered model predictions.

## Why validate?

Many integrated models required changes to run in a standardised benchmark
setting: fixing broken imports, adapting input/output formats, updating
deprecated API calls. The validation stage provides empirical evidence that
these changes did not affect predictions.

## How it works

Each model stores reference data in `models/{model}/validation/`:
- `reference_input.fasta`: a set of sequences (ideally 100-200)
- `reference_output.tsv`: the model's predictions on those sequences from
  the original, pre-modification code

The pipeline re-runs the adapted model on `reference_input.fasta` and
compares the new predictions against `reference_output.tsv`. For
classifiers, it checks that predicted labels match exactly and probability
scores agree within tolerance. For regressors, it checks that MIC values
agree within relative tolerance. A summary report lists which models passed
and which failed.

## Running validation

```bash
snakemake --profile profile/ validate
cat results/validation/validation_summary.tsv
cat results/validation/my-model/validation_report.json
```

## Generating reference data

Reference data must be generated **before** modifying any model code.

### From existing cached predictions

```bash
python workflow/scripts/generate_reference.py \
    --predictions /path/to/cached_predictions.tsv \
    --fasta /path/to/input_sequences.fasta \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

If the cached predictions contain sequence IDs rather than sequences,
provide the ID column name:

```bash
python workflow/scripts/generate_reference.py \
    --predictions /path/to/cached_predictions.tsv \
    --fasta /path/to/input_sequences.fasta \
    --seq-id-col seq_id \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

### From a fresh inference run

If you have the original unmodified model available:

```bash
conda activate my-model-original

python workflow/scripts/generate_reference.py \
    --run-inference \
    --fasta /path/to/benchmark_sequences.fasta \
    --model-dir models/my-model \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

For the strongest coverage, use sequences from both training and test sets,
and ensure both positive and negative examples are represented. The script
handles stratified sampling automatically.

## Adjusting tolerance

Some models produce slightly different float values across environments
(different library versions, CPU vs GPU). The tolerance is configurable in
`model.yaml`:

```yaml
validation_tolerance: 1e-4   # absolute for probabilities, relative for MIC
```

If validation fails only due to small float differences, increase the
tolerance. If predicted labels or classifications change, something is wrong.

## Reference file formats

`reference_input.fasta` -- standard FASTA; matching is done by sequence
content, not headers.

`reference_output.tsv` (classifier):

```
sequence	Prediction	Probability_score
FLGKVFKLASKVFPAVFGKV	AMP	0.9667
GRFKRFRKKFKKLFKKLSPVIPLLHLG	AMP	0.9031
```

`reference_output.tsv` (regressor):

```
sequence	MIC	MIC_unit
FLGKVFKLASKVFPAVFGKV	4.2	ug/ml
GRFKRFRKKFKKLFKKLSPVIPLLHLG	16.0	uM
```
