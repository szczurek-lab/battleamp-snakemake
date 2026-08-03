# Reference inputs and outputs

Files for building the web front end. All were generated from real model
predictions, not invented values.

## Scoring

```bash
snakemake --profile slurm/ score \
    --config fasta="examples/example_input.fasta" \
             output="/path/to/results.csv"
```

This writes `results.csv` and `results.report.json`. Compare against
`reference_output/`, which contains the same two inputs scored by all 22 model
variants:

| File | Contents |
|---|---|
| `example_input.fasta` | 20 valid peptides |
| `example_input_with_errors.fasta` | The same, plus a duplicate, a D-amino-acid peptide, an unknown residue, a gap character, and a 300-mer |
| `reference_output/example_input.csv` | 20 rows, 22 model columns |
| `reference_output/example_input.report.json` | Status per model, no rejections |
| `reference_output/example_input_with_errors.csv` | 21 rows, 22 model columns |
| `reference_output/example_input_with_errors.report.json` | 25 records, 21 scorable, 4 rejected |

A correct run does not score every peptide with every model. Length limits
differ per model, so `hydramp` (1-25 aa) scores 7 of 20 while `ampeppy` (no
limit) scores all 20. If every model scores every peptide, something is wrong.

## Column naming

`{variant}_prob` for classifiers, `{variant}_MIC_{unit}` for regressors.

The prefix is the variant, not the model. One model can contribute several
columns: `apex` produces six (`apex-ecoli`, `apex-saureus` and so on) from a
single run, and `sensexamp` produces one classifier and two regressors. Variant
names are unique, so columns never collide. Use the `columns` array in the
report to label the interface rather than parsing column names.

## Empty cells

An empty cell is `null` in JSON and blank in CSV. `report.json` gives the
reason under `models[variant]`:

| `status` | Meaning |
|---|---|
| `ok` | Model scored every peptide |
| `partial` | Some peptides fall outside the model's length range |
| `failed` | Model produced no output |

A peptide rejected before inference does not appear in `rows` at all. It is
listed in `input.rejected` with a reason.

Snakemake runs with `keep-going`, so a non-zero exit code is expected whenever
any single model fails. Judge a run by `status` in the report, not by the exit
code.

## messages and diagnostics

`messages` are short, human-readable, and safe to display. Each has a
`severity` (`error`, `warning`, `info`) and a `kind` (`input`, `coverage`,
`model`, `run`).

`diagnostics` holds Snakemake log paths, log tails, and the command line. These
contain absolute server paths and are for operators only.

## Caching

The dataset name is `upload_<sha1-of-cleaned-sequences>`, not the uploaded
filename:

- Two users uploading different files both named `peptides.fasta` do not collide.
- Two users submitting identical peptides share the computation.
- Running 3 models today and 5 tomorrow re-runs only the 5. Each response
  contains exactly the models requested in that call.

## Python API

Two functions the interface needs before scoring starts:

```python
import battleamp

battleamp.list_models()                       # JSON: models, variants, length limits
battleamp.validate(fasta_text, models=None)   # JSON: per-peptide validity, model coverage
```

Both return in milliseconds and are safe to call from a request handler. Use
`list_models()` for the model picker and `validate()` to report bad input before
spending GPU time.

`battleamp.score()` returns the same table as JSON, but takes minutes to hours
because it builds conda environments and loads model weights. Call it from a
background worker, never from a request handler.

| File | Shows |
|---|---|
| `example_list_models.json` | Model catalogue for the picker |
| `example_validate.json` | Validation of the clean input |
| `example_validate_with_errors.json` | Every rejection reason and warning state |
| `example_score.json` / `.tsv` | Full `score()` response, MIC in µg/ml |
| `example_score_uM.json` / `.tsv` | The same with `unit="uM"`; column names and values both change |

## Server setup

`score()` invokes Snakemake, which invokes `conda`. Neither is on the PATH of a
typical service worker:

```bash
export BATTLEAMP_SNAKEMAKE=$HOME/.venvs/battleamp-snakemake/bin/snakemake
source $HOME/miniforge3/etc/profile.d/conda.sh
export BATTLEAMP_PROFILE=slurm/     # or profile/ for a single machine
```

## Regenerating

```bash
python scripts/generate_api_examples.py
python tests/test_parity.py
```
