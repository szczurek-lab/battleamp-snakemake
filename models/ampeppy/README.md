# ampeppy adapter for battleamp-snakemake

This directory should contain the BattleAMP-amPEPpy repository contents
alongside these adapter files.

## Setup

```bash
# From the battleamp-snakemake root directory:

# Option A: Clone the repo directly into models/ampeppy/
git clone https://github.com/szczurek-lab/BattleAMP-amPEPpy.git models/ampeppy

# Then copy the adapter files into the directory (they replace originals):
# model.yaml, environment.yaml, setup.sh, inference.sh are the adapters.

# Option B: Add as a git submodule (for the paper repo)
git submodule add https://github.com/szczurek-lab/BattleAMP-amPEPpy.git models/ampeppy
```

## What the adapter files do

- **model.yaml** -- declares ampeppy as a CPU classifier with no length constraints
- **environment.yaml** -- conda env with python, scikit-learn, numpy, pandas, biopython, tqdm
- **setup.sh** -- runs `pip install -e .` to install the `ampep` CLI (no conda env creation)
- **inference.sh** -- calls `ampep predict`, converts raw output to pipeline format
  (lowercase `sequence` column, `non-AMP` with hyphen)

## Expected directory structure after setup

```
models/ampeppy/
    model.yaml             # adapter (this repo)
    environment.yaml       # adapter (this repo)
    setup.sh               # adapter (replaces original)
    inference.sh           # adapter (replaces original)
    amPEPpy/               # from BattleAMP-amPEPpy
    pretrained_models/     # from BattleAMP-amPEPpy (contains amPEP.model, 17MB)
    inference.py           # from BattleAMP-amPEPpy (not used by adapter)
    benchmark_utils.py     # from BattleAMP-amPEPpy (not used by adapter)
    setup.py               # from BattleAMP-amPEPpy (used by setup.sh)
    ...
```

## Notes

- The original setup.sh installed cuml-cu11 (a GPU library). This is unnecessary
  for a random forest model and has been removed.
- The pretrained model file (pretrained_models/amPEP.model) is included in the
  repo, so no weight download is needed.
- ampeppy has no sequence length constraints.

## Validation reference data

The `validation/` directory contains a minimal reference (3 sequences from
`sample.fasta`). You should generate a larger reference set using the
benchmark data:

```bash
# From the battleamp-snakemake root, using cached predictions:
python workflow/scripts/generate_reference.py \
    --predictions models/ampeppy/eval_results/eval_data/dbaasp/activity/all_32_128.tsv \
    --fasta /path/to/benchmark_sequences.fasta \
    --seq-id-col seq_id \
    --model-dir models/ampeppy \
    --model-type classifier \
    --n-samples 200 \
    --seed 42
```

Then verify with: `snakemake validate --use-conda --cores 1`
