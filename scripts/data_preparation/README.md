# Data Preparation Scripts

This directory contains all scripts used to create the BATTLE-AMP benchmark
datasets from raw sources. The pipeline starts from the
[DBAASP](https://dbaasp.org/) database and produces the label files consumed by
the Snakemake evaluation workflow.

> **Provenance.** These scripts were originally developed in the
> [`AMP-data`](https://github.com/szczurek-lab/AMP-data) repository and are
> included here for reproducibility.

## Directory layout

```
scripts/data_preparation/
  config.py                    # Filtering thresholds (MIC, CFU, media, etc.)
  dbaasp_utils.py              # Shared helpers (FASTA I/O, activity parsing)
  moi.json                     # Microbes-of-interest registry
  aa_mass.json                 # Amino-acid molecular masses
  requirements.in              # Python dependencies

  # Numbered pipeline steps (run in order)
  01_dbaasp_download.sh        # Download DBAASP entries as JSON
  02_dbaasp_jsons_to_csv.py    # Convert JSON collection to a single CSV
  03_clean_dbaasp.py           # Sequence and measurement quality filtering
  04_dbaasp_taxonomy.py        # Fetch UniProt taxonomy for target species
  05_get_activity_datasets.py  # Binary classification + MIC regression sets
  06_get_hemolytic_dataset.py  # Hemolytic/non-hemolytic classification set
  07_cluster_cdhit.py          # CD-HIT clustering at 80% identity
  08_activity_cliffs.py        # Extract activity-cliff pairs from clusters
  09_get_slay.py               # Download and process SLAY library from GEO

  # Auxiliary generators
  generate_shuffled.py         # Shuffled negatives (permuted AMPs)
  generate_synthetic.py        # Random and realistic synthetic negatives

  # Task generators
  build_benchmark.py           # Build tasks + unified FASTA from pipeline outputs
  generate_length_tasks.py     # Length-stratified evaluation bins
  generate_homology_tasks.py   # Homology-reduced evaluation subsets
```

## Quick start (from committed intermediates)

The DBAASP intermediates are committed in `data/dbaasp/`, so you can skip the
download and cleaning steps (01-04) and start from step 05:

```bash
# Step 05: build activity datasets
python scripts/data_preparation/05_get_activity_datasets.py \
    data/dbaasp/dbaasp_activity.csv \
    data/dbaasp/dbaasp_sequences.csv \
    data/dbaasp/dbaasp_taxonomy.csv \
    scripts/data_preparation/moi.json \
    data/activity

# Steps 07-08: clustering and activity cliffs
python scripts/data_preparation/07_cluster_cdhit.py \
    data/activity/broad.fasta data/activity/broad.csv \
    data/syntax/clustered broad 5

python scripts/data_preparation/08_activity_cliffs.py \
    data/syntax/clustered/broad_clustered.csv \
    data/activity/broad.csv \
    data/syntax/clustered/broad_activitycliffs.csv

# Step 09: SLAY dataset (downloads from GEO)
python scripts/data_preparation/09_get_slay.py data/slay/

# Synthetic negatives
python scripts/data_preparation/generate_shuffled.py \
    data/activity/broad_positive.fasta \
    data/syntax/broad_positive_shuffled.fasta -n 10

python scripts/data_preparation/generate_synthetic.py \
    data/activity/broad_positive.fasta \
    data/syntax/synthetic_random.fasta -n 10000 --mode random

python scripts/data_preparation/generate_synthetic.py \
    data/activity/broad_positive.fasta \
    data/syntax/synthetic_realistic.fasta -n 10000 --mode realistic

# Final step: build all task labels + unified FASTA
python scripts/data_preparation/build_benchmark.py \
    --amp-data-root . --repo-root .
```

## Environment setup

```bash
pip install pandas numpy requests modlamp pycdhit biopython
```

For the DBAASP download step (01), `curl` is required.

## Full pipeline reference

All commands are run from the repository root. Intermediate outputs go into
`data/`.

### Steps 01-04: DBAASP download and cleaning

These steps require network access and produce the three intermediates committed
in `data/dbaasp/`. They do not need to be rerun unless updating the DBAASP
snapshot.

**Step 1: Download DBAASP**

```bash
bash scripts/data_preparation/01_dbaasp_download.sh \
    --output-path data/dbaasp --end 22878
```

Downloads one JSON per peptide entry from the DBAASP REST API (22,878 entries
as of November 2024).

**Step 2: Convert to CSV**

```bash
python scripts/data_preparation/02_dbaasp_jsons_to_csv.py \
    data/dbaasp/json/ data/dbaasp
```

Produces `data/dbaasp/dbaasp.csv`.

**Step 3: Clean sequences and measurements**

```bash
python scripts/data_preparation/03_clean_dbaasp.py \
    data/dbaasp/dbaasp.csv data/dbaasp/
```

Applies the filters defined in `config.py`:

| Filter | Value |
|--------|-------|
| Standard amino acids only | ACDEFGHIKLMNPRSTWYQV |
| C-terminus | Standard or amidated (AMD) |
| N-terminus | Standard only |
| CFU range | 1E5 to 1E6 |
| Activity measure | MIC only |
| Media | MHB, MHA, CAMHB, M7H10A, M7H11A, M7H9B |
| pH | Unspecified (standard) |
| Salt | None (no added salt) |

Produces `dbaasp_sequences.csv` (filtered peptides) and `dbaasp_activity.csv`
(clean MIC measurements).

**Step 4: Fetch taxonomy**

```bash
python scripts/data_preparation/04_dbaasp_taxonomy.py \
    data/dbaasp/dbaasp_activity.csv data/dbaasp/
```

Queries the UniProt taxonomy API to assign species, genus, family, order, class,
and phylum to each target organism. Produces `dbaasp_taxonomy.csv`.

### Step 5: Build activity datasets

```bash
python scripts/data_preparation/05_get_activity_datasets.py \
    data/dbaasp/dbaasp_activity.csv \
    data/dbaasp/dbaasp_sequences.csv \
    data/dbaasp/dbaasp_taxonomy.csv \
    scripts/data_preparation/moi.json \
    data/activity
```

Creates binary classification datasets (active: MIC <= 32 ug/ml, inactive:
MIC >= 128 ug/ml) and MIC regression datasets for each level defined in
`moi.json`:

- **Broad**: all five species pooled
- **Gram+/Gram-**: grouped by Gram stain
- **Species**: per-species datasets (5 species)
- **Strain**: per-strain datasets (7 strains)
- **MIC**: regression sets for E. coli ATCC 25922, S. aureus ATCC 25923

### Step 6: Hemolytic dataset (optional)

```bash
python scripts/data_preparation/06_get_hemolytic_dataset.py \
    data/dbaasp/dbaasp.csv \
    data/dbaasp/dbaasp_sequences.csv \
    data/toxicity
```

### Step 7: CD-HIT clustering

```bash
python scripts/data_preparation/07_cluster_cdhit.py \
    data/activity/broad.fasta \
    data/activity/broad.csv \
    data/syntax/clustered broad 5
```

Clusters sequences at 80% identity using CD-HIT. Used as input for
activity-cliff extraction (Step 8).

### Step 8: Activity cliffs

```bash
python scripts/data_preparation/08_activity_cliffs.py \
    data/syntax/clustered/broad_clustered.csv \
    data/activity/broad.csv \
    data/syntax/clustered/broad_activitycliffs.csv
```

From each CD-HIT cluster, extracts the most active and least active peptide to
form pairs with high sequence similarity but opposing activity labels.

### Step 9: SLAY dataset

```bash
python scripts/data_preparation/09_get_slay.py data/slay/
```

Downloads the SLAY screening library from GEO (GSE94529). Peptides with
log2 fold change <= -1 are classified as active (7,968 peptides, ~1.7%
of the library).

### Synthetic negatives

```bash
# Shuffled: 10 permutations per active peptide + original
python scripts/data_preparation/generate_shuffled.py \
    data/activity/broad_positive.fasta \
    data/syntax/broad_positive_shuffled.fasta -n 10

# Random: uniform amino-acid sampling
python scripts/data_preparation/generate_synthetic.py \
    data/activity/broad_positive.fasta \
    data/syntax/synthetic_random.fasta -n 10000 --mode random

# Realistic: empirical amino-acid frequency sampling
python scripts/data_preparation/generate_synthetic.py \
    data/activity/broad_positive.fasta \
    data/syntax/synthetic_realistic.fasta -n 10000 --mode realistic
```

**Note:** `generate_shuffled.py` outputs the original sequence plus `n`
shuffled variants per peptide, so `-n 10` produces 11 records per input
sequence.

### Length-stratified and homology-reduced tasks

```bash
# Length-stratified bins (1-10, 11-20, 21-30, 31-50 aa)
python scripts/data_preparation/generate_length_tasks.py

# Homology-reduced subsets (CD-HIT at 80%, 60%, 40%)
python scripts/data_preparation/generate_homology_tasks.py
```

### Build benchmark (final step)

```bash
python scripts/data_preparation/build_benchmark.py \
    --amp-data-root . --repo-root .
```

1. Reads all FASTAs and CSVs produced by the steps above
2. Generates `tasks/*/labels.tsv` for all classification, regression, and
   paired tasks
3. Assembles the de-duplicated `datasets/battleamp-all/sequences.fasta`

Any missing inputs (e.g. SLAY not downloaded, paired_safe not curated) are
skipped with warnings rather than errors.

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--mic-column` | `activity` | CSV column for regression MIC values. `activity` uses the per-row value from the first measurement. `MIC` uses min-aggregated values (more conservative). |
| `--convert-units` | off | Convert uM values to ug/ml using per-peptide MW. Off by default; values are written as-is with the target-unit label regardless of source unit. |
| `--target-unit` | `ug/ml` | Unit label written to regression outputs. |
| `--skip-slay` | off | Skip SLAY task (requires pre-downloaded CSV from GEO). |
| `--skip-fasta` | off | Skip unified FASTA assembly. |
| `--dry-run` | off | Print plan without writing files. |

**Tasks not generated automatically:**

- `paired_safe`: Requires hemolytic toxicity data from step 06; must be curated
  separately.
- `slay`: Requires running `09_get_slay.py` first to download from GEO.
- `datasets/deep-amp/`: External dataset, not derived from the DBAASP pipeline.

**Regression label note:** By default the script uses the `activity` column from
the first measurement per peptide (after deduplication) and labels all values as
ug/ml regardless of the source unit. About 60% of E. coli regression values and
a similar fraction of S. aureus values are originally in uM. Use
`--convert-units` for proper MW-based conversion.

## Configuration reference

### `config.py`

Central filtering parameters for DBAASP cleaning:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_MIC` | 32 | Active threshold (ug/ml) |
| `MAX_MIC` | 128 | Inactive threshold (ug/ml) |
| `MIN_CFU` | 1e5 | Minimum colony-forming units |
| `MAX_CFU` | 1e6 | Maximum colony-forming units |
| `ACT_MEASURE` | `'MIC'` | Activity measure type |
| `MEDIUM` | MHB, MHA, CAMHB, ... | Accepted growth media |

### `moi.json`

Defines the microbes of interest at each taxonomic level (broad, Gram+/-,
species, strain, MIC regression).

## Output structure

After running the full pipeline, the label files in `tasks/*/labels.tsv` and
the unified FASTA in `datasets/battleamp-all/sequences.fasta` are ready for the
Snakemake evaluation workflow.