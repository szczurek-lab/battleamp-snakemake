#!/bin/bash
# ampeppy inference wrapper for battleamp-snakemake
#
# Interface contract:
#   $1 = path to input FASTA file (absolute)
#   $2 = path to output TSV file (absolute)
#
# Output columns: sequence  Prediction  Probability_score
#
# The conda env is already active (managed by Snakemake).
# Do NOT run conda activate here.

set -euo pipefail

INPUT_FASTA="$1"
OUTPUT_TSV="$2"

if [ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_TSV" ]; then
    echo "Usage: inference.sh <input.fasta> <output.tsv>" >&2
    exit 1
fi

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: input FASTA not found: $INPUT_FASTA" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_FILE="$SCRIPT_DIR/pretrained_models/amPEP.model"
TMP_OUTPUT="${OUTPUT_TSV}.raw"

# Step 1: Run ampep predict
echo "Running ampep predict..." >&2
ampep predict -m "$MODEL_FILE" -i "$INPUT_FASTA" -o "$TMP_OUTPUT" --seed 2012

# Step 2: Convert to standard pipeline format
# Raw ampep output columns: probability_nonAMP, probability_AMP, predicted, seq_id, seq
# Pipeline standard columns: sequence, Prediction, Probability_score
python3 - "$TMP_OUTPUT" "$OUTPUT_TSV" << 'PYEOF'
import sys
import pandas as pd

raw_path = sys.argv[1]
out_path = sys.argv[2]

df = pd.read_csv(raw_path, sep="\t")

out = pd.DataFrame({
    "sequence": df["seq"],
    "Prediction": df["predicted"].replace({"nonAMP": "non-AMP"}),
    "Probability_score": df["probability_AMP"],
})

out.to_csv(out_path, sep="\t", index=False)
print(f"Converted {len(out)} predictions to {out_path}", file=sys.stderr)
PYEOF

# Step 3: Clean up
rm -f "$TMP_OUTPUT"

echo "ampeppy inference complete: $(wc -l < "$OUTPUT_TSV") lines" >&2
