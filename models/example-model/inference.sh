#!/bin/bash
# Example model: random baseline classifier
# Reads a FASTA file and produces random AMP/non-AMP predictions.
#
# Arguments:
#   $1 = input FASTA path
#   $2 = output TSV path
#   $3+ = extra args (ignored for this model)

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

python3 - "$INPUT_FASTA" "$OUTPUT_TSV" << 'PYEOF'
import sys
import random

random.seed(42)

def parse_fasta(path):
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            elif line:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)

input_fasta = sys.argv[1]
output_tsv = sys.argv[2]

with open(output_tsv, "w") as out:
    out.write("sequence\tPrediction\tProbability_score\n")
    for header, seq in parse_fasta(input_fasta):
        prob = random.random()
        pred = "AMP" if prob > 0.5 else "non-AMP"
        out.write(f"{seq}\t{pred}\t{prob:.4f}\n")

print(f"Example model: wrote predictions to {output_tsv}", file=sys.stderr)
PYEOF
