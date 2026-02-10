"""
Pre-filter sequences by length constraints from model.yaml.

Produces:
  - A filtered FASTA containing only sequences within [length_min, length_max]
  - A manifest TSV listing all sequences and whether they were included or skipped

Usage: called by Snakemake via script: directive
"""

import sys
from pathlib import Path


def parse_fasta(fasta_path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(fasta_path) as f:
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


def write_fasta(fasta_path, records):
    """Write a list of (header, sequence) tuples to a FASTA file."""
    with open(fasta_path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def main(snakemake):
    input_fasta = snakemake.input.fasta
    output_fasta = snakemake.output.fasta
    manifest_path = snakemake.output.manifest

    length_min = snakemake.params.length_min
    length_max = snakemake.params.length_max

    # Convert None / "None" to actual None
    if length_min is not None and str(length_min).lower() == "none":
        length_min = None
    if length_max is not None and str(length_max).lower() == "none":
        length_max = None

    if length_min is not None:
        length_min = int(length_min)
    if length_max is not None:
        length_max = int(length_max)

    kept = []
    manifest_rows = []
    total = 0
    skipped = 0

    for header, seq in parse_fasta(input_fasta):
        total += 1
        seq_len = len(seq)
        include = True

        if length_min is not None and seq_len < length_min:
            include = False
        if length_max is not None and seq_len > length_max:
            include = False

        if include:
            kept.append((header, seq))
            manifest_rows.append((header, seq, seq_len, "included"))
        else:
            skipped += 1
            manifest_rows.append((header, seq, seq_len, "skipped_length"))

    write_fasta(output_fasta, kept)

    with open(manifest_path, "w") as f:
        f.write("header\tsequence\tlength\tstatus\n")
        for header, seq, seq_len, status in manifest_rows:
            f.write(f"{header}\t{seq}\t{seq_len}\t{status}\n")

    log_msg = (
        f"Pre-filter: {total} total, {total - skipped} kept, {skipped} skipped "
        f"(length_min={length_min}, length_max={length_max})"
    )
    print(log_msg, file=sys.stderr)


main(snakemake)
