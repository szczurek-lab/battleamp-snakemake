#!/usr/bin/env python3
"""
Clean a FASTA dataset before model-specific prefiltering.

Runs once per dataset. Filters out sequences that would cause errors or
ambiguity in downstream models:

  1. Sequences containing lowercase letters (d-amino acids)
  2. Sequences with non-standard amino acids (anything outside ACDEFGHIKLMNPQRSTVWY)
  3. Empty sequences
  4. Sequences with whitespace or gap characters
  5. Duplicate sequences (keeps first occurrence)

Writes a cleaned FASTA and a JSON report summarizing what was removed.
"""

import argparse
import json
import sys
from collections import OrderedDict
from pathlib import Path

from Bio import SeqIO

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_dataset(input_fasta, output_fasta, report_path):
    seen_seqs = set()
    kept = []

    counts = {
        "total": 0,
        "kept": 0,
        "dropped_empty": 0,
        "dropped_lowercase": 0,
        "dropped_nonstandard_aa": 0,
        "dropped_whitespace": 0,
        "dropped_duplicate": 0,
    }

    # Track a few example IDs per drop reason for the report
    examples = {
        "dropped_lowercase": [],
        "dropped_nonstandard_aa": [],
        "dropped_whitespace": [],
        "dropped_duplicate": [],
    }
    max_examples = 5

    for record in SeqIO.parse(input_fasta, "fasta"):
        counts["total"] += 1
        seq_str = str(record.seq)

        # 1. Empty
        if len(seq_str) == 0:
            counts["dropped_empty"] += 1
            continue

        # 2. Whitespace or gap characters (-, .)
        if any(c in seq_str for c in (" ", "\t", "-", ".")):
            counts["dropped_whitespace"] += 1
            if len(examples["dropped_whitespace"]) < max_examples:
                examples["dropped_whitespace"].append(record.id)
            continue

        # 3. Lowercase letters (d-amino acids)
        if any(c.islower() for c in seq_str):
            counts["dropped_lowercase"] += 1
            if len(examples["dropped_lowercase"]) < max_examples:
                examples["dropped_lowercase"].append(record.id)
            continue

        # 4. Non-standard amino acids (sequence is already all-uppercase here)
        if not all(c in STANDARD_AA for c in seq_str):
            bad_chars = sorted(set(seq_str) - STANDARD_AA)
            counts["dropped_nonstandard_aa"] += 1
            if len(examples["dropped_nonstandard_aa"]) < max_examples:
                examples["dropped_nonstandard_aa"].append(
                    f"{record.id} (chars: {','.join(bad_chars)})"
                )
            continue

        # 5. Duplicate (by sequence, not by ID)
        if seq_str in seen_seqs:
            counts["dropped_duplicate"] += 1
            if len(examples["dropped_duplicate"]) < max_examples:
                examples["dropped_duplicate"].append(record.id)
            continue

        seen_seqs.add(seq_str)
        kept.append(record)

    counts["kept"] = len(kept)
    counts["dropped_total"] = counts["total"] - counts["kept"]

    # Write cleaned FASTA
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(kept, output_fasta, "fasta")

    # Write report
    report = {
        "input": str(input_fasta),
        "output": str(output_fasta),
        "counts": counts,
        "examples": examples,
    }
    Path(report_path).parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    # Log summary to stderr
    print(f"Dataset cleaning: {input_fasta}", file=sys.stderr)
    print(f"  Total input:        {counts['total']}", file=sys.stderr)
    print(f"  Kept:               {counts['kept']}", file=sys.stderr)
    print(f"  Dropped total:      {counts['dropped_total']}", file=sys.stderr)
    print(f"    empty:            {counts['dropped_empty']}", file=sys.stderr)
    print(f"    lowercase (d-aa): {counts['dropped_lowercase']}", file=sys.stderr)
    print(f"    non-standard aa:  {counts['dropped_nonstandard_aa']}", file=sys.stderr)
    print(f"    whitespace/gaps:  {counts['dropped_whitespace']}", file=sys.stderr)
    print(f"    duplicate:        {counts['dropped_duplicate']}", file=sys.stderr)


def main():
    # When called from Snakemake, read params from snakemake object
    try:
        clean_dataset(
            input_fasta=snakemake.input.sequences,
            output_fasta=snakemake.output.sequences,
            report_path=snakemake.output.report,
        )
    except NameError:
        # Called from command line
        parser = argparse.ArgumentParser(
            description="Clean a FASTA dataset for the BattleAMP pipeline"
        )
        parser.add_argument("--input", required=True, help="Input FASTA")
        parser.add_argument("--output", required=True, help="Output cleaned FASTA")
        parser.add_argument("--report", required=True, help="Output JSON report")
        args = parser.parse_args()
        clean_dataset(args.input, args.output, args.report)


if __name__ == "__main__":
    main()
