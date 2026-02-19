#!/usr/bin/env python3
"""Generate length-binned classification tasks from an existing task's labels.

Each output task contains only the sequences whose length falls within a given
bin.  The bins are defined as command-line arguments or fall back to the
defaults used in the BattleAMP benchmark.

Usage
-----
    python scripts/generate_length_tasks.py

    # Custom bins and base task
    python scripts/generate_length_tasks.py \
        --base-labels tasks/broad_activity/labels.tsv \
        --bins 1-10 11-20 21-30 31-50 \
        --output-dir tasks

The script writes one ``labels.tsv`` per bin into
``<output-dir>/length_<lo>_<hi>/labels.tsv``.
"""

import argparse
import csv
import os
import sys
from collections import Counter


DEFAULT_BASE = "tasks/broad_activity/labels.tsv"
DEFAULT_BINS = ["1-10", "11-20", "21-30", "31-50"]
DEFAULT_OUTPUT_DIR = "tasks"


def parse_bin(spec: str) -> tuple[int, int]:
    lo, hi = spec.split("-", 1)
    return int(lo), int(hi)


def read_labels(path: str) -> list[tuple[str, str]]:
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for row in reader:
            seq = row[0].strip()
            label = row[1].strip()
            rows.append((seq, label))
    return rows


def write_labels(rows: list[tuple[str, str]], path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("sequence\tlabel\n")
        for seq, label in rows:
            fh.write(f"{seq}\t{label}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate length-binned classification tasks."
    )
    parser.add_argument(
        "--base-labels",
        default=DEFAULT_BASE,
        help="Path to the base task labels.tsv (default: %(default)s)",
    )
    parser.add_argument(
        "--bins",
        nargs="+",
        default=DEFAULT_BINS,
        help="Length bins as LO-HI pairs (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Root directory for output task folders (default: %(default)s)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.base_labels):
        sys.exit(
            f"Error: {args.base_labels} not found.\n"
            f"Run this script from the repository root, or pass "
            f"--base-labels with the correct path."
        )

    bins = [parse_bin(b) for b in args.bins]
    labels = read_labels(args.base_labels)
    print(f"Loaded {len(labels)} sequences from {args.base_labels}\n")

    for lo, hi in bins:
        task_name = f"length_{lo:02d}_{hi:02d}"
        subset = [(seq, lab) for seq, lab in labels if lo <= len(seq) <= hi]
        counts = Counter(lab for _, lab in subset)
        outpath = os.path.join(args.output_dir, task_name, "labels.tsv")
        write_labels(subset, outpath)
        print(
            f"  {task_name}: {len(subset)} sequences "
            f"(AMP={counts.get('AMP', 0)}, non-AMP={counts.get('non-AMP', 0)}) "
            f"-> {outpath}"
        )


if __name__ == "__main__":
    main()
