#!/usr/bin/env python3
"""Generate homology-split classification tasks using CD-HIT clustering.

For each identity threshold the script clusters the sequences from a base
task with CD-HIT and keeps only the cluster representatives.  The result is
a reduced label file where no two sequences share more than the given
identity, producing progressively harder evaluation subsets as the threshold
decreases.

Prerequisites
-------------
    cd-hit must be installed and available on $PATH.
    Install via conda (``conda install -c bioconda cd-hit``) or your system
    package manager.

Usage
-----
    python scripts/generate_homology_tasks.py

    # Custom thresholds and base task
    python scripts/generate_homology_tasks.py \
        --base-labels tasks/broad_activity/labels.tsv \
        --thresholds 0.4 0.6 0.8 \
        --output-dir tasks

The script writes one ``labels.tsv`` per threshold into
``<output-dir>/homology_<pct>/labels.tsv``.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter


DEFAULT_BASE = "tasks/broad_activity/labels.tsv"
DEFAULT_THRESHOLDS = [0.4, 0.6, 0.8]
DEFAULT_OUTPUT_DIR = "tasks"


def read_labels(path: str) -> list[tuple[str, str]]:
    rows = []
    with open(path, newline="") as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                rows.append((parts[0].strip(), parts[1].strip()))
    return rows


def write_labels(rows: list[tuple[str, str]], path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write("sequence\tlabel\n")
        for seq, label in rows:
            fh.write(f"{seq}\t{label}\n")


def parse_fasta(path: str) -> dict[str, str]:
    """Return {header: sequence} from a FASTA file."""
    seqs = {}
    header = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if header is not None:
        seqs[header] = "".join(buf)
    return seqs


def word_size_for_threshold(threshold: float) -> int:
    """Pick the CD-HIT word size appropriate for a given identity threshold."""
    if threshold >= 0.7:
        return 5
    if threshold >= 0.6:
        return 4
    if threshold >= 0.5:
        return 3
    return 2


def run_cdhit(
    fasta_in: str, fasta_out: str, threshold: float, threads: int = 4
) -> None:
    word_size = word_size_for_threshold(threshold)
    cmd = [
        "cd-hit",
        "-i", fasta_in,
        "-o", fasta_out,
        "-c", str(threshold),
        "-n", str(word_size),
        "-M", "0",
        "-T", str(threads),
        "-d", "0",
    ]
    print(f"    {' '.join(cmd)}")
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate homology-split classification tasks via CD-HIT."
    )
    parser.add_argument(
        "--base-labels",
        default=DEFAULT_BASE,
        help="Path to the base task labels.tsv (default: %(default)s)",
    )
    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=DEFAULT_THRESHOLDS,
        help="Sequence identity thresholds for CD-HIT (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Root directory for output task folders (default: %(default)s)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of CD-HIT threads (default: %(default)s)",
    )
    args = parser.parse_args()

    if not shutil.which("cd-hit"):
        sys.exit(
            "Error: cd-hit not found on $PATH. "
            "Install with: conda install -c bioconda cd-hit"
        )

    if not os.path.isfile(args.base_labels):
        sys.exit(
            f"Error: {args.base_labels} not found.\n"
            f"Run this script from the repository root, or pass "
            f"--base-labels with the correct path."
        )

    labels = read_labels(args.base_labels)
    print(f"Loaded {len(labels)} sequences from {args.base_labels}\n")

    tmpdir = tempfile.mkdtemp(prefix="battleamp_homology_")
    try:
        # Write sequences as FASTA for CD-HIT
        input_fasta = os.path.join(tmpdir, "input.fasta")
        idx_to_seq = {}
        with open(input_fasta, "w") as fh:
            for i, (seq, _) in enumerate(labels):
                tag = f"seq_{i}"
                idx_to_seq[tag] = seq
                fh.write(f">{tag}\n{seq}\n")

        label_map = {seq: lab for seq, lab in labels}

        for threshold in sorted(args.thresholds, reverse=True):
            pct = int(threshold * 100)
            task_name = f"homology_{pct}"
            out_fasta = os.path.join(tmpdir, f"cdhit_{pct}")

            print(f"  {task_name} (identity threshold {threshold}):")
            run_cdhit(input_fasta, out_fasta, threshold, args.threads)

            # Representatives are the sequences in the CD-HIT output FASTA
            rep_seqs = set(parse_fasta(out_fasta).values())

            subset = [(seq, lab) for seq, lab in labels if seq in rep_seqs]
            counts = Counter(lab for _, lab in subset)
            outpath = os.path.join(args.output_dir, task_name, "labels.tsv")
            write_labels(subset, outpath)
            print(
                f"    {len(subset)} representative sequences "
                f"(AMP={counts.get('AMP', 0)}, non-AMP={counts.get('non-AMP', 0)}) "
                f"-> {outpath}\n"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
