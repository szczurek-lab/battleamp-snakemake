#!/usr/bin/env python3
"""
Generate reference validation data for a model.

This script creates a validation/ directory inside a model's directory
containing:
  - reference_input.fasta: subset of sequences for validation
  - reference_output.tsv: the model's predictions on those sequences

Run this BEFORE making any code changes to capture the original model's
behavior. Then after adapting the model, the pipeline's `validate` target
compares new outputs against these references.

Usage:
    # From an existing prediction TSV that has a sequence column:
    python generate_reference.py \
        --predictions /path/to/existing_predictions.tsv \
        --fasta /path/to/original_input.fasta \
        --model-dir /path/to/models/ampeppy \
        --model-type classifier \
        --n-samples 200 \
        --seed 42

    # Or: run inference fresh using the model's inference.sh
    python generate_reference.py \
        --run-inference \
        --fasta /path/to/all_sequences.fasta \
        --model-dir /path/to/models/ampeppy \
        --model-type classifier \
        --n-samples 200 \
        --seed 42

The script samples sequences stratified by prediction class (for classifiers)
or uniformly (for regressors) to ensure both positive and negative cases are
represented.
"""

import argparse
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


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
    """Write (header, sequence) tuples to FASTA."""
    with open(fasta_path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def sample_sequences_from_fasta(fasta_path, n_samples, seed):
    """Randomly sample n sequences from a FASTA file."""
    random.seed(seed)
    records = list(parse_fasta(fasta_path))
    if len(records) <= n_samples:
        return records
    return random.sample(records, n_samples)


def sample_from_predictions(pred_df, seq_col, model_type, n_samples, seed):
    """Sample sequences from predictions, stratified by class for classifiers."""
    random.seed(seed)

    if model_type == "classifier" and "Prediction" in pred_df.columns:
        # Stratified: half from each class
        n_per_class = n_samples // 2
        positives = pred_df[pred_df["Prediction"] == "AMP"]
        negatives = pred_df[pred_df["Prediction"] != "AMP"]

        n_pos = min(len(positives), n_per_class)
        n_neg = min(len(negatives), n_samples - n_pos)

        sampled = pd.concat([
            positives.sample(n=n_pos, random_state=seed),
            negatives.sample(n=n_neg, random_state=seed),
        ])
    else:
        n = min(len(pred_df), n_samples)
        sampled = pred_df.sample(n=n, random_state=seed)

    return sampled


def main():
    parser = argparse.ArgumentParser(
        description="Generate reference validation data for a model.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--model-dir",
        required=True,
        help="Path to the model directory (e.g., models/ampeppy)",
    )
    parser.add_argument(
        "--model-type",
        required=True,
        choices=["classifier", "regressor"],
        help="Model type",
    )
    parser.add_argument(
        "--predictions",
        help="Path to existing prediction TSV (must have sequence column). "
        "If provided, reference data is sampled from this file.",
    )
    parser.add_argument(
        "--fasta",
        help="Path to input FASTA file. Required if --run-inference is set, "
        "or to extract sequences when --predictions lacks them.",
    )
    parser.add_argument(
        "--seq-col",
        default="sequence",
        help="Name of the sequence column in TSV files (default: sequence)",
    )
    parser.add_argument(
        "--seq-id-col",
        default=None,
        help="Name of the sequence ID column (for joining FASTA headers to "
        "predictions when predictions lack sequences)",
    )
    parser.add_argument(
        "--run-inference",
        action="store_true",
        help="Run inference.sh on --fasta to generate predictions. "
        "The model's conda env must be active.",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=200,
        help="Number of sequences to include in reference (default: 200)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for sampling (default: 42)",
    )
    args = parser.parse_args()

    model_dir = Path(args.model_dir)
    validation_dir = model_dir / "validation"
    validation_dir.mkdir(exist_ok=True)

    ref_input_path = validation_dir / "reference_input.fasta"
    ref_output_path = validation_dir / "reference_output.tsv"

    # -------------------------------------------------------------------
    # Mode 1: Generate from existing predictions TSV
    # -------------------------------------------------------------------
    if args.predictions and not args.run_inference:
        print(f"Loading predictions from {args.predictions}")
        pred_df = pd.read_csv(args.predictions, sep="\t")

        # Check if predictions have sequences or just IDs
        has_sequences = args.seq_col in pred_df.columns
        has_ids = args.seq_id_col and args.seq_id_col in pred_df.columns

        if not has_sequences and has_ids and args.fasta:
            # Join sequences from FASTA via ID
            print(f"Predictions lack '{args.seq_col}' column, joining from FASTA via '{args.seq_id_col}'")
            fasta_records = {h: s for h, s in parse_fasta(args.fasta)}
            pred_df[args.seq_col] = pred_df[args.seq_id_col].map(fasta_records)
            pred_df = pred_df.dropna(subset=[args.seq_col])
            print(f"  Matched {len(pred_df)} sequences")
        elif not has_sequences:
            print(
                f"ERROR: predictions lack '{args.seq_col}' column. "
                f"Provide --fasta and --seq-id-col to join sequences.",
                file=sys.stderr,
            )
            sys.exit(1)

        # Sample
        sampled = sample_from_predictions(
            pred_df, args.seq_col, args.model_type, args.n_samples, args.seed
        )
        print(f"Sampled {len(sampled)} sequences for reference")

        # Write reference FASTA
        fasta_records = []
        for i, (_, row) in enumerate(sampled.iterrows()):
            seq = row[args.seq_col]
            seq_id = row.get(args.seq_id_col or "seq_id", f"ref_{i}")
            fasta_records.append((str(seq_id), seq))
        write_fasta(ref_input_path, fasta_records)

        # Write reference output TSV (in pipeline standard format)
        if args.model_type == "classifier":
            ref_out = sampled[[args.seq_col, "Prediction", "Probability_score"]].copy()
            ref_out = ref_out.rename(columns={args.seq_col: "sequence"})
        elif args.model_type == "regressor":
            cols = [args.seq_col, "MIC"]
            if "MIC_unit" in sampled.columns:
                cols.append("MIC_unit")
            ref_out = sampled[cols].copy()
            ref_out = ref_out.rename(columns={args.seq_col: "sequence"})

        ref_out.to_csv(ref_output_path, sep="\t", index=False)

    # -------------------------------------------------------------------
    # Mode 2: Run inference fresh
    # -------------------------------------------------------------------
    elif args.run_inference:
        if not args.fasta:
            print("ERROR: --fasta is required with --run-inference", file=sys.stderr)
            sys.exit(1)

        # Sample sequences from input FASTA
        sampled_records = sample_sequences_from_fasta(
            args.fasta, args.n_samples, args.seed
        )
        print(f"Sampled {len(sampled_records)} sequences from {args.fasta}")

        # Write sampled FASTA
        write_fasta(ref_input_path, sampled_records)

        # Run inference
        inference_sh = model_dir / "inference.sh"
        if not inference_sh.exists():
            print(f"ERROR: {inference_sh} not found", file=sys.stderr)
            sys.exit(1)

        print(f"Running inference.sh on {len(sampled_records)} sequences...")
        result = subprocess.run(
            ["bash", str(inference_sh), str(ref_input_path), str(ref_output_path)],
            cwd=str(model_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"ERROR: inference.sh failed:\n{result.stderr}", file=sys.stderr)
            sys.exit(1)

        print(f"Inference complete: {result.stderr.strip()}")

    else:
        print(
            "ERROR: provide either --predictions or --run-inference",
            file=sys.stderr,
        )
        sys.exit(1)

    # Verify outputs
    if ref_input_path.exists() and ref_output_path.exists():
        ref_df = pd.read_csv(ref_output_path, sep="\t")
        n_fasta = sum(1 for _ in parse_fasta(ref_input_path))
        print(f"\nReference data generated:")
        print(f"  {ref_input_path}: {n_fasta} sequences")
        print(f"  {ref_output_path}: {len(ref_df)} predictions")
        print(f"  Columns: {list(ref_df.columns)}")
    else:
        print("ERROR: reference files not created", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
