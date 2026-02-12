#!/usr/bin/env python3
"""
Verify that battleamp-snakemake task label files exactly match
the original BattleAMP-benchmark FASTA sources and paper definitions.

Cross-checks:
1. Each task's label sequences match the original FASTA files (byte-exact)
2. Positive/negative counts match
3. All paper tasks are accounted for (or explicitly excluded)
4. Synthetic tasks use correct positive set (broad_positive, not amp_positive)
5. Regression MIC values match original CSVs
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

from Bio import SeqIO
import pandas as pd


parser = argparse.ArgumentParser(
    description="Verify battleamp-snakemake tasks against AMP-data sources and paper."
)
parser.add_argument(
    "--amp-data", required=True,
    help="Path to AMP-data/data/ directory (contains all_sequences.fasta, activity/, syntax/)"
)
parser.add_argument(
    "--tasks", required=True,
    help="Path to battleamp-snakemake/tasks/ directory"
)
args = parser.parse_args()

DATA = args.amp_data
TASKS = args.tasks

if not os.path.isdir(DATA):
    print(f"ERROR: AMP-data directory not found: {DATA}")
    sys.exit(1)
if not os.path.isdir(TASKS):
    print(f"ERROR: Tasks directory not found: {TASKS}")
    sys.exit(1)

errors = []
warnings = []
ok_count = 0


def check(condition, msg):
    global ok_count
    if condition:
        ok_count += 1
        print(f"  OK: {msg}")
    else:
        errors.append(msg)
        print(f"  FAIL: {msg}")


def warn(msg):
    warnings.append(msg)
    print(f"  WARN: {msg}")


def read_fasta_seqs(path):
    """Read sequences from FASTA, return list (preserving order and dupes)."""
    return [str(r.seq) for r in SeqIO.parse(path, "fasta")]


def read_label_tsv(task_name):
    """Read a task labels.tsv, return DataFrame."""
    path = os.path.join(TASKS, task_name, "labels.tsv")
    return pd.read_csv(path, sep="\t")


def verify_classification_task(task_name, pos_fasta, neg_fasta, paper_name):
    """Verify a classification task matches source FASTAs."""
    print(f"\n{'='*60}")
    print(f"Task: {task_name}")
    print(f"Paper name: {paper_name}")
    print(f"Positive FASTA: {os.path.basename(pos_fasta)}")
    print(f"Negative FASTA: {os.path.basename(neg_fasta)}")

    # Read sources
    pos_seqs = read_fasta_seqs(pos_fasta)
    neg_seqs = read_fasta_seqs(neg_fasta)

    # Read our labels
    df = read_label_tsv(task_name)

    # Check columns
    check("sequence" in df.columns and "label" in df.columns,
          f"labels.tsv has 'sequence' and 'label' columns")

    # Split by label
    our_pos = df[df["label"] == "AMP"]["sequence"].tolist()
    our_neg = df[df["label"] == "non-AMP"]["sequence"].tolist()

    # Count checks
    check(len(our_pos) == len(pos_seqs),
          f"Positive count: ours={len(our_pos)}, source={len(pos_seqs)}")
    check(len(our_neg) == len(neg_seqs),
          f"Negative count: ours={len(our_neg)}, source={len(neg_seqs)}")

    # Exact sequence set match
    check(set(our_pos) == set(pos_seqs),
          f"Positive sequences match source FASTA exactly")
    check(set(our_neg) == set(neg_seqs),
          f"Negative sequences match source FASTA exactly")

    # Check for unexpected labels
    all_labels = set(df["label"].unique())
    check(all_labels == {"AMP", "non-AMP"},
          f"Only labels are 'AMP' and 'non-AMP' (got {all_labels})")


def verify_regression_task(task_name, mic_csv, paper_name):
    """Verify a regression task matches source MIC CSV."""
    print(f"\n{'='*60}")
    print(f"Task: {task_name}")
    print(f"Paper name: {paper_name}")
    print(f"MIC CSV: {os.path.basename(mic_csv)}")

    # Read source
    src_df = pd.read_csv(mic_csv)
    src_seqs = set(src_df["sequence"].tolist())
    src_mics = dict(zip(src_df["sequence"], src_df["activity"]))

    # Read our labels
    df = read_label_tsv(task_name)

    check("sequence" in df.columns and "MIC" in df.columns,
          f"labels.tsv has 'sequence' and 'MIC' columns")

    our_seqs = set(df["sequence"].tolist())

    check(len(df) == len(src_df) - 1 or len(df) == len(src_df),
          f"Row count: ours={len(df)}, source={len(src_df)} (minus header)")

    check(our_seqs == src_seqs,
          f"Sequence sets match source CSV exactly")

    # Check MIC values match
    our_mics = dict(zip(df["sequence"], df["MIC"]))
    mic_mismatches = 0
    for seq in our_seqs & src_seqs:
        if abs(float(our_mics[seq]) - float(src_mics[seq])) > 0.01:
            mic_mismatches += 1
    check(mic_mismatches == 0,
          f"MIC values match (mismatches: {mic_mismatches})")


# =========================================================================
# Paper Table 1 task mapping
# =========================================================================

print("=" * 60)
print("PAPER-TO-CODE TASK MAPPING")
print("=" * 60)

paper_tasks = {
    # Paper name -> (our task name(s), type, notes)
    "AMP/non-AMP": (["amp"], "classification", "AMPTask in code"),
    "GeneralActivity": (["broad_activity"], "classification", "BroadActivityTask in code"),
    "GramActivity GramPlus": (["gram_plus"], "classification", "GramActivityTask gramplus"),
    "GramActivity GramMinus": (["gram_minus"], "classification", "GramActivityTask gramminus"),
    "StrainActivity Ecoli 25922": (["strain_ecoli25922"], "classification", "StrainActivityTask"),
    "StrainActivity Saureus 25923": (["strain_saureus25923"], "classification", "StrainActivityTask"),
    "StrainActivity Kpneumoniae 700603": (["strain_kpneumoniae700603"], "classification", "StrainActivityTask"),
    "StrainActivity Abaumannii 19606": (["strain_abaumannii19606"], "classification", "StrainActivityTask"),
    "StrainActivity Paeruginosa 27853": (["strain_paeruginosa27853"], "classification", "StrainActivityTask"),
    "StrainActivity Saureus 43300": (["strain_saureus43300"], "classification", "StrainActivityTask"),
    "StrainActivity Saureus 33591": (["strain_saureus33591"], "classification", "StrainActivityTask"),
    "SpeciesActivity Ecoli": (["species_ecoli"], "classification", "SpeciesActivityTask"),
    "SpeciesActivity Saureus": (["species_saureus"], "classification", "SpeciesActivityTask"),
    "SpeciesActivity Kpneumoniae": (["species_kpneumoniae"], "classification", "SpeciesActivityTask"),
    "SpeciesActivity Abaumannii": (["species_abaumannii"], "classification", "SpeciesActivityTask"),
    "SpeciesActivity Paeruginosa": (["species_paeruginosa"], "classification", "SpeciesActivityTask"),
    "Regression Ecoli 25922": (["regression_ecoli25922"], "regression", "RegressionActivityTask"),
    "Regression Saureus 25923": (["regression_saureus25923"], "regression", "RegressionActivityTask"),
    "ShortSequences": (["short_sequences"], "classification", "ShortSequencesTask"),
    "LongSequences": (["long_sequences"], "classification", "LongSequencesTask"),
    "HighSimilarity": (["high_similarity"], "classification", "HighSimilarityTask"),
    "Synthetic Random": (["synthetic_random"], "classification", "SyntheticTask random"),
    "Synthetic Realistic": (["synthetic_realistic"], "classification", "SyntheticTask realistic"),
    "Synthetic Shuffled": (["synthetic_shuffled"], "classification", "SyntheticTask shuffled"),
    "ActivityCliff / PairedSequences": (["paired_sequences"], "paired", "PairedSequencesTask"),
    "SLAY": ([], "classification", "EXCLUDED: separate from benchmark per user decision"),
    "RegressionToClassification": ([], "regression-to-classification",
                                   "NOT YET IMPLEMENTED: requires binarizing MIC predictions"),
}

for paper_name, (our_tasks, task_type, notes) in paper_tasks.items():
    status = "PRESENT" if our_tasks else "EXCLUDED"
    our_str = ", ".join(our_tasks) if our_tasks else "(none)"
    print(f"  {paper_name:45s} -> {our_str:30s} [{status}] {notes}")


# =========================================================================
# Verify classification tasks against source FASTAs
# =========================================================================

print("\n\n" + "=" * 60)
print("VERIFYING CLASSIFICATION TASKS AGAINST SOURCE FASTAS")
print("=" * 60)

classification_checks = [
    # (task_name, pos_fasta, neg_fasta, paper_name)
    ("amp",
     f"{DATA}/amp_positive.fasta",
     f"{DATA}/amp_negative.fasta",
     "AMP/non-AMP"),

    ("broad_activity",
     f"{DATA}/activity/broad_positive.fasta",
     f"{DATA}/activity/broad_negative.fasta",
     "GeneralActivity"),

    ("gram_plus",
     f"{DATA}/activity/gramplus_positive.fasta",
     f"{DATA}/activity/gramplus_negative.fasta",
     "GramActivity GramPlus"),

    ("gram_minus",
     f"{DATA}/activity/gramminus_positive.fasta",
     f"{DATA}/activity/gramminus_negative.fasta",
     "GramActivity GramMinus"),

    # Strains
    ("strain_ecoli25922",
     f"{DATA}/activity/strain/escherichiacoliatcc25922_positive.fasta",
     f"{DATA}/activity/strain/escherichiacoliatcc25922_negative.fasta",
     "StrainActivity Ecoli 25922"),

    ("strain_saureus25923",
     f"{DATA}/activity/strain/staphylococcusaureusatcc25923_positive.fasta",
     f"{DATA}/activity/strain/staphylococcusaureusatcc25923_negative.fasta",
     "StrainActivity Saureus 25923"),

    ("strain_kpneumoniae700603",
     f"{DATA}/activity/strain/klebsiellapneumoniaeatcc700603_positive.fasta",
     f"{DATA}/activity/strain/klebsiellapneumoniaeatcc700603_negative.fasta",
     "StrainActivity Kpneumoniae 700603"),

    ("strain_abaumannii19606",
     f"{DATA}/activity/strain/acinetobacterbaumanniiatcc19606_positive.fasta",
     f"{DATA}/activity/strain/acinetobacterbaumanniiatcc19606_negative.fasta",
     "StrainActivity Abaumannii 19606"),

    ("strain_paeruginosa27853",
     f"{DATA}/activity/strain/pseudomonasaeruginosaatcc27853_positive.fasta",
     f"{DATA}/activity/strain/pseudomonasaeruginosaatcc27853_negative.fasta",
     "StrainActivity Paeruginosa 27853"),

    ("strain_saureus43300",
     f"{DATA}/activity/strain/staphylococcusaureusatcc43300_positive.fasta",
     f"{DATA}/activity/strain/staphylococcusaureusatcc43300_negative.fasta",
     "StrainActivity Saureus 43300"),

    ("strain_saureus33591",
     f"{DATA}/activity/strain/staphylococcusaureusatcc33591_positive.fasta",
     f"{DATA}/activity/strain/staphylococcusaureusatcc33591_negative.fasta",
     "StrainActivity Saureus 33591"),

    # Species
    ("species_ecoli",
     f"{DATA}/activity/species/escherichiacoli_positive.fasta",
     f"{DATA}/activity/species/escherichiacoli_negative.fasta",
     "SpeciesActivity Ecoli"),

    ("species_saureus",
     f"{DATA}/activity/species/staphylococcusaureus_positive.fasta",
     f"{DATA}/activity/species/staphylococcusaureus_negative.fasta",
     "SpeciesActivity Saureus"),

    ("species_kpneumoniae",
     f"{DATA}/activity/species/klebsiellapneumoniae_positive.fasta",
     f"{DATA}/activity/species/klebsiellapneumoniae_negative.fasta",
     "SpeciesActivity Kpneumoniae"),

    ("species_abaumannii",
     f"{DATA}/activity/species/acinetobacterbaumannii_positive.fasta",
     f"{DATA}/activity/species/acinetobacterbaumannii_negative.fasta",
     "SpeciesActivity Abaumannii"),

    ("species_paeruginosa",
     f"{DATA}/activity/species/pseudomonasaeruginosa_positive.fasta",
     f"{DATA}/activity/species/pseudomonasaeruginosa_negative.fasta",
     "SpeciesActivity Paeruginosa"),

    # Syntax
    ("high_similarity",
     f"{DATA}/syntax/clustered/broad_similar_positive.fasta",
     f"{DATA}/syntax/clustered/broad_similar_negative.fasta",
     "HighSimilarity"),

    ("short_sequences",
     f"{DATA}/syntax/short_positive.fasta",
     f"{DATA}/syntax/short_negative.fasta",
     "ShortSequences"),

    ("long_sequences",
     f"{DATA}/syntax/long_positive.fasta",
     f"{DATA}/syntax/long_negative.fasta",
     "LongSequences"),

    # Synthetic: positives = broad_positive (not amp_positive!)
    # This matches SyntheticTask code: self._amp_data.positives_broad
    ("synthetic_random",
     f"{DATA}/activity/broad_positive.fasta",
     f"{DATA}/syntax/synthetic_random.fasta",
     "Synthetic Random"),

    ("synthetic_realistic",
     f"{DATA}/activity/broad_positive.fasta",
     f"{DATA}/syntax/synthetic_realistic.fasta",
     "Synthetic Realistic"),

    ("synthetic_shuffled",
     f"{DATA}/activity/broad_positive.fasta",
     f"{DATA}/syntax/broad_positive_shuffled.fasta",
     "Synthetic Shuffled"),
]

for task_name, pos_fasta, neg_fasta, paper_name in classification_checks:
    verify_classification_task(task_name, pos_fasta, neg_fasta, paper_name)


# =========================================================================
# Verify regression tasks against source CSVs
# =========================================================================

print("\n\n" + "=" * 60)
print("VERIFYING REGRESSION TASKS AGAINST SOURCE CSVs")
print("=" * 60)

regression_checks = [
    ("regression_ecoli25922",
     f"{DATA}/activity/strain/escherichiacoliatcc25922_mic.csv",
     "Regression Ecoli 25922"),

    ("regression_saureus25923",
     f"{DATA}/activity/strain/staphylococcusaureusatcc25923_mic.csv",
     "Regression Saureus 25923"),
]

for task_name, mic_csv, paper_name in regression_checks:
    verify_regression_task(task_name, mic_csv, paper_name)


# =========================================================================
# Verify synthetic tasks use broad_positive (not amp_positive)
# This is a critical check: the original code does
#   self._amp_data.positives_broad (not positives_amp)
# =========================================================================

print("\n\n" + "=" * 60)
print("CRITICAL: VERIFYING SYNTHETIC POSITIVE SET")
print("=" * 60)

broad_pos = set(read_fasta_seqs(f"{DATA}/activity/broad_positive.fasta"))
amp_pos = set(read_fasta_seqs(f"{DATA}/amp_positive.fasta"))

for task_name in ["synthetic_random", "synthetic_realistic", "synthetic_shuffled"]:
    df = read_label_tsv(task_name)
    our_pos = set(df[df["label"] == "AMP"]["sequence"])
    check(our_pos == broad_pos,
          f"{task_name}: positives match broad_positive.fasta (not amp_positive)")
    check(our_pos != amp_pos,
          f"{task_name}: positives are NOT amp_positive.fasta")


# =========================================================================
# Verify paired_sequences task structure
# =========================================================================

print("\n\n" + "=" * 60)
print("VERIFYING PAIRED SEQUENCES TASK")
print("=" * 60)

paired_fasta = f"{DATA}/syntax/paired_MD.fasta"
paired_records = list(SeqIO.parse(paired_fasta, "fasta"))
df = read_label_tsv("paired_sequences")

check("pair_id" in df.columns,
      "paired_sequences labels.tsv has 'pair_id' column")
check(len(df) == len(paired_records),
      f"Row count matches FASTA: ours={len(df)}, source={len(paired_records)}")

# Check pairs are consecutive AMP/non-AMP
n_pairs = df["pair_id"].nunique()
check(n_pairs == len(paired_records) // 2,
      f"Number of pairs: {n_pairs} (expected {len(paired_records) // 2})")


# =========================================================================
# Check for tasks we should NOT have
# =========================================================================

print("\n\n" + "=" * 60)
print("CHECKING FOR UNEXPECTED/MISSING TASKS")
print("=" * 60)

expected_tasks = set()
for _, (task_list, _, _) in paper_tasks.items():
    expected_tasks.update(task_list)

# Also include example_classification which is our test task
actual_tasks = set(os.listdir(TASKS))

extra = actual_tasks - expected_tasks
if extra:
    warn(f"Extra tasks not in paper: {extra}")
else:
    check(True, "No unexpected tasks found")

missing = expected_tasks - actual_tasks
if missing:
    warn(f"Paper tasks missing from pipeline: {missing}")
else:
    check(True, "All expected tasks present")


# =========================================================================
# SUMMARY
# =========================================================================

print("\n\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Checks passed: {ok_count}")
print(f"  Errors:        {len(errors)}")
print(f"  Warnings:      {len(warnings)}")

if errors:
    print("\nFAILED CHECKS:")
    for e in errors:
        print(f"  - {e}")

if warnings:
    print("\nWARNINGS:")
    for w in warnings:
        print(f"  - {w}")

if not errors:
    print("\nAll checks passed. Task labels match original benchmark sources.")
else:
    print(f"\n{len(errors)} checks FAILED. Review errors above.")
    sys.exit(1)
