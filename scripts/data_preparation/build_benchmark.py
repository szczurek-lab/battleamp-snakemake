#!/usr/bin/env python3
"""
build_benchmark.py
==================
Convert AMP-data pipeline outputs into the task labels and unified FASTA
consumed by the BattleAMP Snakemake workflow.

Reads from:
    AMP_DATA_ROOT/data/          (FASTAs, CSVs produced by 01-09 scripts)

Writes to:
    REPO_ROOT/tasks/*/labels.tsv
    REPO_ROOT/datasets/battleamp-all/sequences.fasta

Usage:
    python scripts/data_preparation/build_benchmark.py \\
        --amp-data-root  path/to/AMP-data   \\
        --repo-root      path/to/battleamp-snakemake

    # Dry-run (print what would be created, write nothing):
    python scripts/data_preparation/build_benchmark.py \\
        --amp-data-root path/to/AMP-data --repo-root . --dry-run
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict
from pathlib import Path

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_fasta_sequences(fasta_path):
    """Return list of (header, sequence) tuples from a FASTA file."""
    records = []
    header, seq_parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def fasta_to_seqs(fasta_path):
    """Return plain sequence list from a FASTA."""
    return [seq for _, seq in read_fasta_sequences(fasta_path)]


def write_classification_labels(out_path, pos_seqs, neg_seqs,
                                pos_label="AMP", neg_label="non-AMP"):
    """Write a two-column labels.tsv (sequence, label)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sequence", "label"])
        for s in pos_seqs:
            writer.writerow([s, pos_label])
        for s in neg_seqs:
            writer.writerow([s, neg_label])
    return len(pos_seqs) + len(neg_seqs)


def write_regression_labels(out_path, rows):
    """Write a three-column labels.tsv (sequence, MIC, MIC_unit)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sequence", "MIC", "MIC_unit"])
        for seq, mic, unit in rows:
            writer.writerow([seq, mic, unit])
    return len(rows)


def write_paired_labels(out_path, pairs):
    """Write a three-column labels.tsv (sequence, label, pair_id)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sequence", "label", "pair_id"])
        for seq, label, pair_id in pairs:
            writer.writerow([seq, label, pair_id])
    return len(pairs)


def write_fasta(out_path, seq_dict):
    """Write an OrderedDict {header: sequence} to FASTA."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        for header, seq in seq_dict.items():
            fh.write(f">{header}\n{seq}\n")
    return len(seq_dict)


# ---------------------------------------------------------------------------
# Task builders
# ---------------------------------------------------------------------------

def build_binary_from_fastas(task_name, pos_fasta, neg_fasta, repo_root,
                             dry_run=False):
    """Standard classification task from positive/negative FASTAs."""
    pos = fasta_to_seqs(pos_fasta)
    neg = fasta_to_seqs(neg_fasta)
    out = repo_root / "tasks" / task_name / "labels.tsv"
    if dry_run:
        print(f"  {task_name}: {len(pos)} AMP + {len(neg)} non-AMP -> {out}")
        return pos, neg
    n = write_classification_labels(str(out), pos, neg)
    print(f"  {task_name}: {n} rows written")
    return pos, neg


def build_synthetic_task(task_name, pos_fasta, neg_fasta, repo_root,
                         dry_run=False):
    """Synthetic tasks: broad_positive = AMP, synthetic file = non-AMP."""
    pos = fasta_to_seqs(pos_fasta)
    neg = fasta_to_seqs(neg_fasta)
    out = repo_root / "tasks" / task_name / "labels.tsv"
    if dry_run:
        print(f"  {task_name}: {len(pos)} AMP + {len(neg)} non-AMP -> {out}")
        return pos, neg
    n = write_classification_labels(str(out), pos, neg)
    print(f"  {task_name}: {n} rows written")
    return pos, neg


def build_slay_task(slay_csv, repo_root, dry_run=False):
    """SLAY task from the CSV produced by 09_get_slay.py."""
    import csv as csv_mod
    pos, neg = [], []
    with open(slay_csv) as fh:
        reader = csv_mod.DictReader(fh)
        for row in reader:
            seq = row["Sequence"].strip() if row["Sequence"] else ""
            if not seq:
                continue
            lfc = float(row["lfcMLE"])
            if lfc <= -1:
                pos.append(seq)
            else:
                neg.append(seq)

    out = repo_root / "tasks" / "slay" / "labels.tsv"
    if dry_run:
        print(f"  slay: {len(pos)} AMP + {len(neg)} non-AMP -> {out}")
        return pos, neg
    n = write_classification_labels(str(out), pos, neg)
    print(f"  slay: {n} rows written")
    return pos, neg


def build_regression_task(task_name, mic_csv, repo_root, dry_run=False,
                         mic_column="activity", target_unit="ug/ml",
                         convert_units=False, aa_mass=None):
    """MIC regression task from the CSV produced by step 05.

    Parameters
    ----------
    mic_column : str
        Which column to read as the MIC value.  The CSV contains both:
        - "MIC": min-aggregated across duplicate measurements per peptide
        - "activity": per-row value from the first surviving measurement
        Default is "activity".  "MIC" is the scientifically more
        conservative choice (min across multiple measurements); it differs
        for ~0.5% of sequences.
    target_unit : str
        Unit label to write in the output.  Default: "ug/ml".
    convert_units : bool
        If True, convert values whose source unit differs from target_unit
        using molecular weight.  If False (default), just label all values
        as target_unit regardless of source unit.
    aa_mass : dict
        Amino acid residue masses for MW calculation (only used when
        convert_units=True).
    """
    if aa_mass is None:
        aa_mass = {
            "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
            "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
            "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
            "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
            "T": 101.1051, "W": 186.2132, "Y": 163.176, "V": 99.1326,
        }
    water_mass = 18.01528

    def compute_mw(sequence):
        return sum(aa_mass.get(aa, 110.0) for aa in sequence) + water_mass

    def parse_unit(raw_unit):
        """Classify a raw unit string as 'ug/ml' or 'uM'."""
        raw = raw_unit.lower().replace("\u00b5", "u")
        if "g" in raw and "ml" in raw:
            return "ug/ml"
        if "um" in raw or raw.endswith("m"):
            return "uM"
        return "ug/ml"

    import csv as csv_mod
    seen = set()
    rows = []
    with open(mic_csv) as fh:
        reader = csv_mod.DictReader(fh)
        for row in reader:
            seq = row["sequence"].strip()
            if seq in seen:
                continue
            seen.add(seq)
            mic = float(row[mic_column])

            if convert_units:
                src_unit = parse_unit(row.get("unit", "ug/ml"))
                if src_unit != target_unit:
                    mw = compute_mw(seq)
                    if src_unit == "uM" and target_unit == "ug/ml":
                        mic = mic * mw / 1000.0
                    elif src_unit == "ug/ml" and target_unit == "uM":
                        mic = mic * 1000.0 / mw

            rows.append((seq, mic, target_unit))

    out = repo_root / "tasks" / task_name / "labels.tsv"
    if dry_run:
        print(f"  {task_name}: {len(rows)} rows -> {out}")
        return [r[0] for r in rows]
    n = write_regression_labels(str(out), rows)
    print(f"  {task_name}: {n} rows written")
    return [r[0] for r in rows]


def build_paired_task(task_name, tsv_path, repo_root, dry_run=False):
    """Paired task from paired_MD.tsv (alternating AMP/non-AMP rows)."""
    import csv as csv_mod
    records = []
    with open(tsv_path) as fh:
        reader = csv_mod.DictReader(fh, delimiter="\t")
        for row in reader:
            records.append(row["sequence"].strip())

    pairs = []
    for i in range(0, len(records), 2):
        pair_id = i // 2
        pairs.append((records[i], "AMP", pair_id))
        if i + 1 < len(records):
            pairs.append((records[i + 1], "non-AMP", pair_id))

    out = repo_root / "tasks" / task_name / "labels.tsv"
    if dry_run:
        print(f"  {task_name}: {len(pairs)} rows ({len(pairs)//2} pairs) "
              f"-> {out}")
        return [p[0] for p in pairs]
    n = write_paired_labels(str(out), pairs)
    print(f"  {task_name}: {n} rows written")
    return [p[0] for p in pairs]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Task registry: maps task_name -> (builder, *args relative to AMP_DATA_ROOT)
# Args use forward slashes; resolved at runtime.

BINARY_TASKS = {
    # AMP identification (from external databases, not DBAASP pipeline)
    "amp": ("data/amp_positive.fasta", "data/amp_negative.fasta"),
    # Broad activity
    "broad_activity": (
        "data/activity/broad_positive.fasta",
        "data/activity/broad_negative.fasta",
    ),
    # Gram-type
    "gram_minus": (
        "data/activity/gramminus_positive.fasta",
        "data/activity/gramminus_negative.fasta",
    ),
    "gram_plus": (
        "data/activity/gramplus_positive.fasta",
        "data/activity/gramplus_negative.fasta",
    ),
    # Species-level
    "species_ecoli": (
        "data/activity/species/escherichiacoli_positive.fasta",
        "data/activity/species/escherichiacoli_negative.fasta",
    ),
    "species_saureus": (
        "data/activity/species/staphylococcusaureus_positive.fasta",
        "data/activity/species/staphylococcusaureus_negative.fasta",
    ),
    "species_kpneumoniae": (
        "data/activity/species/klebsiellapneumoniae_positive.fasta",
        "data/activity/species/klebsiellapneumoniae_negative.fasta",
    ),
    "species_abaumannii": (
        "data/activity/species/acinetobacterbaumannii_positive.fasta",
        "data/activity/species/acinetobacterbaumannii_negative.fasta",
    ),
    "species_paeruginosa": (
        "data/activity/species/pseudomonasaeruginosa_positive.fasta",
        "data/activity/species/pseudomonasaeruginosa_negative.fasta",
    ),
    # Strain-level
    "strain_ecoli25922": (
        "data/activity/strain/escherichiacoliatcc25922_positive.fasta",
        "data/activity/strain/escherichiacoliatcc25922_negative.fasta",
    ),
    "strain_saureus25923": (
        "data/activity/strain/staphylococcusaureusatcc25923_positive.fasta",
        "data/activity/strain/staphylococcusaureusatcc25923_negative.fasta",
    ),
    "strain_kpneumoniae700603": (
        "data/activity/strain/klebsiellapneumoniaeatcc700603_positive.fasta",
        "data/activity/strain/klebsiellapneumoniaeatcc700603_negative.fasta",
    ),
    "strain_abaumannii19606": (
        "data/activity/strain/acinetobacterbaumanniiatcc19606_positive.fasta",
        "data/activity/strain/acinetobacterbaumanniiatcc19606_negative.fasta",
    ),
    "strain_paeruginosa27853": (
        "data/activity/strain/pseudomonasaeruginosaatcc27853_positive.fasta",
        "data/activity/strain/pseudomonasaeruginosaatcc27853_negative.fasta",
    ),
    "strain_saureus43300": (
        "data/activity/strain/staphylococcusaureusatcc43300_positive.fasta",
        "data/activity/strain/staphylococcusaureusatcc43300_negative.fasta",
    ),
    "strain_saureus33591": (
        "data/activity/strain/staphylococcusaureusatcc33591_positive.fasta",
        "data/activity/strain/staphylococcusaureusatcc33591_negative.fasta",
    ),

}

SYNTHETIC_TASKS = {
    "synthetic_random": (
        "data/activity/broad_positive.fasta",
        "data/syntax/synthetic_random.fasta",
    ),
    "synthetic_realistic": (
        "data/activity/broad_positive.fasta",
        "data/syntax/synthetic_realistic.fasta",
    ),
    "synthetic_shuffled": (
        "data/activity/broad_positive.fasta",
        "data/syntax/broad_positive_shuffled.fasta",
    ),
}

REGRESSION_TASKS = {
    "regression_ecoli25922":
        "data/activity/strain/escherichiacoliatcc25922_mic.csv",
    "regression_saureus25923":
        "data/activity/strain/staphylococcusaureusatcc25923_mic.csv",
}

PAIRED_TASKS = {
    # paired_sequences and paired_activity_cliffs are identical in the
    # shipped benchmark; both derive from paired_MD.tsv
    "paired_sequences": "data/syntax/paired_MD.tsv",
    "paired_activity_cliffs": "data/syntax/paired_MD.tsv",
}

# paired_safe is NOT generated by the AMP-data pipeline; it requires
# hemolytic data that must be curated separately.  If the file exists
# it will be copied; otherwise a warning is printed.
MANUAL_TASKS = {
    "paired_safe": "tasks/paired_safe/labels.tsv",
}


def resolve(root, relpath):
    return root / relpath


def main():
    parser = argparse.ArgumentParser(
        description="Build BattleAMP benchmark tasks from AMP-data outputs."
    )
    parser.add_argument(
        "--amp-data-root", required=True, type=Path,
        help="Root of the AMP-data repository (contains data/ directory)."
    )
    parser.add_argument(
        "--repo-root", required=True, type=Path,
        help="Root of the battleamp-snakemake repository."
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print what would be created without writing files."
    )
    parser.add_argument(
        "--skip-fasta", action="store_true",
        help="Skip building the unified sequences.fasta."
    )
    parser.add_argument(
        "--skip-slay", action="store_true",
        help="Skip the SLAY task (requires network or pre-downloaded CSV)."
    )
    parser.add_argument(
        "--mic-column", default="activity", choices=["MIC", "activity"],
        help="Which CSV column to use for regression MIC values. "
             "'activity' uses the per-row value from the first measurement "
             "(default). 'MIC' uses the min-aggregated value "
             "(scientifically more conservative; differs for ~0.5%% of "
             "sequences)."
    )
    parser.add_argument(
        "--target-unit", default="ug/ml", choices=["ug/ml", "uM"],
        help="Unit label for regression MIC values. Default: ug/ml."
    )
    parser.add_argument(
        "--convert-units", action="store_true",
        help="Actually convert uM values to ug/ml using per-peptide MW. "
             "By default, values are written as-is with the target-unit "
             "label regardless of source unit."
    )
    args = parser.parse_args()

    amp = args.amp_data_root
    repo = args.repo_root
    dry = args.dry_run

    # Collect all sequences for the unified FASTA
    # key = sequence, value = first header seen
    all_sequences = OrderedDict()

    def collect_fasta(fasta_path):
        """Add sequences from a FASTA to the global pool."""
        for header, seq in read_fasta_sequences(fasta_path):
            if seq not in all_sequences:
                all_sequences[seq] = header

    def collect_seqs(seq_list, prefix="seq"):
        """Add plain sequences (no headers) to the global pool."""
        for i, seq in enumerate(seq_list):
            if seq not in all_sequences:
                all_sequences[seq] = f"{prefix}_{i}"

    # ---------------------------------------------------------------
    print("Building classification tasks...")
    errors = []
    for task_name, (pos_rel, neg_rel) in BINARY_TASKS.items():
        pos_path = resolve(amp, pos_rel)
        neg_path = resolve(amp, neg_rel)
        if not pos_path.exists() or not neg_path.exists():
            errors.append(f"  SKIP {task_name}: missing input "
                          f"({pos_rel} or {neg_rel})")
            continue
        pos, neg = build_binary_from_fastas(
            task_name, str(pos_path), str(neg_path), repo, dry
        )
        collect_fasta(str(pos_path))
        collect_fasta(str(neg_path))

    # ---------------------------------------------------------------
    print("\nBuilding synthetic tasks...")
    for task_name, (pos_rel, neg_rel) in SYNTHETIC_TASKS.items():
        pos_path = resolve(amp, pos_rel)
        neg_path = resolve(amp, neg_rel)
        if not pos_path.exists() or not neg_path.exists():
            errors.append(f"  SKIP {task_name}: missing input")
            continue
        pos, neg = build_synthetic_task(
            task_name, str(pos_path), str(neg_path), repo, dry
        )
        collect_fasta(str(pos_path))
        collect_fasta(str(neg_path))

    # ---------------------------------------------------------------
    # Load amino acid masses for unit conversion
    aa_mass = None
    aa_mass_path = amp / "data" / "aa_mass.json"
    if aa_mass_path.exists():
        import json
        with open(aa_mass_path) as fh:
            aa_mass = json.load(fh)

    print("\nBuilding regression tasks...")
    for task_name, csv_rel in REGRESSION_TASKS.items():
        csv_path = resolve(amp, csv_rel)
        if not csv_path.exists():
            errors.append(f"  SKIP {task_name}: missing {csv_rel}")
            continue
        seqs = build_regression_task(
            task_name, str(csv_path), repo, dry,
            mic_column=args.mic_column,
            target_unit=args.target_unit,
            convert_units=args.convert_units,
            aa_mass=aa_mass,
        )
        # Regression MIC CSVs also have FASTA counterparts
        fasta_rel = csv_rel.replace(".csv", ".fasta")
        fasta_path = resolve(amp, fasta_rel)
        if fasta_path.exists():
            collect_fasta(str(fasta_path))
        else:
            collect_seqs(seqs, prefix=task_name)

    # ---------------------------------------------------------------
    print("\nBuilding paired tasks...")
    for task_name, tsv_rel in PAIRED_TASKS.items():
        tsv_path = resolve(amp, tsv_rel)
        if not tsv_path.exists():
            errors.append(f"  SKIP {task_name}: missing {tsv_rel}")
            continue
        seqs = build_paired_task(task_name, str(tsv_path), repo, dry)
        # Paired FASTA
        fasta_rel = tsv_rel.replace(".tsv", ".fasta")
        fasta_path = resolve(amp, fasta_rel)
        if fasta_path.exists():
            collect_fasta(str(fasta_path))
        else:
            collect_seqs(seqs, prefix=task_name)
        # Also collect the shuffled paired FASTA if present
        shuf_path = resolve(amp, "data/syntax/pared_MD_shuffled.fasta")
        if shuf_path.exists():
            collect_fasta(str(shuf_path))

    # ---------------------------------------------------------------
    if not args.skip_slay:
        print("\nBuilding SLAY task...")
        slay_csv = resolve(amp, "data/slay/slay_all.csv")
        if slay_csv.exists():
            pos, neg = build_slay_task(str(slay_csv), repo, dry)
            # Collect SLAY sequences
            slay_fasta = resolve(amp, "data/slay/slay_all.fasta")
            if slay_fasta.exists():
                collect_fasta(str(slay_fasta))
            else:
                collect_seqs(pos, prefix="slay_pos")
                collect_seqs(neg, prefix="slay_neg")
        else:
            errors.append("  SKIP slay: data/slay/slay_all.csv not found "
                          "(run 09_get_slay.py first)")

    # ---------------------------------------------------------------
    print("\nManual/external tasks...")
    for task_name, label_rel in MANUAL_TASKS.items():
        src = resolve(repo, label_rel)
        if src.exists():
            print(f"  {task_name}: already present at {src}")
            # Collect its sequences into the FASTA pool
            import csv as csv_mod
            with open(src) as fh:
                reader = csv_mod.DictReader(fh, delimiter="\t")
                for row in reader:
                    seq = row["sequence"].strip()
                    if seq and seq not in all_sequences:
                        all_sequences[seq] = f"{task_name}_{len(all_sequences)}"
        else:
            errors.append(
                f"  SKIP {task_name}: {label_rel} not found. "
                f"This task requires manual curation."
            )

    # ---------------------------------------------------------------
    # Also collect any remaining FASTAs listed in fastas_list.json
    fastas_json = resolve(amp, "data/fastas_list.json")
    if fastas_json.exists():
        import json
        with open(fastas_json) as fh:
            flist = json.load(fh)["fastas"]
        for rel in flist:
            # Paths in fastas_list.json are relative to AMP-data root
            fpath = resolve(amp, rel.lstrip("./"))
            if fpath.exists():
                collect_fasta(str(fpath))

    # ---------------------------------------------------------------
    if not args.skip_fasta:
        print(f"\nAssembling unified FASTA: {len(all_sequences)} unique "
              f"sequences...")
        fasta_out = repo / "datasets" / "battleamp-all" / "sequences.fasta"
        if dry:
            print(f"  Would write {len(all_sequences)} sequences to "
                  f"{fasta_out}")
        else:
            # Write as header -> sequence
            fasta_dict = OrderedDict()
            for seq, header in all_sequences.items():
                fasta_dict[header] = seq
            n = write_fasta(str(fasta_out), fasta_dict)
            print(f"  Wrote {n} sequences to {fasta_out}")

    # ---------------------------------------------------------------
    if errors:
        print(f"\nWarnings ({len(errors)}):")
        for e in errors:
            print(e)

    print("\nDone.")


if __name__ == "__main__":
    main()