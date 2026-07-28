#!/usr/bin/env python3
"""Generate the 50-peptide reference dataset used to smoke test every variant.

The point of this dataset is coverage, not biology. Sequences are sampled from
datasets/battleamp-all so they are real peptides, but the length bins are chosen
to straddle every model's declared limits, so that a single run exercises:

  - peptides every model accepts,
  - peptides some models must skip for length (status "partial"),
  - peptides no model accepts,
  - records rejected before any model runs (D-amino acids, non-standard
    residues, gap characters, duplicates).

Model length limits this is built against:

    mbc-attention    5-60      sensexamp     6-25     hydramp    1-25
    ampscanner      10-200     apex          1-50     deep-amp   1-48
    amplify          2-200     ampredictor   1-65     ampeppy    any

Deterministic: same output every run.

Usage:
    python scripts/generate_smoke_dataset.py
"""

import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from battleamp import sequences  # noqa: E402

SOURCE = REPO_ROOT / "datasets" / "battleamp-all" / "sequences.fasta"
OUTPUT = REPO_ROOT / "examples" / "smoke_reference.fasta"
SEED = 20260728

# (label, min_len, max_len, count, what it tests)
BINS = [
    ("tiny",   1,   5,  7, "below mbc-attention(5), sensexamp(6), ampscanner(10)"),
    ("short",  6,   9,  7, "sensexamp/hydramp accept; ampscanner too short"),
    ("core",  10,  25,  8, "every model accepts"),
    ("mid",   26,  50,  7, "apex/deep-amp ceiling; sensexamp/hydramp too long"),
    ("long",  51,  65,  6, "ampredictor ceiling; apex too long"),
    ("xlong", 66, 200,  6, "only amplify/ampscanner/ampeppy accept"),
    ("huge", 201, 999,  5, "above every declared maximum"),
]

# Records that must be rejected before any model runs.
INVALID = [
    ("invalid_lowercase_dAA",   "giGKFLHSAKKFGKAF",  "lowercase = D-amino acids"),
    ("invalid_nonstandard_X",   "GIGKFLHSAXKFGKAF",  "X is not a standard residue"),
    ("invalid_gap_character",   "GIGK-FLHSAKKFGKAF", "gap character"),
]


def main():
    if not SOURCE.exists():
        print(f"ERROR: source dataset not found: {SOURCE}", file=sys.stderr)
        return 1

    records = sequences.parse_fasta(SOURCE.read_text())
    by_length = {}
    for _id, seq in records:
        by_length.setdefault(len(seq), []).append(seq)

    rng = random.Random(SEED)
    chosen = []
    used = set()

    for label, lo, hi, count, _why in BINS:
        pool = sorted(
            {s for L, seqs in by_length.items() if lo <= L <= hi for s in seqs}
        )
        pool = [s for s in pool if s not in used]
        if len(pool) < count:
            print(f"WARNING: only {len(pool)} sequences for bin {label} "
                  f"({lo}-{hi}), wanted {count}", file=sys.stderr)
            count = len(pool)
        picked = rng.sample(pool, count)
        for i, seq in enumerate(picked, start=1):
            used.add(seq)
            chosen.append((f"smoke_{label}_{len(seq):03d}aa_{i:02d}", seq))

    lines = []
    for name, seq in chosen:
        lines.append(f">{name}")
        for i in range(0, len(seq), 80):
            lines.append(seq[i:i + 80])

    # A duplicate of the first core peptide: must collapse onto one row and
    # keep both headers.
    first_core = next(s for n, s in chosen if "_core_" in n)
    lines.append(">invalid_duplicate_of_core")
    lines.append(first_core)

    for name, seq, _why in INVALID:
        lines.append(f">{name}")
        lines.append(seq)

    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text("\n".join(lines) + "\n")

    n_records = len(chosen) + 1 + len(INVALID)
    print(f"Wrote {OUTPUT.relative_to(REPO_ROOT)}: {n_records} records")
    for label, lo, hi, count, why in BINS:
        actual = sum(1 for n, _ in chosen if f"_{label}_" in n)
        print(f"  {label:<6} {lo:>3}-{hi:<3} aa  n={actual}   {why}")
    print(f"  {'invalid':<6} {'':>3} {'':<3}     n={len(INVALID) + 1}   "
          f"rejected before any model runs (incl. 1 duplicate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
