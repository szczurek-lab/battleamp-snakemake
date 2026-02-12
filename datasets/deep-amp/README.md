# Deep-AMP Training Dataset -- BattleAMP Benchmark

Training data from the Deep-AMP model (Pandi et al., Nature Communications 2023),
reformatted as a BattleAMP benchmark dataset.

## Source

- Paper: "Cell-free biosynthesis combined with deep learning accelerates
  de novo-development of antimicrobial peptides"
- DOI: 10.1038/s41467-023-42434-9
- GitHub: https://github.com/amirpandi/Deep_AMP

## Data Origin

The dataset aggregates MIC measurements from public databases. MIC values
are provided per gram stain category (gram-negative and gram-positive),
not per individual species. Non-AMPs are sourced separately and have no
MIC annotation; they are assigned label=0 in classification tasks.

36 sequences appeared in both the AMP and non-AMP sets in the original
data. These were removed from the non-AMP set to avoid conflicts.

## Dataset Statistics

| Subset                | Sequences |
|-----------------------|-----------|
| Unique AMPs (total)   | 5,064     |
| -- with gram-neg MIC  | 4,619     |
| -- with gram-pos MIC  | 4,175     |
| -- with both          | 3,730     |
| Non-AMPs              | 7,668     |
| Total                 | 12,732    |

MIC range (uM):
- Gram-negative: 0.02 to 3,200  (median 13.5)
- Gram-positive: 0.002 to 7,251 (median 12.5)

Sequence lengths: 1 to 190 aa (median 19 for AMPs, 34 for non-AMPs).
139 AMPs exceed 48 aa (the Deep-AMP model's max input length).

## Task Files

### all_sequences.fasta
All 12,732 sequences in FASTA format.

### task_classification.tsv
General AMP activity classification (binary).
Active = minimum MIC across gram types <= 32 uM.
Inactive = all non-AMPs + AMPs with minimum MIC >= 128 uM.
Grey zone (32 < MIC < 128 for all measurements) excluded.
- 11,727 sequences: 3,730 active, 7,997 inactive

### task_species_gramneg.tsv
Gram-negative-specific classification.
Active = gram-neg MIC <= 32 uM.
Inactive = non-AMPs + AMPs with gram-neg MIC >= 128 uM.
- 11,111 sequences: 3,062 active, 8,049 inactive

### task_species_grampos.tsv
Gram-positive-specific classification.
Active = gram-pos MIC <= 32 uM.
Inactive = non-AMPs + AMPs with gram-pos MIC >= 128 uM.
- 10,826 sequences: 2,789 active, 8,037 inactive

### task_regression_gramneg.tsv
Gram-negative MIC regression. Columns: sequence, MIC, MIC_unit.
- 4,619 sequences with MIC in uM

### task_regression_grampos.tsv
Gram-positive MIC regression. Columns: sequence, MIC, MIC_unit.
- 4,175 sequences with MIC in uM

## Data Leakage Warning

Deep-AMP models (all 4 variants) were trained on this exact dataset.
Results for deep-amp-lstm-gramneg, deep-amp-lstm-grampos,
deep-amp-cnn-gramneg, and deep-amp-cnn-grampos on these tasks reflect
in-distribution / memorization performance, not generalization. This is
documented intentionally: comparing Deep-AMP's in-distribution
performance against other models' out-of-distribution performance
reveals whether other models can generalize to this data distribution
without having seen it during training.

For all other models in BattleAMP, this dataset represents a held-out
evaluation set, assuming no sequence overlap with their training data.
Overlap should be checked and reported (e.g., via the QMAP-style
max-identity analysis).

## Thresholds

Classification thresholds follow BattleAMP conventions:
- Active: MIC <= 32 uM (equivalent to <= 32 ug/mL for a typical
  ~1000 Da peptide)
- Inactive: MIC >= 128 uM
- Grey zone excluded from classification tasks

## Citation

Pandi A, Adam D, Zare A, et al. Cell-free biosynthesis combined with
deep learning accelerates de novo-development of antimicrobial
peptides. Nat Commun 14, 7197 (2023).
