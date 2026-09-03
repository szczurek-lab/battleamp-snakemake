# holdout-2026

A temporal holdout: 600 peptides curated by DBAASP after the January 2025 snapshot that
every other BattleAMP dataset derives from. No benchmarked model can have trained on them,
so this dataset measures generalization rather than recall.

## Provenance

A refreshed 2026 DBAASP download, restricted to peptides whose sequence appears in no earlier
BattleAMP dataset, then processed into activity datasets using the same MIC thresholds as the
main benchmark (active at 32 ug/ml, inactive at 128 ug/ml).

The two DBAASP id ranges are disjoint. The January 2025 snapshot that every other BattleAMP
dataset derives from spans ids 10 to 22,878; this dataset spans 22,881 to 25,325. Every
peptide here was curated after that snapshot, and every benchmarked model was published
before it, so none could have been trained on these sequences.

DBAASP records no submission date, so the id range is the motivation rather than the proof.
Novelty is enforced by sequence, and the overlap table below is computed by exact matching
against each model's training set.

## Contents

600 unique sequences, 6 to 85 residues, standard 20 amino acids only.

| Task | Labels file | AMP | non-AMP |
|---|---|---|---|
| `holdout2026_broad_activity` | `labels_broad_activity.tsv` | 430 | 91 |
| `holdout2026_gram_minus` | `labels_gram_minus.tsv` | 325 | 91 |
| `holdout2026_gram_plus` | `labels_gram_plus.tsv` | 264 | 109 |
| `holdout2026_species_ecoli` | `labels_species_ecoli.tsv` | 213 | 96 |
| `holdout2026_species_saureus` | `labels_species_saureus.tsv` | 264 | 109 |
| `holdout2026_species_paeruginosa` | `labels_species_paeruginosa.tsv` | 212 | 80 |
| `holdout2026_species_kpneumoniae` | `labels_species_kpneumoniae.tsv` | 112 | 28 |
| `holdout2026_species_abaumannii` | `labels_species_abaumannii.tsv` | 163 | 34 |

| Task | Labels file | Peptides |
|---|---|---|
| `holdout2026_regression_ecoli` | `labels_mic_ecoli.tsv` | 360 |
| `holdout2026_regression_saureus` | `labels_mic_saureus.tsv` | 451 |
| `holdout2026_regression_paeruginosa` | `labels_mic_paeruginosa.tsv` | 372 |
| `holdout2026_regression_kpneumoniae` | `labels_mic_kpneumoniae.tsv` | 186 |
| `holdout2026_regression_abaumannii` | `labels_mic_abaumannii.tsv` | 237 |

MIC values are in ug/ml, taken as the minimum across measurements for that peptide and
species. The source `unit` column describes the original measurement and is not used.

Task names mirror their `battleamp-all` counterparts so results can be compared directly.
There is no strain-level data, so the mirror stops at species level.

## Overlap with model training data

Checked against every benchmarked model's training set by exact sequence matching. Counts are
distinct real peptides the released checkpoint was fitted on, derived from each model's own
training script rather than from row counts, since most apply length filters or subsampling.

| Model | Training sequences | Overlap | Basis |
|---|---|---|---|
| apex | 1,388 | 2 | CV group; Ind split held out |
| ampscanner | 2,132 | 0 | train and validation merged, as the shipped `FULL_MODEL` weights imply |
| ampredictor | 2,952 | 0 | train plus validation, the latter driving early stopping |
| hydramp-mic-classifier | 3,444 | 0 | positives up to 25 residues |
| mbc-attention | 3,929 | 0 | `whole_train.mdl` covers all of `EC.csv` |
| ampeppy | 6,536 | 0 | 3,268 positives and a 3,268 negative subsample |
| amplify | 6,676 | 0 | the balanced model, which is the one deployed |
| hydramp-amp-classifier | 11,131 | 0 | positives up to 25 residues |
| ampredmfa | 11,320 | 0 | full pool; the shipped checkpoint records no trial or fold |
| deep-amp | 12,593 | 0 | up to 48 residues, union of the four variants |
| sensexamp | 20,587 | 0 | training splits of the three variants |

Counted per benchmarked variant: data belonging to other components of the same project, such
as HydrAMP's VAE generator or Deep-AMP's toxicity head, is out of scope.

HydrAMP trains on random fragments of UniProt sequences rather than the sequences themselves,
so its negative class cannot be matched exactly. It was checked against the full source pool,
which is conservative.

The two APEX matches, `PMAKLLPRIKKKILAAAFK` and `PMARNKPKILKRILAKIFK`, are in that model's
cross-validation group and are labelled active here.

One sequence, `RRWWRRWW`, also appears in `battleamp-all` as the shuffled decoy
`19949_shuffled_7`. That is a collision between a random shuffle and a real octapeptide, not
shared provenance.
