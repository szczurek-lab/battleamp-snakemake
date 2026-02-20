"""
Extract predictions for MD validation datasets from battleamp-all inference.

Reads:
  - control-md.tsv        (control peptides: Sequence, true_class, ...)
  - activity-cliffs.tsv    (paired peptides: sequence, MIC [ug/ml], ...)
  - inference/{variant}/battleamp-all/predictions.tsv  (per-model)

Produces:
  - control-md_predictions.tsv        (wide: one row per peptide, columns per model)
  - activity-cliffs_predictions.tsv   (long: one row per model x peptide)
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
)
log = logging.getLogger(__name__)

# ── Parameters ───────────────────────────────────────────────────────
inference_dir = Path(snakemake.params.inference_dir)
dataset = snakemake.params.dataset

control = pd.read_csv(snakemake.input.control_labels, sep="\t")
cliffs = pd.read_csv(snakemake.input.cliff_labels, sep="\t")

control_seqs = set(control["Sequence"].str.upper())
cliff_seqs = set(cliffs["sequence"].str.upper())

log.info("Control peptides: %d sequences", len(control_seqs))
log.info("Activity cliff peptides: %d sequences", len(cliff_seqs))

SKIP_VARIANTS = {"example-model"}

# ── Helpers ──────────────────────────────────────────────────────────


def find_sequence_col(df):
    """Find the sequence column in a predictions dataframe."""
    for candidate in ("sequence", "Sequence", "seq", "Seq", "peptide"):
        if candidate in df.columns:
            return candidate
    return None


def classify_model(df):
    """Determine if a model is a classifier or regressor."""
    has_prob = (
        "Probability_score" in df.columns
        and df["Probability_score"].notna().any()
    )
    has_mic = "MIC" in df.columns and df["MIC"].notna().any()
    if has_prob:
        return "classifier"
    if has_mic:
        return "regressor"
    return None


# ── Walk inference results ───────────────────────────────────────────
control_wide = {}  # seq_upper -> {col: val}
cliff_long = []    # list of dicts

variant_dirs = sorted(inference_dir.iterdir()) if inference_dir.exists() else []

for vdir in variant_dirs:
    pred_file = vdir / dataset / "predictions.tsv"
    if not pred_file.exists():
        continue

    variant = vdir.name
    if variant in SKIP_VARIANTS:
        continue

    try:
        preds = pd.read_csv(pred_file, sep="\t")
    except Exception as e:
        log.warning("Failed to read %s: %s", pred_file, e)
        continue

    seq_col = find_sequence_col(preds)
    if seq_col is None:
        log.info(
            "%s: no sequence column found (columns: %s), skipping",
            variant,
            list(preds.columns)[:5],
        )
        continue

    model_type = classify_model(preds)
    if model_type is None:
        log.info("%s: could not determine model type, skipping", variant)
        continue

    preds["_seq_upper"] = preds[seq_col].astype(str).str.upper()

    # ── Control peptides (wide format) ───────────────────────────
    ctrl_match = preds[preds["_seq_upper"].isin(control_seqs)]
    n_ctrl = len(ctrl_match)

    for _, row in ctrl_match.iterrows():
        seq = row["_seq_upper"]
        if seq not in control_wide:
            control_wide[seq] = {}

        if model_type == "classifier":
            if "Prediction" in preds.columns:
                control_wide[seq][f"{variant}_Prediction"] = row.get(
                    "Prediction"
                )
            if "Probability_score" in preds.columns:
                control_wide[seq][f"{variant}_Probability_score"] = row.get(
                    "Probability_score"
                )
        else:
            if "MIC" in preds.columns:
                control_wide[seq][f"{variant}_MIC"] = row.get("MIC")

    # ── Activity cliff peptides (long format) ────────────────────
    cliff_match = preds[preds["_seq_upper"].isin(cliff_seqs)]
    n_cliff = len(cliff_match)

    for _, row in cliff_match.iterrows():
        entry = {"variant": variant, "sequence": row[seq_col]}
        if model_type == "classifier":
            entry["Probability_score"] = row.get("Probability_score")
            entry["Prediction"] = row.get("Prediction")
            entry["MIC"] = np.nan
        else:
            entry["Probability_score"] = np.nan
            entry["Prediction"] = np.nan
            entry["MIC"] = row.get("MIC")
        cliff_long.append(entry)

    log.info(
        "%s (%s): control=%d/%d  cliffs=%d/%d",
        variant,
        model_type,
        n_ctrl,
        len(control_seqs),
        n_cliff,
        len(cliff_seqs),
    )

# ── Build control peptides wide table ────────────────────────────────
ctrl_preds_df = pd.DataFrame.from_dict(control_wide, orient="index")
ctrl_preds_df.index.name = "_seq_upper"
ctrl_preds_df.reset_index(inplace=True)

# Merge with labels (keep metadata columns from the input file)
control["_seq_upper"] = control["Sequence"].str.upper()

# Strip any existing prediction columns from the input so we replace them
meta_cols = ["sim_round", "id", "name", "Sequence", "true_class"]
extra_meta = ["classifier_rank", "regressor_rank"]
keep_cols = [c for c in meta_cols + extra_meta if c in control.columns]
ctrl_labels = control[keep_cols + ["_seq_upper"]].copy()

ctrl_out = ctrl_labels.merge(ctrl_preds_df, on="_seq_upper", how="left")
ctrl_out.drop(columns=["_seq_upper"], inplace=True)

# Sort columns: metadata first, then predictions alphabetically
pred_cols = sorted(c for c in ctrl_out.columns if c not in keep_cols)
ctrl_out = ctrl_out[keep_cols + pred_cols]

ctrl_out.to_csv(snakemake.output.control, sep="\t", index=False)
log.info(
    "Saved control predictions: %d rows x %d columns -> %s",
    len(ctrl_out),
    len(ctrl_out.columns),
    snakemake.output.control,
)

# ── Build activity cliff predictions long table ──────────────────────
cliff_out = pd.DataFrame(cliff_long)
col_order = ["variant", "sequence", "Prediction", "Probability_score", "MIC"]
cliff_out = cliff_out[[c for c in col_order if c in cliff_out.columns]]
cliff_out.sort_values(["variant", "sequence"], inplace=True)

cliff_out.to_csv(snakemake.output.cliffs, sep="\t", index=False)
log.info(
    "Saved activity cliff predictions: %d rows -> %s",
    len(cliff_out),
    snakemake.output.cliffs,
)