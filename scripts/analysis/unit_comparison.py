#!/usr/bin/env python3
"""Compare regression results reported in ug/ml against uM.

MIC in ug/ml depends on molecular weight, so converting to uM reorders
peptides rather than rescaling them.  This script quantifies what that does
to the regression metrics and to the activity threshold.

Only the reporting unit changes.  Predictions are reused as they are: the
benchmark unit is a configuration setting and evaluate.py converts every
predicted and measured MIC per peptide from a molecular weight computed from
its sequence, so no re-inference is needed.  mic_clamp stays configured in
ug/ml and its bounds are converted per peptide, so clamping is physically
identical in both runs and the comparison isolates the unit itself.

Classification labels were assigned in ug/ml when the datasets were built, so
changing the unit would redefine those tasks rather than re-score them.  The
regression tasks are therefore scored in both units, and the threshold effect
is reported as a reclassification count.  The uM analogue of the 32 ug/ml
activity threshold is not unique: 32 uM is a round-number convention used
elsewhere in the field, including Deep-AMP's own training labels, while 16 uM
is closer to 32 ug/ml back-converted at a typical 2.5 kDa peptide.  Both are
reported.

Prerequisites
-------------
    results/inference/ must be populated, which the pipeline produces with
    ``snakemake --profile profile/ score``.

Usage
-----
    python scripts/analysis/unit_comparison.py

Outputs
-------
    unit_metrics.tsv      every regressor and regression task, scored in both units
    classification_metrics.tsv
                          the same peptides scored as a classification task with
                          activity defined in each unit
    threshold_effect.tsv  per task: molecular weight range, peptide-level rank
                          correlation between units, and how many peptides change
                          activity class at each uM threshold
    peptide_units.tsv     per peptide: sequence, length, molecular weight, MIC in both units
"""
import argparse
import os
import sys

import pandas as pd
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "workflow", "scripts"))

from evaluate import compute_mw_array, convert_mic_ugml_to_um, evaluate_task

# The two MIC regression tasks of the benchmark.  The deep-amp tasks run on
# that model's own training data and the holdout tasks are not part of the
# published benchmark, so neither belongs in this comparison.
REGRESSION_TASKS = ("regression_ecoli25922", "regression_saureus25923")
UM_THRESHOLDS = (8, 16, 32)
UNITS = ("ug/ml", "uM")


def read_config():
    with open(os.path.join(REPO, "config/config.yaml")) as handle:
        return yaml.safe_load(handle)


def variant_index(config):
    """Map variant name to (model, type), mirroring the Snakefile's VARIANTS block."""
    index = {}
    for model in config["models"]:
        with open(os.path.join(REPO, "models", model, "model.yaml")) as handle:
            meta = yaml.safe_load(handle)
        for variant in meta.get("variants") or [{"name": model}]:
            index[variant["name"]] = (model, variant.get("type", meta["type"]))
    return index


def score_both_units(config, variants, tasks, work_dir):
    """Score every regressor on every regression task, once per unit."""
    rows = []
    for task in tasks:
        dataset = config["tasks"][task]["dataset"]
        for variant, (model, variant_type) in variants.items():
            predictions = f"results/inference/{variant}/{dataset}/predictions.tsv"
            if variant_type != "regressor" or not os.path.exists(predictions):
                continue
            for unit in UNITS:
                out_path = os.path.join(work_dir, variant, task,
                                        unit.replace("/", "") + ".json")
                os.makedirs(os.path.dirname(out_path), exist_ok=True)
                report = evaluate_task(
                    predictions_path=predictions,
                    labels_path=config["tasks"][task]["labels"],
                    output_path=out_path,
                    task_config=config["tasks"][task],
                    variant_type=variant_type,
                    sequence_column=config.get("sequence_column", "sequence"),
                    benchmark_unit=unit,
                    activity_thresholds=config["activity_thresholds"],
                    mic_clamp=config.get("mic_clamp"),
                )
                rows.append(dict(variant=variant, model=model, task=task, unit=unit,
                                 **report["metrics"]))
        print(f"  scored {task}")
    return pd.DataFrame(rows)


def threshold_effect(config, tasks):
    """How the activity threshold and the peptide ordering change with the unit."""
    active_ugml = config["activity_thresholds"]["active"]
    summary, peptides = [], []
    for task in tasks:
        labels = pd.read_csv(config["tasks"][task]["labels"], sep="\t")
        if "MIC" not in labels.columns:
            continue
        lengths = labels["sequence"].str.len()
        mw = compute_mw_array(labels["sequence"].tolist())
        mic_um = convert_mic_ugml_to_um(labels["MIC"].values, mw)
        is_active = labels["MIC"].values <= active_ugml

        row = {
            "task": task,
            "n": len(labels),
            "mw_min": round(float(mw.min())),
            "mw_max": round(float(mw.max())),
            "len_min": int(lengths.min()),
            "len_max": int(lengths.max()),
            "peptide_spearman": round(float(pd.Series(labels["MIC"].values).corr(
                pd.Series(mic_um), method="spearman")), 4),
        }
        for threshold in UM_THRESHOLDS:
            changed = int((is_active != (mic_um <= threshold)).sum())
            row[f"reclassified_at_{threshold}uM"] = changed
            row[f"reclassified_at_{threshold}uM_pct"] = round(100 * changed / len(labels), 2)
        summary.append(row)
        peptides.append(pd.DataFrame({
            "task": task, "sequence": labels["sequence"], "length": lengths, "mw": mw,
            "mic_ugml": labels["MIC"].values, "mic_uM": mic_um,
        }))
    return pd.DataFrame(summary), pd.concat(peptides, ignore_index=True)


def classification_from_mic(config, variants, tasks, work_dir):
    """Score the same peptides as a classification task under each unit.

    Classification labels in the benchmark were assigned in ug/ml when the
    datasets were built, so they cannot simply be re-scored in uM.  The
    regression tasks carry a MIC per peptide, so an equivalent classification
    task can be derived from them under either unit, applying the benchmark's
    own active and inactive thresholds and dropping the grey zone between.
    The peptides and the predictions are identical across the two runs; only
    the unit that defines activity differs.
    """
    active = config["activity_thresholds"]["active"]
    inactive = config["activity_thresholds"]["inactive"]
    rows = []
    for task in tasks:
        source = config["tasks"][task]
        labels = pd.read_csv(source["labels"], sep="\t")
        mw = compute_mw_array(labels["sequence"].tolist())
        mic = {"ug/ml": labels["MIC"].values,
               "uM": convert_mic_ugml_to_um(labels["MIC"].values, mw)}

        for unit in UNITS:
            values = mic[unit]
            called = pd.DataFrame({"sequence": labels["sequence"], "label": None})
            called.loc[values <= active, "label"] = "AMP"
            called.loc[values >= inactive, "label"] = "non-AMP"
            called = called.dropna(subset=["label"])

            labels_path = os.path.join(work_dir, "derived",
                                       f"{task}_{unit.replace('/', '')}_labels.tsv")
            os.makedirs(os.path.dirname(labels_path), exist_ok=True)
            called.to_csv(labels_path, sep="\t", index=False)

            task_config = {"type": "classification", "dataset": source["dataset"],
                           "labels": labels_path}
            thresholds = {"active": active, "inactive": inactive, "unit": unit}
            for variant, (model, variant_type) in variants.items():
                predictions = f"results/inference/{variant}/{source['dataset']}/predictions.tsv"
                if not os.path.exists(predictions):
                    continue
                out_path = os.path.join(work_dir, variant, task,
                                        f"clf_{unit.replace('/', '')}.json")
                os.makedirs(os.path.dirname(out_path), exist_ok=True)
                report = evaluate_task(
                    predictions_path=predictions,
                    labels_path=labels_path,
                    output_path=out_path,
                    task_config=task_config,
                    variant_type=variant_type,
                    sequence_column=config.get("sequence_column", "sequence"),
                    benchmark_unit=unit,
                    activity_thresholds=thresholds,
                    mic_clamp=config.get("mic_clamp"),
                )
                rows.append(dict(variant=variant, model=model, model_type=variant_type,
                                 task=task, unit=unit, n_labelled=len(called),
                                 **report["metrics"]))
        print(f"  derived classification {task}")
    return pd.DataFrame(rows)


def main(out_dir, work_dir):
    os.chdir(REPO)
    config = read_config()
    tasks = [name for name in REGRESSION_TASKS if name in config["tasks"]]
    os.makedirs(out_dir, exist_ok=True)

    metrics = score_both_units(config, variant_index(config), tasks, work_dir)
    metrics.to_csv(os.path.join(out_dir, "unit_metrics.tsv"), sep="\t", index=False)

    classification = classification_from_mic(config, variant_index(config), tasks, work_dir)
    classification.to_csv(os.path.join(out_dir, "classification_metrics.tsv"),
                          sep="\t", index=False)

    summary, peptides = threshold_effect(config, tasks)
    summary.to_csv(os.path.join(out_dir, "threshold_effect.tsv"), sep="\t", index=False)
    peptides.to_csv(os.path.join(out_dir, "peptide_units.tsv"), sep="\t", index=False)

    print(f"\n{len(metrics)} metric rows over {len(tasks)} tasks and "
          f"{metrics['variant'].nunique()} regressors -> {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out-dir", default="results/unit_comparison")
    parser.add_argument("--work-dir", default="results/unit_comparison/_work")
    args = parser.parse_args()
    main(args.out_dir, args.work_dir)
