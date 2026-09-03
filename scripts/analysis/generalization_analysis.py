"""Generalization analysis: training-data overlap and its effect on every task.

Evaluates every benchmarked variant on three versions of each task:

  battleamp    the task as published
  minustrain   the same task with that model's own training sequences removed
  holdout      the holdout-2026 counterpart, curated after every model was
               published and therefore seen by none of them

Metrics come from workflow/scripts/evaluate.py, the same code path that produces
the published numbers, so the three sets are directly comparable.

Training-set sizes are read from each model's training script rather than from
row counts: most apply a length filter or subsample, so the file row count
overstates what the released checkpoint saw. The derivation for each model is
recorded in TRAINING_SETS below and echoed into the output.

Usage:
    python scripts/analysis/generalization_analysis.py
"""
import argparse
import glob
import json
import os
import sys
import types

import pandas as pd
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
M = os.path.join(REPO, "models")


# --------------------------------------------------------------------------
# Training sets
# --------------------------------------------------------------------------

def _fasta(path, maxlen=None):
    seqs, chunk = set(), []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if chunk: seqs.add("".join(chunk).upper()); chunk = []
            elif line:
                chunk.append(line)
    if chunk: seqs.add("".join(chunk).upper())
    return {s for s in seqs if maxlen is None or len(s) <= maxlen}


def _column(path, name, maxlen=None):
    df = pd.read_csv(path, low_memory=False)
    col = next(c for c in df.columns if c.lower() == name.lower())
    seqs = {str(v).strip().upper() for v in df[col].dropna()}
    return {s for s in seqs if maxlen is None or len(s) <= maxlen}


def training_sets():
    """model -> (n_train, sequences_to_match_against, derivation note).

    n_train counts the distinct real peptides the released checkpoint was fitted
    on. The match set can be larger where the exact subset is unrecorded, which
    is the conservative choice for a leakage check.
    """
    s = {}

    pos = _fasta(f"{M}/ampeppy/training_data/M_model_train_AMP_sequence.numbered.fasta")
    sub1 = _fasta(f"{M}/ampeppy/training_data/M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta")
    sub2 = _fasta(f"{M}/ampeppy/training_data/M_model_train_nonAMP_sequence.numbered.randomsubsample.fasta")
    s["ampeppy"] = (len(pos) + len(sub1), pos | sub1 | sub2,
                    "3,268 positives and a 3,268 negative subsample; which of the two "
                    "shipped subsamples produced the checkpoint is unrecorded, so both are matched")

    pos = _fasta(f"{M}/amplify/data/AMPlify_AMP_train_common.fa")
    neg = _fasta(f"{M}/amplify/data/AMPlify_non_AMP_train_balanced.fa")
    s["amplify"] = (len(pos | neg), pos | neg, "balanced model, which is the deployed one")

    tr = _fasta(f"{M}/ampscanner/training_data/AMP.tr.fa") | _fasta(f"{M}/ampscanner/training_data/DECOY.tr.fa")
    ev = _fasta(f"{M}/ampscanner/training_data/AMP.eval.fa") | _fasta(f"{M}/ampscanner/training_data/DECOY.eval.fa")
    s["ampscanner"] = (len(tr | ev), tr | ev,
                       "train and validation merged, as the shipped FULL_MODEL weights imply; test held out")

    pos = _fasta(f"{M}/ampredmfa/dataset/our_dataset/amps.fasta")
    neg = _fasta(f"{M}/ampredmfa/dataset/our_dataset/non_amps.fasta")
    s["ampredmfa"] = (len(pos | neg), pos | neg,
                      "full pool; the shipped checkpoint records neither trial nor fold")

    seqs = set()
    for name in ("train.txt", "valid.txt"):
        with open(f"{M}/ampredictor/data/{name}") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    seqs.add(parts[1].strip().upper())
    s["ampredictor"] = (len(seqs), seqs, "train plus validation, the latter driving early stopping")

    gn = _column(f"{M}/deep-amp/data/train_data/AMP/regressor_train_AMP_gr_neg.csv", "sequence", 48)
    gp = _column(f"{M}/deep-amp/data/train_data/AMP/regressor_train_AMP_gr_pos.csv", "sequence", 48)
    nn = _column(f"{M}/deep-amp/data/train_data/nonAMP/regressor_train_nonAMP.csv", "sequence", 48)
    s["deep-amp"] = (len(gn | gp | nn), gn | gp | nn,
                     "up to 48 residues (code/utils.py chopping); union of the four variants; "
                     "generator and toxicity data excluded")

    seqs = _column(f"{M}/mbc-attention/data/EC.csv", "SEQUENCE")
    s["mbc-attention"] = (len(seqs), seqs, "whole_train.mdl covers all of EC.csv")

    seqs = set()
    for name in ("cls_benchmark_balanced/train.csv", "regression_benchmark/E.coli/train.csv",
                 "regression_benchmark/S.aureus/train.csv"):
        seqs |= _column(f"{M}/sensexamp/datasets/ori_datasets/{name}", "Sequence")
    s["sensexamp"] = (len(seqs), seqs, "training splits of the three benchmarked variants")

    negpool = _column(f"{M}/hydramp-amp-classifier/data/unlabelled_negative.csv", "Sequence")
    amp_pos = _column(f"{M}/hydramp-amp-classifier/data/unlabelled_positive.csv", "Sequence", 25)
    mic_pos = _column(f"{M}/hydramp-mic-classifier/data/mic_data.csv", "sequence", 25)
    note = ("positives up to 25 residues; negatives are random fragments of UniProt sequences "
            "drawn at training time, so the full source pool is matched instead")
    s["hydramp-amp-classifier"] = (len(amp_pos), amp_pos | negpool, note)
    s["hydramp-mic-classifier"] = (len(mic_pos), mic_pos | negpool, note)

    apex = pd.read_csv(f"{M}/apex/training_data/APEXDB.csv")
    cv = {str(x).strip().upper() for x in apex[apex["Group (APEX training)"] == "CV"]["Sequence"].dropna()}
    s["apex"] = (len(cv), cv, "cross-validation group; the independent split is held out")

    return s


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------

def load_config():
    with open(os.path.join(REPO, "config/config.yaml")) as fh:
        return yaml.safe_load(fh)


def variant_index(cfg):
    """variant -> (model, type), mirroring the Snakefile's VARIANTS construction."""
    out = {}
    for model in cfg["models"]:
        with open(f"{M}/{model}/model.yaml") as fh:
            meta = yaml.safe_load(fh)
        if meta.get("variants"):
            for v in meta["variants"]:
                out[v["name"]] = (model, v.get("type", meta["type"]))
        else:
            out[model] = (model, meta["type"])
    return out


def evaluate(cfg, variant, vtype, task, labels_path, out_path):
    """Run workflow/scripts/evaluate.py with rule evaluate_task's parameters."""
    src_path = os.path.join(REPO, "workflow/scripts/evaluate.py")
    with open(src_path) as fh:
        src = fh.read()
    tc = dict(cfg["tasks"][task])
    tc["labels"] = labels_path
    ns = types.SimpleNamespace
    sm = ns(
        input=ns(predictions=f"results/inference/{variant}/{tc['dataset']}/predictions.tsv",
                 labels=labels_path, validation=None),
        output=ns(metrics=out_path),
        params=ns(task_config=tc, variant_type=vtype,
                  sequence_column=cfg.get("sequence_column", "sequence"),
                  benchmark_unit=cfg.get("benchmark_unit", "ug/ml"),
                  activity_thresholds=cfg["activity_thresholds"],
                  mic_clamp=cfg.get("mic_clamp")),
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    exec(compile(src, src_path, "exec"), {"__name__": "__evaluate__", "snakemake": sm})
    with open(out_path) as fh:
        return json.load(fh)["metrics"]


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(out_dir, work_dir, only_tasks=None):
    os.chdir(REPO)
    cfg = load_config()
    variants = variant_index(cfg)
    tsets = training_sets()

    holdout_of = {}          # battleamp task -> holdout counterpart
    for task in cfg["tasks"]:
        if task.startswith("holdout2026_"):
            holdout_of[task.replace("holdout2026_", "", 1)] = task

    base_tasks = [t for t in cfg["tasks"]
                  if not t.startswith("holdout2026_") and t != "example_classification"]
    if only_tasks:
        base_tasks = [t for t in base_tasks if t in only_tasks]

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    # Per-model training set summary.
    pd.DataFrame([
        {"model": m, "n_train": n, "n_matched_against": len(seqs), "derivation": note}
        for m, (n, seqs, note) in sorted(tsets.items())
    ]).to_csv(os.path.join(out_dir, "training_sets.tsv"), sep="\t", index=False)

    rows = []
    for task in base_tasks:
        tc = cfg["tasks"][task]
        req = tc.get("model_type")
        req = [req] if isinstance(req, str) else req
        labels = pd.read_csv(tc["labels"], sep="\t")
        labels["_u"] = labels["sequence"].astype(str).str.upper()
        keep = [c for c in labels.columns if c != "_u"]

        for variant, (model, vtype) in variants.items():
            if req is not None and vtype not in req:
                continue
            pred = f"results/inference/{variant}/{tc['dataset']}/predictions.tsv"
            if not os.path.exists(pred):
                continue

            n_train, train_seqs, _ = tsets.get(model, (None, set(), ""))
            overlap = labels["_u"].isin(train_seqs)

            base = dict(variant=variant, model=model, model_type=vtype, task=task,
                        n_train=n_train, overlap_n=int(overlap.sum()),
                        overlap_pct=round(100.0 * overlap.mean(), 2))

            full = evaluate(cfg, variant, vtype, task, tc["labels"],
                            os.path.join(work_dir, variant, task, "battleamp.json"))
            rows.append({**base, "evaluation_set": "battleamp", **full})

            lab_path = os.path.join(work_dir, variant, task, "minustrain_labels.tsv")
            os.makedirs(os.path.dirname(lab_path), exist_ok=True)
            labels.loc[~overlap, keep].to_csv(lab_path, sep="\t", index=False)
            mt = evaluate(cfg, variant, vtype, task, lab_path,
                          os.path.join(work_dir, variant, task, "minustrain.json"))
            rows.append({**base, "evaluation_set": "minustrain", **mt})

            htask = holdout_of.get(task)
            if htask:
                htc = cfg["tasks"][htask]
                hpred = f"results/inference/{variant}/{htc['dataset']}/predictions.tsv"
                if os.path.exists(hpred):
                    hm = evaluate(cfg, variant, vtype, htask, htc["labels"],
                                  os.path.join(work_dir, variant, task, "holdout.json"))
                    rows.append({**base, "evaluation_set": "holdout", **hm})
        print(f"  done {task}", file=sys.stderr)

    df = pd.DataFrame(rows)
    lead = ["variant", "model", "model_type", "task", "evaluation_set",
            "n_train", "overlap_n", "overlap_pct"]
    df = df[lead + [c for c in df.columns if c not in lead]]
    out = os.path.join(out_dir, "generalization_metrics.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}: {len(df)} rows, {df.task.nunique()} tasks, "
          f"{df.variant.nunique()} variants", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="results/generalization")
    ap.add_argument("--work-dir", default="results/generalization/_work")
    ap.add_argument("--tasks", nargs="*", default=None,
                    help="restrict to these task names (default: all non-holdout tasks)")
    a = ap.parse_args()
    main(a.out_dir, a.work_dir, a.tasks)
