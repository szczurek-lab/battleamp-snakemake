"""Public API for the BattleAMP web service.

Every function takes and returns strings, so a web layer never has to touch the
filesystem or import pandas. Input is FASTA text; output is JSON text.

    >>> import battleamp
    >>> battleamp.list_models()                      # fast
    >>> battleamp.validate(fasta_text)               # fast
    >>> battleamp.score(fasta_text, models=["ampeppy"], unit="uM")   # SLOW

Timing matters for where you call these from:

    list_models()  milliseconds  -- safe in an HTTP request handler
    validate()     milliseconds  -- safe in an HTTP request handler
    score()        minutes-hours -- background worker ONLY

score() launches Snakemake, which builds conda environments and loads model
weights onto a GPU. Calling it from a request handler will exhaust any normal
gateway timeout. Run it from a task queue and poll for the result.

Two tiers of error reporting are returned:

    "messages"    short, already human-readable, safe to render to end users
    "diagnostics" Snakemake log paths and tails -- operators only, since they
                  contain absolute server paths and environment details
"""

import hashlib
import json
import os
import subprocess
from pathlib import Path

from . import aggregate, registry, sequences, units

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
UPLOADS_DIR = RESULTS_DIR / "uploads"

# Snakemake usually lives in a virtualenv that is not on the PATH of a service
# worker (on bury it is ~/.venvs/battleamp-snakemake/bin/snakemake). Point
# BATTLEAMP_SNAKEMAKE at the executable rather than relying on PATH.
SNAKEMAKE_BIN = os.environ.get("BATTLEAMP_SNAKEMAKE", "snakemake")

# Which Snakemake profile to run with: profile/ for a single machine,
# slurm/ to submit jobs to the cluster queue.
SNAKEMAKE_PROFILE = os.environ.get("BATTLEAMP_PROFILE", "profile/")


# ---------------------------------------------------------------------------
# Model catalogue
# ---------------------------------------------------------------------------

def list_models():
    """JSON describing the models a user can choose from.

    Use this to render the model-selection UI: it gives each model's type,
    accepted peptide length range, and the columns it will contribute.
    """
    variants = registry.load_variants()

    by_model = {}
    for name in sorted(variants):
        info = variants[name]
        model = info["model"]
        entry = by_model.setdefault(model, {
            "name": model,
            "type": info["type"],
            "framework": info.get("framework"),
            "gpu_required": info.get("gpu_required", False),
            "length_min": info.get("length_min"),
            "length_max": info.get("length_max"),
            "variants": [],
        })
        entry["variants"].append({"name": name, "type": info["type"]})

    # A single-variant model reports itself as its own variant; that is noise
    # in a UI, so flatten it away.
    for entry in by_model.values():
        if len(entry["variants"]) == 1 and entry["variants"][0]["name"] == entry["name"]:
            entry["variants"] = []

    ordered = [by_model[m] for m in registry.enabled_models() if m in by_model]
    return _dumps({"models": ordered, "units": list(units.VALID_UNITS)})


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def validate(fasta_text, models=None):
    """JSON report on whether a FASTA is usable, without running any model.

    Fast and side-effect free. Call this when the user uploads or pastes their
    sequences, before offering them a "run" button.

    Reports, per peptide, whether it survives cleaning, and per model, how many
    of the peptides fall inside that model's supported length range.
    """
    try:
        records = sequences.parse_fasta(fasta_text)
    except ValueError as e:
        return _dumps({
            "valid": False,
            "n_records": 0,
            "n_sequences": 0,
            "sequences": [],
            "rejected": [],
            "model_coverage": {},
            "messages": [_message("error", str(e))],
        })

    if not records:
        return _dumps({
            "valid": False,
            "n_records": 0,
            "n_sequences": 0,
            "sequences": [],
            "rejected": [],
            "model_coverage": {},
            "messages": [_message(
                "error",
                "No sequences found. Expected FASTA format, for example: "
                ">peptide_1 on one line, then GIGKFLHSAKKFGKAFVGEIMNS on the next.",
            )],
        })

    accepted, rejected = sequences.classify_records(records)

    for item in rejected:
        item["message"] = sequences.reason_text(item["reason"], item.get("detail"))

    messages = []
    by_reason = {}
    for item in rejected:
        by_reason.setdefault(item["reason"], []).append(item["id"])
    for reason, ids in by_reason.items():
        severity = "info" if reason == sequences.REASON_DUPLICATE else "warning"
        messages.append(_message(
            severity,
            f"{len(ids)} of {len(records)} sequences: "
            f"{sequences.reason_text(reason)}",
            ids=_sample_ids(ids),
        ))

    coverage = {}
    if accepted:
        variants = registry.load_variants(registry.resolve_models(models))
        for name in sorted(variants):
            info = variants[name]
            lo, hi = info.get("length_min"), info.get("length_max")
            too_short = sum(
                1 for seq, _ in accepted if lo is not None and len(seq) < lo
            )
            too_long = sum(
                1 for seq, _ in accepted if hi is not None and len(seq) > hi
            )
            coverage[name] = {
                "n_eligible": len(accepted) - too_short - too_long,
                "n_too_short": too_short,
                "n_too_long": too_long,
                "length_min": lo,
                "length_max": hi,
            }
            if too_short or too_long:
                messages.append(_message(
                    "info",
                    f"{name} accepts peptides of "
                    f"{_length_range(lo, hi)}; "
                    f"{too_short + too_long} of your {len(accepted)} sequences "
                    f"fall outside that range and will have no score from it.",
                    kind="coverage",
                ))

    if not accepted:
        messages.insert(0, _message(
            "error", "None of the submitted sequences can be scored."
        ))

    return _dumps({
        "valid": bool(accepted),
        "n_records": len(records),
        "n_sequences": len(accepted),
        "sequences": [
            {"ids": ids, "sequence": seq, "length": len(seq)}
            for seq, ids in accepted
        ],
        "rejected": rejected,
        "model_coverage": coverage,
        "messages": messages,
    })


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def score(fasta_text, models=None, unit="ug/ml", timeout=None, cores=None):
    """Run the selected models over the peptides and return a JSON table.

    SLOW -- minutes to hours. Run from a background worker, never from an HTTP
    request handler.

    Parameters
    ----------
    fasta_text : str    FASTA content.
    models : list or comma-separated str or None
        Models to run. None means every model enabled in config/config.yaml.
    unit : str          "ug/ml" or "uM". All MIC columns are converted to this.
    timeout : int or None   Seconds before Snakemake is killed.
    cores : int or None     Override the profile's core count.

    Results are cached: the dataset is keyed by the hash of its cleaned
    sequences, so re-running with more models only runs the new ones, and two
    users submitting identical peptides share the computation. The returned
    table always contains exactly the models requested in this call.
    """
    if unit not in units.VALID_UNITS:
        return _failure(
            f"Unknown unit {unit!r}. Choose one of: {', '.join(units.VALID_UNITS)}."
        )

    try:
        model_names = registry.resolve_models(models)
    except ValueError as e:
        return _failure(str(e))

    validation = json.loads(validate(fasta_text, models=model_names))
    if not validation["valid"]:
        return _dumps({
            "status": "failed",
            "unit": unit,
            "requested_models": model_names,
            "columns": [],
            "rows": [],
            "input": _input_summary(validation),
            "models": {},
            "messages": validation["messages"],
            "diagnostics": {},
        })

    accepted = [(s["sequence"], s["ids"]) for s in validation["sequences"]]
    dataset, fasta_path = stage_fasta(accepted)

    proc = _run_snakemake(fasta_path, model_names, timeout=timeout, cores=cores)

    return _dumps(build_report(
        dataset=dataset,
        model_names=model_names,
        unit=unit,
        validation=validation,
        accepted=accepted,
        proc=proc,
        timeout=timeout,
    ))


def build_report(dataset, model_names, unit, validation, accepted, proc,
                 timeout=None):
    """Assemble the result payload from pipeline output already on disk.

    Returns a dict (score() serialises it; the Snakefile's output_table rule writes
    parts of it to separate files).

    This is the single place where predictions become a user-facing answer.
    Both entry points go through it:

      score()            stages the FASTA, runs Snakemake, then calls this
      --config output=... Snakemake has already run; the output_table rule
                         handler calls this directly

    Keeping one implementation is what stops the two front doors from
    disagreeing about status, messages or column layout.
    """
    variants = registry.load_variants(model_names)
    columns, rows, model_report = aggregate.build_scores(
        dataset=dataset,
        variants=variants,
        results_dir=RESULTS_DIR,
        unit=unit,
        sequences=accepted,
    )

    # validate()'s coverage messages predict what will be missing; the
    # per-model messages below report what actually was. Keep only the latter.
    messages = [m for m in validation["messages"] if m.get("kind") != "coverage"]
    n_failed = sum(1 for r in model_report.values() if r["status"] == "failed")

    if proc["status"] == "timeout":
        messages.append(_message(
            "error",
            f"Scoring exceeded the {timeout} second time limit and was stopped. "
            f"Any models that finished are included below.",
            kind="run",
        ))
    elif proc["status"] == "not_found":
        messages.append(_message(
            "error",
            "Snakemake is not installed or not on PATH on the server.",
            kind="run",
        ))

    for variant, entry in sorted(model_report.items()):
        if entry["status"] == "failed":
            messages.append(_message(
                "warning",
                f"{variant} could not be run and has no scores. "
                f"The other models are unaffected.",
                kind="model",
            ))
        elif entry.get("n_skipped_length"):
            messages.append(_message(
                "info",
                f"{variant} scored {entry['n_scored']} of {len(rows)} peptides; "
                f"{entry['n_skipped_length']} are outside its supported length "
                f"range of {_length_range(entry.get('length_min'), entry.get('length_max'))}.",
                kind="model",
            ))

    if n_failed == len(model_report):
        status = "failed"
        messages.append(_message(
            "error", "No model produced any scores for this submission."
        ))
    elif n_failed or proc["status"] != "ok":
        status = "partial"
    else:
        status = "ok"

    return {
        "status": status,
        "dataset": dataset,
        "unit": unit,
        "requested_models": model_names,
        "columns": columns,
        "rows": rows,
        "input": _input_summary(validation),
        "models": model_report,
        "messages": messages,
        "diagnostics": {
            "snakemake": proc,
            "results_dir": str(RESULTS_DIR),
            "model_logs": {
                v: e.get("log") for v, e in model_report.items() if e.get("log")
            },
            "log_tails": {
                v: e.get("log_tail")
                for v, e in model_report.items() if e.get("log_tail")
            },
        },
    }


def to_tsv(result_json):
    """Convert the JSON returned by score() into TSV text for download."""
    result = json.loads(result_json)
    return aggregate.rows_to_tsv(result["columns"], result["rows"])


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def stage_fasta(accepted):
    """Write cleaned sequences to a content-addressed FASTA. Returns (name, path).

    The dataset name is the hash of the cleaned sequences, not the uploaded
    filename. The pipeline derives the dataset name from the file stem, so
    hashing is what keeps two users who upload different files both called
    'peptides.fasta' from overwriting each other's results -- and lets two users
    who submit identical peptides share the cache.
    """
    body = sequences.to_fasta(accepted)
    digest = hashlib.sha1(body.encode()).hexdigest()[:16]
    dataset = f"upload_{digest}"
    UPLOADS_DIR.mkdir(parents=True, exist_ok=True)
    path = UPLOADS_DIR / f"{dataset}.fasta"
    if not path.exists():
        path.write_text(body)
    return dataset, path


def _run_snakemake(fasta_path, model_names, timeout=None, cores=None):
    """Invoke the pipeline's inference-only target. Never raises."""
    cmd = [
        SNAKEMAKE_BIN,
        "--profile", SNAKEMAKE_PROFILE,
        "score",
        "--config",
        f"fasta={fasta_path}",
        f"run_models={','.join(model_names)}",
    ]
    if cores is not None:
        cmd += ["--cores", str(cores)]

    try:
        completed = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError:
        return {"status": "not_found", "command": " ".join(cmd),
                "returncode": None, "stderr_tail": []}
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "command": " ".join(cmd),
                "returncode": None, "stderr_tail": []}

    return {
        # A non-zero exit is expected whenever any single model fails, because
        # the profile sets keep-going. It is not an error for the run overall.
        "status": "ok" if completed.returncode == 0 else "partial",
        "command": " ".join(cmd),
        "returncode": completed.returncode,
        "stderr_tail": completed.stderr.splitlines()[-aggregate.LOG_TAIL_LINES:],
    }


def _input_summary(validation):
    return {
        "n_records": validation["n_records"],
        "n_sequences": validation["n_sequences"],
        "rejected": validation["rejected"],
    }


def _length_range(lo, hi):
    if lo is None and hi is None:
        return "any length"
    if lo is None:
        return f"up to {hi} aa"
    if hi is None:
        return f"at least {lo} aa"
    return f"{lo}-{hi} aa"


def _sample_ids(ids, limit=10):
    """Keep messages short: show a handful of ids, not thousands."""
    if len(ids) <= limit:
        return ids
    return ids[:limit] + [f"... and {len(ids) - limit} more"]


def _message(severity, text, ids=None, kind="input"):
    """Build a user-facing message.

    ``kind`` lets a front end group messages, and lets score() drop validate()'s
    predictive coverage warnings in favour of its own factual ones.
    Values: input | coverage | model | run.
    """
    msg = {"severity": severity, "kind": kind, "text": text}
    if ids:
        msg["ids"] = ids
    return msg


def _failure(text):
    return _dumps({
        "status": "failed",
        "columns": [], "rows": [], "models": {},
        "input": {"n_records": 0, "n_sequences": 0, "rejected": []},
        "messages": [_message("error", text)],
        "diagnostics": {},
    })


def _dumps(obj):
    return json.dumps(obj, indent=2, sort_keys=False)
