"""Build a per-peptide score table from per-model pipeline outputs.

The pipeline writes one file per model variant:

    results/inference/{variant}/{dataset}/predictions.tsv

keyed on sequence, with a schema that depends on the model type:

    classifier -> sequence, Prediction, Probability_score
    regressor  -> sequence, MIC, MIC_unit

This module transposes those into one row per peptide, one column per variant,
and reports why any cell is empty.

It deliberately has no Snakemake dependency and never raises on a missing or
malformed model output: partial failure is the normal case (the pipeline runs
with keep-going), and the service must still return a table.
"""

import csv
import json
from pathlib import Path

from . import units

# How many lines of a failed model's log to keep for the operator.
LOG_TAIL_LINES = 20


def _read_predictions(path):
    """Read a predictions TSV into {sequence: row_dict}. Never raises.

    Returns (rows, error) where error is None on success or a short string
    describing why the file could not be used.
    """
    try:
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                return {}, "prediction file is empty"
            if "sequence" not in reader.fieldnames:
                return {}, (
                    f"prediction file has no 'sequence' column "
                    f"(found: {', '.join(reader.fieldnames)})"
                )
            rows = {}
            for row in reader:
                seq = (row.get("sequence") or "").strip()
                if seq:
                    rows[seq] = row
            return rows, None
    except OSError as e:
        return {}, f"could not read prediction file: {e}"


def _read_skipped_lengths(prefiltered_dir):
    """Sequences this variant skipped for length, from the prefilter manifest.

    Returns a set of sequences, empty if the manifest is absent.
    """
    manifest = Path(prefiltered_dir) / "manifest.tsv"
    skipped = set()
    if not manifest.exists():
        return skipped
    try:
        with open(manifest, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row.get("status") == "skipped_length":
                    seq = (row.get("sequence") or "").strip()
                    if seq:
                        skipped.add(seq)
    except OSError:
        pass
    return skipped


def _log_tail(log_path):
    """Last few lines of a log file, for operator diagnostics."""
    try:
        with open(log_path, errors="replace") as f:
            lines = f.read().splitlines()
        return lines[-LOG_TAIL_LINES:]
    except OSError:
        return []


def _to_float(raw):
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def build_scores(dataset, variants, results_dir, unit="ug/ml", sequences=None):
    """Assemble the per-peptide table.

    Parameters
    ----------
    dataset : str
        Dataset name used by the pipeline (the FASTA stem it was run under).
    variants : dict
        {variant_name: info} as returned by registry.load_variants.
    results_dir : str or Path
        The pipeline's results/ directory.
    unit : str
        "ug/ml" or "uM". All regressor columns are converted to this.
    sequences : list of (sequence, [ids]) or None
        The peptides to emit, in order. Rows are produced for every entry even
        if no model scored it. If None, the union of sequences seen across
        model outputs is used instead.

    Returns
    -------
    (columns, rows, model_report)
        columns     -- list of column descriptors (name, variant, kind, unit)
        rows        -- list of dicts, keys are 'id', 'sequence' and column names
        model_report-- {variant: {status, ...}} explaining every variant
    """
    if unit not in units.VALID_UNITS:
        raise ValueError(
            f"Unknown unit {unit!r}. Valid units: {', '.join(units.VALID_UNITS)}."
        )

    results_dir = Path(results_dir)
    suffix = units.unit_suffix(unit)

    columns = []
    model_report = {}
    per_variant_values = {}

    for variant in sorted(variants):
        info = variants[variant]
        pred_path = results_dir / "inference" / variant / dataset / "predictions.tsv"
        log_path = (
            results_dir / "logs" / "inference"
            / info.get("log_name", variant) / f"{dataset}.log"
        )
        prefiltered = (
            results_dir / "prefiltered"
            / info.get("prefilter_variant", variant) / dataset
        )

        is_regressor = info["type"] == "regressor"
        column = (
            f"{variant}_MIC_{suffix}" if is_regressor else f"{variant}_prob"
        )
        columns.append({
            "name": column,
            "variant": variant,
            "model": info["model"],
            "kind": "mic" if is_regressor else "probability",
            "unit": unit if is_regressor else None,
            "length_min": info.get("length_min"),
            "length_max": info.get("length_max"),
        })

        if not pred_path.exists():
            model_report[variant] = {
                "model": info["model"],
                "type": info["type"],
                "column": column,
                "status": "failed",
                "reason": "the model did not produce an output file",
                "log": str(log_path) if log_path.exists() else None,
                "log_tail": _log_tail(log_path) if log_path.exists() else [],
                "n_scored": 0,
            }
            per_variant_values[variant] = {}
            continue

        rows, read_error = _read_predictions(pred_path)
        if read_error is not None:
            model_report[variant] = {
                "model": info["model"],
                "type": info["type"],
                "column": column,
                "status": "failed",
                "reason": read_error,
                "log": str(log_path) if log_path.exists() else None,
                "log_tail": _log_tail(log_path) if log_path.exists() else [],
                "n_scored": 0,
            }
            per_variant_values[variant] = {}
            continue

        values = {}
        n_unparseable = 0
        for seq, row in rows.items():
            if is_regressor:
                raw = _to_float(row.get("MIC"))
                if raw is None:
                    n_unparseable += 1
                    continue
                src_unit = (row.get("MIC_unit") or unit).strip()
                if src_unit not in units.VALID_UNITS:
                    # Unknown unit: keep the number rather than silently
                    # mis-converting it, and say so in the report.
                    n_unparseable += 1
                    continue
                values[seq] = units.convert(raw, seq, src_unit, unit)
            else:
                raw = _to_float(row.get("Probability_score"))
                if raw is None:
                    n_unparseable += 1
                    continue
                values[seq] = raw

        per_variant_values[variant] = values
        skipped = _read_skipped_lengths(prefiltered)

        entry = {
            "model": info["model"],
            "type": info["type"],
            "column": column,
            "status": "ok",
            "n_scored": len(values),
            "log": str(log_path) if log_path.exists() else None,
        }
        if skipped:
            entry["status"] = "partial"
            entry["n_skipped_length"] = len(skipped)
            entry["length_min"] = info.get("length_min")
            entry["length_max"] = info.get("length_max")
        if n_unparseable:
            entry["status"] = "partial"
            entry["n_unparseable"] = n_unparseable
        model_report[variant] = entry

    if sequences is None:
        seen = []
        seen_set = set()
        for values in per_variant_values.values():
            for seq in values:
                if seq not in seen_set:
                    seen_set.add(seq)
                    seen.append(seq)
        sequences = [(seq, []) for seq in seen]

    rows_out = []
    for sequence, ids in sequences:
        row = {"id": ";".join(ids), "sequence": sequence}
        for col in columns:
            row[col["name"]] = per_variant_values[col["variant"]].get(sequence)
        rows_out.append(row)

    return columns, rows_out, model_report


def rows_to_tsv(columns, rows):
    """Render the table as TSV text. Empty cells for unscored peptides."""
    header = ["id", "sequence"] + [c["name"] for c in columns]
    lines = ["\t".join(header)]
    for row in rows:
        cells = [row.get("id", ""), row["sequence"]]
        for col in columns:
            value = row.get(col["name"])
            cells.append("" if value is None else _format_number(value))
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


def _format_number(value):
    """Trim float noise without losing precision that matters for MIC."""
    return f"{value:.6g}"


def read_cleaning_report(dataset, results_dir):
    """Counts from the pipeline's own cleaning step, or None if absent."""
    path = Path(results_dir) / "cleaned" / dataset / "cleaning_report.json"
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
