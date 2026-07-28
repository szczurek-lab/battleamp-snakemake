#!/usr/bin/env python3
"""Score a FASTA with the BattleAMP models and write a per-peptide table.

This is a thin command-line front end to the battleamp package. The web service
should import battleamp directly rather than shelling out to this script; the
script exists so the pipeline can be driven and tested from a terminal.

    # check the input without running anything (fast)
    python scripts/score_fasta.py --fasta peptides.fasta --validate-only

    # list the models available for selection (fast)
    python scripts/score_fasta.py --list-models

    # run selected models and write results (slow)
    python scripts/score_fasta.py --fasta peptides.fasta \
        --models ampeppy,amplify --unit uM --out-dir jobs/123

Writes scores.tsv and result.json to --out-dir.

Exit codes:
    0  at least one model produced scores
    1  the input FASTA is unusable
    2  every requested model failed
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import battleamp  # noqa: E402


def main():
    parser = argparse.ArgumentParser(
        description="Score peptides with the BattleAMP models.",
    )
    parser.add_argument("--fasta", help="Input FASTA file")
    parser.add_argument(
        "--models",
        help="Comma-separated model names. Default: all models in config.yaml.",
    )
    parser.add_argument(
        "--unit", default="ug/ml", choices=["ug/ml", "uM"],
        help="Unit for MIC columns (default: ug/ml)",
    )
    parser.add_argument("--out-dir", help="Directory for scores.tsv and result.json")
    parser.add_argument(
        "--validate-only", action="store_true",
        help="Check the input and report model coverage without running models",
    )
    parser.add_argument(
        "--list-models", action="store_true",
        help="Print the available models as JSON and exit",
    )
    parser.add_argument("--timeout", type=int, help="Seconds before giving up")
    parser.add_argument("--cores", type=int, help="Override the profile core count")
    args = parser.parse_args()

    if args.list_models:
        print(battleamp.list_models())
        return 0

    if not args.fasta:
        parser.error("--fasta is required unless --list-models is given")

    fasta_text = Path(args.fasta).read_text()

    if args.validate_only:
        report = json.loads(battleamp.validate(fasta_text, models=args.models))
        print(battleamp.validate(fasta_text, models=args.models))
        _print_messages(report["messages"])
        return 0 if report["valid"] else 1

    result_json = battleamp.score(
        fasta_text,
        models=args.models,
        unit=args.unit,
        timeout=args.timeout,
        cores=args.cores,
    )
    result = json.loads(result_json)

    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "result.json").write_text(result_json)
        (out_dir / "scores.tsv").write_text(battleamp.to_tsv(result_json))
        print(f"Wrote {out_dir / 'scores.tsv'} and {out_dir / 'result.json'}",
              file=sys.stderr)
    else:
        print(battleamp.to_tsv(result_json))

    _print_summary(result)

    if result["status"] == "failed":
        return 2 if result["rows"] else 1
    return 0


def _print_messages(messages):
    for m in messages:
        print(f"  [{m['severity']}] {m['text']}", file=sys.stderr)


def _print_summary(result):
    print(
        f"\nStatus: {result['status']}  "
        f"({len(result['rows'])} peptides, {len(result['columns'])} model columns)",
        file=sys.stderr,
    )
    for variant, entry in sorted(result.get("models", {}).items()):
        print(f"  {entry['status']:>8}  {variant}  "
              f"(scored {entry.get('n_scored', 0)})", file=sys.stderr)
    _print_messages(result["messages"])


if __name__ == "__main__":
    sys.exit(main())
