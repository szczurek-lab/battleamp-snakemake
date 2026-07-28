#!/usr/bin/env python3
"""Regenerate the reference request/response examples in examples/.

These files are the contract handed to whoever builds the web front end: they
show the exact JSON shape each API function returns, so the UI can be built and
tested before the cluster is wired up.

The score example is built from real model predictions already on disk
(results/inference/.../example-dataset/) rather than from invented numbers, and
it is assembled by the same code path score() uses, so it cannot drift from the
real response shape.

Usage:
    python scripts/generate_api_examples.py
"""

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import battleamp  # noqa: E402
from battleamp import api, sequences  # noqa: E402

EXAMPLES = REPO_ROOT / "examples"

# A mix that exercises every column kind and both failure modes:
#   example-model  classifier, always succeeds
#   apex           regressor, multi-output (6 columns from one run)
#   sensexamp      6-25 aa limit -> partial coverage
#   deep-amp       reports MIC in uM -> exercises unit conversion
#   mbc-attention  has no predictions for this dataset -> demonstrates a
#                  failed model, so the UI has a failure state to render
DEMO_MODELS = ["example-model", "apex", "sensexamp", "deep-amp", "mbc-attention"]
DEMO_DATASET = "example-dataset"


def write(name, text):
    path = EXAMPLES / name
    path.write_text(text)
    print(f"  wrote {path.relative_to(REPO_ROOT)}  ({len(text):,} bytes)")
    return path


def main():
    EXAMPLES.mkdir(exist_ok=True)
    print("Generating API examples...")

    # --- input -------------------------------------------------------------
    base = (REPO_ROOT / "datasets" / DEMO_DATASET / "sequences.fasta").read_text()
    write("example_input.fasta", base)

    # A second input that deliberately contains bad records, so the front end
    # has something to render its warning states against.
    messy = base.rstrip("\n") + "\n" + (
        ">duplicate_of_seq_001\nGIGKFLHSAKKFGKAFVGEIMNS\n"
        ">has_d_amino_acids\nGIGkFLHSAKKFGKAF\n"
        ">has_unknown_residue\nGIGKFLHSAXKFGKAF\n"
        ">has_a_gap\nGIGK-FLHSAKKFGKAF\n"
        ">far_too_long\n" + "A" * 300 + "\n"
    )
    write("example_input_with_errors.fasta", messy)

    # --- list_models -------------------------------------------------------
    write("example_list_models.json", battleamp.list_models())

    # --- validate ----------------------------------------------------------
    write("example_validate.json",
          battleamp.validate(base, models=DEMO_MODELS))
    write("example_validate_with_errors.json",
          battleamp.validate(messy, models=DEMO_MODELS))

    # --- score -------------------------------------------------------------
    # Assembled from predictions already on disk for DEMO_DATASET.
    validation = json.loads(battleamp.validate(base, models=DEMO_MODELS))
    accepted = [(s["sequence"], s["ids"]) for s in validation["sequences"]]
    proc = {
        "status": "partial",
        "command": (
            "snakemake --profile profile/ score "
            "--config fasta=results/uploads/upload_<hash>.fasta "
            f"run_models={','.join(DEMO_MODELS)}"
        ),
        "returncode": 1,
        "stderr_tail": ["(example: one model failed, keep-going returned 1)"],
    }

    for unit, suffix in (("ug/ml", ""), ("uM", "_uM")):
        result = api._assemble_result(
            dataset=DEMO_DATASET,
            model_names=DEMO_MODELS,
            unit=unit,
            validation=validation,
            accepted=accepted,
            proc=proc,
        )
        write(f"example_score{suffix}.json", result)
        write(f"example_score{suffix}.tsv", battleamp.to_tsv(result))

    print("\nDone. See examples/README.md for what each file demonstrates.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
