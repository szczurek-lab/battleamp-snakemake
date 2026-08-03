#!/usr/bin/env python3
"""Guard the two deliberate duplications in the battleamp package.

The web API re-implements two things the pipeline already does, for good
reasons (speed, and avoiding a scikit-learn dependency in the web server). Both
duplications are only safe if they cannot drift:

  1. Input validation. If battleamp.sequences accepts a peptide that
     clean_dataset.py later drops, the user is told their input is fine and
     then silently gets no row for it.
  2. Molecular weights. If the two AA_MW tables diverge, MIC values shown to
     users disagree with the values in the published benchmark.

Run standalone (no pytest needed):

    python3 tests/test_parity.py
"""

import ast
import json
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from battleamp import sequences, units  # noqa: E402

# Sequences chosen to hit every branch of both implementations.
FIXTURE = [
    ("ok_plain",       "GIGKFLHSAKKFGKAFVGEIMNS"),
    ("ok_short",       "KWK"),
    ("ok_long",        "A" * 120),
    ("ok_all_aa",      "ACDEFGHIKLMNPQRSTVWY"),
    ("dup_of_plain",   "GIGKFLHSAKKFGKAFVGEIMNS"),
    ("bad_lowercase",  "GIGkFLHSAKKF"),
    ("bad_all_lower",  "gigkflhsakkf"),
    ("bad_x",          "GIGKFLHSAXKF"),
    ("bad_b",          "GIGKFLHSABKF"),
    ("bad_star",       "GIGKFLHSA*KF"),
    ("bad_gap",        "GIGK-FLHSAKKF"),
    ("bad_dot",        "GIGK.FLHSAKKF"),
    ("bad_empty",      ""),
]


def _fixture_fasta():
    return "".join(f">{name}\n{seq}\n" for name, seq in FIXTURE)


class Skip(Exception):
    """Raised when a test cannot run in this environment."""


def test_validation_matches_pipeline():
    """battleamp.sequences must keep exactly what clean_dataset.py keeps."""
    try:
        import clean_dataset
    except ModuleNotFoundError as e:
        # clean_dataset.py needs biopython, which lives in the pipeline's
        # virtualenv rather than the system interpreter on some machines.
        raise Skip(
            f"{e.name} is not installed; run this with the interpreter that "
            f"runs the pipeline (e.g. ~/.venvs/battleamp-snakemake/bin/python)"
        ) from None

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        src = tmp / "in.fasta"
        src.write_text(_fixture_fasta())
        out = tmp / "clean.fasta"
        clean_dataset.clean_dataset(str(src), str(out), str(tmp / "report.json"))

        pipeline_kept = [
            str(r.seq) for r in __import__("Bio.SeqIO", fromlist=["SeqIO"]).parse(
                str(out), "fasta"
            )
        ]

    accepted, _ = sequences.classify_records(sequences.parse_fasta(_fixture_fasta()))
    api_kept = [seq for seq, _ids in accepted]

    assert api_kept == pipeline_kept, (
        f"validation drift.\n"
        f"  web API keeps: {api_kept}\n"
        f"  pipeline keeps: {pipeline_kept}"
    )


def test_molecular_weights_match_evaluate():
    """units.AA_MW must equal the table evaluate.py uses for the benchmark.

    evaluate.py calls main(snakemake) at module level so it cannot be imported;
    the constants are read out of its source instead.
    """
    source = (REPO_ROOT / "workflow" / "scripts" / "evaluate.py").read_text()
    tree = ast.parse(source)

    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in ("AA_MW", "WATER_MW"):
                found[target.id] = ast.literal_eval(node.value)

    assert "AA_MW" in found, "could not find AA_MW in evaluate.py"
    assert "WATER_MW" in found, "could not find WATER_MW in evaluate.py"
    assert found["AA_MW"] == units.AA_MW, "AA_MW drifted from evaluate.py"
    assert found["WATER_MW"] == units.WATER_MW, "WATER_MW drifted from evaluate.py"


def test_unit_round_trip():
    """Converting out and back must return the original value."""
    for seq in ("GIGKFLHSAKKFGKAFVGEIMNS", "KWK", "ACDEFGHIKLMNPQRSTVWY"):
        for value in (0.25, 32.0, 512.0):
            there = units.convert(value, seq, "ug/ml", "uM")
            back = units.convert(there, seq, "uM", "ug/ml")
            assert abs(back - value) < 1e-9, f"{seq} {value} -> {there} -> {back}"


def test_duplicate_ids_are_preserved():
    """A sequence submitted twice keeps both headers on its single row."""
    accepted, rejected = sequences.classify_records(
        sequences.parse_fasta(_fixture_fasta())
    )
    by_seq = {seq: ids for seq, ids in accepted}
    assert by_seq["GIGKFLHSAKKFGKAFVGEIMNS"] == ["ok_plain", "dup_of_plain"]
    assert any(r["reason"] == sequences.REASON_DUPLICATE for r in rejected)


def test_multioutput_artefact_paths_match_the_pipeline():
    """Multi-output models write two artefacts under a name that is not the variant.

    workflow/rules/inference.smk writes the prefilter manifest under the FIRST
    variant (whose FASTA the others reuse) and the inference log under the MODEL
    name. Reading either from the variant path yields nothing: length-exclusion
    counts vanish and a failed model reports no log to its operator.
    """
    from battleamp import registry

    variants = registry.load_variants(["apex", "sensexamp", "ampeppy"])

    # apex: 6 variants, one inference run, log under "apex", manifest under "apex-min"
    for name in ("apex-ecoli", "apex-saureus", "apex-min"):
        assert variants[name]["log_name"] == "apex", name
        assert variants[name]["prefilter_variant"] == "apex-min", name

    # sensexamp: mixed classifier/regressor variants, same rule
    for name in ("sensexamp-classifier", "sensexamp-ecoli", "sensexamp-saureus"):
        assert variants[name]["log_name"] == "sensexamp", name
        assert variants[name]["prefilter_variant"] == "sensexamp-classifier", name

    # single-variant model: everything lives under its own name
    assert variants["ampeppy"]["log_name"] == "ampeppy"
    assert variants["ampeppy"]["prefilter_variant"] == "ampeppy"


def test_delimiter_and_report_path_follow_the_filename():
    """output=out.csv must produce commas; anything else stays tab-separated."""
    from battleamp import aggregate

    assert aggregate.delimiter_for("/tmp/out.csv") == ","
    assert aggregate.delimiter_for("/tmp/out.CSV") == ","
    assert aggregate.delimiter_for("/tmp/out.tsv") == "\t"
    assert aggregate.delimiter_for("/tmp/out.txt") == "\t"
    assert aggregate.delimiter_for("/tmp/out") == "\t"

    assert str(aggregate.report_path_for("/tmp/out.csv")) == "/tmp/out.report.json"
    assert str(aggregate.report_path_for("/tmp/out.tsv")) == "/tmp/out.report.json"
    assert str(aggregate.report_path_for("/tmp/out")) == "/tmp/out.report.json"


def test_csv_and_tsv_render_the_same_table():
    """Only the separator may differ between the two renderings."""
    from battleamp import aggregate

    columns = [{"name": "m_prob"}, {"name": "r_MIC_ugml"}]
    rows = [
        {"id": "a;b", "sequence": "KWK", "m_prob": 0.5, "r_MIC_ugml": None},
        {"id": "c", "sequence": "GIGK", "m_prob": None, "r_MIC_ugml": 12.25},
    ]
    csv = aggregate.rows_to_delimited(columns, rows, delimiter=",")
    tsv = aggregate.rows_to_tsv(columns, rows)

    assert csv.replace(",", "\t") == tsv
    # Unscored cells must be empty, not "None" or "nan".
    assert csv.splitlines()[1] == "a;b,KWK,0.5,"
    assert csv.splitlines()[2] == "c,GIGK,,12.25"


def test_report_counts_rejected_records_from_the_raw_input():
    """The report must describe the user's file, not the cleaned staging copy.

    output= replaces config["fasta"] with a content-addressed copy of the *accepted*
    sequences. Building the report from that copy silently reports zero
    rejections, hiding the peptides the user most needs to hear about.
    """
    from battleamp import api

    raw = (
        ">good\nGIGKFLHSAKKFGKAFVGEIMNS\n"
        ">dup\nGIGKFLHSAKKFGKAFVGEIMNS\n"
        ">bad_daa\nGIGkFLHSAKKF\n"
        ">bad_x\nGIGKFLHSAXKF\n"
    )
    validation = json.loads(api.validate(raw, models=["example-model"]))
    assert validation["n_records"] == 4, validation["n_records"]
    assert validation["n_sequences"] == 1, validation["n_sequences"]
    assert len(validation["rejected"]) == 3

    # The duplicate keeps its header on the surviving row.
    assert validation["sequences"][0]["ids"] == ["good", "dup"]

    # Staging deliberately drops all of that -- which is why the report must be
    # built from the raw text.
    accepted = [(s["sequence"], s["ids"]) for s in validation["sequences"]]
    staged = sequences.to_fasta(accepted)
    staged_validation = json.loads(api.validate(staged, models=["example-model"]))
    assert staged_validation["n_records"] == 1
    assert staged_validation["rejected"] == []


def test_headerless_input_is_rejected_clearly():
    """A pasted bare peptide list must produce an actionable error."""
    try:
        sequences.parse_fasta("GIGKFLHSAKKFGKAFVGEIMNS\nKWKLFKKIEK\n")
    except ValueError as e:
        assert "header" in str(e).lower()
    else:
        raise AssertionError("headerless input should raise ValueError")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failures = skipped = 0
    for test in tests:
        try:
            test()
        except Skip as e:
            skipped += 1
            print(f"SKIP  {test.__name__}: {e}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL  {test.__name__}\n      {e}")
        except Exception as e:  # noqa: BLE001
            failures += 1
            print(f"ERROR {test.__name__}: {type(e).__name__}: {e}")
        else:
            print(f"ok    {test.__name__}")
    summary = f"\n{len(tests) - failures - skipped}/{len(tests)} passed"
    if skipped:
        summary += f", {skipped} skipped"
    print(summary)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
