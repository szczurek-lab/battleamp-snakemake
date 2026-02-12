#!/usr/bin/env python3
"""Generate the model table from models/registry.yaml and inject it into README.md.

Usage:
    python scripts/generate_model_table.py          # updates README.md in place
    python scripts/generate_model_table.py --check   # exit 1 if README is out of date (for CI)
    python scripts/generate_model_table.py --stdout   # print the table to stdout instead

The script looks for two marker comments in README.md:

    <!-- MODEL_TABLE_START -->
    <!-- MODEL_TABLE_END -->

Everything between them is replaced with the generated table.
"""

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "models" / "registry.yaml"
README = ROOT / "README.md"

START_MARKER = "<!-- MODEL_TABLE_START -->"
END_MARKER = "<!-- MODEL_TABLE_END -->"


def length_str(min_len, max_len):
    """Format the accepted length range as a human-readable string."""
    if min_len is None and max_len is None:
        return "Unlimited"
    if min_len is None:
        return f"<= {max_len}"
    if max_len is None:
        return f">= {min_len}"
    return f"{min_len} to {max_len}"


def generate_table(models):
    lines = []
    lines.append("| Model | Variants | Type | Framework | Accepted lengths | GPU |")
    lines.append("|-------|----------|------|-----------|------------------|-----|")

    classifiers = [m for m in models if m["type"] == "classifier"]
    regressors = [m for m in models if m["type"] == "regressor"]

    for group in [classifiers, regressors]:
        for m in group:
            name = m["name"]
            variants = ", ".join(m.get("variants", [])) if m.get("variants") else "(single)"
            mtype = m["type"]
            framework = m.get("framework", "")
            lengths = length_str(m.get("min_length"), m.get("max_length"))
            gpu = "yes" if m.get("gpu_required") else "no"
            lines.append(f"| {name} | {variants} | {mtype} | {framework} | {lengths} | {gpu} |")

    n_classifiers = len(classifiers)
    n_regressors = len(regressors)
    n_variants = sum(max(len(m.get("variants", [])), 1) for m in models)
    lines.append("")
    lines.append(f"Total: {n_classifiers} classifiers, {n_regressors} regressors, "
                 f"{n_variants} variants")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="Check if README is up to date; exit 1 if not")
    parser.add_argument("--stdout", action="store_true",
                        help="Print the table to stdout instead of updating README")
    args = parser.parse_args()

    if not REGISTRY.exists():
        print(f"Error: {REGISTRY} not found", file=sys.stderr)
        sys.exit(1)

    with open(REGISTRY) as f:
        data = yaml.safe_load(f)

    table = generate_table(data.get("models", []))

    if args.stdout:
        print(table)
        return

    if not README.exists():
        print(f"Error: {README} not found", file=sys.stderr)
        sys.exit(1)

    readme_text = README.read_text()

    if START_MARKER not in readme_text or END_MARKER not in readme_text:
        print(f"Error: could not find {START_MARKER} and {END_MARKER} in {README}",
              file=sys.stderr)
        sys.exit(1)

    before = readme_text.split(START_MARKER)[0]
    after = readme_text.split(END_MARKER)[1]
    new_text = f"{before}{START_MARKER}\n{table}\n{END_MARKER}{after}"

    if args.check:
        if new_text != readme_text:
            print("README.md model table is out of date. Run:")
            print("  python scripts/generate_model_table.py")
            sys.exit(1)
        else:
            print("README.md model table is up to date.")
            return

    README.write_text(new_text)
    print(f"Updated {README}")


if __name__ == "__main__":
    main()
