"""
Rule: clean_dataset
Runs once per dataset before any model-specific prefiltering.
Removes sequences with d-amino acids (lowercase), non-standard amino acids,
whitespace/gaps, empty sequences, and duplicates.
"""


rule clean_dataset_fasta:
    """Clean a raw dataset FASTA for use in the pipeline."""
    input:
        sequences=lambda wc: get_dataset_fasta(wc.dataset),
    output:
        sequences=f"{OUTPUT_DIR}/cleaned/{{dataset}}/sequences.fasta",
        report=f"{OUTPUT_DIR}/cleaned/{{dataset}}/cleaning_report.json",
    log:
        f"{OUTPUT_DIR}/logs/cleaning/{{dataset}}.log",
    script:
        "../scripts/clean_dataset.py"
