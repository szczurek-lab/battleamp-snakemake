"""
Rule: sequence pre-filtering
Creates a model-specific FASTA file containing only sequences within the
model's declared length range (length_min, length_max from model.yaml).
"""


rule prefilter_sequences:
    """Filter dataset FASTA to only sequences the model can handle."""
    input:
        fasta = lambda wc: get_dataset_fasta(wc.dataset),
    output:
        fasta = f"{OUTPUT_DIR}/prefiltered/{{variant}}/{{dataset}}/sequences.fasta",
        manifest = f"{OUTPUT_DIR}/prefiltered/{{variant}}/{{dataset}}/manifest.tsv",
    params:
        length_min = lambda wc: VARIANTS[wc.variant].get("length_min"),
        length_max = lambda wc: VARIANTS[wc.variant].get("length_max"),
    log:
        f"{OUTPUT_DIR}/logs/prefilter/{{variant}}/{{dataset}}.log",
    script:
        "../scripts/prefilter_sequences.py"
