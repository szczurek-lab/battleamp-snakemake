"""FASTA parsing and per-record input validation.

The rules here mirror ``workflow/scripts/clean_dataset.py``, which is what the
pipeline itself applies before any model sees a sequence. Keeping a second
implementation lets the web service tell a user *up front* which of their
peptides will be rejected, without spinning up Snakemake.

The two implementations must not drift: if the web service accepts a peptide
that the pipeline later drops, the user gets a silently missing row.
``tests/test_validation_parity.py`` asserts they agree.
"""

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Reason codes are part of the public JSON payload. Keep them stable.
REASON_EMPTY = "empty"
REASON_WHITESPACE = "whitespace_or_gap"
REASON_LOWERCASE = "lowercase_d_amino_acid"
REASON_NONSTANDARD = "non_standard_amino_acid"
REASON_DUPLICATE = "duplicate_sequence"

_REASON_TEXT = {
    REASON_EMPTY: "Sequence is empty.",
    REASON_WHITESPACE: "Sequence contains whitespace or gap characters "
                       "(space, tab, '-' or '.').",
    REASON_LOWERCASE: "Sequence contains lowercase letters, which denote "
                      "D-amino acids. No integrated model supports them.",
    REASON_NONSTANDARD: "Sequence contains characters outside the 20 standard "
                        "amino acids.",
    REASON_DUPLICATE: "Sequence is identical to an earlier record and is "
                      "scored only once.",
}


def reason_text(reason, detail=None):
    """Human-readable explanation for a reason code, safe to show a user."""
    text = _REASON_TEXT.get(reason, reason)
    if detail:
        return f"{text} ({detail})"
    return text


def parse_fasta(text):
    """Parse FASTA text into a list of (record_id, sequence) tuples.

    Sequence lines are concatenated and stripped of surrounding whitespace, but
    internal whitespace is preserved so that ``classify_sequence`` can report it
    rather than silently repairing it.

    Raises ValueError if the text contains non-comment content before the first
    header, which is the usual symptom of a user pasting a bare peptide list.
    """
    records = []
    header = None
    seq_lines = []

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        else:
            if header is None:
                raise ValueError(
                    f"Line {lineno} contains sequence data before any '>' header. "
                    f"Every sequence must be preceded by a header line, "
                    f"e.g. '>peptide_1'."
                )
            seq_lines.append(line)

    if header is not None:
        records.append((header, "".join(seq_lines)))

    return records


def classify_sequence(sequence, seen=None):
    """Decide whether a single sequence survives cleaning.

    Returns ``(ok, reason, detail)``. ``reason`` is None when ok is True.
    ``seen`` is a set of sequences already accepted; pass it to detect
    duplicates, or omit it to check a sequence in isolation.

    Order matches clean_dataset.py exactly: empty, whitespace, lowercase,
    non-standard, duplicate.
    """
    if len(sequence) == 0:
        return False, REASON_EMPTY, None

    if any(c in sequence for c in (" ", "\t", "-", ".")):
        return False, REASON_WHITESPACE, None

    if any(c.islower() for c in sequence):
        return False, REASON_LOWERCASE, None

    if not all(c in STANDARD_AA for c in sequence):
        bad = sorted(set(sequence) - STANDARD_AA)
        return False, REASON_NONSTANDARD, "found: " + ", ".join(bad)

    if seen is not None and sequence in seen:
        return False, REASON_DUPLICATE, None

    return True, None, None


def classify_records(records):
    """Apply cleaning rules to parsed records, in file order.

    Returns ``(accepted, rejected)`` where:
      accepted -- list of (sequence, [record_id, ...]) in first-seen order.
                  A sequence submitted under several headers appears once, with
                  all of its ids, matching the pipeline's dedup-by-sequence.
      rejected -- list of dicts with keys id, sequence, reason, detail.
    """
    seen = {}
    accepted = []
    rejected = []

    for record_id, sequence in records:
        ok, reason, detail = classify_sequence(sequence, seen=seen)
        if ok:
            seen[sequence] = len(accepted)
            accepted.append((sequence, [record_id]))
        elif reason == REASON_DUPLICATE:
            # Not an error: attach this id to the sequence already accepted so
            # the user still sees their header in the results.
            accepted[seen[sequence]][1].append(record_id)
            rejected.append({
                "id": record_id,
                "sequence": sequence,
                "reason": reason,
                "detail": detail,
            })
        else:
            rejected.append({
                "id": record_id,
                "sequence": sequence,
                "reason": reason,
                "detail": detail,
            })

    return accepted, rejected


def to_fasta(accepted):
    """Render accepted (sequence, ids) pairs back to FASTA text.

    The first id is used as the header, so the FASTA written to disk stays
    traceable to the user's input.
    """
    lines = []
    for sequence, ids in accepted:
        lines.append(f">{ids[0]}")
        for i in range(0, len(sequence), 80):
            lines.append(sequence[i:i + 80])
    return "\n".join(lines) + "\n"
