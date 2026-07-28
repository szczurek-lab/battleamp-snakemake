"""BattleAMP web-service API.

    import battleamp
    battleamp.list_models()                 -> JSON str   (fast)
    battleamp.validate(fasta_text)          -> JSON str   (fast)
    battleamp.score(fasta_text, models=...) -> JSON str   (SLOW, worker only)
    battleamp.to_tsv(result_json)           -> TSV str

See battleamp/api.py for the timing constraints that decide where each one may
be called from.
"""

from .api import list_models, score, to_tsv, validate

__all__ = ["list_models", "validate", "score", "to_tsv"]
