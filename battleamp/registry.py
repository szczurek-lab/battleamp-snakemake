"""Enumerate the model variants the pipeline can run.

A "variant" is the unit the pipeline actually produces predictions for. Most
models have exactly one (variant name == model name), but multi-output models
like APEX emit several from a single inference pass, each with its own name,
type and output file.

This mirrors the variant expansion in workflow/Snakefile so that the web UI
offers exactly the columns the pipeline will produce.
"""

from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
MODELS_DIR = REPO_ROOT / "models"
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"


def _load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def enabled_models():
    """Model names enabled in config/config.yaml, in declared order."""
    return list(_load_config()["models"])


def load_variants(model_names=None):
    """Return {variant_name: info} for the given models (default: all enabled).

    Each info dict carries: model, variant, type, length_min, length_max,
    gpu_required, framework, multioutput.
    """
    if model_names is None:
        model_names = enabled_models()

    variants = {}
    for model_name in model_names:
        meta_path = MODELS_DIR / model_name / "model.yaml"
        if not meta_path.exists():
            raise FileNotFoundError(
                f"Model metadata not found: {meta_path}. "
                f"Each model directory must contain a model.yaml file."
            )
        with open(meta_path) as f:
            meta = yaml.safe_load(f)

        base = {
            "model": model_name,
            "type": meta["type"],
            "length_min": meta.get("length_min"),
            "length_max": meta.get("length_max"),
            "gpu_required": meta.get("gpu_required", False),
            "framework": meta.get("framework"),
            "multioutput": bool(meta.get("multioutput")),
        }

        if meta.get("variants"):
            # A multi-output model runs once and writes files for every variant,
            # so two of its artefacts do NOT live under the variant name (see
            # workflow/rules/inference.smk):
            #   - the prefilter manifest is written under the FIRST variant,
            #     whose FASTA all the others reuse;
            #   - the inference log is written under the MODEL name.
            # Getting either wrong silently loses length-exclusion counts or
            # leaves operators with no log for a failed model.
            primary = meta["variants"][0]["name"]
            for v in meta["variants"]:
                name = v["name"]
                variants[name] = {
                    **base,
                    "variant": name,
                    # Multi-output models may override type per variant, e.g.
                    # SenseXAMP ships one classifier and two regressors.
                    "type": v.get("type", base["type"]),
                    "prefilter_variant": primary if base["multioutput"] else name,
                    "log_name": model_name if base["multioutput"] else name,
                }
        else:
            variants[model_name] = {
                **base,
                "variant": model_name,
                "prefilter_variant": model_name,
                "log_name": model_name,
            }

    return variants


def resolve_models(requested):
    """Validate a user's model selection against config/config.yaml.

    ``requested`` may be None (meaning all enabled models), a comma-separated
    string, or a list. Returns the list of model names.

    Raises ValueError naming the unknown entries and the valid options, so the
    message can be shown to the user directly.
    """
    known = enabled_models()
    if requested is None:
        return known

    if isinstance(requested, str):
        names = [m.strip() for m in requested.split(",") if m.strip()]
    else:
        names = [str(m).strip() for m in requested if str(m).strip()]

    if not names:
        return known

    unknown = [m for m in names if m not in known]
    if unknown:
        raise ValueError(
            f"Unknown model(s): {', '.join(unknown)}. "
            f"Available: {', '.join(known)}."
        )
    return names


def variants_for_models(model_names):
    """Variant names produced by the given models, sorted for stable output."""
    return sorted(load_variants(model_names).keys())
