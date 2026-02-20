"""
Aggregate Snakemake benchmark logs into resource usage tables.

Reads benchmark TSVs produced by Snakemake's ``benchmark:`` directive
from the directory tree under ``snakemake.params.benchmark_dir``:

    benchmarks/inference/{variant}/{dataset}.tsv
    benchmarks/setup/{model}.tsv

Produces two output files:

    resource_usage.tsv    -- one row per (stage, variant, dataset)
    resource_summary.tsv  -- per-variant totals (inference only)
"""

import sys
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def read_benchmark(path: Path) -> dict | None:
    """Read a single Snakemake benchmark TSV and return key columns."""
    try:
        df = pd.read_csv(path, sep="\t")
        row = df.iloc[0]
        result = {
            "s": row["s"],
            "max_rss": row.get("max_rss"),
            "max_vms": row.get("max_vms"),
            "max_uss": row.get("max_uss"),
            "cpu_time": row.get("cpu_time"),
        }
        return result
    except Exception as e:
        log.warning("Failed to read %s: %s", path, e)
        return None


def main():
    benchmark_dir = Path(snakemake.params.benchmark_dir)
    rows = []

    # Inference benchmarks: benchmarks/inference/{variant}/{dataset}.tsv
    inference_dir = benchmark_dir / "inference"
    if inference_dir.exists():
        for tsv in sorted(inference_dir.glob("*/*.tsv")):
            variant = tsv.parent.name
            dataset = tsv.stem
            data = read_benchmark(tsv)
            if data is not None:
                rows.append({"stage": "inference", "variant": variant,
                             "dataset": dataset, **data})
    else:
        log.warning("No inference benchmark directory: %s", inference_dir)

    # Setup benchmarks: benchmarks/setup/{model}.tsv
    setup_dir = benchmark_dir / "setup"
    if setup_dir.exists():
        for tsv in sorted(setup_dir.glob("*.tsv")):
            model = tsv.stem
            data = read_benchmark(tsv)
            if data is not None:
                rows.append({"stage": "setup", "variant": model,
                             "dataset": "-", **data})
    else:
        log.warning("No setup benchmark directory: %s", setup_dir)

    if not rows:
        log.error("No benchmark files found under %s", benchmark_dir)
        sys.exit(1)

    # Full per-run table
    usage = pd.DataFrame(rows)
    usage.to_csv(snakemake.output.usage, sep="\t", index=False)
    log.info("Wrote %d rows to %s", len(usage), snakemake.output.usage)

    # Per-variant summary (inference only)
    inf = usage[usage["stage"] == "inference"]
    summary = (
        inf.groupby("variant")
        .agg(
            total_wall_s=("s", "sum"),
            max_peak_rss_mb=("max_rss", "max"),
            n_datasets=("dataset", "count"),
            mean_wall_s=("s", "mean"),
        )
        .sort_values("total_wall_s", ascending=False)
        .reset_index()
    )
    summary.to_csv(snakemake.output.summary, sep="\t", index=False)
    log.info("Wrote %d variants to %s", len(summary), snakemake.output.summary)


main()
