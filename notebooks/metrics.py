#!/usr/bin/env python3
"""
metrics.py
==========

Standalone module for computing sequence diversity and similarity metrics
used across BATTLE-AMP figures. Wraps CD-HIT (coverage) and MMseqs2
(alignment-based similarity) with a consistent interface.

All functions accept plain lists of sequences and return either a single
scalar (coverage) or a numpy array (per-sequence scores).

Dependencies:
    pip install numpy pandas biopython
    External: cd-hit, mmseqs (both on PATH)

Usage:
    from metrics import compute_coverage, compute_self_similarity, compute_similarity_to_reference

    cov = compute_coverage(sequences)
    self_sim = compute_self_similarity(sequences, n=1000)
    ref_sim = compute_similarity_to_reference(query_sequences, reference_sequences, n=1000)
"""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# ── Defaults (can be overridden per call) ─────────────────────────────────

CDHIT_EXEC = "cd-hit"
MMSEQS_EXEC = "mmseqs"

CDHIT_IDENTITY = 0.6
CDHIT_WORDSIZE = 5
MMSEQS_EVALUE = 1000


# ══════════════════════════════════════════════════════════════════════════
# INTERNAL HELPERS
# ══════════════════════════════════════════════════════════════════════════

def _write_fasta(sequences: List[str], path: Path, prefix: str = "seq"):
    """Write sequences to a FASTA file."""
    with open(path, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">{prefix}_{i}\n{seq}\n")


def _resolve_exec(name: str) -> Optional[str]:
    """Check whether an executable is available on PATH."""
    return shutil.which(name)


def _subsample(sequences: List[str], n: Optional[int], replace: bool = False) -> List[str]:
    """Subsample sequences if n is set and smaller than the total."""
    if n is None or n >= len(sequences):
        return sequences
    return list(np.random.choice(sequences, n, replace=replace))


# ══════════════════════════════════════════════════════════════════════════
# COVERAGE (CD-HIT)
# ══════════════════════════════════════════════════════════════════════════

def compute_coverage(
    sequences: List[str],
    identity: float = CDHIT_IDENTITY,
    wordsize: int = CDHIT_WORDSIZE,
    cd_hit_exec: str = CDHIT_EXEC,
) -> float:
    """
    Compute sequence coverage using CD-HIT clustering.

    Coverage = n_clusters / n_sequences

    A value near 1.0 indicates high internal diversity (each sequence in
    its own cluster). Values near 0 indicate high redundancy.

    Parameters
    ----------
    sequences : list of str
        Amino acid sequences.
    identity : float
        CD-HIT global sequence identity threshold (default 0.6).
    wordsize : int
        CD-HIT word length parameter ``-n`` (default 5; valid for identity >= 0.6).
    cd_hit_exec : str
        Path or name of the cd-hit executable.

    Returns
    -------
    float
        Coverage in [0, 1], or NaN on failure.
    """
    if len(sequences) == 0:
        return float("nan")

    if _resolve_exec(cd_hit_exec) is None:
        log.error(f"'{cd_hit_exec}' not found on PATH. Install CD-HIT.")
        return float("nan")

    with TemporaryDirectory(prefix="cdhit_") as tmp:
        tmpdir = Path(tmp)
        inp = tmpdir / "input.fasta"
        out = tmpdir / "clustered.fasta"

        _write_fasta(sequences, inp)

        cmd = [
            cd_hit_exec,
            "-i", str(inp),
            "-o", str(out),
            "-c", str(identity),
            "-n", str(wordsize),
            "-T", "0",      # all threads
            "-M", "0",      # unlimited memory
            "-d", "0",      # unlimited description length
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            log.error(f"cd-hit executable not found: {cd_hit_exec}")
            return float("nan")
        except subprocess.CalledProcessError as e:
            log.error(f"cd-hit failed (exit {e.returncode}): {e.stderr[:500]}")
            return float("nan")

        clstr_path = Path(str(out) + ".clstr")
        if not clstr_path.exists():
            log.error("cd-hit .clstr file not found after run")
            return float("nan")

        n_clusters = 0
        with open(clstr_path) as f:
            for line in f:
                if line.startswith(">Cluster"):
                    n_clusters += 1

        coverage = n_clusters / len(sequences)
        log.info(
            f"CD-HIT coverage: {n_clusters} clusters / {len(sequences)} sequences "
            f"= {coverage:.4f} (identity={identity})"
        )
        return coverage


# ══════════════════════════════════════════════════════════════════════════
# SIMILARITY (MMseqs2)
# ══════════════════════════════════════════════════════════════════════════

def _run_mmseqs_search(
    query_fasta: Path,
    target_fasta: Path,
    out_tsv: Path,
    tmpdir: Path,
    evalue: float = MMSEQS_EVALUE,
    mmseqs_exec: str = MMSEQS_EXEC,
    extra_args: Optional[List[str]] = None,
) -> bool:
    """Run mmseqs easy-search. Returns True on success."""
    cmd = [
        mmseqs_exec, "easy-search",
        str(query_fasta), str(target_fasta), str(out_tsv),
        str(tmpdir / "mmseqs_tmp"),
        "-v", "0",
        "-e", str(evalue),
    ]
    if extra_args:
        cmd += extra_args

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except FileNotFoundError:
        log.error(f"mmseqs not found: {mmseqs_exec}")
        return False
    except subprocess.CalledProcessError as e:
        log.error(f"mmseqs failed (exit {e.returncode}): {e.stderr[:500]}")
        return False


def _safe_load_tsv(path: Path) -> pd.DataFrame:
    """Load a TSV, returning an empty DataFrame if the file is missing or empty."""
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", header=None)


def _compute_normalized_bitscores(
    query_sequences: List[str],
    reference_sequences: List[str],
    evalue: float = MMSEQS_EVALUE,
    mmseqs_exec: str = MMSEQS_EXEC,
) -> np.ndarray:
    """
    Core MMseqs2 computation: for each query, compute normalised bitscore
    against a reference set.

    Normalised bitscore = best_hit_bitscore / self_bitscore

    Sequences identical to a reference entry receive score 1.0 without
    running alignment. Sequences with no hit receive score 0.0.

    Returns
    -------
    np.ndarray of shape (len(query_sequences),) with values in [0, 1].
    """
    from Bio import SeqIO
    from Bio.Seq import Seq

    n = len(query_sequences)
    if n == 0:
        return np.array([])

    results = np.zeros(n)

    with TemporaryDirectory(prefix="mmseqs_") as tmp:
        tmpdir = Path(tmp)

        # Write FASTAs
        query_fasta = tmpdir / "query.fasta"
        records = [
            SeqIO.SeqRecord(seq=Seq(s), id=f"q{i}", description="")
            for i, s in enumerate(query_sequences)
        ]
        SeqIO.write(records, query_fasta, "fasta")

        ref_fasta = tmpdir / "reference.fasta"
        ref_records = [
            SeqIO.SeqRecord(seq=Seq(s), id=f"r{i}", description="")
            for i, s in enumerate(reference_sequences)
        ]
        SeqIO.write(ref_records, ref_fasta, "fasta")

        # Exact matches get 1.0 immediately
        ref_set = set(reference_sequences)
        qid_to_idx = {}
        filtered_records = []
        for i, s in enumerate(query_sequences):
            qid = f"q{i}"
            qid_to_idx[qid] = i
            if s in ref_set:
                results[i] = 1.0
            else:
                filtered_records.append(records[i])

        if not filtered_records:
            return results

        filtered_fasta = tmpdir / "filtered.fasta"
        SeqIO.write(filtered_records, filtered_fasta, "fasta")

        # Query vs reference
        db_tsv = tmpdir / "vs_ref.tsv"
        if not _run_mmseqs_search(
            filtered_fasta, ref_fasta, db_tsv, tmpdir,
            evalue=evalue, mmseqs_exec=mmseqs_exec,
        ):
            return results

        df_db = _safe_load_tsv(db_tsv)
        if df_db.empty:
            return results

        bit_col = df_db.columns[-1]
        idx = df_db.groupby(0)[bit_col].idxmax()
        best_db = df_db.loc[idx, [0, bit_col]].copy()
        best_db.columns = ["qid", "bitscore"]

        # Self-alignment for normalisation
        mmseqs_tmp = tmpdir / "mmseqs_tmp"
        if mmseqs_tmp.exists():
            shutil.rmtree(mmseqs_tmp)

        self_tsv = tmpdir / "self.tsv"
        if not _run_mmseqs_search(
            filtered_fasta, filtered_fasta, self_tsv, tmpdir,
            evalue=evalue, mmseqs_exec=mmseqs_exec,
            extra_args=["--add-self-matches", "--max-seqs", "1"],
        ):
            return results

        df_self = _safe_load_tsv(self_tsv)
        if df_self.empty:
            return results

        idx_s = df_self.groupby(0)[df_self.columns[-1]].idxmax()
        best_self = df_self.loc[idx_s, [0, df_self.columns[-1]]].copy()
        best_self.columns = ["qid", "self_bitscore"]

        merged = best_db.merge(best_self, on="qid")
        merged["norm"] = (merged["bitscore"] / merged["self_bitscore"]).clip(0, 1)

        for _, row in merged.iterrows():
            if row["qid"] in qid_to_idx:
                results[qid_to_idx[row["qid"]]] = row["norm"]

    return results


def compute_self_similarity(
    sequences: List[str],
    n: Optional[int] = 1000,
    evalue: float = MMSEQS_EVALUE,
    mmseqs_exec: str = MMSEQS_EXEC,
) -> np.ndarray:
    """
    Compute within-dataset self-similarity via MMseqs2.

    Each (subsampled) sequence is aligned against the full dataset.
    Returns per-sequence normalised bitscores.

    Parameters
    ----------
    sequences : list of str
        The full dataset.
    n : int or None
        Number of query sequences to subsample (for speed). The full
        dataset is always used as the reference. None = no subsampling.
    evalue : float
        MMseqs2 e-value threshold.
    mmseqs_exec : str
        Path or name of the mmseqs executable.

    Returns
    -------
    np.ndarray of shape (n_subsampled,) with values in [0, 1].
    """
    if _resolve_exec(mmseqs_exec) is None:
        log.error(f"'{mmseqs_exec}' not found on PATH. Install MMseqs2.")
        return np.array([])

    sample = _subsample(sequences, n)
    log.info(f"Self-similarity: {len(sample)} queries vs {len(sequences)} reference sequences")
    return _compute_normalized_bitscores(sample, sequences, evalue=evalue, mmseqs_exec=mmseqs_exec)


def compute_similarity_to_reference(
    query_sequences: List[str],
    reference_sequences: List[str],
    n: Optional[int] = 1000,
    evalue: float = MMSEQS_EVALUE,
    mmseqs_exec: str = MMSEQS_EXEC,
) -> np.ndarray:
    """
    Compute similarity of query sequences to a reference dataset via MMseqs2.

    Parameters
    ----------
    query_sequences : list of str
        Sequences to evaluate.
    reference_sequences : list of str
        Reference dataset (e.g. GeneralActivity positives).
    n : int or None
        Subsample query sequences for speed. None = use all.
    evalue : float
        MMseqs2 e-value threshold.
    mmseqs_exec : str
        Path or name of the mmseqs executable.

    Returns
    -------
    np.ndarray of shape (n_subsampled,) with values in [0, 1].
    """
    if _resolve_exec(mmseqs_exec) is None:
        log.error(f"'{mmseqs_exec}' not found on PATH. Install MMseqs2.")
        return np.array([])

    sample = _subsample(query_sequences, n)
    log.info(
        f"Similarity to reference: {len(sample)} queries vs "
        f"{len(reference_sequences)} reference sequences"
    )
    return _compute_normalized_bitscores(
        sample, reference_sequences, evalue=evalue, mmseqs_exec=mmseqs_exec,
    )


# ══════════════════════════════════════════════════════════════════════════
# CONVENIENCE: compute all metrics for a collection of datasets
# ══════════════════════════════════════════════════════════════════════════

def compute_all_metrics(
    datasets: dict,
    reference_key: str = "GeneralActivity",
    similarity_n: int = 1000,
    cdhit_identity: float = CDHIT_IDENTITY,
    cdhit_wordsize: int = CDHIT_WORDSIZE,
) -> Tuple[dict, dict, dict]:
    """
    Compute coverage, self-similarity, and similarity-to-reference for
    all datasets in a single call.

    Parameters
    ----------
    datasets : dict[str, list[str]]
        Mapping of dataset name to list of sequences.
    reference_key : str
        Which dataset to use as the reference for similarity-to-reference.
    similarity_n : int
        Subsample size for MMseqs2 queries.
    cdhit_identity : float
        CD-HIT identity threshold.
    cdhit_wordsize : int
        CD-HIT word size.

    Returns
    -------
    coverage : dict[str, float]
    self_similarity : dict[str, list[float]]
    sim_to_reference : dict[str, list[float]]
    """
    coverage = {}
    self_similarity = {}
    sim_to_reference = {}

    ref_seqs = datasets.get(reference_key)
    if ref_seqs is None:
        log.warning(f"Reference dataset '{reference_key}' not found. "
                     "Similarity-to-reference will be empty.")

    for name, seqs in datasets.items():
        log.info(f"--- {name} ({len(seqs):,} sequences) ---")

        # Coverage
        coverage[name] = compute_coverage(seqs, identity=cdhit_identity, wordsize=cdhit_wordsize)

        # Self-similarity
        scores = compute_self_similarity(seqs, n=similarity_n)
        self_similarity[name] = scores.tolist()

        # Similarity to reference
        if ref_seqs is not None:
            scores = compute_similarity_to_reference(seqs, ref_seqs, n=similarity_n)
            sim_to_reference[name] = scores.tolist()

    return coverage, self_similarity, sim_to_reference


# ══════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(description="Compute sequence diversity metrics")
    parser.add_argument("fasta", help="FASTA file with sequences")
    parser.add_argument("--reference", help="Reference FASTA for similarity-to-reference")
    parser.add_argument("--n", type=int, default=1000, help="Subsample size for similarity")
    parser.add_argument("--identity", type=float, default=0.6, help="CD-HIT identity threshold")
    args = parser.parse_args()

    from Bio import SeqIO

    seqs = [str(r.seq) for r in SeqIO.parse(args.fasta, "fasta")]
    print(f"Loaded {len(seqs)} sequences from {args.fasta}")

    cov = compute_coverage(seqs, identity=args.identity)
    print(f"Coverage: {cov:.4f}")

    self_sim = compute_self_similarity(seqs, n=args.n)
    print(f"Self-similarity: mean={np.mean(self_sim):.4f}, std={np.std(self_sim):.4f}")

    if args.reference:
        ref = [str(r.seq) for r in SeqIO.parse(args.reference, "fasta")]
        print(f"Reference: {len(ref)} sequences from {args.reference}")
        ref_sim = compute_similarity_to_reference(seqs, ref, n=args.n)
        print(f"Similarity to reference: mean={np.mean(ref_sim):.4f}, std={np.std(ref_sim):.4f}")
