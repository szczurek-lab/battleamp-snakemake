#!/usr/bin/env bash
# Download and unpack the cached BattleAMP predictions from Zenodo.
# The benchmark's figures and tables regenerate without this archive: they read
# the aggregated tables and figure caches tracked in the repository. Fetch this
# only to recompute metrics from raw per-model predictions.
#
#
# Zenodo record: https://doi.org/10.5281/zenodo.22661234 (all versions)
#
# Usage:  bash scripts/fetch_predictions.sh
set -euo pipefail

ZENODO_RECORD="${ZENODO_RECORD:-22661235}"
ARCHIVE="battleamp-predictions-v1.tar.gz"
SHA256="b3faa3e9221596ecf70713f80981599acea2234801ded2694f84bdf8ee5d9b64"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ "$ZENODO_RECORD" == REPLACE_* ]]; then
    echo "ERROR: Zenodo record ID not set." >&2
    echo "Edit ZENODO_RECORD in this script, or run:" >&2
    echo "  ZENODO_RECORD=<id> bash scripts/fetch_predictions.sh" >&2
    exit 1
fi

URL="https://zenodo.org/records/${ZENODO_RECORD}/files/${ARCHIVE}?download=1"

if [ ! -f "$ARCHIVE" ]; then
    echo "Downloading ${ARCHIVE} (309 MB) from Zenodo record ${ZENODO_RECORD}..."
    curl -L --fail --progress-bar -o "${ARCHIVE}.part" "$URL"
    mv "${ARCHIVE}.part" "$ARCHIVE"
else
    echo "${ARCHIVE} already present, skipping download."
fi

echo "Verifying checksum..."
echo "${SHA256}  ${ARCHIVE}" | sha256sum -c -

echo "Unpacking into ${REPO_ROOT}/results/inference/ ..."
tar -xzf "$ARCHIVE" -C "$REPO_ROOT"

cat <<'MSG'

Done. 292 prediction files unpacked under results/inference/.

To recompute metrics without rerunning inference:

    snakemake --profile profile/ --touch
    snakemake --profile profile/ --forcerun evaluate_task

The --touch is required. Unpacked files carry archive timestamps, so Snakemake
would otherwise treat them as stale and rerun inference.
MSG
