#!/bin/bash
# ampeppy setup: install the ampep CLI tool.
# The conda env is already active (managed by Snakemake).
# Do NOT create or activate any conda env here.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Install amPEPpy package (provides the 'ampep' CLI command)
pip install -e . --quiet

# Verify the CLI is available
ampep -h > /dev/null 2>&1 || { echo "ERROR: ampep CLI not found after install" >&2; exit 1; }

# The pretrained model (pretrained_models/amPEP.model) is already in the repo.
# No weight download needed.
echo "ampeppy setup complete" >&2
