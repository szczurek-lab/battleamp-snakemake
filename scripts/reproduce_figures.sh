#!/usr/bin/env bash
# Regenerate every figure in the paper from the result tables tracked in this
# repository. Model submodules, conda environments and a GPU are not required.
#
#   pip install -r requirements-figures.txt
#   bash scripts/reproduce_figures.sh
#
# Figures are written to figures/ as .png and .pdf.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Render to file rather than to a display, so the scripts also run headless.
export MPLBACKEND=Agg

FIGURES=(
    figure2
    figure3
    figure4
    figure5
    figure_supp_mic_distributions
    figure_supp_mic_scatter
    figure_supp_overlap
    figure_supp_sankey
    figure_supp_synthetic
    figure_supp_units
)

pass=0; fail=0; failed_names=()
for f in "${FIGURES[@]}"; do
    printf '%-34s ' "$f"
    if out=$(python3 "notebooks/${f}.py" 2>&1); then
        printf 'ok\n'; pass=$((pass + 1))
    else
        printf 'FAILED\n'; fail=$((fail + 1)); failed_names+=("$f")
        printf '%s\n' "$out" | tail -5 | sed 's/^/    /'
    fi
done

echo
echo "${pass} figures regenerated, ${fail} failed."

if [ "$fail" -gt 0 ]; then
    echo "Failed: ${failed_names[*]}" >&2
    echo "Check that requirements-figures.txt is installed." >&2
    exit 1
fi

cat <<'MSG'

Figures written to figures/ as .png and .pdf.
The tables they are drawn from are in results/aggregated/ and
results/generalization/. See the README for their contents.
MSG
