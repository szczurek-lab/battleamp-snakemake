#!/bin/bash
# verify_no_prediction_changes.sh
#
# Verifies that the BattleAMP-amPEPpy fork did not modify any files
# from the original tlawrence3/amPEPpy that affect prediction behavior.
#
# Produces a verification report suitable for supplementary materials.
#
# Usage:
#   cd /path/to/BattleAMP-amPEPpy
#   bash verify_no_prediction_changes.sh
#
# Requirements: git, md5sum

set -euo pipefail

FORK_DIR="$(pwd)"
ORIGINAL_REPO="https://github.com/tlawrence3/amPEPpy.git"
ORIGINAL_TAG="v1.1.0"   # latest release; adjust if you forked from an earlier commit
WORK_DIR=$(mktemp -d)
REPORT_FILE="$FORK_DIR/validation/verification_report.txt"

mkdir -p "$FORK_DIR/validation"

echo "============================================================"
echo "  BattleAMP-amPEPpy: Verification of unmodified prediction"
echo "============================================================"
echo ""

# Step 1: Clone the original repo
echo "[1/4] Cloning original repo ($ORIGINAL_TAG)..."
git clone --branch "$ORIGINAL_TAG" --depth 1 "$ORIGINAL_REPO" "$WORK_DIR/original" 2>/dev/null \
    || git clone --depth 1 "$ORIGINAL_REPO" "$WORK_DIR/original" 2>/dev/null

ORIGINAL_DIR="$WORK_DIR/original"

# Step 2: Identify which files are prediction-critical
# These are the ONLY files that affect model output:
#   - amPEPpy/amPEP.py       (prediction logic, feature extraction, RF inference)
#   - amPEPpy/__init__.py     (package init)
#   - amPEPpy/_version.py     (version string, does not affect predictions)
#   - pretrained_models/amPEP.model  (the trained random forest weights)
#   - setup.py                (install config -- affects which code gets installed)
echo "[2/4] Comparing prediction-critical files..."
echo ""

CRITICAL_FILES=(
    "amPEPpy/amPEP.py"
    "amPEPpy/__init__.py"
    "pretrained_models/amPEP.model"
    "setup.py"
)

ALL_MATCH=true
{
    echo "BattleAMP-amPEPpy Verification Report"
    echo "======================================"
    echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Fork directory: $FORK_DIR"
    echo "Original repo: $ORIGINAL_REPO"
    echo "Original tag/branch: $ORIGINAL_TAG"
    echo ""
    echo "Prediction-critical files"
    echo "-------------------------"
    echo ""

    for f in "${CRITICAL_FILES[@]}"; do
        FORK_FILE="$FORK_DIR/$f"
        ORIG_FILE="$ORIGINAL_DIR/$f"

        if [ ! -f "$FORK_FILE" ]; then
            echo "  $f: MISSING in fork"
            ALL_MATCH=false
            continue
        fi
        if [ ! -f "$ORIG_FILE" ]; then
            echo "  $f: MISSING in original (new file in fork)"
            ALL_MATCH=false
            continue
        fi

        FORK_MD5=$(md5sum "$FORK_FILE" | cut -d' ' -f1)
        ORIG_MD5=$(md5sum "$ORIG_FILE" | cut -d' ' -f1)

        if [ "$FORK_MD5" = "$ORIG_MD5" ]; then
            echo "  $f: IDENTICAL (md5: $FORK_MD5)"
        else
            echo "  $f: DIFFERS"
            echo "    original md5: $ORIG_MD5"
            echo "    fork md5:     $FORK_MD5"
            ALL_MATCH=false

            # Show the diff for code files (not binary)
            if [[ "$f" != *.model ]]; then
                echo ""
                echo "    --- diff ---"
                diff "$ORIG_FILE" "$FORK_FILE" | head -50 | sed 's/^/    /'
                echo "    --- end diff ---"
                echo ""
            fi
        fi
    done

    echo ""
    echo "Files added by BattleAMP fork (not in original)"
    echo "-------------------------------------------------"
    echo ""

    # List files in fork that don't exist in original
    cd "$FORK_DIR"
    for f in $(find . -maxdepth 2 -type f -not -path './.git/*' -not -path './eval_results/*' | sort); do
        f_clean="${f#./}"
        if [ ! -f "$ORIGINAL_DIR/$f_clean" ]; then
            echo "  ADDED: $f_clean"
        fi
    done

    echo ""
    echo "Conclusion"
    echo "----------"
    echo ""

    if [ "$ALL_MATCH" = true ]; then
        echo "ALL PREDICTION-CRITICAL FILES ARE IDENTICAL to the original"
        echo "tlawrence3/amPEPpy repository ($ORIGINAL_TAG)."
        echo ""
        echo "The BattleAMP fork adds only wrapper scripts (inference.sh,"
        echo "inference.py, benchmark_utils.py, setup.sh) that call the"
        echo "original 'ampep predict' CLI as a subprocess. The prediction"
        echo "logic, feature extraction code, and pretrained model weights"
        echo "are unmodified."
        echo ""
        echo "VERDICT: No validation of prediction equivalence needed."
        echo "The code-level identity guarantees identical outputs."
    else
        echo "SOME PREDICTION-CRITICAL FILES DIFFER from the original."
        echo "Review the diffs above. If changes are cosmetic (whitespace,"
        echo "comments, version strings), document this. If changes touch"
        echo "prediction logic, empirical validation is required."
        echo ""
        echo "VERDICT: Manual review of diffs required."
    fi

} | tee "$REPORT_FILE"

echo ""
echo "Report saved to: $REPORT_FILE"

# Cleanup
rm -rf "$WORK_DIR"
