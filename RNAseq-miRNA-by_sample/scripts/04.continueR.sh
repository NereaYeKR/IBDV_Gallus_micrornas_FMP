#!/bin/bash
# ==============================================================================
# 04 - R analysis launcher (downstream analysis)
#
# Description:
#   Loads the variables from config.sh (exporting them so that the R script
#   can read them via Sys.getenv), checks that the sample metadata table
#   (SampleInfo.tab) exists, and runs the main R script (04.continueR.R),
#   aborting execution if it fails.
#
# Usage:
#   bash scripts/04.continueR.sh
# ==============================================================================

# Load the configuration and export the variables so R can read them
set -a
[ -f config.sh ] && source config.sh
set +a

echo "=========================================="
echo "STEP 4 - Running R Downstream Analysis"
echo "=========================================="

if [ ! -f "$SAMPLE_INFO" ]; then
    echo "Error: SampleInfo.tab not found in $SAMPLE_INFO"
    exit 1
fi

# Run the R script and abort if it fails
Rscript scripts/04.continueR.R || { echo "Error: R script failed. Aborting."; exit 1; }

echo "R Analysis completed. All files are in $RESULTS"
