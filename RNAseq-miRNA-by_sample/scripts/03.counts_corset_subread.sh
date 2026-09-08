#!/bin/bash

# Corset analysis for mature miRNA
# Usage:
#   bash 03.counts_corset_subread.sh [input_bam_dir] [output_dir]
#
# Example (individual):
#   bash 03.counts_corset_subread.sh results/alig_subread results/corset_mature

# Load configuration for environment variables
[ -f config.sh ] && source config.sh

# Assign arguments with defaults
BAMS_DIR=${1:-$ALIGN_DIR}
OUT_DIR=${2:-"results/corset_mature"}

# Ensure Corset is available 
CORSET_BIN="/home/jr/contrib/corset/corset"

mkdir -p "$OUT_DIR"

# Collect BAM files
BAMS=$(ls -1 "$BAMS_DIR"/*.bam 2>/dev/null)

if [ -z "$BAMS" ]; then
    echo "Error: No BAM files found in $BAMS_DIR"
    exit 1
fi

echo "Running Corset on mature miRNAs..."

$CORSET_BIN \
    -p "$OUT_DIR/mature" \
    -m 1 \
    $BAMS

echo ">>> Corset analysis for mature miRNAs completed."
echo ">>> Results: $OUT_DIR"
