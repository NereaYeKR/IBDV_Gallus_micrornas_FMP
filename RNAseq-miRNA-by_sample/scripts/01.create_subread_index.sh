#!/bin/bash
# ==============================================================================
# 01 - Build the Subread index from the reference genome
#
# Description:
#   Generates the Subread index (subread-buildindex) required to align the
#   miRNA-seq reads against the given reference genome.
#
# Requirements:
#   - Subread installed and available in the PATH.
#     Download: http://subread.sourceforge.net/
#   - If no arguments are given, GENOME and SUBREAD_INDEX_DIR are read from
#     config.sh (if present in the working directory).
#
# Usage (standalone):
#   bash scripts/01.create_subread_index.sh /path/to/genome.fa [output_dir]
#
# Example (standalone):
#   bash scripts/01.create_subread_index.sh \
#       data/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa
# ==============================================================================

source config.sh 2>/dev/null

genome=${1:-$GENOME}
outdir=${2:-$SUBREAD_INDEX_DIR}

mkdir -p "$outdir"

filename=$(basename "$genome")
index_name=${filename%.*}

echo "Building Subread index..."
echo "Genome: $genome"
echo "Output: $outdir"

subread-buildindex -F -B -o "$outdir/$index_name" "$genome"

echo "Done!"
