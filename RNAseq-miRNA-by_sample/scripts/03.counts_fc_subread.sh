#!/bin/bash
# ==============================================================================
# 03 - miRNA quantification with featureCounts (POSIX)
#
# Description:
#   Runs featureCounts on the BAM files aligned in step 02, using the
#   mature miRNA GFF annotation and allowing fractional assignment of
#   multi-mapped reads (-M --fraction). Produces the original
#   featureCounts count table and a "clean" version (without comment
#   headers and with simplified sample names) ready for the R analysis.
#
# Usage:
#   bash scripts/03.counts_fc_subread.sh [alignment_dir] [gff] [output_file]
# ==============================================================================

# Load the configuration for the environment variables (ALIGN_DIR, MATURE_GFF, etc.)
[ -f config.sh ] && source config.sh

# Assign arguments, with defaults taken from config.sh
ALN_DIR=${1:-$ALIGN_DIR}
GTF_FILE=${2:-$MATURE_GFF}
OUTPUT=${3:-$COUNTS_DIR/counts.txt}

CLEAN_OUTPUT="$(dirname "$OUTPUT")/mat_count_fc_clean.txt"

# Make sure the output directory exists
mkdir -p "$(dirname "$OUTPUT")"

echo "Running featureCounts with 32 threads..."

# 1. Run featureCounts
featureCounts -T 32 \
              -a "$GTF_FILE" \
              -t miRNA \
              -g Name \
              -M --fraction \
              -o "$OUTPUT" \
              "$ALN_DIR"/*.bam

echo "Cleaning the output..."

TEMP_FILE="$(dirname "$OUTPUT")/temp_clean.txt"

# 2. Remove the featureCounts comment line and keep the miRNA ID and the
#    per-sample count columns
grep -v "^#" "$OUTPUT" | cut -f 1,7- > "$TEMP_FILE"

# 3. Simplify the header (remove path and suffix from the BAM names) and
#    append it back to the count rows
head -n 1 "$TEMP_FILE" | sed 's|results/alig_subread/||g' | sed 's|_sorted.bam||g' > "$CLEAN_OUTPUT"
tail -n +2 "$TEMP_FILE" >> "$CLEAN_OUTPUT"

rm "$TEMP_FILE"

echo ">>> Process completed."
echo ">>> Original file: $OUTPUT"
echo ">>> Cleaned file: $CLEAN_OUTPUT"
