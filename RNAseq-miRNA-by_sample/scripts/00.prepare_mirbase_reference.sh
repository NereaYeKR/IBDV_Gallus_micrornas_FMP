#!/bin/bash
# ==============================================================================
# 00 - Prepare a miRBase mature reference for a given species
#
# Description:
#   Downloads (if not already present) the miRBase mature-sequence file,
#   cleans up the headers and non-standard bases, filters the sequences by
#   the given organism prefix (e.g. "gga" for Gallus gallus), and produces
#   both the RNA version and its DNA conversion (U -> T), ready to be used
#   as an alignment reference.
#
# Requirements:
#   - seqkit and wget available in the PATH.
#
# Usage:
#   bash scripts/00.prepare_mirbase_reference.sh <species_code> [output_dir]
#
# Arguments:
#   $1 = miRBase 3-letter organism code (e.g. gga, hsa, mmu)
#   $2 = output directory (default: "data")
#
# Example:
#   bash scripts/00.prepare_mirbase_reference.sh gga data
# ==============================================================================

set -e

SPECIES=$1
OUTDIR=${2:-data}

if [ -z "$SPECIES" ]; then
    echo "ERROR: no organism code provided"
    echo "Usage: $0 <species_code> [output_dir]"
    echo "Example: $0 gga data"
    exit 1
fi

mkdir -p "$OUTDIR"

echo "Preparing miRNA reference for: $SPECIES"

# Temporary files
TMP_MATURE="$OUTDIR/mature.fa"
TMP_CLEAN="$OUTDIR/mature_clean.fa"

# Download the miRBase mature sequences if they do not exist yet
if [ ! -f "$TMP_MATURE" ]; then
    echo "Downloading mature.fa from miRBase..."
    wget --no-check-certificate \
        https://www.mirbase.org/download/mature.fa \
        -O "$TMP_MATURE"
else
    echo "Using existing file $TMP_MATURE"
fi

echo "Cleaning sequences..."

# Replace any character other than A/U/G/C with N
sed '/^[^>]/ s/[^AUGCaugc]/N/g' \
    "$TMP_MATURE" > "$TMP_CLEAN"

# Keep only the miRNA identifier in the FASTA headers
sed -i 's/\s.*//' "$TMP_CLEAN"

echo "Extracting ${SPECIES} miRNAs..."

# Extract the sequences specific to the given organism
seqkit grep \
    -r \
    --pattern "^${SPECIES}-" \
    "$TMP_CLEAN" \
    > "$OUTDIR/mature-${SPECIES}.fa"

echo "Converting RNA to DNA..."

# Convert U -> T
seqkit seq \
    --rna2dna \
    "$OUTDIR/mature-${SPECIES}.fa" \
    > "$OUTDIR/mature-${SPECIES}-DNA.fa"

echo "Reference created successfully"
echo "RNA reference:"
echo "  $OUTDIR/mature-${SPECIES}.fa"
echo "DNA reference:"
echo "  $OUTDIR/mature-${SPECIES}-DNA.fa"
