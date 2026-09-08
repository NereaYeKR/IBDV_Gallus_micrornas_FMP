#!/bin/bash

# Create Bowtie2 index from a reference genome
# Usage: ./01.create_bowtie2_index.sh /path/to/genome.fa
# Example: ./01.create_bowtie2_index.sh Gallus_gallus.fa
#
# USAGE (individual use):
#       bash 01.create_bowtie2_index.sh /path/to/genome.fa
#
# EXAMPLE (individual use):
#       bash 01.create_bowtie2_index.sh data/genome.fa
#
# NOTE: This script requires Bowtie2 installed
# Install: conda install -c bioconda bowtie2


source config.sh 2>/dev/null

genome=${1:-$GENOME}
outdir=${2:-$SUBREAD_INDEX_DIR}

if [ -z "$outdir" ]; then
    outdir="index_bowtie2"
fi

mkdir -p "$outdir"

filename=$(basename "$genome")
index_name=${filename%.*}

echo "Building Bowtie2 index..."
echo "Genome: $genome"
echo "Output: $outdir"

bowtie2-build "$genome" "$outdir/$index_name"

echo "Done!"
