#!/bin/bash
# ==============================================================================
# 02 - Convert FASTA to GFF with organism filtering
#
# Description:
#   Converts a FASTA file (e.g. the mature miRBase sequences) into a simple
#   GFF file, keeping only the sequences whose identifier matches the given
#   organism prefix (e.g. gga, hsa, mmu). This GFF is subsequently used as
#   the annotation for subread-align and featureCounts.
#
# Usage (standalone):
#   bash scripts/02.fasta_to_gff.sh /path/to/input.fa organism_prefix [output.gff]
#
# Example (standalone):
#   bash scripts/02.fasta_to_gff.sh data/mature.fa gga
# ==============================================================================

file=$1
org=$2
out=$3

# Check the required arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 input.fa organism_prefix [output.gff]"
    exit 1
fi

# Check that the input file exists
if [ ! -f "$file" ]; then
    echo "Error: input file not found: $file"
    exit 1
fi

# If no output path is given, save it next to the input file
if [ -z "$out" ]; then
    dir=$(dirname "$file")
    base=$(basename "$file")
    base="${base%.*}"
    out="$dir/${base}_${org}.gff"
fi

echo "Input file: $file"
echo "Organism filter: $org"
echo "Output file: $out"

# For every FASTA sequence whose ID starts with the organism prefix, write
# a GFF line of type "miRNA" spanning the full length of the mature sequence.
awk -v org="$org" '
/^>/ {
    id=$1
    sub(">", "", id)

    getline seq

    if (id ~ "^" org) {
        print id "\t.\tmiRNA\t1\t" length(seq) "\t.\t+\t.\tID=" id ";Name=" id
    }
}
' "$file" > "$out"

echo "Done! Filtered GFF created: $out"
